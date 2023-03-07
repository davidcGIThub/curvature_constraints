"""
mavsim_python: drawing tools
    - Beard & McLain, PUP, 2012
    - Update history:
        4/15/2019 - RWB
        3/30/2022 - RWB
"""

import numpy as np
import sys
sys.path.append('..')
from planning.dubins_parameters import DubinsParameters
from message_types.msg_path import MsgPath


class PathManager:
    def __init__(self):
        # message sent to path follower
        self.path = MsgPath()
        # pointers to previous, current, and next waypoints
        self.ptr_previous = 0
        self.ptr_current = 1
        self.ptr_next = 2
        self.num_waypoints = 0
        self.halfspace_n = np.inf * np.ones((3,1))
        self.halfspace_r = np.inf * np.ones((3,1))
        # state of the manager state machine
        self.manager_state = 1
        self.manager_requests_waypoints = True
        self.dubins_path = DubinsParameters()

    def update(self, waypoints, radius, state):
        if waypoints.num_waypoints == 0:
            self.manager_requests_waypoints = True
        if self.manager_requests_waypoints is True \
                and waypoints.flag_waypoints_changed is True:
            self.manager_requests_waypoints = False
        if waypoints.type == 'straight_line':
            self.line_manager(waypoints, state)
        elif waypoints.type == 'fillet':
            self.fillet_manager(waypoints, radius, state)
        elif waypoints.type == 'dubins':
            self.dubins_manager(waypoints, radius, state)
        else:
            print('Error in Path Manager: Undefined waypoint type.')
        return self.path

    def line_manager(self, waypoints, state):
        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        # if the waypoints have changed, update the waypoint pointer
        if waypoints.flag_waypoints_changed is True:
            waypoints.flag_manager_requests_waypoints = False
            waypoints.flag_waypoints_changed = False
            self.num_waypoints = waypoints.num_waypoints
            self.initialize_pointers()
            self.construct_line(waypoints)

        # entered into the half plane separating waypoint segments
        if self.inHalfSpace(mav_pos):
            self.increment_pointers()
            self.construct_line(waypoints)
            # requests new waypoints when reach end of current list
            if self.ptr_current == 0:
                self.manager_requests_waypoints = True

    def fillet_manager(self, waypoints, radius, state):
        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        # if the waypoints have changed, update the waypoint pointer
        if waypoints.flag_waypoints_changed is True:
            #waypoints.flag_manager_requests_waypoints = False
            waypoints.flag_waypoints_changed = False
            self.num_waypoints = waypoints.num_waypoints
            self.initialize_pointers()
            self.construct_fillet_line(waypoints, radius)
            self.manager_state = 1
        # state machine for fillet path
        if self.manager_state == 1:
            # follow straight line path from previous to current
            if self.inHalfSpace(mav_pos):
                # entered into the half plane H1±
                self.construct_fillet_circle(waypoints, radius)
                self.manager_state = 2
        elif self.manager_state == 2:
            # follow start orbit until out of H2
            if not self.inHalfSpace(mav_pos):
                self.manager_state = 3
        elif self.manager_state == 3:
            # follow orbit from previous->current to current->next
             if self.inHalfSpace(mav_pos):
                # entered into the half plane H2
                self.increment_pointers()
                self.construct_fillet_line(waypoints, radius)
                self.manager_state = 1
                # requests new waypoints when reach end of current list
                if self.ptr_current == 0:
                    self.manager_requests_waypoints = True

    def dubins_manager(self, waypoints, radius, state):
        mav_pos = np.array([[state.north, state.east, -state.altitude]]).T
        close_distance = 10
        # if the waypoints have changed, update the waypoint pointer
        if waypoints.flag_waypoints_changed is True:
            waypoints.flag_waypoints_changed = False
            self.num_waypoints = waypoints.num_waypoints
            self.initialize_pointers()
            # dubins path parameters
            self.dubins_path.update(
                ps=waypoints.ned[:, self.ptr_previous:self.ptr_previous+1],
                chis=waypoints.course.item(self.ptr_previous),
                pe=waypoints.ned[:, self.ptr_current:self.ptr_current+1],
                chie=waypoints.course.item(self.ptr_current),
                R=radius)
            self.construct_dubins_circle_start(waypoints, self.dubins_path)
            if np.linalg.norm(self.dubins_path.p_s - self.dubins_path.r1) < close_distance:
                self.construct_dubins_line(waypoints, self.dubins_path)
                self.manager_state = 3
            elif self.inHalfSpace(mav_pos):
                self.manager_state = 1
            else:
                self.manager_state = 2
        # state machine for dubins path
        if self.manager_state == 1:
            # skip the first circle if distance along circle is small
            if ((np.linalg.norm(self.dubins_path.p_s - self.dubins_path.r1) < close_distance)
                    # follow start orbit until out of H1
                    or not self.inHalfSpace(mav_pos)):
                self.manager_state=2
        elif self.manager_state == 2:
            # skip the first circle if distance along circle is small
            if ((np.linalg.norm(self.dubins_path.p_s - self.dubins_path.r1) < close_distance)
                    # follow start orbit until cross into H1
                    or self.inHalfSpace(mav_pos)):
                self.construct_dubins_line(waypoints, self.dubins_path)
                self.manager_state = 3
        elif self.manager_state == 3:
            # skip line if it is short
            if ((np.linalg.norm(self.dubins_path.r1 - self.dubins_path.r2) < close_distance)
                    or self.inHalfSpace(mav_pos)):
                self.construct_dubins_circle_end(waypoints, self.dubins_path)
                if self.inHalfSpace(mav_pos):
                    self.manager_state = 4
                else:
                    self.manager_state = 5
        elif self.manager_state == 4:
            # distance along end orbit is small
            if ((np.linalg.norm(self.dubins_path.r2 - self.dubins_path.p_e) < close_distance)
                    # follow start orbit until out of H3
                    or not self.inHalfSpace(mav_pos)):
                self.manager_state = 5
        elif self.manager_state == 5:
            # skip circle if small
            if ((np.linalg.norm(self.dubins_path.r2 - self.dubins_path.p_e) < close_distance)
                    # follow start orbit until cross into H3
                    or self.inHalfSpace(mav_pos)):
                self.increment_pointers()
                self.dubins_path.update(
                    waypoints.ned[:, self.ptr_previous:self.ptr_previous+1],
                    waypoints.course.item(self.ptr_previous),
                    waypoints.ned[:, self.ptr_current:self.ptr_current+1],
                    waypoints.course.item(self.ptr_current),
                    radius)
                self.construct_dubins_circle_start(waypoints, self.dubins_path)
                self.manager_state = 1
                # requests new waypoints when reach end of current list
                if self.ptr_current == 0:
                    self.manager_requests_waypoints = True

    def initialize_pointers(self):
        if self.num_waypoints >= 3:
            self.ptr_previous = 0
            self.ptr_current = 1
            self.ptr_next = 2
        else:
            print('Error Path Manager: need at least three waypoints')

    def increment_pointers(self):
        self.ptr_previous = self.ptr_current
        self.ptr_current = self.ptr_next
        self.ptr_next = self.ptr_next + 1
        if self.ptr_next > self.num_waypoints-1:
            #self.ptr_next = 9999
            self.ptr_next = 0
        if self.ptr_current > self.num_waypoints-1:
            #self.ptr_current = 9999
            self.ptr_current = 0

    def construct_line(self, waypoints):
        previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous+1]
        if self.ptr_current == 9999:
            current = previous + 100*self.path.line_direction
        else:
            current = waypoints.ned[:, self.ptr_current:self.ptr_current+1]
        if self.ptr_next == 9999:
            next = previous + 200*self.path.line_direction
        else:
            next = waypoints.ned[:, self.ptr_next:self.ptr_next+1]
        self.path.type = 'line'
        self.path.airspeed = waypoints.airspeed.item(self.ptr_current)
        self.path.line_origin = previous
        q_previous = (current - previous)/ np.linalg.norm(current - previous)
        self.path.line_direction = q_previous
        q_next = (next - current) / np.linalg.norm(next - current)
        self.halfspace_n = (q_previous + q_next) / 2
        self.halfspace_n = self.halfspace_n / np.linalg.norm(self.halfspace_n)
        self.halfspace_r = current
        self.path.plot_updated = False

    def construct_fillet_line(self, waypoints, radius):
        previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous+1]
        if self.ptr_current == 9999:
            current = previous + 100*self.path.line_direction
        else:
            current = waypoints.ned[:, self.ptr_current:self.ptr_current+1]
        if self.ptr_next == 9999:
            next = previous + 200*self.path.line_direction
        else:
            next = waypoints.ned[:, self.ptr_next:self.ptr_next+1]
        self.path.type = 'line'
        self.path.airspeed = waypoints.airspeed.item(self.ptr_current)
        self.path.line_origin = previous
        q_previous = (current - previous) / np.linalg.norm(current - previous)
        self.path.line_direction = q_previous
        q_next = (next - current) / np.linalg.norm(next - current)
        beta = np.arccos(-q_previous.T @ q_next)
        self.halfspace_n = q_previous
        self.halfspace_r = current - radius / np.tan(beta/2) * q_previous
        self.path.plot_updated = False

    def construct_fillet_circle(self, waypoints, radius):
        previous = waypoints.ned[:, self.ptr_previous:self.ptr_previous+1]
        if self.ptr_current == 9999:
            current = previous + 100*self.path.line_direction
        else:
            current = waypoints.ned[:, self.ptr_current:self.ptr_current+1]
        if self.ptr_next == 9999:
            next = previous + 200*self.path.line_direction
        else:
            next = waypoints.ned[:, self.ptr_next:self.ptr_next+1]
        self.path.type = 'orbit'
        self.path.airspeed = waypoints.airspeed.item(self.ptr_current)
        q_previous = (current - previous) / np.linalg.norm(current - previous)
        q_next = (next - current) / np.linalg.norm(next - current)
        varrho = np.arccos(-q_previous.T @ q_next)
        q_tmp = (q_previous - q_next) / np.linalg.norm(q_previous - q_next)
        self.path.orbit_center = current - radius / np.sin(varrho/2.0) * q_tmp
        self.path.orbit_radius = radius
        if np.sign(q_previous.item(0)*q_next.item(1) - q_previous.item(1)*q_next.item(0)) > 0:
            self.path.orbit_direction = 'CW'
        else:
            self.path.orbit_direction = 'CCW'
        self.halfspace_n = q_next
        self.halfspace_r = current + radius / np.tan(varrho/2.0) * q_next
        self.path.plot_updated = False

    def construct_dubins_circle_start(self, waypoints, dubins_path):
        self.path.type = 'orbit'
        self.path.airspeed = waypoints.airspeed.item(self.ptr_current)
        self.path.orbit_radius = dubins_path.radius
        self.path.orbit_center = dubins_path.center_s
        if dubins_path.dir_s == 1:
            self.path.orbit_direction = 'CW'
        else:
            self.path.orbit_direction = 'CCW'
        self.halfspace_n = dubins_path.n1
        self.halfspace_r = dubins_path.r1
        self.path.plot_updated = False

    def construct_dubins_line(self, waypoints, dubins_path):
        self.path.type = 'line'
        self.path.airspeed = waypoints.airspeed.item(self.ptr_current)
        self.path.line_origin = dubins_path.r1
        self.path.line_direction = dubins_path.n1
        self.halfspace_n = dubins_path.n1
        self.halfspace_r = dubins_path.r2
        self.path.plot_updated = False

    def construct_dubins_circle_end(self, waypoints, dubins_path):
        self.path.type = 'orbit'
        self.path.airspeed = waypoints.airspeed.item(self.ptr_current)
        self.path.orbit_radius = dubins_path.radius
        self.path.orbit_center = dubins_path.center_e
        if dubins_path.dir_e == 1:
            self.path.orbit_direction = 'CW'
        else:
            self.path.orbit_direction = 'CCW'
        self.halfspace_n = dubins_path.n3
        self.halfspace_r = dubins_path.r3
        self.path.plot_updated = False

    def inHalfSpace(self, pos):
        if (pos-self.halfspace_r).T @ self.halfspace_n >= 0:
            return True
        else:
            return False

