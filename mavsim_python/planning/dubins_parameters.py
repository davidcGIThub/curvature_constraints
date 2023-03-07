# dubins_parameters
#   - Dubins parameters that define path between two configurations
#
# mavsim_matlab 
#     - Beard & McLain, PUP, 2012
#     - Update history:  
#         3/26/2019 - RWB
#         4/2/2020 - RWB
#         3/30/2022 - RWB

import numpy as np
import sys
sys.path.append('..')


class DubinsParameters:

    def update(self, ps, chis, pe, chie, R):
         self.p_s = ps
         self.chi_s = chis
         self.p_e = pe
         self.chi_e = chie
         self.radius = R
         self.compute_parameters()

    def compute_parameters(self):
        ps = self.p_s
        pe = self.p_e
        chis = self.chi_s
        chie = self.chi_e
        R = self.radius
        ell = np.linalg.norm(ps[0:2] - pe[0:2])
        if ell < 2 * R:
            print('Error in Dubins Parameters: The distance between nodes must be larger than 2R.')
        else:
            # compute start and end circles
            crs = ps + R * rotz(np.pi / 2) @ np.array([[np.cos(chis)], [np.sin(chis)], [0]])
            cls = ps + R * rotz(-np.pi / 2) @ np.array([[np.cos(chis)], [np.sin(chis)], [0]])
            cre = pe + R * rotz(np.pi / 2) @ np.array([[np.cos(chie)], [np.sin(chie)], [0]])
            cle = pe + R * rotz(-np.pi / 2) @ np.array([[np.cos(chie)], [np.sin(chie)], [0]])

            # compute L1
            theta = np.arctan2(cre.item(1) - crs.item(1), cre.item(0) - crs.item(0))
            L1 = np.linalg.norm(crs - cre) \
                 + R * mod(2 * np.pi + mod(theta - np.pi / 2) - mod(chis - np.pi / 2)) \
                 + R * mod(2 * np.pi + mod(chie - np.pi / 2) - mod(theta - np.pi / 2))
            # compute L2
            ell = np.linalg.norm(cle - crs)
            theta = np.arctan2(cle.item(1) - crs.item(1), cle.item(0) - crs.item(0))
            if 2*R > ell:
                L2 = np.inf
            else:
                theta2 = theta - np.pi / 2 + np.arcsin(2 * R / ell)
                L2 = np.sqrt(ell ** 2 - 4 * R ** 2) \
                     + R * mod(2 * np.pi + mod(theta2) - mod(chis - np.pi / 2)) \
                     + R * mod(2 * np.pi + mod(theta2 + np.pi) - mod(chie + np.pi / 2))

            # compute L3
            ell = np.linalg.norm(cre - cls)
            theta = np.arctan2(cre.item(1) - cls.item(1), cre.item(0) - cls.item(0))
            if (2*R > ell):
                L3 = np.inf
            else:
                theta2 = np.arccos(2 * R / ell)
                L3 = np.sqrt(ell ** 2 - 4 * R ** 2) \
                     + R * mod(2 * np.pi - mod(theta + theta2) + mod(chis + np.pi / 2)) \
                     + R * mod(2 * np.pi - mod(theta + theta2 - np.pi) + mod(chie - np.pi / 2))
            # compute L4
            theta = np.arctan2(cle.item(1) - cls.item(1), cle.item(0) - cls.item(0))
            L4 = np.linalg.norm(cls - cle) \
                 + R * mod(2 * np.pi - mod(theta + np.pi / 2) + mod(chis + np.pi / 2)) \
                 + R * mod(2 * np.pi - mod(chie + np.pi / 2) + mod(theta + np.pi / 2))
            # L is the minimum distance
            L = np.min([L1, L2, L3, L4])
            idx = np.argmin([L1, L2, L3, L4])
            e1 = np.array([[1, 0, 0]]).T
            if idx == 0:
                cs = crs
                lams = 1
                ce = cre
                lame = 1
                q1 = (ce - cs) / np.linalg.norm(ce - cs)
                w1 = cs + R * rotz(-np.pi / 2) @ q1
                w2 = ce + R * rotz(-np.pi / 2) @ q1
            elif idx == 1:
                cs = crs
                lams = 1
                ce = cle
                lame = -1
                ell = np.linalg.norm(ce - cs)
                theta = np.arctan2(ce.item(1) - cs.item(1), ce.item(0) - cs.item(0))
                theta2 = theta - np.pi / 2 + np.arcsin(2 * R / ell)
                q1 = rotz(theta2 + np.pi / 2) @ e1
                w1 = cs + R * rotz(theta2) @ e1
                w2 = ce + R * rotz(theta2 + np.pi) @ e1
            elif idx == 2:
                cs = cls
                lams = -1
                ce = cre
                lame = 1
                ell = np.linalg.norm(ce - cs)
                theta = np.arctan2(ce.item(1) - cs.item(1), ce.item(0) - cs.item(0))
                theta2 = np.arccos(2 * R / ell)
                q1 = rotz(theta + theta2 - np.pi / 2) @ e1
                w1 = cs + R * rotz(theta + theta2) @ e1
                w2 = ce + R * rotz(-np.pi + theta + theta2) @ e1
            elif idx == 3:
                cs = cls
                lams = -1
                ce = cle
                lame = -1
                q1 = (ce - cs) / np.linalg.norm(ce - cs)
                w1 = cs + R * rotz(np.pi / 2) @ q1
                w2 = ce + R * rotz(np.pi / 2) @ q1
            w3 = pe
            q3 = rotz(chie) @ e1
            self.length = L
            self.center_s = cs
            self.dir_s = lams
            self.center_e = ce
            self.dir_e = lame
            self.r1 = w1
            self.n1 = q1
            self.r2 = w2
            self.r3 = w3
            self.n3 = q3

    def compute_points(self):
        Del = 0.1  # distance between point

        # points along start circle
        th1 = np.arctan2(self.p_s.item(1) - self.center_s.item(1),
                         self.p_s.item(0) - self.center_s.item(0))
        th1 = mod(th1)
        th2 = np.arctan2(self.r1.item(1) - self.center_s.item(1),
                         self.r1.item(0) - self.center_s.item(0))
        th2 = mod(th2)
        th = th1
        theta_list = [th]
        if self.dir_s > 0:
            if th1 >= th2:
                while th < th2 + 2*np.pi - Del:
                    th += Del
                    theta_list.append(th)
            else:
                while th < th2 - Del:
                    th += Del
                    theta_list.append(th)
        else:
            if th1 <= th2:
                while th > th2 - 2*np.pi + Del:
                    th -= Del
                    theta_list.append(th)
            else:
                while th > th2 + Del:
                    th -= Del
                    theta_list.append(th)

        points = np.array([[self.center_s.item(0) + self.radius * np.cos(theta_list[0]),
                            self.center_s.item(1) + self.radius * np.sin(theta_list[0]),
                            self.center_s.item(2)]])
        for angle in theta_list:
            new_point = np.array([[self.center_s.item(0) + self.radius * np.cos(angle),
                                   self.center_s.item(1) + self.radius * np.sin(angle),
                                   self.center_s.item(2)]])
            points = np.concatenate((points, new_point), axis=0)

        # points along straight line
        sig = 0
        while sig <= 1:
            new_point = np.array([[(1 - sig) * self.r1.item(0) + sig * self.r2.item(0),
                                   (1 - sig) * self.r1.item(1) + sig * self.r2.item(1),
                                   (1 - sig) * self.r1.item(2) + sig * self.r2.item(2)]])
            points = np.concatenate((points, new_point), axis=0)
            sig += Del

        # points along end circle
        th2 = np.arctan2(self.p_e.item(1) - self.center_e.item(1),
                         self.p_e.item(0) - self.center_e.item(0))
        th2 = mod(th2)
        th1 = np.arctan2(self.r2.item(1) - self.center_e.item(1),
                         self.r2.item(0) - self.center_e.item(0))
        th1 = mod(th1)
        th = th1
        theta_list = [th]
        if self.dir_e > 0:
            if th1 >= th2:
                while th < th2 + 2 * np.pi - Del:
                    th += Del
                    theta_list.append(th)
            else:
                while th < th2 - Del:
                    th += Del
                    theta_list.append(th)
        else:
            if th1 <= th2:
                while th > th2 - 2 * np.pi + Del:
                    th -= Del
                    theta_list.append(th)
            else:
                while th > th2 + Del:
                    th -= Del
                    theta_list.append(th)
        for angle in theta_list:
            new_point = np.array([[self.center_e.item(0) + self.radius * np.cos(angle),
                                   self.center_e.item(1) + self.radius * np.sin(angle),
                                   self.center_e.item(2)]])
            points = np.concatenate((points, new_point), axis=0)
        return points


def rotz(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]])


def mod(x):
    while x < 0:
        x += 2*np.pi
    while x > 2*np.pi:
        x -= 2*np.pi
    return x


