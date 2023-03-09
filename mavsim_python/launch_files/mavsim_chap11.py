"""
mavsim_python
    - Chapter 11 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        3/26/2019 - RWB
        2/27/2020 - RWB
"""
import sys
sys.path.append('..')
import numpy as np
import pyqtgraph as pg
import parameters.simulation_parameters as SIM
import parameters.planner_parameters as PLAN
from models.mav_dynamics_control import MavDynamics
from models.wind_simulation import WindSimulation
from control.autopilot import Autopilot
from planning.spline_path_follower import BsplinePathFollower
from viewers.path_data_viewer import PathDataViewer
from viewers.mav_spline_viewer import MAVAndSplineViewer
from tools.quit_listener import QuitListener
from message_types.msg_autopilot import MsgAutopilot


quitter = QuitListener()

VIDEO = False
DATA_PLOTS = False
ANIMATION = True
SAVE_PLOT_IMAGE = False

# video initialization
if VIDEO is True:
    from viewers.video_writer import VideoWriter
    video = VideoWriter(video_name="chap11_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

#initialize the visualization
if ANIMATION or DATA_PLOTS:
    app = pg.QtWidgets.QApplication([]) # use the same main process for Qt applications
if ANIMATION:
    mav_path_view = MAVAndSplineViewer(app=app)  # initialize the mav waypoint viewer
if DATA_PLOTS:
    data_view = PathDataViewer(app=app,dt=SIM.ts_simulation, plot_period=SIM.ts_plot_refresh, 
                           data_recording_period=SIM.ts_plot_record_data, time_window_length=30)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
autopilot = Autopilot(SIM.ts_simulation)
order = 3
spline_path_follower = BsplinePathFollower(order)

# initialize the simulation time
sim_time = SIM.start_time
end_time = 200
path_update_flag = True
control_points = np.array([[-295.79238267,   17.39955685,  226.19415528,  330.59143755,  330.59143456,
   226.19415857,   17.3995552,  -295.79237938],
   [  27.55578793,  -13.77789397,   27.55578793,  134.36206838,  265.63791371,
   372.44420061,  413.7778997,   372.44420061],
 [ 103.44446534,   98.27776733,  103.44446534,  116.79526483,  133.20474561,
   146.5555204,   151.7222398,   146.5555204 ]])

commands = MsgAutopilot()

# Va_command = Signals(dc_offset=25.0,
#                      amplitude=3.0,
#                      start_time=2.0,
#                      frequency=0.01)
# altitude_rate_command = Signals(dc_offset=0.0,
#                            amplitude=5,
#                            start_time=0.0,
#                            frequency=0.02)
# course_command = Signals(dc_offset=np.radians(180),
#                          amplitude=np.radians(45),
#                          start_time=5.0,
#                          frequency=0.015)

# main simulation loop
print("Press 'Esc' to exit...")
while sim_time < end_time:
    # -------observer-------------

    # -------path manager-------------
    # path = path_manager.update(waypoints, PLAN.R_min, mav.true_state)

    # -------path follower-------------
    # autopilot_commands = path_follower.update(path, mav.true_state)
    position = np.array([[mav.true_state.north],[mav.true_state.east],[mav.true_state.altitude]])
    desired_airspeed = 30
    [commands.course_command, commands.altitude_rate_command, commands.airspeed_command] = \
        spline_path_follower.get_commands(control_points, position, desired_airspeed)

    # -------autopilot-------------
    delta, commanded_state = autopilot.update(commands, mav.true_state)

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    # -------update viewer-------------
    if ANIMATION:
        mav_path_view.update_mav(mav.true_state)  # plot path and MAV
        if path_update_flag == True:
            mav_path_view.update_path(control_points, order)
            path_update_flag = False
    
    if DATA_PLOTS:
        plot_time = sim_time
        data_view.update(mav.true_state,  # true states
                         commanded_state)  # commanded states
    if ANIMATION or DATA_PLOTS:
        app.processEvents()

    # -------Check to Quit the Loop-------
    if quitter.check_quit():
        break

    # -------increment time-------------
    sim_time += SIM.ts_simulation

# Save an Image of the Plot
if SAVE_PLOT_IMAGE and DATA_PLOTS:
        data_view.save_plot_image("ch11_data_plot")

if VIDEO is True:
    video.close()




