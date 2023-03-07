"""
mavsim_python
    - Chapter 6 assignment for Beard & McLain, PUP, 2012
    - Last Update:
        2/5/2019 - RWB
        2/24/2020 - RWB
"""
import sys
sys.path.append('..')
import numpy as np
import parameters.simulation_parameters as SIM

from tools.signals import Signals
from models.mav_dynamics_control import MavDynamics
from models.wind_simulation import WindSimulation
from control.autopilot import Autopilot
# from control.autopilot_lqr import Autopilot
# from control.autopilot_tecs import Autopilot
from viewers.mav_viewer import MavViewer
from viewers.path_data_viewer import PathDataViewer
import pyqtgraph as pg
from tools.quit_listener import QuitListener
from tools.rotations import Euler2Rotation

quitter = QuitListener()

VIDEO = False
PATH_PLOTS = True
ANIMATION = True
SAVE_PLOT_IMAGE = False

# video initialization
if VIDEO is True:
    from viewers.video_writer import VideoWriter
    video = VideoWriter(video_name="chap6_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

#initialize the visualization
if ANIMATION or PATH_PLOTS:
    app = pg.QtWidgets.QApplication([]) # use the same main process for Qt applications
if ANIMATION:
    mav_view = MavViewer(app=app)  # initialize the mav viewer
if PATH_PLOTS:
    # initialize view of data plots
    path_data_view = PathDataViewer(app=app,dt=SIM.ts_simulation, plot_period=SIM.ts_plot_refresh, 
                           data_recording_period=SIM.ts_plot_record_data, time_window_length=30)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
autopilot = Autopilot(SIM.ts_simulation)

# autopilot commands
from message_types.msg_autopilot import MsgAutopilot
commands = MsgAutopilot()
Va_command = Signals(dc_offset=25.0,
                     amplitude=3.0,
                     start_time=2.0,
                     frequency=0.01)
altitude_rate_command = Signals(dc_offset=0.0,
                           amplitude=5,
                           start_time=0.0,
                           frequency=0.02)
course_command = Signals(dc_offset=np.radians(180),
                         amplitude=np.radians(45),
                         start_time=5.0,
                         frequency=0.015)

# initialize the simulation time
sim_time = SIM.start_time
end_time = 200

# main simulation loop
print("Press 'Esc' to exit...")
while sim_time < end_time:

    # -------autopilot commands-------------
    commands.airspeed_command = Va_command.square(sim_time)
    commands.course_command = course_command.square(sim_time)
    commands.altitude_rate_command = altitude_rate_command.square(sim_time)

    # -------autopilot-------------
    delta, commanded_state = autopilot.update(commands, mav.true_state)

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    # ------- animation -------
    if ANIMATION:
        mav_view.update(mav.true_state)  # plot body of MAV
    if PATH_PLOTS:
        plot_time = sim_time
        # hdot = mav.true_state.Va*np.sin(mav.true_state.theta)
        # mav.true_state.altitude =  hdot
        path_data_view.update(mav.true_state,  # true states
                            commanded_state)  # commanded states
    if ANIMATION or PATH_PLOTS:
        app.processEvents()
    if VIDEO is True:
        video.update(sim_time)
        
    # -------Check to Quit the Loop-------
    if quitter.check_quit():
        break

    # -------increment time-------------
    sim_time += SIM.ts_simulation

if SAVE_PLOT_IMAGE:
    path_data_view.save_plot_image("ch6_plot")

if VIDEO is True:
    video.close()




