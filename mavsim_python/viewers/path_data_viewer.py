from plotter.plotter import Plotter
import numpy as np
from tools.wrap import wrap

class PathDataViewer:
    def __init__(self, app,  dt = 0.01,
                 time_window_length = 30, # number of data points plotted at a time
                 plot_period = 0.2, # time interval between a plot update
                 data_recording_period = 0.1): # time interval between recording a data update
        self._dt = dt
        self._data_window_length= time_window_length/data_recording_period
        self._update_counter = 0
        self._plots_per_row = 3
        self._plotter = Plotter(app=app, plots_per_row=self._plots_per_row)  # plot last time_window seconds of data
        self._plot_period = plot_period
        self._data_recording_period = data_recording_period
        self._plot_delay = 0
        self._data_recording_delay = 0
        self._time = 0

        #define colors
        truth_color = (0,255,0)
        control_color = (0,0,255)

        # define first row
        self._plotter.create_plot_widget(plot_id='pn', xlabel='Time (s)', ylabel='pn (m)',
                                        window_length=self._data_window_length)
        self._plotter.create_plot_widget(plot_id='pe', xlabel='Time (s)', ylabel='pe (m)',
                                       window_length=self._data_window_length)
        self._plotter.create_plot_widget(plot_id='h', xlabel='Time (s)', ylabel='h (m)',
                                       window_length=self._data_window_length)
        self._plotter.create_data_set(plot_id="pn", data_label="pn", data_color=truth_color)
        self._plotter.create_data_set(plot_id="pn", data_label="pn_d", data_color=control_color)
        self._plotter.create_data_set(plot_id="pe", data_label="pe", data_color=truth_color)
        self._plotter.create_data_set(plot_id="pe", data_label="pe_d", data_color=control_color)
        self._plotter.create_data_set(plot_id="h", data_label="h", data_color=truth_color)
        self._plotter.create_data_set(plot_id="h", data_label="h_d", data_color=control_color)
        
        # define second row
        self._plotter.create_plot_widget(plot_id='phi', xlabel='Time (s)', ylabel='phi (deg)',
                                       window_length=self._data_window_length)
        self._plotter.create_plot_widget(plot_id='theta', xlabel='Time (s)', ylabel='theta (deg)',
                                       window_length=self._data_window_length)
        self._plotter.create_plot_widget(plot_id='psi', xlabel='Time (s)', ylabel='psi (deg)',
                                       window_length=self._data_window_length)
        self._plotter.create_data_set(plot_id="phi", data_label="phi", data_color=truth_color)
        self._plotter.create_data_set(plot_id="phi", data_label="phi_c", data_color=control_color)
        self._plotter.create_data_set(plot_id="theta", data_label="theta", data_color=truth_color)
        self._plotter.create_data_set(plot_id="theta", data_label="theta_c", data_color=control_color)
        self._plotter.create_data_set(plot_id="psi", data_label="psi", data_color=truth_color)
        self._plotter.create_data_set(plot_id="psi", data_label="psi_c", data_color=control_color)

        # define third row
        self._plotter.create_plot_widget(plot_id='Va', xlabel='Time (s)', ylabel='Va (m/s)',
                                       window_length=self._data_window_length)
        self._plotter.create_plot_widget(plot_id='chi', xlabel='Time (s)', ylabel='chi (deg)',
                                       window_length=self._data_window_length)
        self._plotter.create_plot_widget(plot_id='h_dot', xlabel='Time (s)', ylabel='h_dot (m/s)',
                                       window_length=self._data_window_length)
        self._plotter.create_data_set(plot_id="Va", data_label="Va", data_color=truth_color)
        self._plotter.create_data_set(plot_id="Va", data_label="Va_c", data_color=control_color)
        self._plotter.create_data_set(plot_id="chi", data_label="chi", data_color=truth_color)
        self._plotter.create_data_set(plot_id="chi", data_label="chi_c", data_color=control_color)
        self._plotter.create_data_set(plot_id="h_dot", data_label="h_dot", data_color=truth_color)
        self._plotter.create_data_set(plot_id="h_dot", data_label="h_dot_c", data_color=control_color)
    
        self._plotter.show_window()

    def update(self, true_state, commanded_state):
        if self._data_recording_delay >= self._data_recording_period:
            self.__update_data(true_state, commanded_state, self._time)
            self._data_recording_delay = 0
        if self._plot_delay >= self._plot_period:
            self.__update_plot()
            self._plot_delay = 0
        self._plot_delay += self._dt
        self._data_recording_delay += self._dt
        self._time += self._dt
        
    def __update_data(self, true_state, commanded_state, t):
        #add the commanded state data
        if commanded_state != None:
            self._plotter.add_data_point(plot_id='pn', data_label='pn_d', xvalue=t, yvalue=commanded_state.north)
            self._plotter.add_data_point(plot_id='pe', data_label='pe_d', xvalue=t, yvalue=commanded_state.east)
            self._plotter.add_data_point(plot_id='h', data_label='h_d', xvalue=t, yvalue=commanded_state.altitude)
            self._plotter.add_data_point(plot_id='phi', data_label='phi_c', xvalue=t, yvalue=self.__rad_to_deg(commanded_state.phi))
            self._plotter.add_data_point(plot_id='theta', data_label='theta_c', xvalue=t, yvalue=self.__rad_to_deg(commanded_state.theta))
            self._plotter.add_data_point(plot_id='psi', data_label='psi_c', xvalue=t, yvalue=self.__rad_to_deg(commanded_state.psi))
            self._plotter.add_data_point(plot_id='Va', data_label='Va_c', xvalue=t, yvalue=commanded_state.Va)
            self._plotter.add_data_point(plot_id='chi', data_label='chi_c', xvalue=t, yvalue=self.__rad_to_deg(commanded_state.chi))
            self._plotter.add_data_point(plot_id='h_dot', data_label='h_dot_c', xvalue=t, yvalue=commanded_state.h_dot)

        #add the true state data
        if true_state != None:
            self._plotter.add_data_point(plot_id='pn', data_label='pn', xvalue=t, yvalue=true_state.north)
            self._plotter.add_data_point(plot_id='pe', data_label='pe', xvalue=t, yvalue=true_state.east)
            self._plotter.add_data_point(plot_id='h', data_label='h', xvalue=t, yvalue=true_state.altitude)
            self._plotter.add_data_point(plot_id='phi', data_label='phi', xvalue=t, yvalue=self.__rad_to_deg(true_state.phi))
            self._plotter.add_data_point(plot_id='theta', data_label='theta', xvalue=t, yvalue=self.__rad_to_deg(true_state.theta))
            self._plotter.add_data_point(plot_id='psi', data_label='psi', xvalue=t, yvalue=self.__rad_to_deg(true_state.psi))
            self._plotter.add_data_point(plot_id='Va', data_label='Va', xvalue=t, yvalue=true_state.Va)
            self._plotter.add_data_point(plot_id='chi', data_label='chi', xvalue=t, yvalue=self.__rad_to_deg(true_state.chi))
            self._plotter.add_data_point(plot_id='h_dot', data_label='h_dot', xvalue=t, yvalue=true_state.h_dot)

    def process_app(self):
        self._plotter.process_app(0)

    def __update_plot(self):
        self._plotter.update_plots()

    def close_data_viewer(self):
        self._plotter.close_window()

    def save_plot_image(self, plot_name):
        self._plotter.save_image(plot_name)

    def __rad_to_deg(self, radians):
        rad = wrap(radians,0)
        return rad*180/np.pi


