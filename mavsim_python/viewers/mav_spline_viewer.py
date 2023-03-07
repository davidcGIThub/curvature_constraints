"""
mavsim_python: waypoint viewer (for chapter 11)
    - Beard & McLain, PUP, 2012
    - Update history:
        4/15/2019 - BGM
"""
import sys
sys.path.append("..")
import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from viewers.draw_mav import DrawMav
from viewers.draw_spline_path import DrawSplinePath


class MAVAndSplineViewer:
    def __init__(self, app):
        self.scale = 1000
        # initialize Qt gui application and window
        self.app = app  # initialize QT
        self.window = gl.GLViewWidget()  # initialize the view object
        self.window.setWindowTitle('World Viewer')
        grid = gl.GLGridItem() # make a grid to represent the ground
        grid.scale(self.scale/20, self.scale/20, self.scale/20) # set the size of the grid (distance between each line)
        self.window.addItem(grid) # add grid to viewer
        self.window.setCameraPosition(distance=self.scale, elevation=90, azimuth=-90)
        self.window.setBackgroundColor('k')  # set background color to black
        self.window.setGeometry(0, 0, 750, 750)  # args: upper_left_x, upper_right_y, width, height
        center = self.window.cameraPosition()
        center.setX(250)
        center.setY(250)
        center.setZ(0)
        self.window.setCameraPosition(pos=center, distance=self.scale, elevation=50, azimuth=-90)
        self.window.show()  # display configured window
        self.window.raise_()  # bring window to the front
        self.mav_plot_initialized = False  # has the mav been plotted yet?
        self.path_plot_initialized = False  # has the mav been plotted yet?
        self.mav_plot = []
        self.path_plot = []

    def update_mav(self, state):
        # initialize the drawing the first time update() is called
        if not self.mav_plot_initialized:
            self.mav_plot = DrawMav(state, self.window)
            self.mav_plot_initialized = True
        # else update drawing on all other calls to update()
        else:
            self.mav_plot.update(state)


    def update_path(self, control_points, order):
        if not self.path_plot_initialized:
            self.path_plot = DrawSplinePath(control_points, order, self.window)
            self.path_plot_initialized = True
        else:
            self.path_plot.update(control_points)

    def process_app(self):
        self.app.processEvents()

    def clear_viewer(self):
        self.window.clear()
