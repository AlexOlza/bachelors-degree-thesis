#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script is part of the end-of-degree project by Alexander Olza Rodriguez (UPV/EHU).

It goes with the file amuse_tutorial.py. This file plots the data obtained in amuse-tutorial

"""
def plot_error(dE):
    import matplotlib.pyplot as pyplot
    import numpy
    fig=pyplot.figure()
    ax=fig.add_subplot(111)
    pyplot.title('Energy error integrating with ph4')
    pyplot.xlabel('Integration steps')
    pyplot.ylabel('Relative error')
    ax.plot(numpy.absolute(dE))
    
def plot_trajectories(sun_position,mercury_position,venus_position,earth_position):
    import matplotlib.pyplot as pyplot
    from mpl_toolkits.mplot3d import Axes3D
    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')
    #we need to reorganize the data in separate arrays to make a scatter plot
    sun_x,sun_y,sun_z=[],[],[]
    for i in range(len(sun_position)):
        sun_x.append(sun_position[i][0])
        sun_y.append(sun_position[i][1])
        sun_z.append(sun_position[i][2])
        
    merc_x,merc_y,merc_z=[],[],[]
    for i in range(len(mercury_position)):
        merc_x.append(mercury_position[i][0])
        merc_y.append(mercury_position[i][1])
        merc_z.append(mercury_position[i][2])
        
    ven_x,ven_y,ven_z=[],[],[]
    for i in range(len(venus_position)):
        ven_x.append(venus_position[i][0])
        ven_y.append(venus_position[i][1])
        ven_z.append(venus_position[i][2])
        
    ear_x,ear_y,ear_z=[],[],[]
    for i in range(len(earth_position)):
        ear_x.append(earth_position[i][0])
        ear_y.append(earth_position[i][1])
        ear_z.append(earth_position[i][2])
    
    ax.scatter(sun_x,sun_y,sun_z,s=2)
    ax.scatter(merc_x,merc_y,merc_z,s=2)
    ax.scatter(ven_x,ven_y,ven_z,s=2)
    ax.scatter(ear_x,ear_y,ear_z,s=2)
    
    pyplot.locator_params(axis='y', nbins=6)
    pyplot.locator_params(axis='x', nbins=6)
    pyplot.title('Trajectories of some planets')
    pyplot.xlabel(' x (AU)')
    pyplot.ylabel(' y (AU)')
    pyplot.legend(('Sun','Mercury','Venus','Earth'))
    pyplot.show()
    
