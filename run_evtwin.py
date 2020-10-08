#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 19:06:09 2020

@author: alex
"""

from amuse.lab import EVtwin 
from amuse.units import units
from functions import uniform_star,cycle,save_results
import os

def run_simulation(Z,H,composition):
    if not os.path.exists('evtwin_{0}_results'.format(composition)):
        os.makedirs('evtwin_{0}_results'.format(composition))
    
    mass = 1.0 | units.MSun

    evtwin=EVtwin() 
    evtwin.set_RGB_wind_setting(0.0)
    evtwin.set_Ostar_wind_setting(0.0)
    #IMPORTANT: Currently AGB wind cannot be switched off with AMUSE/EVtwin. It is set to the default value 1, corresponding to Watcher's fitting formula. To see other parameters, type 'evtwin=EVtwin()' and 'print evtwin.parameters' in your python console.
    #%%
    """##################################################################
                    CREATING A STAR WITH UNIFORM COMPOSITION
        #################################################################"""
    composition_file='{0}.csv'.format(composition)
    sun=uniform_star(evtwin,'evtwin',mass,Z,H,composition_file)
    #%%
    """##################################################################
                EVOLVING THE STAR THROUGH DIFFERENT STELLAR TYPES
        #################################################################"""
    i=0 #un contador de los pasos de avance de t
    t=0 |units.Gyr
    path='evtwin_{0}_results/evtwin_{0}_results_'.format(composition)
    filename=path+str(i)+'.csv'
    
    print 'Stellar type is {0} and age is {1} '.format(sun.stellar_type,sun.age)
    save_results(sun, filename,composition,'evtwin')
    
    print 'Evolving'
    i,t=cycle(i,100,t,sun,path,composition,evtwin,'evtwin') #MAIN SEQUENCE
    i,t=cycle(i,500,t,sun,path,composition,evtwin,'evtwin') #HERTZSPRUNG GAP (Beginning of the FGB)

    
    evtwin.stop()
    print 'Your results are in evtwin_{0}_results'.format(composition)
    print 'end'
    
"""##################################################################
            MAIN: RUN SIMULATION WITH DESIRED COMPOSITION
    #################################################################"""
    
which_simulation='Invalid answer'
while which_simulation=='Invalid answer':
    which_simulation=raw_input('Choose a composition. Type AGS or GS98\n')
    if which_simulation=='AGS':
        Z=1.220445E-02
        H=0.75830
        run_simulation(Z,H,'AGS')
    elif which_simulation=='GS98':
        Z=1.695595E-02
        H=0.74045
        run_simulation(Z,H,'GS98')
    else:
        print 'Invalid answer'
        which_simulation='Invalid answer'
 
