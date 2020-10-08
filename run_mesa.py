#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script is part of the end-of-degree project by Alexander Olza Rodriguez (UPV/EHU).

Created on Sat Mar 28 11:36:40 2020

@author: alex

"""

from amuse.lab import MESA 
from amuse.units import units
from functions import uniform_star,cycle,save_results
import os

def run_simulation(Z,H,composition):
    if not os.path.exists('mesa_{0}_results'.format(composition)):
        os.makedirs('mesa_{0}_results'.format(composition))
    
    mass = 1.0 | units.MSun

    mesa=MESA() 
    # We switch off the wind
    mesa.parameters.AGB_wind_scheme=0
    mesa.parameters.RGB_wind_scheme=0
    mesa.parameters.blocker_wind_efficiency=0
    mesa.parameters.de_jager_wind_efficiency=0
    mesa.parameters.dutch_wind_efficiency=0
    
    """##################################################################
                    CREATING A STAR WITH UNIFORM COMPOSITION
        #################################################################"""
    composition_file='{0}.csv'.format(composition)
    sun=uniform_star(mesa,'mesa',mass,Z,H,composition_file)
    
    """##################################################################
                EVOLVING THE STAR THROUGH DIFFERENT STELLAR TYPES
        #################################################################"""
    i=0 #a counter for the time intervals
    t=0 |units.Gyr
    path='mesa_{0}_results/mesa_{0}_results_'.format(composition)
    filename=path+str(i)+'.csv'
    composition=composition
    
    print 'Stellar type is {0} and age is {1} '.format(sun.stellar_type,sun.age)
    save_results(sun, filename,composition,'mesa')
    
    print 'Evolving'
    i,t=cycle(i,100,t,sun,path,composition,mesa,'mesa') #MAIN SEQUENCE
    i,t=cycle(i,500,t,sun,path,composition,mesa,'mesa') #FIRST GIANT BRANCH 
    mesa.stop()
    print 'Your results are in mesa_{0}_results'.format(composition)
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
        

