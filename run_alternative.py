#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 11:36:40 2020

@author: alex

"""

from amuse.lab import MESA,Particle
from amuse.units import units
from functions import cycle,save_results
import os
import numpy as np
#%%
  
if not os.path.exists('mesa_alternative_results'):
    os.makedirs('mesa_alternative_results')
  
mass = 1.0 | units.MSun

mesa=MESA() 
# We switch off the wind.
mesa.parameters.AGB_wind_scheme=0
mesa.parameters.RGB_wind_scheme=0
mesa.parameters.blocker_wind_efficiency=0
mesa.parameters.de_jager_wind_efficiency=0
mesa.parameters.dutch_wind_efficiency=0

"""##################################################################
            CREATING A STAR WITH UNIFORM COMPOSITION
#################################################################"""
estrella_default=mesa.particles.add_particle(Particle(mass=mass)) 

zones=estrella_default.get_number_of_zones()

composition=np.ones(zones)
sun=mesa.new_particle_from_model(dict(
    mass=(estrella_default.get_cumulative_mass_profile(number_of_zones=zones)*estrella_default.mass),
    radius=estrella_default.get_radius_profile(number_of_zones=zones),
    rho=estrella_default.get_density_profile(number_of_zones=zones),
    temperature=estrella_default.get_temperature_profile(number_of_zones=zones),
    luminosity=estrella_default.get_luminosity_profile(number_of_zones=zones),
    X_H=composition*0.706,
    X_He=composition*0.275,
    X_C=composition*3.03e-3,
    X_N=composition*1.10e-3,
    X_O=composition*9.59e-3,
    X_Ne=composition*1.62e-3,
    X_Mg=composition*5.15e-4,
    X_He3=composition*0.0,
    X_Si=composition*6.53e-4,
    X_Fe=composition*1.17e-3), 0.0 | units . Myr )

"""##################################################################
        EVOLVING THE STAR THROUGH DIFFERENT STELLAR TYPES
#################################################################"""
i=0 #un contador de los pasos de avance de t
t=0 |units.Gyr
composition='alternative'
path='mesa_{0}_results/mesa_{0}_results_'.format(composition)
filename=path+str(i)+'.csv'


print 'Stellar type is {0} and age is {1} '.format(sun.stellar_type,sun.age)
save_results(sun, filename,composition,'mesa')

print 'Evolving'
i,t=cycle(i,100,t,sun,path,composition,mesa,'mesa') #MAIN SEQUENCE
i,t=cycle(i,1000,t,sun,path,composition,mesa,'mesa') #FIRST GIANT BRANCH
mesa.stop()
print 'Your results are in mesa_{0}_results'.format(composition)
print 'end'


    

        

