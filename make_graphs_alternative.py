#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 12:49:09 2020

@author: alex
"""

import os
#from amuse.lab import * 
#from amuse.units import units
from matplotlib import pyplot 
from matplotlib import style
import numpy as np
#import csv
from matplotlib import rc
rc('text', usetex=True)

#from properties import plot_track
def diagrama_hr(luminosity,temperature):
    style.use('dark_background')
    figure = pyplot.figure(figsize=(6, 6))
    plot = figure.add_subplot(1, 1, 1)
    plot.set_title('Hertzsprung-Russell diagram', fontsize=14)
    
    plot.loglog(
        temperature,
        luminosity,
        marker="*",color='yellow')  # first few points show transient

    plot.set_xlabel('Effective temperature [$K$]')
    plot.set_ylabel('Luminosity [$L_\odot$]')
    plot.set_xlim(10**4.6, 10**3.5)
    plot.set_ylim(1.0e-2, 1.e5)
    pyplot.savefig('hr-diagram.png')
    print 'Saved hr-diagram.png'
    pyplot.close()

def structure(radius,mass,density,luminosity,temperature,age_label,filename):
    style.use('ggplot')
    masa_cumulativa=np.cumsum(mass)
    figure=pyplot.figure(figsize=(11, 11))
    pyplot.suptitle('Structure at present with alternative composition')
    fig_masa=figure.add_subplot(221)
#    fig_masa.title.set_text('Distribucion cumulativa de masa')
    fig_masa.set_xlabel('r ($R_\odot$)')
    fig_masa.set_ylabel('MASS ($M_\odot$)')
    fig_masa.plot(radius,masa_cumulativa,'r-',label=age_label)
    
    fig_dens=figure.add_subplot(222)
#    fig_dens.title.set_text('Densidad')
    fig_dens.set_xlabel('r ($R_\odot$)')
    fig_dens.set_ylabel('DENSITY ($g/cm^3$)')
    fig_dens.plot(radius,density,'r-',label=age_label)

    
    fig_lumin=figure.add_subplot(223)
#    fig_lumin.set_title('Luminosidad')
    fig_lumin.set_xlabel('r ($R_\odot$)')
    fig_lumin.set_ylabel('LUMINOSITY ($L_\odot$)')
    fig_lumin.plot(radius,luminosity,'r-',label=age_label)

    
    fig_temp=figure.add_subplot(224)
#    fig_temp.set_title('Temperatura')
    fig_temp.set_xlabel('r ($R_\odot$)')
    fig_temp.set_ylabel('TEMPERATURE ($K$)')
    fig_temp.plot(radius,temperature,'r-',label=age_label)

    pyplot.savefig(filename,facecolor='white')
    print 'Saved ',filename
    pyplot.close()

#    
def plot_evolution(radius_history,core_radius_history,pressure_history,temperature_history,
                   luminosity_history,central_density_history,central_temperature_history,time):
    style.use('ggplot')
    figure=pyplot.figure()
    figure.patch.set_facecolor('white')
    plot = figure.add_subplot(1, 1, 1)
    pyplot.title('Evolution with alternative composition using MESA ',color='black',fontsize=12)

    pyplot.axvline(x=10.58,linestyle='dashed',color='gray')
    pyplot.axvline(x=11.444,linestyle='dashed',color='gray')
    pressure_history=pressure_history/pressure_history[0]
    radius_history=radius_history/radius_history[0]
    core_radius_history=np.array(core_radius_history)
    core_radius_history=core_radius_history/core_radius_history[0]
    temperature_history=temperature_history/temperature_history[0]
    central_temperature_history=central_temperature_history/central_temperature_history[0]
    luminosity_history=luminosity_history/luminosity_history[0]
    central_density_history=central_density_history/central_density_history[0]
    
 
    plot.semilogy(time,temperature_history,'y',label='Effective temperature')
    plot.semilogy(time,luminosity_history,'b',label='Luminosity')
    plot.semilogy(time, central_temperature_history,'r',label='Central temperature')
    plot.semilogy(time, central_density_history,color='purple',label='Central density')
    plot.semilogy(time,radius_history,color='cyan',label='Radius')
    plot.semilogy(time,core_radius_history,color='pink',label='Core radius')
    
    pyplot.xlim(9.5,11.6)
    pyplot.xlabel('$t (Gyr)$')
    pyplot.ylabel('$F(t)/F(0)$')
    l=pyplot.legend()

    for text in l.get_texts():
        text.set_color("black")
    pyplot.savefig('evol-comp-alternative-fgb.png',facecolor='white')
    pyplot.show()
    print 'Saved ','evol-comp-alternative-fgb.png'
    
def simple_composition(file_mesa1,file_mesa2,file_mesa3,file_mesa4,j,composition,xlim=1):
    figure = pyplot.figure(figsize=(10, 8))
        
    with open(file_mesa1,'r') as f:
        for i in range(4):
            line=f.readline()
        age_mesa=float(f.readline())
        line=f.readline()

    fig1=figure.add_subplot(221)
#    pyplot.subplot(1.5, 1, 1)
    fig1.set_xlim(left=0,right=xlim)
    fig1.set_ylim(0.0, 1.0)
    
 
    
    radius=np.loadtxt(file_mesa1,skiprows=7,unpack=True, usecols=1,delimiter=',')
    H,He,C,N,O,Ne,Mg=np.loadtxt(file_mesa1,skiprows=7,unpack=True, usecols=(range(6,13)),delimiter=',')
    metals= C+ N+ O+ Ne+ Mg
    fig1.plot(radius,H,'r-',label='H')
    fig1.plot(radius,He,'g-',label='He')
    fig1.plot(radius,metals,'b-',label='Metals')
    

    fig1.set_xlabel('r ($R_\odot$), for age  $t_{MESA}=$'+str(round(age_mesa,3))+str(' $Gyr$'))
    fig1.set_ylabel('Mass fraction')

        
    with open(file_mesa2,'r') as f:
        for i in range(4):
            line=f.readline()
        age_mesa=float(f.readline())
        line=f.readline()

    fig2=figure.add_subplot(222)
#    pyplot.subplot(1.5, 1, 1)
    fig2.set_xlim(left=0,right=xlim)
    fig2.set_ylim(0.0, 1.0)    
    
    radius=np.loadtxt(file_mesa2,skiprows=7,unpack=True, usecols=1,delimiter=',')
    H,He,C,N,O,Ne,Mg=np.loadtxt(file_mesa2,skiprows=7,unpack=True, usecols=(range(6,13)),delimiter=',')
    metals= C+ N+ O+ Ne+ Mg
    fig2.plot(radius,H,'r-',label='H')
    fig2.plot(radius,He,'g-',label='He')
    fig2.plot(radius,metals,'b-',label='Metals')
    

    fig2.set_xlabel('r ($R_\odot$), for age $t_{MESA}=$'+str(round(age_mesa,3))+str(' $Gyr$'))
    fig2.set_ylabel('Mass fraction')
    
    with open(file_mesa3,'r') as f:
        for i in range(4):
            line=f.readline()
        age_mesa=float(f.readline())
        line=f.readline()

    fig3=figure.add_subplot(223)
#    pyplot.subplot(1.5, 1, 1)
    fig3.set_xlim(left=0,right=xlim)
    fig3.set_ylim(0.0, 1.0)    
    
    radius=np.loadtxt(file_mesa3,skiprows=7,unpack=True, usecols=1,delimiter=',')
    H,He,C,N,O,Ne,Mg=np.loadtxt(file_mesa3,skiprows=7,unpack=True, usecols=(range(6,13)),delimiter=',')
    metals= C+ N+ O+ Ne+ Mg
    fig3.plot(radius,H,'r-',label='H')
    fig3.plot(radius,He,'g-',label='He')
    fig3.plot(radius,metals,'b-',label='Metals')
    

    fig3.set_xlabel('r ($R_odot$), for age $t_{MESA}=$'+str(round(age_mesa,3))+str(' $Gyr$'))
    fig3.set_ylabel('Mass fraction')
    
    with open(file_mesa4,'r') as f:
        for i in range(4):
            line=f.readline()
        age_mesa=float(f.readline())
        line=f.readline()

    fig4=figure.add_subplot(224)
#    pyplot.subplot(1.5, 1, 1)
    fig4.set_xlim(left=0,right=xlim)
    fig4.set_ylim(0.0, 1.0)    
    
    radius=np.loadtxt(file_mesa4,skiprows=7,unpack=True, usecols=1,delimiter=',')
    H,He,C,N,O,Ne,Mg=np.loadtxt(file_mesa4,skiprows=7,unpack=True, usecols=(range(6,13)),delimiter=',')
    metals= C+ N+ O+ Ne+ Mg
    fig4.plot(radius,H,'r-',label='H')
    fig4.plot(radius,He,'g-',label='He')
    fig4.plot(radius,metals,'b-',label='Metals')
    

    fig4.set_xlabel('r ($R_odot$), for age $t_{MESA}=$'+str(round(age_mesa,3))+str(' $Gyr$'))
    fig4.set_ylabel('Mass fraction')
    l=pyplot.legend()    
    for text in l.get_texts():
        text.set_color("black")
    pyplot.suptitle('Abundances for MESA and EVtwin from {0} composition'.format(composition), color='black',size=16)
    pyplot.savefig('simplecomp_{0}_{1}.png'.format(composition,str(j)),facecolor='white')
    print 'Saved ','simplecomp_{0}_{1}.png'.format(composition,str(j))
    pyplot.close(figure)
#%%
#FILE STRUCTURE: 'dm (MSun), R (RSun), T (K), rho (g/cm^3), P(Pa), L (LSun), h1, he4, he3, c12, n14, o16, ne20, mg24
#modelo inicial    
#mass,radius,temperature,density,pressure,luminosity,h1, he4, he3, c12, n14, o16, ne20, mg24=np.loadtxt('mesa_results/mesa_results_1.csv',skiprows=7,unpack=True, usecols=range(14),delimiter=',')
#modelo_inicial(radius,mass,density,luminosity,temperature)
#mass,radius,temperature,density,pressure,luminosity,h1, he4, he3, c12, n14, o16, ne20, mg24=np.loadtxt('mesa_results/mesa_results_0.csv',skiprows=7,unpack=True, usecols=range(14),delimiter=',')
#composicion_inicial(radius,h1, he4, he3, c12, n14, o16, ne20, mg24)
#%% PRESENTE

#Leemos los datos de presente
mass,radius,temperature,density,pressure,luminosity=np.loadtxt('mesa_alternative_results/mesa_alternative_results_100.csv',skiprows=7,unpack=True, usecols=range(6),delimiter=',')
#Leemos los datos de referencia

structure(radius,mass,density,luminosity,temperature,'Presente','present_alternative.png')
composition='alternative'
file_mesa1='mesa_{0}_results/mesa_{0}_results_300.csv'.format(composition)
file_mesa2='mesa_{0}_results/mesa_{0}_results_360.csv'.format(composition)
file_mesa3='mesa_{0}_results/mesa_{0}_results_424.csv'.format(composition)
file_mesa4='mesa_{0}_results/mesa_{0}_results_448.csv'.format(composition)
simple_composition(file_mesa1,file_mesa2,file_mesa3,file_mesa4,2,'alternative',xlim=1)

#%%
radius_history=[]
core_radius_history=[]
pressure_history=[]
temperature_history=[]
luminosity_history=[]
central_density_history=[]
central_temperature_history=[]
time=[]

if not os.path.exists('estructura'):
    os.makedirs('estructura')
j=1
nfiles=input('How many files do you have in mesa_alternative_results?')
for i in range(1,nfiles+1):
    filename='mesa_alternative_results/mesa_alternative_results_'+str(i)+'.csv'
    with open(filename,'r') as f:
        for i in range(4):
            line=f.readline()
        age=float(f.readline())
        stellar_type=f.readline()
        core_radius=float(f.readline())
    mass,radius,temperature,density,pressure,luminosity,h1, he4, he3, c12, n14, o16, ne20, mg24=np.loadtxt(filename,skiprows=7,unpack=True, usecols=range(14),delimiter=',')
    temperature_history.append(temperature[-1])
    luminosity_history.append(luminosity[-1])
    central_density_history.append(density[0])
    central_temperature_history.append(temperature[0])
    radius_history.append(radius[-1])
    pressure_history.append(pressure[0])
    time.append(age)
    core_radius_history.append(core_radius)

for i in range(len(time)):
    if central_temperature_history[i]>10**8:
        index=i
        break
print 'Tc>1e8 at time, index ',time[i],index
    
#    plot_temperature(radius,temperature,age,i)
#    plot_composition(h1, he4, he3, c12, n14, o16, ne20, mg24,age,i)
#    estructura(radius_ini/radius_ini[-1],mass_ini,density_ini,luminosity_ini,temperature_ini,radius/radius[-1],mass,density,luminosity,temperature,mass_bench,radius_bench/radius_bench[-1],temp_bench,dens_bench,lumin_bench,str(age)+' Gyr','estructura/estructura_'+str(j)+'.png')

#%%
plot_evolution(radius_history,core_radius_history,pressure_history,temperature_history,luminosity_history,central_density_history,central_temperature_history,time)
diagrama_hr(luminosity_history,temperature_history)


