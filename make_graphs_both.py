#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 24 17:08:38 2020

@author: alex
"""

import os
from matplotlib import pyplot 
from matplotlib import style
import numpy as np
from csv import reader
from matplotlib import rc
rc('text', usetex=True)

    
def structure(age,radius_GS,mass_GS,density_GS,luminosity_GS,temperature_GS,
              radius_AGS,mass_AGS,density_AGS,luminosity_AGS,temperature_AGS):
    style.use('ggplot')
    masa_cumulativa_GS=np.cumsum(mass_GS)
    masa_cumulativa_AGS=np.cumsum(mass_AGS)
    figure=pyplot.figure(figsize=(8, 8))
    fig_masa=figure.add_subplot(221)
    fig_masa.title.set_text('Cumulative mass distribution')
    fig_masa.set_xlabel('r ($R_\odot$)')
    fig_masa.set_ylabel('MASS ($M_\odot$)')
    fig_masa.plot(radius_GS,masa_cumulativa_GS,'b--')
    fig_masa.plot(radius_AGS,masa_cumulativa_AGS,'g--')

    fig_dens=figure.add_subplot(222)
    fig_dens.title.set_text('Density')
    fig_dens.set_xlabel('r ($R_\odot$)')
    fig_dens.set_ylabel('DENSITY ($g/cm^3$)')
    fig_dens.plot(radius_GS,density_GS,'b--')
    fig_dens.plot(radius_AGS,density_AGS,'g--')
    
    fig_lumin=figure.add_subplot(223)
    fig_lumin.set_title('Luminosity')
    fig_lumin.set_xlabel('r ($R_\odot$)')
    fig_lumin.set_ylabel('LUMINOSITY ($L_\odot$)')
    fig_lumin.plot(radius_GS,luminosity_GS,'b--')
    fig_lumin.plot(radius_AGS,luminosity_AGS,'g--')
    
    fig_temp=figure.add_subplot(224)
    fig_temp.set_title('Temperature')
    fig_temp.set_xlabel('r ($R_\odot$)')
    fig_temp.set_ylabel('TEMPERATURE ($K$)')
    fig_temp.plot(radius_GS,temperature_GS,'b--',label='GS98')
    fig_temp.plot(radius_AGS,temperature_AGS,'g--',label='AGS')
    
    pyplot.legend()
    pyplot.suptitle('Structure at age {0}'.format(age),size=16)
    pyplot.savefig('figures/{0}'.format(filename))
    pyplot.show()
    print 'Saved ',filename
    pyplot.close()    

def present(mass_mesa,radius_mesa,temperature_mesa,density_mesa,luminosity_mesa,
        mass,radius,temperature,density,luminosity,
        mass_bench,radius_bench,temp_bench,dens_bench,lumin_bench,composition,filename):
    style.use('ggplot')
    masa_cumulativa_mesa=np.cumsum(mass_mesa)
    masa_cumulativa=np.cumsum(mass)
    figure=pyplot.figure(figsize=(10, 10))
    figure.patch.set_facecolor('white')
    fig_masa=figure.add_subplot(221)
#    fig_masa.title.set_text('Distribucion cumulativa de masa')
    fig_masa.set_xlabel('r ($R_\odot$)')
    fig_masa.set_ylabel('MASS ($M_\odot$)')
    fig_masa.plot(radius,masa_cumulativa,'r-',label='EVtwin')
    fig_masa.plot(radius_mesa,masa_cumulativa_mesa,'g-',label='EVtwin')
    fig_masa.plot(radius_bench,mass_bench,'k--',label='EVtwin')

    fig_dens=figure.add_subplot(222)
#    fig_dens.title.set_text('Densidad')
    fig_dens.set_xlabel('r ($R_\odot$)')
    fig_dens.set_ylabel('DENSITY ($g/cm^3$)')
    fig_dens.plot(radius,density,'r-',label='EVtwin')
    fig_dens.plot(radius_mesa,density_mesa,'g-')
    fig_dens.plot(radius_bench,dens_bench,'k--')
    
    fig_lumin=figure.add_subplot(223)
#    fig_lumin.set_title('Luminosidad')
    luminosity=luminosity*3.8418/3.839
    fig_lumin.set_xlabel('r ($R_\odot$)')
    fig_lumin.set_ylabel('LUMINOSITY ($L_\odot$)')
    fig_lumin.plot(radius,luminosity,'r-',label='EVtwin')
    fig_lumin.plot(radius_mesa,luminosity_mesa,'g-')
    fig_lumin.plot(radius_bench,lumin_bench,'k--')
    
    fig_temp=figure.add_subplot(224)
#    fig_temp.set_title('Temperatura')
    fig_temp.set_xlabel('r ($R_\odot$)')
    fig_temp.set_ylabel('TEMPERATURE ($K$)')
    fig_temp.plot(radius,temperature,'r-',label='EVtwin')
    fig_temp.plot(radius_mesa,temperature_mesa,'g-',label='MESA')
    fig_temp.plot(radius_bench,temp_bench,'k--',label='Benchmark')
    pyplot.legend()
    pyplot.suptitle('Structure at present with composition {0}'.format(composition),size=16)
    pyplot.savefig('figures/{0}'.format(filename),facecolor='white')
    pyplot.show()
    
    print 'Saved ',filename
    pyplot.close()    
    
def plot_evolution(radius_history,core_radius_history,pressure_history,temperature_history,
                   luminosity_history,central_density_history,central_temperature_history,time,
                   radius_history_evtwin,core_radius_history_evtwin,pressure_history_evtwin,temperature_history_evtwin,
                   luminosity_history_evtwin,central_density_history_evtwin,central_temperature_history_evtwin,time_evtwin,composition):
    style.use('ggplot')
    figure=pyplot.figure()
    figure.patch.set_facecolor('white')
    plot = figure.add_subplot(1, 1, 1)
    pyplot.title('Evolution with {0} using MESA (solid) and EVtwin (dashed)'.format(composition),color='black',fontsize=12)
    if composition=='AGS':
        pyplot.axvline(x=11.716,linestyle='dashed',color='gray')
    pressure_history=pressure_history/pressure_history[0]
    radius_history=radius_history/radius_history[0]
    core_radius_history=np.array(core_radius_history)
    core_radius_history=core_radius_history/core_radius_history[0]
    temperature_history=temperature_history/temperature_history[0]
    central_temperature_history=central_temperature_history/central_temperature_history[0]
    luminosity_history=luminosity_history/luminosity_history[0]
    central_density_history=central_density_history/central_density_history[0]
    
    pressure_history_evtwin=pressure_history_evtwin/pressure_history_evtwin[0]
    radius_history_evtwin=radius_history_evtwin/radius_history_evtwin[0]
    core_radius_history_evtwin=np.array(core_radius_history_evtwin)
    core_radius_history_evtwin=core_radius_history_evtwin/core_radius_history_evtwin[0]
    temperature_history_evtwin=temperature_history_evtwin/temperature_history_evtwin[0]
    central_temperature_history_evtwin=central_temperature_history_evtwin/central_temperature_history_evtwin[0]
    luminosity_history_evtwin=luminosity_history_evtwin/luminosity_history_evtwin[0]
    central_density_history_evtwin=central_density_history_evtwin/central_density_history_evtwin[0]
    
    plot.semilogy(time,temperature_history,'y',label='Effective temperature')
    plot.semilogy(time,luminosity_history,'b',label='Luminosity')
    plot.semilogy(time, central_temperature_history,'r',label='Central temperature')
    plot.semilogy(time, central_density_history,color='purple',label='Central density')
    plot.semilogy(time,radius_history,color='cyan',label='Radius')
    plot.semilogy(time,core_radius_history,color='pink',label='Core radius')
    
    plot.semilogy(time_evtwin,temperature_history_evtwin,'y--')
    plot.semilogy(time_evtwin,luminosity_history_evtwin,'b--')
    plot.semilogy(time_evtwin, central_temperature_history_evtwin,'r--')
    plot.semilogy(time_evtwin, central_density_history_evtwin,color='purple',linestyle='dashed')
    plot.semilogy(time_evtwin,radius_history_evtwin,color='cyan',linestyle='dashed')
    plot.semilogy(time_evtwin,core_radius_history_evtwin,color='pink',linestyle='dashed')

    pyplot.xlabel('$t (Gyr)$')
    pyplot.ylabel('$F(t)/F(0)$')
    l=pyplot.legend()

    for text in l.get_texts():
        text.set_color("black")
    pyplot.savefig('figures/evol-comp-both-{0}.png'.format(composition),facecolor='white')
    pyplot.show()
    print 'Saved ','figures/evol-comp-both-{0}.png'.format(composition)
    
#def evolution_data(radius_history,core_radius_history,pressure_history,temperature_history,
#                   luminosity_history,central_density_history,central_temperature_history,time,composition):
#    pressure_history=pressure_history/pressure_history[0]
#    print 'max radio', max(radius_history)
#    radius_history=radius_history/radius_history[0]
#    core_radius_history=np.array(core_radius_history)
#    core_radius_history=core_radius_history/core_radius_history[0]
#    print 'max radio relativo', max(radius_history)
#    print 'cambio Teff', temperature_history[0], max(temperature_history)
#    print 'cambio Teff relativo',max(temperature_history)/temperature_history[0]
#    temperature_history=temperature_history/temperature_history[0]
#    print 'max t central', max(central_temperature_history)
#    central_temperature_history=central_temperature_history/central_temperature_history[0]
#    print 'max t central relativa', max(central_temperature_history)   
#    print 'cambio T central relativo', central_temperature_history[-1]
#    print 'max lum', max(luminosity_history)
#    luminosity_history=luminosity_history/luminosity_history[0]
#    print 'max lum relativa', max(luminosity_history), np.argmax(luminosity_history)
#    print 'Lugar llamativo en fgb', time[int(np.argmax(central_density_history))],time[int(np.argmax(luminosity_history))]
#    print 'max densidad', max(central_density_history)
#    central_density_history=central_density_history/central_density_history[0]
#    print 'max densidad relativa', max(central_density_history)
def relative_error(x,y):
    return (x-y)/y
def selected_values(file_mesa,file_evtwin):
    with open(file_mesa,'r') as f:
        for i in range(6):
            line=f.readline()
        core_radius_mesa=float(f.readline())

    with open(file_evtwin,'r') as f:
        for i in range(6):
            line=f.readline()
        core_radius_evtwin=float(f.readline())        
    mass_evtwin,radius_evtwin,temperature_evtwin,density_evtwin,pressure_evtwin,luminosity_evtwin=np.loadtxt(file_evtwin,skiprows=7,unpack=True, usecols=range(6),delimiter=',')
    mass_mesa,radius_mesa,temperature_mesa,density_mesa,pressure_mesa,luminosity_mesa=np.loadtxt(file_mesa,skiprows=7,unpack=True, usecols=range(6),delimiter=',')
    luminosity_evtwin=luminosity_evtwin*3.8418/3.839
    luminosity_mesa=luminosity_mesa*3.8418/3.839
    r_b,rc_b, pres_b, tc_b, teff_b, lumin_b=benchmark_selected_values(composition)
    print 'Selected values for mesa, error, evtwin, error, benchmark'
    print 'Radius',radius_mesa[-1],relative_error(radius_mesa[-1],r_b),radius_evtwin[-1],relative_error(radius_evtwin[-1],r_b),r_b
    print 'Core radius',core_radius_mesa,relative_error(core_radius_mesa,rc_b),core_radius_evtwin,relative_error(core_radius_evtwin,rc_b),rc_b
    print 'Central pressure', pressure_mesa[0],relative_error(pressure_mesa[0],pres_b/10.0),pressure_evtwin[0],relative_error(pressure_evtwin[0],pres_b/10.0),pres_b/10.0
    print 'Central Temperature', temperature_mesa[0],relative_error(temperature_mesa[0],tc_b),temperature_evtwin[0],relative_error(temperature_evtwin[0],tc_b),tc_b
    print 'Effective Temperature', temperature_mesa[-1],relative_error(temperature_mesa[-1],teff_b),temperature_evtwin[-1],relative_error(temperature_evtwin[-1],teff_b),teff_b
    print 'Luminosity', luminosity_mesa[-1],relative_error(luminosity_mesa[-1],lumin_b),luminosity_evtwin[-1],relative_error(luminosity_evtwin[-1],lumin_b),lumin_b
def benchmark_core_radius(composition):
    mass=[]
    radius=[]
    with open('struct_{0}.dat'.format(composition), 'r') as read_obj:
        csv_reader = reader(read_obj)
        for i in range(9):
            next(csv_reader)
        for line in csv_reader:
            row=line[0].split('  ')
    #        print row[0],row[1]
            mass.append(float(row[0]))
            radius.append(float(row[1]))
    for i in range(len(mass)):
        if mass[i]>0.3:
            index=i
            break
    return radius[index]    
def benchmark_selected_values(composition,verbose=False):       
    mass=[]
    radius=[]
    temperature=[]
    density=[]
    pressure=[]
    luminosity=[]
    with open('struct_{0}.dat'.format(composition), 'r') as read_obj:
        csv_reader = reader(read_obj)
        for i in range(9):
            next(csv_reader)
        for line in csv_reader:
            row=line[0].split('  ')
    #        print row[0],row[1]
            mass.append(float(row[0]))
            radius.append(float(row[1]))
            temperature.append(float(row[2]))
            density.append(float(row[3]))
            pressure.append(float(row[4]))
            luminosity.append(float(row[5]))
    if verbose==True:
        print 'Selected values for benchmark-{0}'.format(composition)
        print 'Radius',radius[-1]
        print 'Core radius',benchmark_core_radius(composition)
        print 'Central pressure', pressure[0]
        print 'Central Temperature', temperature[0]
        print 'Effective Temperature', temperature[-1]
        print 'Luminosity', luminosity[-1]
        print 'WARNING: Pressure is in dyn/cm^2'
        print 'WARNING: LSun=3.8418e33 erg/s, whereas the AMUSE unit module uses LSun=3.839e33 erg/s'

    return radius[-1],benchmark_core_radius(composition), pressure[0], temperature[0], temperature[-1], luminosity[-1]

def simple_composition_evtwin_vs_mesa(file_evtwin1,file_evtwin2,file_mesa1,file_mesa2,j,composition,xlim=1):
    figure = pyplot.figure(figsize=(16, 8))
    with open(file_evtwin1,'r') as f:
        for i in range(4):
            line=f.readline()
        age_evtwin=float(f.readline())
        line=f.readline()
        
    with open(file_mesa1,'r') as f:
        for i in range(4):
            line=f.readline()
        age_mesa=float(f.readline())
        line=f.readline()
    
    radius=np.loadtxt(file_evtwin1,skiprows=7,unpack=True, usecols=1,delimiter=',')
    H,He,C,N,O,Ne,Mg,Si,Fe=np.loadtxt(file_evtwin1,skiprows=7,unpack=True, usecols=(range(6,15)),delimiter=',')
    metals= C+ N+ O+ Ne+ Mg
    fig1=figure.add_subplot(121)
#    pyplot.subplot(1.5, 1, 1)
    fig1.set_xlim(left=0,right=xlim)
    fig1.set_ylim(0.0, 1.0)
    
    fig1.plot(radius,H,'r--')
    fig1.plot(radius,He,'g--')
    fig1.plot(radius,metals,'b--')
    
    
    radius=np.loadtxt(file_mesa1,skiprows=7,unpack=True, usecols=1,delimiter=',')
    H,He,C,N,O,Ne,Mg=np.loadtxt(file_mesa1,skiprows=7,unpack=True, usecols=(range(6,13)),delimiter=',')
    metals= C+ N+ O+ Ne+ Mg
    fig1.plot(radius,H,'r-',label='H')
    fig1.plot(radius,He,'g-',label='He')
    fig1.plot(radius,metals,'b-',label='Metals')
    

    fig1.set_xlabel('r ($R_\odot$), for age $t_{EVtwin}=$'+str(round(age_evtwin,3))+' $t_{MESA}=$'+str(round(age_mesa,3))+str(' $Gyr$'))
    fig1.set_ylabel('Mass fraction')
    
    with open(file_evtwin2,'r') as f:
        for i in range(4):
            line=f.readline()
        age_evtwin=float(f.readline())
        line=f.readline()
        
    with open(file_mesa2,'r') as f:
        for i in range(4):
            line=f.readline()
        age_mesa=float(f.readline())
        line=f.readline()

#    fig_masa.title.set_text('Distribucion cumulativa de masa')
    
    radius=np.loadtxt(file_evtwin2,skiprows=7,unpack=True, usecols=1,delimiter=',')
    H,He,C,N,O,Ne,Mg,Si,Fe=np.loadtxt(file_evtwin2,skiprows=7,unpack=True, usecols=(range(6,15)),delimiter=',')
    metals= C+ N+ O+ Ne+ Mg
    fig2=figure.add_subplot(122)
#    pyplot.subplot(1.5, 1, 1)
    fig2.set_xlim(left=0,right=xlim)
    fig2.set_ylim(0.0, 1.0)
    
    fig2.plot(radius,H,'r--')
    fig2.plot(radius,He,'g--')
    fig2.plot(radius,metals,'b--')
    
    
    radius=np.loadtxt(file_mesa2,skiprows=7,unpack=True, usecols=1,delimiter=',')
    H,He,C,N,O,Ne,Mg=np.loadtxt(file_mesa2,skiprows=7,unpack=True, usecols=(range(6,13)),delimiter=',')
    metals= C+ N+ O+ Ne+ Mg
    fig2.plot(radius,H,'r-',label='H')
    fig2.plot(radius,He,'g-',label='He')
    fig2.plot(radius,metals,'b-',label='Metals')
    

    fig2.set_xlabel('r ($R_\odot$), for age $t_{EVtwin}=$'+str(round(age_evtwin,3))+' $t_{MESA}=$'+str(round(age_mesa,3))+str(' $Gyr$'))
    fig2.set_ylabel('Mass fraction')
    l=pyplot.legend()

    for text in l.get_texts():
        text.set_color("black")
    pyplot.suptitle('Abundances for MESA and EVtwin from {0}'.format(composition), color='black',size=16)
    pyplot.savefig('figures/simplecomp_{0}_{1}.png'.format(composition,str(j)),facecolor='white')
    print 'Saved figures/simplecomp_{0}_{1}.png'.format(composition,str(j))
    
    if not os.path.exists('figures/simple_composition_both_{0}'.format(composition)):
        os.makedirs('figures/simple_composition_both_{0}'.format(composition))
    pyplot.savefig('figures/simple_composition_both_{0}/simplecomp_{1}.png'.format(composition,str(j)))  
    pyplot.close()
    
def simple_composition_only_one(file_evtwin1,file_mesa1,j,composition,xlim=1):
    figure = pyplot.figure(figsize=(12, 6))
    with open(file_evtwin1,'r') as f:
        for i in range(4):
            line=f.readline()
        age_evtwin=float(f.readline())
        line=f.readline()
        
    with open(file_mesa1,'r') as f:
        for i in range(4):
            line=f.readline()
        age_mesa=float(f.readline())
        line=f.readline()

#    fig_masa.title.set_text('Distribucion cumulativa de masa')
    
    radius=np.loadtxt(file_evtwin1,skiprows=7,unpack=True, usecols=1,delimiter=',')
    H,He,C,N,O,Ne,Mg,Si,Fe=np.loadtxt(file_evtwin1,skiprows=7,unpack=True, usecols=(range(6,15)),delimiter=',')
    metals= C+ N+ O+ Ne+ Mg
    fig1=figure.add_subplot(111)
#    pyplot.subplot(1.5, 1, 1)
    fig1.set_xlim(left=0,right=xlim)
    fig1.set_ylim(0.0, 1.0)
    
    fig1.plot(radius,H,'r--')
    fig1.plot(radius,He,'g--')
    fig1.plot(radius,metals,'b--')
    
    
    radius=np.loadtxt(file_mesa1,skiprows=7,unpack=True, usecols=1,delimiter=',')
    H,He,C,N,O,Ne,Mg=np.loadtxt(file_mesa1,skiprows=7,unpack=True, usecols=(range(6,13)),delimiter=',')
    metals= C+ N+ O+ Ne+ Mg
    fig1.plot(radius,H,'r-',label='H')
    fig1.plot(radius,He,'g-',label='He')
    fig1.plot(radius,metals,'b-',label='Metals')
    

    fig1.set_xlabel('r ($R_\odot$), for age $t_{EVtwin}=$'+str(round(age_evtwin,3))+' $t_{MESA}=$'+str(round(age_mesa,3))+str(' $Gyr$'))
    fig1.set_ylabel('Mass fraction')
    
    l=pyplot.legend()
    for text in l.get_texts():
        text.set_color("black")
    pyplot.suptitle('Abundances for MESA and EVtwin from {0}'.format(composition), color='black',size=16)
    pyplot.savefig('figures/simplecomp_{0}_{1}.png'.format(composition,str(j)),facecolor='white') 
    print 'Saved figures/simplecomp_{0}_{1}.png'.format(composition,str(j))

if not os.path.exists('figures'):
    os.makedirs('figures')
    
composition='Invalid answer'
while composition=='Invalid answer':
    composition=raw_input('Choose a composition. Type AGS or GS98\n')
    if composition!='AGS' and composition!='GS98':
        print 'Invalid answer'
        composition='Invalid answer'
    elif composition=='GS98':
        nfiles_mesa=input('How many files do you have in mesa_GS98_results?')
        nfiles_evtwin=input('How many files do you have in evtwin_GS98_results?')
    elif composition=='AGS':
        nfiles_mesa=input('How many files do you have in mesa_AGS_results?')
        nfiles_evtwin=input('How many files do you have in evtwin_AGS_results?')
        
file_evtwin='evtwin_{0}_results/evtwin_{0}_results_1.csv'.format(composition)
file_mesa='mesa_{0}_results/mesa_{0}_results_1.csv'.format(composition)
print 'INITIAL VALUES'
selected_values(file_mesa,file_evtwin)
print
print 'VALUES AT PRESENT'
file_evtwin='evtwin_{0}_results/evtwin_{0}_results_100.csv'.format(composition)
file_mesa='mesa_{0}_results/mesa_{0}_results_100.csv'.format(composition)
selected_values(file_mesa,file_evtwin)

#File sructure dm (MSun), R (RSun), T (K), rho (g/cm^3), P(Pa), L (LSun), h1, he4, he3, c12, n14, o16, ne20, mg24
mass_evtwin,radius_evtwin,temperature_evtwin,density_evtwin,pressure_evtwin,luminosity_evtwin=np.loadtxt(file_evtwin,skiprows=7,unpack=True, usecols=range(6),delimiter=',')
mass_mesa,radius_mesa,temperature_mesa,density_mesa,pressure_mesa,luminosity_mesa=np.loadtxt(file_mesa,skiprows=7,unpack=True, usecols=range(6),delimiter=',')

mass_bench=[]
radius_bench=[]
temp_bench=[]
dens_bench=[]
pres_bench=[]
lumin_bench=[]

# file structure  Mass     Radius     Temp      Rho       Pres       Lumi   
with open('struct_{0}.dat'.format(composition), 'r') as read_obj:
    csv_reader = reader(read_obj)
    for i in range(9):
        next(csv_reader)
    for line in csv_reader:
        row=line[0].split('  ')
#        print row[0],row[1]
        mass_bench.append(float(row[0]))
        radius_bench.append(float(row[1]))
        temp_bench.append(float(row[2]))
        dens_bench.append(float(row[3]))
        pres_bench.append(float(row[4]))
        lumin_bench.append(float(row[5]))

print 'csv was read'
present(mass_mesa,radius_mesa,temperature_mesa,density_mesa,luminosity_mesa,
        mass_evtwin,radius_evtwin,temperature_evtwin,density_evtwin,luminosity_evtwin,
        mass_bench,radius_bench,temp_bench,dens_bench,lumin_bench,'{0}'.format(composition),'present_new_{0}.png'.format(composition))

#%%
radius_history=[]
pressure_history=[]
temperature_history=[]
luminosity_history=[]
central_density_history=[]
central_temperature_history=[]
core_radius_history=[]
time=[]
radius_history_evtwin=[]
pressure_history_evtwin=[]
temperature_history_evtwin=[]
luminosity_history_evtwin=[]
central_density_history_evtwin=[]
central_temperature_history_evtwin=[]
core_radius_history_evtwin=[]
time_evtwin=[]

file_mesa1='mesa_{0}_results/mesa_{0}_results_99.csv'.format(composition)
file_mesa2='mesa_{0}_results/mesa_{0}_results_249.csv'.format(composition)
file_evtwin1='evtwin_{0}_results/evtwin_{0}_results_100.csv'.format(composition)
file_evtwin2='evtwin_{0}_results/evtwin_{0}_results_250.csv'.format(composition)
simple_composition_evtwin_vs_mesa(file_evtwin1,file_evtwin2,file_mesa1,file_mesa2,0,composition,xlim=0.6)

file_mesa1='mesa_{0}_results/mesa_{0}_results_{1}.csv'.format(composition,str(nfiles_mesa-50))
file_mesa2='mesa_{0}_results/mesa_{0}_results_{1}.csv'.format(composition,str(nfiles_mesa-5))
file_evtwin1='evtwin_{0}_results/evtwin_{0}_results_{1}.csv'.format(composition,str(nfiles_evtwin-49))
file_evtwin2='evtwin_{0}_results/evtwin_{0}_results_{1}.csv'.format(composition,str(nfiles_evtwin))
simple_composition_evtwin_vs_mesa(file_evtwin1,file_evtwin2,file_mesa1,file_mesa2,1,composition,xlim=1.0)
simple_composition_only_one(file_evtwin2,file_mesa2,0,composition)
#%%
for j in range(1,nfiles_mesa+1):
    file_mesa='mesa_{0}_results/mesa_{0}_results_{1}.csv'.format(composition,str(j))
    with open(file_mesa,'r') as f:
        for i in range(4):
            line=f.readline()
        age=float(f.readline())
        stellar_type=f.readline()
        core_radius=float(f.readline())

#dm (MSun), R (RSun), T (K), rho (g/cm^3), P(Pa), L (LSun), H, He, C, N, O, Ne, Mg, Si, Fe
    mass,radius,temperature,density,pressure,luminosity,H,He,C,N,O,Ne,Mg=np.loadtxt(file_mesa,skiprows=7,unpack=True, usecols=range(13),delimiter=',')

    
#    estructura(radius_ini,mass_ini,density_ini,luminosity_ini,temperature_ini,radius,mass,density,luminosity,temperature,mass_bench,radius_bench,temp_bench,dens_bench,lumin_bench,str(age)+' Gyr','estructura/estructura_'+str(j)+'.png')
    temperature_history.append(temperature[-1])
    luminosity_history.append(luminosity[-1])
    central_density_history.append(density[0])
    central_temperature_history.append(temperature[0])
    radius_history.append(radius[-1])
    pressure_history.append(pressure[0])
    time.append(age)
    core_radius_history.append(core_radius)

for j in range(1,nfiles_evtwin+1):
    file_evtwin='evtwin_{0}_results/evtwin_{0}_results_{1}.csv'.format(composition,str(j))
    with open(file_evtwin,'r') as g:  
        for i in range(4):
            line=g.readline()
        age_evtwin=float(g.readline())
        stellar_type=g.readline()
        core_radius_evtwin=float(g.readline())
    mass_evtwin,radius_evtwin,temperature_evtwin,density_evtwin,pressure_evtwin,luminosity_evtwin=np.loadtxt(file_evtwin,skiprows=7,unpack=True, usecols=range(6),delimiter=',')
    
    temperature_history_evtwin.append(temperature_evtwin[-1])
    luminosity_history_evtwin.append(luminosity_evtwin[-1])
    central_density_history_evtwin.append(density_evtwin[0])
    central_temperature_history_evtwin.append(temperature_evtwin[0])
    radius_history_evtwin.append(radius_evtwin[-1])
    pressure_history_evtwin.append(pressure_evtwin[0])
    time_evtwin.append(age_evtwin)
    core_radius_history_evtwin.append(core_radius_evtwin)

plot_evolution(radius_history,core_radius_history,pressure_history,temperature_history,
                   luminosity_history,central_density_history,central_temperature_history,time,
                   radius_history_evtwin,core_radius_history_evtwin,pressure_history_evtwin,
                   temperature_history_evtwin,luminosity_history_evtwin,
                   central_density_history_evtwin,central_temperature_history_evtwin,time_evtwin,composition)    
