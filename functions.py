#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 18:20:47 2020

@author: alex
"""
from amuse.lab import * 
from amuse.units import units
import numpy as np
import csv


def uniform_star(code,code_name,mass,Z,H,composition_file):
    
    """----------------READING THE COMPOSITION FILE------------------"""
    with open(composition_file, mode='r') as infile:
        reader = csv.reader(infile,delimiter='\t')
        for i in range(3):
            next(reader)
        mass_fraction = {rows[1]:(float(rows[3])*Z) for rows in reader}
    metals_keys=['C ','N ','O ','Ne ','Mg ','Si ','S ','Fe ']#greater than 1e-4
    #notice the keys have a blank space character because of the file structure
    metal_sum=0
    print 'Creating a star with uniform composition:'
    for key in metals_keys:
        print key,mass_fraction[key]
        metal_sum+=mass_fraction[key]
    He=(1.0-H-metal_sum)
    print 'Metal sum=',metal_sum
    print 'H=',H
    print 'He=',(1.0-H-metal_sum)
    if code_name=='evtwin':
        """---------------CREATING A STAR WITH DEFAULT STRUCTURE-------------"""
        estrella_default=code.particles.add_particle(Particle(mass=mass)) 
        structure=estrella_default.get_internal_structure()
        zones = estrella_default.get_number_of_zones()
        """----------------CREATING A STAR WITH OUR COMPOSITION--------------"""
        new_structure=Grid(199)
         
        new_structure.X_Fe=np.ones(199)*  mass_fraction['Fe ']
        new_structure.X_C= np.ones(199)*  mass_fraction['C ']
        new_structure.X_Mg=np.ones(199)*  mass_fraction['Mg ']
        new_structure.X_N= np.ones(199)*  mass_fraction['N ']
        new_structure.X_Ne=np.ones(199)*  mass_fraction['Ne ']
        new_structure.X_O= np.ones(199)*  mass_fraction['O ']
        new_structure.X_Si=np.ones(199)*  mass_fraction['Si ']
        new_structure.X_H= np.ones(199)*  H
        new_structure.X_He=np.ones(199)*  He
        
        new_structure.d_mass=structure.d_mass
        new_structure.entropy=structure.entropy
        new_structure.luminosity=structure.luminosity
        new_structure.mass=structure.mass
        new_structure.molecular_weight=structure.molecular_weight
        new_structure.pressure=structure.pressure
        new_structure.radius=structure.radius
        new_structure.rho=structure.rho
        new_structure.temperature=structure.temperature
        
        sun=code.new_particle_from_model(new_structure)
        
        print 'Created an EVtwin star with uniform composition'
        return sun
    
    elif code_name=='mesa':
        """---------------CREATING A STAR WITH DEFAULT STRUCTURE-------------"""
        estrella_default=code.particles.add_particle(Particle(mass=mass)) 
    
        #%%
        """----------------CREATING A STAR WITH OUR COMPOSITION--------------"""
        zones = estrella_default.get_number_of_zones()
        composition = np.ones(zones)
        sun = code.new_particle_from_model(dict(
            mass=(estrella_default.get_cumulative_mass_profile(
                number_of_zones=zones) * estrella_default.mass),
            radius=estrella_default.get_radius_profile(
                number_of_zones=zones),
            rho=estrella_default.get_density_profile(
                number_of_zones=zones),
            temperature=estrella_default.get_temperature_profile(
                number_of_zones=zones),
            luminosity=estrella_default.get_luminosity_profile(
                number_of_zones=zones),
            X_H=composition * H,
            X_He=composition* He,
            X_C=composition * mass_fraction['C '],
            X_N=composition * mass_fraction['N '],
            X_O=composition * mass_fraction['O '],
            X_Ne=composition* mass_fraction['Ne '],
            X_Mg=composition* mass_fraction['Mg '],
            X_He3=composition * 0.0,
            X_Si=composition * mass_fraction['Si '],
            X_Fe=composition * mass_fraction['Fe ']), 0.0 | units.Myr)
            
        print 'Created a MESA star with uniform composition'
        return sun
    else:
        print 'Error in uniform_star: Invalid code name.'
        print 'Valid names: "mesa", "evtwin"'

def core_radius(sun,code):
    if code=='mesa':
        number_of_zones=sun.get_number_of_zones()
        radius=sun.get_radius_profile(number_of_zones)
        dm=sun.get_mass_profile(number_of_zones)
        mass=np.cumsum(dm)
        for i in range(len(mass)):
            if mass[i]>0.1:
                index=i
                break
        return radius[index]  
    elif code=='evtwin':
        structure=sun.get_internal_structure()
#        number_of_zones=sun.get_number_of_zones()
        radius=structure.radius
        dm=structure.d_mass.value_in(units.MSun)
        mass=np.cumsum(dm)
        for i in range(len(mass)):
            if mass[i]>0.1:
                index=i
                break
        return radius[index]
    else:
        print 'Error in core_radius: Invalid code name'
        return -1

def structure_from_star(star):
#    star.evolve_for(age)

    radius_profile = star.get_radius_profile()
    density_profile = star.get_density_profile()
    if hasattr(star, "get_mass_profile"):
        mass_profile = star.get_mass_profile() * star.mass
    else:
            radii_cubed = radius_profile**3
            radii_cubed.prepend(0 | units.m**3)
            mass_profile = (
                (4.0 / 3.0 * np.pi) * density_profile
                * (radii_cubed[1:] - radii_cubed[:-1])
            )
            print("Derived mass profile from density and radius.")

    return dict(
        radius=radius_profile.as_quantity_in(units.RSun),
        density=density_profile,
        mass=mass_profile,
        luminosity=star.get_luminosity_profile(),
        temperature=star.get_temperature_profile(),
        pressure=star.get_pressure_profile(),
        composition=star.get_chemical_abundance_profiles(),
        species_names=star.get_names_of_species()
    )

def write_header(f):
    f.write('This file is part of the end-of-degree project by Alexander Olza Rodriguez (UPV/EHU).\n')
    
def file_description(f,composition,code_name):
    f.write('Results from script run_{1}.py obtained with {1} with composition {0}\n'.format(composition,code_name))
    f.write('Row 1: Age (Gyr).Row 2: Stellar type. Row 3: coreradius (RSun). Rest of the rows are divided in columns as follows:\n')
    if code_name=='mesa':
        f.write('dm (MSun), R (RSun), T (K), rho (g/cm^3), P(Pa), L (LSun), h1, he4, he3, c12, n14, o16, ne20, mg24\n')
    elif code_name=='evtwin':
        f.write('dm (MSun), R (RSun), T (K), rho (g/cm^3), P(Pa), L (LSun), H, He, C, N, O, Ne, Mg, Si, Fe\n')
    else:
        print 'Error writing file_description: Invalid code name.'
        print 'Valid options "mesa", "evtwin"'
        
def save_results(sun,filename,composition,code_name):
    if code_name=='evtwin':
        with open(filename,'w') as file:
            writer = csv.writer(file)
            write_header(file)
            file_description(file,composition,code_name)
            writer.writerow([sun.age.value_in(units.Gyr)])
            writer.writerow([sun.stellar_type])
            writer.writerow([core_radius(sun,'evtwin').value_in(units.RSun)])
            structure=sun.get_internal_structure()
            for i in range(len(structure)):
                writer.writerow([
                        structure.d_mass.value_in(units.MSun)[i],   structure.radius[i].value_in(units.RSun),
                        structure.temperature[i].value_in(units.K), structure.rho[i].value_in(units.g/(units.cm)**3),
                        structure.pressure[i].value_in(units.Pa),   structure.luminosity[i].value_in(units.LSun),
                        structure.X_H[i],structure.X_He[i],structure.X_C[i],structure.X_N[i],structure.X_O[i],
                        structure.X_Ne[i],structure.X_Mg[i],structure.X_Si[i],structure.X_Fe[i]
                        ])
    elif code_name=='mesa':
        data=structure_from_star(sun)
        with open(filename,'w') as file:
            writer = csv.writer(file)
            write_header(file)
            file_description(file,composition,code_name)
            writer.writerow([sun.age.value_in(units.Gyr)])
            writer.writerow([sun.stellar_type])
            writer.writerow([core_radius(sun,'mesa').value_in(units.RSun)])

            for i in range(len(data['radius'])-1):
                writer.writerow([
                        data['mass'][i].value_in(units.MSun), data['radius'][i].value_in(units.RSun),
                        data['temperature'][i].value_in(units.K), data['density'][i].value_in(units.g/(units.cm)**3),
                        data['pressure'][i].value_in(units.Pa), data['luminosity'][i].value_in(units.LSun),
                        data['composition'][data['species_names'].index('h1')][i],
                        data['composition'][data['species_names'].index('he4')][i],
                        data['composition'][data['species_names'].index('he3')][i],
                        data['composition'][data['species_names'].index('c12')][i],
                        data['composition'][data['species_names'].index('n14')][i],
                        data['composition'][data['species_names'].index('o16')][i],
                        data['composition'][data['species_names'].index('ne20')][i],
                        data['composition'][data['species_names'].index('mg24')][i]])
            writer.writerow([
                data['mass'][-1].value_in(units.MSun), data['radius'][-1].value_in(units.RSun),
                data['temperature'][-1].value_in(units.K), data['density'][-1].value_in(units.g/(units.cm)**3),
                data['pressure'][-1].value_in(units.Pa), data['luminosity'][-1].value_in(units.LSun),
                data['composition'][data['species_names'].index('h1')][-1],
                data['composition'][data['species_names'].index('he4')][-1],
                data['composition'][data['species_names'].index('he3')][-1],
                data['composition'][data['species_names'].index('c12')][-1],
                data['composition'][data['species_names'].index('n14')][-1],
                data['composition'][data['species_names'].index('o16')][-1],
                data['composition'][data['species_names'].index('ne20')][-1],
                data['composition'][data['species_names'].index('mg24')][-1]])
                    
    else:
        print 'Error in save_results: Invalid code name'
        return -1

def cycle(i,n,t,sun,path,composition,code,code_name):
    tipo=sun.stellar_type
    present=4.5395 |units.Gyr
    ask=1
    answer='Invalid answer'
    while sun.stellar_type==tipo:
        i+=1
        filename=path+str(i)+'.csv'
        t+=present/n #Timestep for saving the results (NOT internal timestep of mesa/evtwin)
        code.evolve_model(t)
        
        if sun.age>12.037|units.Gyr and composition=='GS98' and code_name=='mesa':
            print "Calculations are extremely slow at this point. Do you want to continue?"
            while answer=='Invalid answer' and ask==1:
                answer=raw_input("Type yes or no. (Case sensitive)")
                if answer=='no':
                    return i,t
                elif answer==yes:
                    break
                else:
                    print 'Invalid answer'
                    answer='Invalid answer'
            ask=0
        if sun.age>12.069|units.Gyr and composition=='AGS' and code_name=='mesa':
            print "Calculations are extremely slow at this point. Do you want to continue?"
            while answer=='Invalid answer' and ask==1:
                answer=raw_input("Type yes or no. (Case sensitive)")
                if answer=='no':
                    return i,t
                elif answer==yes:
                    break
                else:
                    print 'Invalid answer'
                    answer='Invalid answer'
            ask=0
        if sun.age>12.029|units.Gyr and composition=='GS98' and code_name=='evtwin':
            print "Calculations are extremely slow at this point. Do you want to continue?"
            while answer=='Invalid answer' and ask==1:
                answer=raw_input("Type yes or no. (Case sensitive)")
                if answer=='no':
                    return i,t
                elif answer==yes:
                    break
                else:
                    print 'Invalid answer'
                    answer='Invalid answer'
            ask=0
        if sun.age>12.029|units.Gyr and composition=='AGS' and code_name=='evtwin':
            print "Calculations are extremely slow at this point. Do you want to continue?"
            while answer=='Invalid answer' and ask==1:
                answer=raw_input("Type yes or no. (Case sensitive)")
                if answer=='no':
                    return i,t
                elif answer==yes:
                    break
                else:
                    print 'Invalid answer'
                    answer='Invalid answer'
            ask=0
        save_results(sun, filename,composition,code_name)

        print 'Age=',sun.age.as_quantity_in(units.Gyr)
        
    print 'Now stellar type is {0} and age is {1}'.format(sun.stellar_type,sun.age)
    return i,t
