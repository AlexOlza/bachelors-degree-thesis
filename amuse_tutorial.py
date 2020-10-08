"""
This script is part of the end-of-degree project by Alexander Olza Rodriguez (UPV/EHU).

This is a brief tutorial on how to get started with AMUSE (https://amusecode.github.io/copyright). 

We will work through an exercise similar to one proposed in the book "ASTROPHYSICAL RECIPES: The art of AMUSE" (Portegies Zwart, S. & McMillan, S.L.W., 2018).

The goal is to integrate the motion of the planets in our Solar System.
"""

from amuse.units import units, nbody_system
from amuse.lab import ph4,new_solar_system

import  numpy

import tutorial_plots 

def append_position(planet_positions,bodies,i):
	planet_positions.append([bodies[i].x.value_in(units.AU),bodies[i].y.value_in(units.AU),bodies[i].z.value_in(units.AU)])
	return planet_positions

##############################  MAIN STARTS HERE  ##############################################

bodies=new_solar_system() 
"""
 Now this variable (which is a particle set) contains 9 particles with attributes ID key, name, mass, radius, position (vector attribute x,y,z) and velocity (vx,vy,vz)
 We could have initialized the particles manually, for example:
     
 bodies=Particles(9)

 bodies.mass=[1,1.6e-7,2.4e-6,3e-6,3.2e-7,9.5e-4,2.9e-4,4.4e-5,5.1e-5] |units.MSun
 
 bodies[0].position=[0,0,0]|units.AU
 ...
 bodies[8].position=[...,...,...]|units.AU
 
 bodies[0].velocity=[...,...,...]|units.kms
 ...
 bodies[8].velocity=[...,...,...]|units.kms
 
 The rest of the attributes are not necessary for our task.
 
 Let's take a look at the content of bodies:
"""
 
print(bodies)

"""
Nbody codes work internally in nbody units. Our initial conditions are given in SI units. 
We will ask AMUSE to make the conversion before passing the input to the nbody code. 
The first step to do that is to define a converter between these two systems using the method nbody_to_si from the module nbody_system. 
It takes as input any pair of non equivalent quantities (i.e. mass-lenght, luminosity-velocity...). The choice 1|units.MJupiter,1|units.AU is arbitrary.
"""


converter=nbody_system.nbody_to_si(1|units.MJupiter,1|units.AU)

"""
ph4 is a Nbody solver. There are some others in AMUSE, but this is only a tutorial so the choice has been arbitrary.
With the following line we initialize a "worker" os this code via MPI/Sockets, and tell AMUSE to use the converter before passing input to ph4.
The worker is a completely separate process from this python script, running independently.
It will make a remote copy of all the variables passed as input, and won't modify the local ones (i.e. bodies).
All communications (such as retrieving data) between ph4 and this script must be written explicitly.
"""

gravity=ph4(converter) 
gravity.particles.add_particles(bodies) 

"""
We have passed the initial conditions to ph4. The unit conversion has been done in the background. Gravity is not acting yet.
"""

"""
For validation, we will be checking energy conservation. The energy of the REMOTE COPY of bodies (initially identical to the local copy) is:
"""

Energy_init=gravity.kinetic_energy+gravity.potential_energy
E=[Energy_init] #in nbody units
dE=[0]  # relative energy error (adimensional)

#We initialize time (in years)
t=0.0 |units.yr
dt=10.0|units.day
t_end=10.0 |units.yr

#We will only plot the orbits of these 4 bodies. The following will be arrays of VALUES, not quantities (explained below)
sun_position,mercury_position,venus_position,earth_position=[],[],[],[]
#We append the initial positions
sun_position=append_position(sun_position,bodies,0)
mercury_position=append_position(mercury_position,bodies,1)
venus_position=append_position(venus_position,bodies,2)
earth_position=append_position(earth_position,bodies,3)

"""
Now we define a communication channel between the remote copy, gravity.particles, and our local particle set, bodies.
This is necessary to access the attributes of the particles in ph4, updating the positions and velocities of each body.
However, it's not necessary to get the energy of the system as a whole (which we have already done).
"""
channel_to_bodies=gravity.particles.new_channel_to(bodies)


"""
Finally, we let gravity start to act and we save the positions periodically (each 10 days).
"""

while t<t_end:	
	gravity.evolve_model(t+dt) #This line makes ph4 integrate the equations of the nbody problem
	t+=dt
	channel_to_bodies.copy() #Now the variable bodies contains up-to-date information

	#Retrieving properties from ph4:
	E.append(gravity.kinetic_energy+gravity.potential_energy)
	dE.append((E[0]-E[-1])/E[0])
	#Saving positions:
	sun_position=append_position(sun_position,bodies,0)
	mercury_position=append_position(mercury_position,bodies,1)
	venus_position=append_position(venus_position,bodies,2)
	earth_position=append_position(earth_position,bodies,3)
    
gravity.stop()#This line terminates the ph4 process
"""
Now we save and quick-plot the trajectories and the energy error. 
IMPORTANT: 
Each time channel_to_bodies.copy() is executed, we update the attributes of the object bodies. These are quantities (value+unit).
Notice that, when calling append_position(...), we have used the method .value_in(units.AU). 
This transforms quantities to fixed numerical values. It's necessary for plotting because Gnuplot (or Matplotlib) doesn't know what quantities are! 
"""

numpy.savetxt("sun",sun_position)	
numpy.savetxt("mercury",mercury_position)	
numpy.savetxt("venus",venus_position)	
numpy.savetxt("earth",earth_position)	
numpy.savetxt("energy_error",numpy.abs(dE))

tutorial_plots.plot_error(dE)
tutorial_plots.plot_trajectories(sun_position,mercury_position,venus_position,earth_position)


