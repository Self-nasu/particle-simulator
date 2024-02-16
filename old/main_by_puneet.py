#modules, constants, functions:-
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import functions
from functions import *
import time
import math
K=8.9975*(10**9)
G=6.67430*(10**(-11))

#number of particles
n=int(input("numbers of particles= "))


#time intervals
t=float(input("desired time= "))
t1=int(input("enter t for 10^t wehre 10^t is division of 1 sec= "))
dt=(10**t1)**(-1)
t2=int(t/dt)

#array for parameters of n particles
parameters3d=np.zeros((t2+2,n,5))
velocity=np.zeros((t2+2,n,3)) 

#inputing initial dimension of n particles
parameters=[]
for ab in range(0,n):
    str = input(f"\t\tEnter dimensionas of particle {ab+1}: ")
    list1=eval(str)
    parameters.append(list1)
parameters3d[0]=parameters

#inputing velocities of particles
velocity2d=[]
for ab in range(0,n):
    str=input(f"\t\t\tvelocity of particle {ab+1}:")
    list2=eval(str)
    velocity2d.append(list2)
velocity[0]=velocity2d

# introducing all arrays before initiating loop.
#distance
distance_array_3d=np.zeros((t2+1,n,n))
distance_cube_array_3d=np.zeros((t2+1,n,n))

#forces
electrostatic_force_4d=np.zeros((t2+1,3,n,n))
gravitational_force_4d=np.zeros((t2+1,3,n,n))
force_4d=np.zeros((t2+1,3,n,n))

electrostatic_force=np.zeros((t2+1,n,3))
gravitational_force=np.zeros((t2+1,n,3))
force=np.zeros((t2+1,n,3))

net_electrostatic_force=np.zeros((t2+1,n))
net_gravtational_force=np.zeros((t2+1,n))
net_force=np.zeros((t2+1,n))

#speed and linear momentum
speed=np.zeros((t2+1,n))
linear_momentum=np.zeros((t2+1,n))
total_linear_momuntum=np.zeros((t2+1))

#displacement
displacement_3d=np.zeros((t2+1,n,3))
displacement=np.zeros((n,3))

#energy
KE=np.zeros((t2+1,n))
total_KE=np.zeros(t2+1)
EPE_3d=np.zeros((t2+1,n,n))
EPE_2d=np.zeros((t2+1,n))
EPE=np.zeros((t2+1))
GPE_3d=np.zeros((t2+1,n,n))
GPE_2d=np.zeros((t2+1,n))
GPE=np.zeros((t2+1))
TPE_2d=np.zeros((t2+1,n))
TPE=np.zeros((t2+1))
energy=np.zeros((t2+1,n))
total_energy=np.zeros((t2+1))

start=time.time()

for a in range(0,t2+1):
   
    #distances
    #array storing distances of any particle with each particle
    vary_distance_array_3d=np.zeros((n,n))
    for i in range(0,n):
        for j in range(0,n):
            if i==j: 
                distance_array_3d[a][i][j]=0
            else:
                x=parameters3d[a][i][0]-parameters3d[a][j][0]
                y=parameters3d[a][i][1]-parameters3d[a][j][1]
                z=parameters3d[a][i][2]-parameters3d[a][j][2]
                d=((x**2)+(y**2)+(z**2))**(1/2)
                vary_distance_array_3d[i][j]=vary_distance_array_3d[j][i]=d
                distance_array_3d[a][i][j]=np.array(vary_distance_array_3d[i][j])
    #array storing cube of distances of any particle with each particle
    vary_distance_cube_array_3d=np.zeros((n,n))
    for i in range(0,n):
        for j in range(0,n):
            vary_distance_cube_array_3d[i][j]=(vary_distance_array_3d[i][j])**3
            distance_cube_array_3d[a][i][j]=vary_distance_cube_array_3d[i][j]

    #calculating force
    # arrays of electrostatic, gravitational, total force at each particle by each one of particle.
    vary_electrostatic_force_4d=np.zeros((3,n,n))
    vary_gravitational_force_4d=np.zeros((3,n,n))
    vary_force_4d=np.zeros((3,n,n))
    for j in range(0,3):
        for k in range(0,n):
            for l in range(0,n):
                if (distance_cube_array_3d[a][k][l])!=0:
                    vary_electrostatic_force_4d[j][k][l]=(K*parameters3d[a][k][3]*parameters3d[a][l][3]*(parameters3d[a][l][j]-parameters3d[a][k][j]))/(distance_cube_array_3d[a][k][l])
                    vary_gravitational_force_4d[j][k][l]=(G*parameters3d[a][k][4]*parameters3d[a][l][4]*(parameters3d[a][l][j]-parameters3d[a][k][j]))/(distance_cube_array_3d[a][k][l])
                    vary_force_4d[j][k][l]=vary_electrostatic_force_4d[j][k][l]+vary_gravitational_force_4d[j][k][l]
                else: 
                    vary_electrostatic_force_4d[j][k][l]=0
                    vary_gravitational_force_4d[j][k][l]=0
                    vary_force_4d[j][k][l]=0
    electrostatic_force_4d[a]=np.array(vary_electrostatic_force_4d)
    gravitational_force_4d[a]=np.array(vary_gravitational_force_4d)
    force_4d[a]=np.array(vary_force_4d)
    # making arrays of electrostatic, gravitational, total force at each particle by all particles agrigatedly.
    vary_electrostatic_force=np.zeros((3,n))
    vary_gravitational_force=np.zeros((3,n))
    vary_force=np.zeros((3,n))
    for i in range(0,3):
        vary_electrostatic_force[i]=np.sum(electrostatic_force_4d[a][i],axis=0)
        vary_gravitational_force[i]=np.sum(gravitational_force_4d[a][i],axis=0)
        vary_force[i]=np.sum(force_4d[a][i],axis=0)
    vary_electrostatic_force=np.transpose(vary_electrostatic_force)
    vary_gravitational_force=np.transpose(vary_gravitational_force)
    vary_force=np.transpose(vary_force)
    electrostatic_force[a]=vary_electrostatic_force
    gravitational_force[a]=vary_gravitational_force
    force[a]=vary_force
    #net electrostatic, graitational, total force on any particle.
    for j in range(0,n):
        net_e=0
        net_g=0
        net_t=0
        for k in range(0,3):
            net_e=net_e+((electrostatic_force[a][j][k])**2)
            net_g=net_g+((gravitational_force[a][j][k])**2)
            net_t=net_t+((force[a][j][k])**2)
        net_electrostatic_force[a][j]=(net_e)**(1/2)
        net_gravtational_force[a][j]=(net_g)**(1/2)
        net_force[a][j]=(net_t)**(1/2)

    #calculating velocities
    #creating array for holding velocities
    vary_velocity=np.zeros((n,3))

    #new velocity after force is applied for 10^(-t) seconds.
    if (a<t2):
        for j in range(0,n):
            for k in range (0,3):
                vary_velocity[j][k]=velocity[a][j][k]+((force[a][j][k]*dt)/parameters3d[a][j][4])
        velocity[a+1]=vary_velocity
    else: break

    #displacement determination
    #making an array for storing displacement
    vary_displacement_3d=np.zeros((n,3))
    if (a<t2):
        #evaluating displacement
        for j in range(0,n):
            for k in range(0,3):
                if (parameters3d[a][j][4]!=0): vary_displacement_3d[j][k]=(velocity[a][j][k]*dt)+((force[a][j][k]*dt*dt)/(2*parameters3d[a][j][4]))
                else: vary_displacement_3d[j][k]=0
        #storing evaluated displacement
        displacement_3d[a+1]=vary_displacement_3d
    else: break

    #editing dimensions:-
    if (a<t2):
        for i in range(0,n):
            for j in range(0,3):
                parameters3d[a+1][i][j]=parameters3d[a][i][j]+displacement_3d[a+1][i][j]
    else: break

    #keeping mass and charge constant
    if (a<t2):
        for i in range(0,n):
            parameters3d[a+1][i][3]=parameters3d[a][i][3]
            parameters3d[a+1][i][4]=parameters3d[a][i][4]
    else: break

    #speed
    if (a<t2):
        #initial speed
        for i in range(0,n):
            s=0
            for j in range(0,3):
                s=s+(velocity[0][i][j]**2)
            speed[0][i]=s
        #finding speed of each particle
        vary_speed=np.zeros((n))
        for j in range(0,n):
            net_v=0
            for k in range(0,3):
                net_v=net_v+((velocity[a][j][k])**2)
            vary_speed[j]=(net_v)**(1/2)
        speed[a+1]=vary_speed
    else: break
    
    #linear momentum
    vary_linear_momentum=np.zeros((n))
    if (a<t2):
        for i in range(0,n):linear_momentum[0]=parameters3d[0][i][4]*speed[0][i]
        for j in range(0,n): vary_linear_momentum[j]=parameters3d[a][j][4]*speed[a][j]
        linear_momentum[a+1]=vary_linear_momentum
    else: break
    
    #KE
    if (a<t2):
        for j in range(0,n):KE[0][j]=((parameters3d[0][j][4]*(speed[0][j]**2))/2)
        for j in range(0,n):KE[a+1][j]=((parameters3d[a+1][j][4]*(speed[a+1][j]**2))/2)
    else: break

    #EPE
    #EPE of each particle at any instant duw to any other particle
    for j in range(0,n):
        for k in range(0,n):
            if (distance_array_3d[a][j][k]!=0):
                EPE_3d[a][j][k]=((K*parameters3d[a][j][3]*parameters3d[a][k][3])/distance_array_3d[a][j][k])
            else: 
                EPE_3d[a][j][k]=0
    #EPE of each particle due to all particles agrigately
    EPE_2d[a]=np.sum(EPE_3d[a],axis=0)


    #GPE
    #GPE of each particle at any istant due to any other particle
    for j in range(0,n):
        for k in range(0,n):
            if (distance_array_3d[a][j][k]!=0):
                GPE_3d[a][j][k]=((G*parameters3d[a][j][4]*parameters3d[a][k][4])/distance_array_3d[a][j][k])
            else: GPE_3d[a][j][k]=0
    #GPE of each particle due to all particles agrigately
    GPE_2d[a]=np.sum(GPE_3d[a],axis=0)

    print(f"hello{a}")


#extras
total_linear_momuntum=np.sum(linear_momentum,axis=1)
total_KE=np.sum(KE,axis=1)
EPE=np.sum(EPE_2d,axis=1)
GPE=np.sum(GPE_2d,axis=1)
TPE_2d=GPE_2d+EPE_2d
TPE=GPE+EPE
energy=TPE_2d+KE
total_energy=np.sum(energy,axis=1)

# time
end=time.time()
print("total time = ",end-start)
print("time for 1 iteration = ",(end-start)/(t2+1))


#calling graph drawing functions
force_plot(net_force,t2,n,dt,t)
hello(n,parameters3d,t2)
coordinates_plot(parameters3d,t2,n,t)
speed_plot(speed,dt,n,t2,t)
linear_momentum_plot(linear_momentum,t2,dt,n,t)
KE_plot(KE,t2,dt,t,n)
TPE_plot(dt,n,t,TPE_2d,t2)
energy_plot(energy,t2,n,t,dt)