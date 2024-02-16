#modules, constants, functions:-
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
K=8.9975*(10**9)
G=6.67430*(10**(-11))

def force_plot(net_force,t2,n,dt,t):
    max=0
    min=0
    for i in range(0,t2+1):
        for j in range(0,n):
            if max<net_force[i][j]: max=net_force[i][j]
            if min>net_force[i][j]: min=net_force[i][j]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    forceY=np.transpose(net_force)
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("Force (N)")
    plt.title("Function line plot of each particle.")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    for i in range(0,n):
        plt.plot(timeX,forceY[i],label=f"particle {i}",color=Color[i])
    plt.legend()
    plt.show()

def coordinates_plot(parameters3d,t2,n,t):
    xmax=0
    xmin=0
    ymax=0
    ymin=0
    for i in range(0,t2+1):
        for j in range(0,n):
            if xmax<parameters3d[i][j][0]: xmax=parameters3d[i][j][0]
            if xmin>parameters3d[i][j][0]: xmin=parameters3d[i][j][0]
            if ymax<parameters3d[i][j][1]: ymax=parameters3d[i][j][1]
            if ymin>parameters3d[i][j][1]: ymin=parameters3d[i][j][1]          
    xalpha=((xmax-xmin)*5)/(t2+1)
    yalpha=((ymax-ymin)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    plt.xlabel("X coordinate")
    plt.ylabel("Y coordinate")
    plt.title("coordinates of each particle")
    plt.xlim(xmin-xalpha,xmax+xalpha)
    plt.ylim(ymin-yalpha,ymax+yalpha)
    X=np.zeros((n,t2+1))
    Y=np.zeros((n,t2+1))
    for i in range(0,t2+1):
        for j in range(0,n):
                X[j][i]=parameters3d[i][j][0]
                Y[j][i]=parameters3d[i][j][1]
    for j in range(0,n):
        plt.plot(X[j],Y[j],label=f"particle {j+1}",color=Color[j])
    plt.legend()
    plt.show()

def speed_plot(speed,dt,n,t2,t):
    max=0
    min=0
    for i in range(0,t2+1):
        for j in range(0,n):
            if max<speed[i][j]: max=speed[i][j]
            if min>speed[i][j]: min=speed[i][j]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    speedY=np.transpose(speed)
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("speed (N)")
    plt.title("speed line plot of each particle.")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    for i in range(0,n):
        plt.plot(timeX,speedY[i],label=f"particle {i+1}",color=Color[i])
    plt.legend()
    plt.show()

def linear_momentum_plot(linear_momentum,t2,dt,n,t):
    max=0
    min=0
    for i in range(0,t2+1):
        for j in range(0,n):
            if max<linear_momentum[i][j]: max=linear_momentum[i][j]
            if min>linear_momentum[i][j]: min=linear_momentum[i][j]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    linear_momentumY=np.transpose(linear_momentum)
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("Linear Momentum ")
    plt.title("Linear Momentum line plot of each particle.")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    for i in range(0,n):
        plt.plot(timeX,linear_momentumY[i],label=f"particle {i+1}",color=Color[i])
    plt.legend()
    plt.show()

def total_linear_momentum_plot(total_linear_momuntum,dt,t,t2):
    max=0
    min=0
    for i in range(0,t2+1):
        if max<total_linear_momuntum[i]: max=total_linear_momuntum[i]
        if min>total_linear_momuntum[i]: min=total_linear_momuntum[i]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("Total Linear Momentum ")
    plt.title("Total Linear Momentum line plot")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    plt.plot(timeX,total_linear_momuntum,color="black")
    plt.legend()
    plt.show()

def KE_plot(KE,t2,dt,t,n):
    max=0
    min=0
    for i in range(0,t2+1):
        for j in range(0,n):
            if max<KE[i][j]: max=KE[i][j]
            if min>KE[i][j]: min=KE[i][j]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    KEY=np.transpose(KE)
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("KE (Joules)")
    plt.title("KE line plot of each particle.")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    for i in range(0,n):
        plt.plot(timeX,KEY[i],label=f"particle {i+1}",color=Color[i])
    plt.legend()
    plt.show()

def total_KE_plot(total_KE,t2,dt,t):
    max=0
    min=0
    for i in range(0,t2+1):
        if max<total_KE[i]: max=total_KE[i]
        if min>total_KE[i]: min=total_KE[i]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("Total KE ")
    plt.title("Total KE line plot")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    plt.plot(timeX,total_KE,color="black")
    plt.legend()
    plt.show()

def EPE_plot(t2,n,EPE_2d,dt,t):
    max=0
    min=0
    for i in range(0,t2+1):
        for j in range(0,n):
            if max<EPE_2d[i][j]: max=EPE_2d[i][j]
            if min>EPE_2d[i][j]: min=EPE_2d[i][j]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    EPE_2dY=np.transpose(EPE_2d)
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("EPE (Joules)")
    plt.title("EPE line plot of each particle.")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    for i in range(0,n):
        plt.plot(timeX,EPE_2dY[i],label=f"particle {i+1}",color=Color[i])
    plt.legend()
    plt.show()

def GPE_plot(t2,t,n,dt,GPE_2d):
    max=0
    min=0
    for i in range(0,t2+1):
        for j in range(0,n):
            if max<GPE_2d[i][j]: max=GPE_2d[i][j]
            if min>GPE_2d[i][j]: min=GPE_2d[i][j]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    GPE_2dY=np.transpose(GPE_2d)
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("GPE (Joules)")
    plt.title("GPE line plot of each particle.")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    for i in range(0,n):
        plt.plot(timeX,GPE_2dY[i],label=f"particle {i+1}",color=Color[i])
    plt.legend()
    plt.show()

def TPE_plot(dt,n,t,TPE_2d,t2):
    max=0
    min=0
    for i in range(0,t2+1):
        for j in range(0,n):
            if max<TPE_2d[i][j]: max=TPE_2d[i][j]
            if min>TPE_2d[i][j]: min=TPE_2d[i][j]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    TPE_2dY=np.transpose(TPE_2d)
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("TPE (Joules)")
    plt.title("TPE line plot of each particle.")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    for i in range(0,n):
        plt.plot(timeX,TPE_2dY[i],label=f"particle {i+1}",color=Color[i])
    plt.legend()
    plt.show()

def total_EPE_plot(EPE,t2,dt,t):
    max=0
    min=0
    for i in range(0,t2+1):
        if max<EPE[i]: max=EPE[i]
        if min>EPE[i]: min=EPE[i]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("Total KE ")
    plt.title("Total KE line plot")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    plt.plot(timeX,EPE,color="black")
    plt.legend()
    plt.show()

def total_GPE_plot(GPE,t2,dt,t):
    max=0
    min=0
    for i in range(0,t2+1):
        if max<GPE[i]: max=GPE[i]
        if min>GPE[i]: min=GPE[i]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("Total KE ")
    plt.title("Total KE line plot")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    plt.plot(timeX,GPE,color="black")
    plt.legend()
    plt.show()

def total_TPE_plot(TPE,t2,dt,t):
    max=0
    min=0
    for i in range(0,t2+1):
        if max<TPE[i]: max=TPE[i]
        if min>TPE[i]: min=TPE[i]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("Total KE ")
    plt.title("Total KE line plot")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    plt.plot(timeX,TPE,color="black")
    plt.legend()
    plt.show()

def energy_plot(energy,t2,n,t,dt):
    max=0
    min=0
    for i in range(0,t2+1):
        for j in range(0,n):
            if max<energy[i][j]: max=energy[i][j]
            if min>energy[i][j]: min=energy[i][j]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    energyY=np.transpose(energy)
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel(" (Joules)")
    plt.title("Energy line plot of each particle.")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    for i in range(0,n):
        plt.plot(timeX,energyY[i],label=f"particle {i+1}",color=Color[i])
    plt.legend()
    plt.show()

def total_energy_plot(total_energy,t2,t,dt):
    max=0
    min=0
    for i in range(0,t2+1):
        if max<total_energy[i]: max=total_energy[i]
        if min>total_energy[i]: min=total_energy[i]
    alpha=((max-min)*5)/(t2+1)
    Color=["blue","red","green","black","magenta"]
    timeX=np.zeros(t2+1)
    for i in range(0,t2+1): timeX[i]=i*dt
    plt.xlabel("Time (s)")
    plt.ylabel("Total Energy (Joules) ")
    plt.title("Total energy line plot")
    plt.xlim(0,t)
    plt.ylim(min-alpha,max+alpha)
    plt.plot(timeX,total_energy,color="black")
    plt.legend()
    plt.show()

def simulation_plot(points,xlim,ylim,initial_delay,delay):
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    point, = ax.plot([], [], 'bo', markersize=10)
    print(point)
    def update_point():
        global points
        for i in range(len(points)):
            x, y = points[i]
            point.set_data(x, y)
            plt.pause(delay)
            point.set_data([], [])
            plt.pause(delay)
    plt.pause(initial_delay)
    update_point()

def cros_product(a,b):
        c=np.zeros(3)
        c[0]=((a[1]*b[2])-(a[2]*b[1]))
        c[1]=((a[2]*b[0])-(a[0]*b[2]))
        c[2]=((a[0]*b[1])-(a[1]*b[0]))
        return c

def condition_for_approching_particle(n,t,velocity,area_vector):
    sum=0
    for i in range(0,3):
        sum=sum+(area_vector[i]*velocity[t][n][i])
    magnitude_area_vector=0
    for i in range(0,3):
        magnitude_area_vector=magnitude_area_vector+(area_vector[i]**2)
    magnitude_velocity=0
    for i in range(0,3):
        magnitude_velocity=magnitude_velocity+(velocity[t][n][i]**2)
    product=magnitude_velocity*magnitude_area_vector
    product=product**(1/2)
    angle=sum/product
    if (angle)<0: return True
    else: return False

def new_velocity_after_collision_fo_wall(t,n,velocity,area_vector,e):
    x_component=((area_vector[0]*velocity[t][n][0]+area_vector[1]*velocity[t][n][1]+area_vector[2]*velocity[t][n][2])/((area_vector[0]**2)+(area_vector[1]**2)+(area_vector[2]**2)))+((velocity[t][n][1]**area_vector[2]-velocity[t][n][2]**area_vector[1])/(((area_vector[0]**2)+(area_vector[1]**2)+(area_vector[2]**2))**(1/2)))
    y_component=((area_vector[0]*velocity[t][n][0]+area_vector[1]*velocity[t][n][1]+area_vector[2]*velocity[t][n][2])/((area_vector[0]**2)+(area_vector[1]**2)+(area_vector[2]**2)))+((velocity[t][n][2]**area_vector[0]-velocity[t][n][0]**area_vector[2])/(((area_vector[0]**2)+(area_vector[1]**2)+(area_vector[2]**2))**(1/2)))
    z_component=((area_vector[0]*velocity[t][n][0]+area_vector[1]*velocity[t][n][1]+area_vector[2]*velocity[t][n][2])/((area_vector[0]**2)+(area_vector[1]**2)+(area_vector[2]**2)))+((velocity[t][n][0]**area_vector[1]-velocity[t][n][1]**area_vector[0])/(((area_vector[0]**2)+(area_vector[1]**2)+(area_vector[2]**2))**(1/2)))
    new_instantanous_velocity=[e*x_component,e*y_component,e*z_component]
    return new_instantanous_velocity

def distance_of_particle_from_wall(area_vector,parameters3d,point_on_plane,magnitude_of_area_vector,t,n):
    dot_product=(area_vector[0]*(parameters3d[t][n][0]-point_on_plane[0]))+(area_vector[1]*(parameters3d[t][n][1]-point_on_plane[1]))+(area_vector[2]*(parameters3d[t][n][2]-point_on_plane[2]))
    distance=dot_product*magnitude_of_area_vector
    return distance

def hello(n, parameters3d, t2):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(n):
        x = []
        y = []
        z = []
        for j in range(t2):
            x.append(parameters3d[j][i][0])
            y.append(parameters3d[j][i][1])  # Corrected
            z.append(parameters3d[j][i][2])  # Corrected

        # Plot the points for this iteration
        ax.scatter(x, y, z)

    # Set labels
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()