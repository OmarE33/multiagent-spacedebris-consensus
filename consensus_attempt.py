import spaceObject
import agent
import helper
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

mu = spaceObject.spaceObject.mu
N = 2 #number of space object
agent_num = 1 #number of agents in space


def xdot(xs, t):                     # This function calculates the x_dot term of every agent for all t
    x_vals = np.zeros((6*N,))
    for i in range(N):
        if i < agent_num:
            for j in range(agent_num, N):
                if j != i:
                    if spaceObject.spaceObject.adj_mat[i][j] == 1:
                        # x_vals[3*i] += (x[3*j] - x[3*i]) + object_velocities[3*(j - agent_num)]    # Regular LCP + debris velocity to put into debris' reference frame.
                        # x_vals[3*i + 1] += (x[3*j + 1] - x[3*i + 1]) + object_velocities[3*(j - agent_num) + 1]
                        # x_vals[3*i + 2] += (x[3*j + 2] - x[3*i + 2]) + object_velocities[3*(j - agent_num) + 2]
                        x = xs[6*i]
                        y = xs[6*i + 1]
                        z = xs[6*i + 2]
                        xv = xs[6*i + 3]
                        yv = xs[6*i + 4]
                        zv = xs[6*i + 5]
                        xdeb = xs[6*j]
                        ydeb = xs[6*j + 1]
                        zdeb = xs[6*j + 2]
                        xdebv = xs[6*j + 3]
                        ydebv = xs[6*j + 4]
                        zdebv = xs[6*j + 5]

                        x_vals[6*i] += xv
                        x_vals[6*i + 1] += yv
                        x_vals[6*i + 2] += zv 
                        x_vals[6*i + 3] += -mu*x / (x**2 + y**2 + z**2)**(3 / 2)
                        x_vals[6*i + 4] += -mu*y / (x**2 + y**2 + z**2)**(3 / 2)
                        x_vals[6*i + 5] += -mu*z / (x**2 + y**2 + z**2)**(3 / 2)
                        x_vals[6*j] += xdebv 
                        x_vals[6*j + 1] += ydebv 
                        x_vals[6*j + 2] += zdebv 
                        x_vals[6*j + 3] += -mu*xdeb / (xdeb**2 + ydeb**2 + zdeb**2)**(3 / 2)
                        x_vals[6*j + 4] += -mu*ydeb / (xdeb**2 + ydeb**2 + zdeb**2)**(3 / 2)
                        x_vals[6*j + 5] += -mu*zdeb / (xdeb**2 + ydeb**2 + zdeb**2)**(3 / 2)
    return x_vals

#generate space objects, propogate, and plot them all
d = helper.generate_random_objects(N, agent_num)
print(d)
for deb in d:
    deb.propogate(6*3600)
d[0].plot_orbit_all()
a1 = d[0]
a1.make_adj_mat()
print(spaceObject.spaceObject.adj_mat)

x = np.zeros((6*N,)) #create array of size (3N, number of steps taken in propogate function)
object_velocities = np.zeros((6*N))
for i in range(N):
    x[6*i] = spaceObject.spaceObject.GLOBAL_INSTANCE_LIST[i].X[0]
    x[6*i + 1] = spaceObject.spaceObject.GLOBAL_INSTANCE_LIST[i].Y[0]
    x[6*i + 2] = spaceObject.spaceObject.GLOBAL_INSTANCE_LIST[i].Z[0]
    x[6*i + 3] = spaceObject.spaceObject.GLOBAL_INSTANCE_LIST[i].Vx[0]
    x[6*i + 4] = spaceObject.spaceObject.GLOBAL_INSTANCE_LIST[i].Vy[0]
    x[6*i + 5] = spaceObject.spaceObject.GLOBAL_INSTANCE_LIST[i].Vz[0]


# Time Array
trange = 12*3600
t = np.linspace(0, trange, 200)  # Simulates for a time period of trange hours [s]

# Solving ODE
sol = odeint(xdot, x, t)
Xcon = np.array(sol[:, 0])  # X-coord [km] of satellite over time interval 
Ycon = np.array(sol[:, 1])  # Y-coord [km] of satellite over time interval
Zcon = np.array(sol[:, 2])  # Z-coord [km] of satellite over time interval
Xdeb = np.array(sol[:, 6])
Ydeb = np.array(sol[:, 7])
Zdeb = np.array(sol[:, 8])
print(Xcon.shape)


# Setting up Spherical Earth to Plot
N = 50
phi = np.linspace(0, 2*np.pi, N)
theta = np.linspace(0, np.pi, N)
theta, phi = np.meshgrid(theta, phi)

r_Earth = 6378.14  # Average radius of Earth [km]
X_Earth = r_Earth * np.cos(phi) * np.sin(theta)
Y_Earth = r_Earth * np.sin(phi) * np.sin(theta)
Z_Earth = r_Earth * np.cos(theta)

# Plotting Earth and Orbit
fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
ax.plot(Xcon, Ycon, Zcon)
ax.plot(Xdeb, Ydeb, Zdeb)
ax.plot_surface(X_Earth, Y_Earth, Z_Earth, color='blue', alpha=0.7)
ax.view_init(30, 145)  # Changing viewing angle (adjust as needed)
plt.title('Two-Body Orbit')
ax.set_xlabel('X [km]')
ax.set_ylabel('Y [km]')
ax.set_zlabel('Z [km]')

# Make axes limits
xyzlim = np.array([ax.get_xlim3d(), ax.get_ylim3d(),      
                ax.get_zlim3d()]).T
XYZlim = np.asarray([min(xyzlim[0]), max(xyzlim[1])])
ax.set_xlim3d(XYZlim)
ax.set_ylim3d(XYZlim)
ax.set_zlim3d(XYZlim * 3/4)
plt.show()

