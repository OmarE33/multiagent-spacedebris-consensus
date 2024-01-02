import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class spaceObject:
    mu = 3.986e5  #Earth's gravitational parameter [km^3/s^2]
    adj_mat = np.array([])
    GLOBAL_INSTANCE_LIST = []

    def __init__(self, position, velocity):
        self.position = position #km    
        self.velocity = velocity #km
        self.vmag = np.linalg.norm(velocity) #km/s
        self.rmag = np.linalg.norm(position) #km/s
        self.X = None
        self.Y = None
        self.Z = None
        self.a = None
        self.Vx = None
        self.Vy = None
        self.Vz = None


    def get_orbital_elems(self):
        #semimajor axis
        self.a = 1 / (2/self.rmag - (self.vmag**2)/spaceObject.mu)

        #angular momentum
        h = np.cross(self.position, self.velocity)
        hmag = np.linalg.norm(h)

        #inclination
        i = np.arccos(h[2] / hmag)

        #RAAN
        K = np.array([0, 0, 1])
        N = np.cross(K, h)
        Nmag = np.linalg.norm(N)
        if N[1] < 0:
            RAAN = 2*np.pi - np.arccos(N[0] / Nmag)
        else:
            RAAN = np.arccos(N[0] / Nmag)

        #eccentricity
        e = (1 / spaceObject.mu) * ((self.vmag**2 - (spaceObject.mu / self.rmag))*self.position - np.dot(self.position, self.velocity)*self.velocity)
        emag = np.linalg.norm(e)

        #argument of periapsis
        if e[2] < 0:
            omega = 2*np.pi - np.arccos(np.dot(N, e) / (Nmag * emag))
        else:
            omega = np.arccos(np.dot(N, e) / (Nmag * emag))

        #true anomaly
        if np.dot(self.position, self.velocity) < 0:
            nu = 2*np.pi - np.arccos(np.dot(self.position / self.rmag, e / emag))  
        else:  
            nu = np.arccos(np.dot(self.position / self.rmag, e / emag))

        return {"angular momentum vec": h, 
                "angular momentum mag": hmag, 
                "inclination": np.rad2deg(i), 
                "RAAN": np.rad2deg(RAAN), 
                "eccentricity": e,
                "eccentricity mag": emag, 
                "arugment of periapsis": np.rad2deg(omega), 
                "true anomaly": np.rad2deg(nu)}

    def propogate(self, trange):
        def model_2BP(state, t):
            x = state[0]
            y = state[1]
            z = state[2]
            x_dot = state[3]
            y_dot = state[4]
            z_dot = state[5]
            x_ddot = -spaceObject.mu*x / (x**2 + y**2 + z**2)**(3 / 2)
            y_ddot = -spaceObject.mu*y / (x**2 + y**2 + z**2)**(3 / 2)
            z_ddot = -spaceObject.mu*z / (x**2 + y**2 + z**2)**(3 / 2)
            dstate_dt = [x_dot, y_dot, z_dot, x_ddot, y_ddot, z_ddot]
            return dstate_dt
        
        # Initial Conditions
        X_0 = self.position[0]  # [km]
        Y_0 = self.position[1]  # [km]
        Z_0 = self.position[2]  # [km]
        VX_0 = self.velocity[0]  # [km/s]
        VY_0 = self.velocity[1]  # [km/s]
        VZ_0 = self.velocity[2]  # [km/s]
        state_0 = [X_0, Y_0, Z_0, VX_0, VY_0, VZ_0]

        # Time Array
        t = np.linspace(0, trange, 200)  # Simulates for a time period of trange hours [s]

        # Solving ODE
        sol = odeint(model_2BP, state_0, t)
        self.X = np.array(sol[:, 0])  # X-coord [km] of satellite over time interval 
        self.Y = np.array(sol[:, 1])  # Y-coord [km] of satellite over time interval
        self.Z = np.array(sol[:, 2])  # Z-coord [km] of satellite over time interval
        self.Vx = np.array(sol[:, 3])  # X-coord vel [km/s] of satellite over time interval 
        self.Vy = np.array(sol[:, 4])  # Y-coord vel [km/s] of satellite over time interval
        self.Vz = np.array(sol[:, 5])  # Z-coord vel [km/s] of satellite over time interval

        #update global object list
        spaceObject.GLOBAL_INSTANCE_LIST.append(self)

        return {"x": self.X, "y": self.Y, "z": self.Z}
    
    def plot_orbit(self):
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
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(X_Earth, Y_Earth, Z_Earth, color='blue', alpha=0.7)
        ax.plot3D(self.X, self.Y, self.Z, 'black')
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

    def plot_orbit_all(self):
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
        for i in spaceObject.GLOBAL_INSTANCE_LIST:
            ax.plot(i.X, i.Y, i.Z)
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