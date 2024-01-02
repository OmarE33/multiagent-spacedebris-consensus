import spaceObject
import agent
import numpy as np
import random

mu = 3.986e5  #Earth's gravitational parameter [km^3/s^2]

#orbital velocity calc
def orbital_velocity(rvec, a):
    rmag = np.linalg.norm(rvec)
    return np.sqrt(mu*(2/rmag - 1/a))

#create debris object initial state in MEO
def generate_random_objects(N, agent_num):
    spaceob = []
    agent_count = 0
    for i in range(N):
        #random position vector
        r = random.uniform(9000, 12000)
        theta = random.uniform(0, 2*np.pi)
        alpha = random.uniform(0, 2*np.pi)
        z = r*np.cos(alpha)
        y = r*np.sin(theta)*np.sin(alpha)
        x = r*np.cos(theta)*np.sin(alpha)
        position = np.array([x, y, z])
        rmag = np.linalg.norm(position)

        #random velocity vector
        a = rmag #from position
        v = orbital_velocity(position, a)
        pos_uvec = position / np.linalg.norm(position)
        randomvec = np.array([random.choice([0, 1, -1]), random.choice([0, 1, -1]), random.choice([0, 1, -1])])
        v_uvec = np.cross(position, randomvec) / np.linalg.norm(np.cross(position, randomvec))
        v_vec = v*v_uvec
        vx = v_vec[0]
        vy = v_vec[1]
        vz = v_vec[2]
        velocity = np.array([vx, vy, vz])

        #create object (agents and debris)
        if agent_count < agent_num:
            ag = agent.agent(position, velocity)
            agent_count += 1
            spaceob.append(ag)
        else:
            deb = spaceObject.spaceObject(position, velocity)
            spaceob.append(deb)
    return spaceob