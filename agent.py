import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import spaceObject


class agent(spaceObject.spaceObject): 
    num_agents = 0

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
        agent.num_agents += 1
    
    def make_adj_mat(self):
        tot = len(agent.GLOBAL_INSTANCE_LIST)
        n = agent.num_agents
        m = tot - n
        spaceObject.spaceObject.adj_mat = np.zeros((tot, tot))
        for i in range(n):
            for j in range(m):
                ## MODIFY CONDITIONS FOR CONSENSUS LATER
                if i == j:
                    spaceObject.spaceObject.adj_mat[i][j+n] = 1
