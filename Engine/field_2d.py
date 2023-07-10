# Calculation of the trajectory in a given acceleration field using increments

import numpy as np

class Field2D:
    def __init__(self, r_0=[0, 0], v_0=[0, 0], acc_field=lambda x,y: [0, 0], dt=1, t_0 = 0, t_max=100):
        self.r_0 = r_0
        self.v_0 = v_0
        self.acc_field = acc_field
        self.dt = dt
        self.t_max = t_max

        num_points = int(np.ceil((t_max - t_0) / dt) + 1)    # first and last points
        self.idx_arr = np.arange(start = 0, stop = num_points, step = 1)
        self.time_arr = self.idx_arr * dt

        self.x_ar = np.zeros(num_points)
        self.y_ar = np.zeros(num_points)
        self.vx_ar = np.zeros(num_points)
        self.vy_ar = np.zeros(num_points)

        self.num_points = num_points
        self.x_ar[0], self.y_ar[0] = (self.r_0[0], self.r_0[1])
        self.vx_ar[0], self.vy_ar[0] = (self.v_0[0], self.v_0[1])

        self.orbit_calculated = False

    def calc_next_point(self, 
                  x, y, 
                  vx, vy):
        ax, ay = self.acc_field(x, y)
        dvx, dvy = [ax * self.dt , ay * self.dt ]
        avg_vx, avg_vy = [vx + dvx / 2., vy + dvy / 2.]

        dx, dy = [self.dt * avg_vx, self.dt * avg_vy]
        x_new, y_new = [x + dx, y + dy]
        vx_new, vy_new = [vx + dvx, vy + dvy]

        return (x_new, y_new, vx_new, vy_new)
    
    def calculate_motion(self):
        for idx in range(self.num_points - 1):
            self.x_ar[idx + 1], self.y_ar[idx + 1], self.vx_ar[idx + 1], self.vy_ar[idx + 1] = \
                self.calc_next_point(x=self.x_ar[idx], y=self.y_ar[idx], 
                               vx=self.vx_ar[idx], vy=self.vy_ar[idx])
        
        self.orbit_calculated = True
    
    # calculate the orbital period (primitive)
    # epsilon is the criterion for "sufficiently close"
    def calculate_orbital_period(self, eps=0, skip_points=20, verbose=False):
        if not self.orbit_calculated:
            self.calculate_motion()

        if eps == 0: # set to the half of distanced covered in dt with initial velocity
            eps = (self.v_0[0]**2 + self.v_0[1]**2)**0.5 * self.dt / 2.
        
        # calculate distance to the initial point 
        for idx in range(skip_points, self.num_points):   # skip first few points point
            dist = abs(self.x_ar[idx] - self.x_ar[0]) + abs(self.y_ar[idx] - self.y_ar[0])
            if verbose and dist < 10 * eps: 
                print("Orbital point: ", dist, self.time_arr[idx], idx)
            if dist < eps:
                return (self.time_arr[idx], idx)
            
        return (0.0, 0)
            
        



