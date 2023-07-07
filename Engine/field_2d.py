# Calculation of the trajectory in a given acceleration field using increments

import numpy as np

class Field2D:
    def __init__(self, r_0=[0, 0], v_0=[0, 0], acc_field=lambda x,y: [0, 0], dt=1, t_0 = 0, t_max=100):
        self.r_0 = r_0
        self.v_0 = v_0
        self.acc_field = acc_field
        self.dt = dt
        self.t_max = t_max

        num_points = np.ceil((t_max - t_0) / dt) + 1    # first and last points
        self.idx_arr = np.arange(start = 0, stop = num_points, step = 1)
        self.time_arr = self.idx_arr * dt

        self.x_ar = np.zeros(num_points)
        self.y_ar = np.zeros(num_points)
        self.vx_ar = np.zeros(num_points)
        self.vy_ar = np.zeros(num_points)

        self.num_points = num_points
        self.x_ar[0], self.y_ar[0] = (self.r_0[0], self.r_0[1])
        self.vx_ar[0], self.vy_ar[0] = (self.v_0[0], self.v_0[1])

    def calc_next_point(self, 
                  x, y, 
                  vx, vy):
        ax, ay = self.acc_field(x, y)
        dvx, dvy = [ax, ay] * self.dt
        avg_v = [vx + dvx / 2., vy + dvy / 2.]

        dx, dy = avg_v * self.dt
        x_new, y_new = [x + dx, y + dy]
        vx_new, vy_new = [vx + dvx, vy + dvy]

        return (x_new, y_new, vx_new, vy_new)
    
    def calculate_motion(self):
        for idx in range(self.num_points - 1):
            self.x_ar[idx + 1], self.y_ar[idx + 1], self.vx_ar[idx + 1], self.vy_ar[idx + 1] = \
                self.calc_next_point(x=self.x_ar[idx], y=self.y_ar[idx], 
                               vx=self.vx_ar[idx], vy=self.vy_ar[idx])
    


