# Calculation of the trajectory in a given acceleration field using increments

import numpy as np

class Field2D:
    def __init__(self, r_0=[0, 0], v_0=[0, 0], acc_field=lambda x,y: [0, 0], dt=1, t_0 = 0, t_max=100):
        self.r_0 = r_0
        self.x_0, self.y_0 = self.r_0

        self.v_0 = v_0
        self.v_x_0, self.v_y_0 = self.v_0

        self.acc_field = acc_field

        self.dt = dt
        self.t_max = t_max

        num_points = int(np.ceil((t_max - t_0) / dt) + 1)    # first and last points
        self.idx_arr = np.arange(start = 0, stop = num_points, step = 1)
        self.time_arr = self.idx_arr * dt

        self.x_ar = np.zeros(num_points)
        self.y_ar = np.zeros(num_points)

        self.r_ar = np.zeros(num_points)
        self.th_ar = np.zeros(num_points)

        self.vx_ar = np.zeros(num_points)
        self.vy_ar = np.zeros(num_points)

        self.num_points = num_points
        self.x_ar[0], self.y_ar[0] = (self.x_0, self.y_0)
        self.vx_ar[0], self.vy_ar[0] = (self.v_0[0], self.v_0[1])

        self.orbit_calculated = False
        self.orbital_period = 0

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
        self.calculate_distances()
    
    # calculate the array of euclidean distance to the coordinates origin
    def calculate_distances(self):
        self.r_ar = np.sqrt(self.x_ar * self.x_ar + self.y_ar * self.y_ar)
        return self.r_ar
    
    # calculation of orbital angles from the initial r_0 point
    def calculate_angles(self, verbose=False):
        # np.arctan2(y, x)      <- note the order ( y , x )
        th_0 = np.arctan2(self.y_0, self.x_0) % (2 * np.pi)
        if verbose: print("initial angle: %f" % np.rad2deg(th_0))

        # basic 0 to 2*pi angle calculation
        # self.th_ar = (np.arctan2(self.y_ar, self.x_ar) - th_0) % (2 * np.pi)

        # adjust for the number of turns (positive or negative)
        winding_num = 0
        self.th_ar[0] = 0
        for idx in range(1, self.num_points):  # self.num_points
            self.th_ar[idx] = (np.arctan2(self.y_ar[idx], self.x_ar[idx]) - th_0) % (2 * np.pi) + winding_num * 2 * np.pi
            
            if abs(self.th_ar[idx] - self.th_ar[idx - 1]) > np.pi: # new turn
                winding_diff =  -1 * round((self.th_ar[idx] - self.th_ar[idx - 1]) / (2 * np.pi)) # +1/-1 
                winding_num += winding_diff
                if verbose: print("new winding num: %i | angle: %f " % (winding_num, self.th_ar[idx - 1]))
                self.th_ar[idx] += winding_diff * 2 * np.pi

        return self.th_ar

    # calculate eccentricity of the orbit using max and min distances (for elliptical orbits)
    # for an elliptical orbit:
    # r_max = a + c with a being the major semi-axis, c the focal distance
    # r_min = a - c
    # eccentricity e = c / a = (r_max - r_min) / (r_max + r_min)
    # idx_first and idx_last indicate the range of indices within which 
    # the orbital maxima/minima (i.e. apoapsis, periapsis) are calculated
    def calculate_elliptical_eccentricity(self, idx_first=0, idx_last=0):
        if idx_last == 0: idx_last = self.num_points
        if not self.orbit_calculated: self.calculate_motion()

        idx_rng = range(idx_first, idx_last)
        idx_min, idx_max = [ np.argmin(self.r_ar[idx_rng]),  np.argmax(self.r_ar[idx_rng]) ]

        r_min, r_max = [self.r_ar[idx_min], self.r_ar[idx_max]]
        ecc = (r_max - r_min) / (r_max + r_min)

        return ecc, idx_min, idx_max

    # calculate the orbital period (primitive)
    # epsilon is the criterion for "sufficiently close"
    def calculate_orbital_period(self, eps=0, skip_points=20, verbose=False):
        if not self.orbit_calculated:
            self.calculate_motion()

        if eps == 0: # set to the half of distanced covered in dt with initial velocity
            eps = (self.v_0[0]**2 + self.v_0[1]**2)**0.5 * self.dt / 2.
        
        # calculate distance to the initial point 
        for idx in range(skip_points, self.num_points):   # skip first few points point
            dist = abs(self.x_ar[idx] - self.x_0) + abs(self.y_ar[idx] - self.y_0)
            if verbose and dist < 10 * eps: 
                print("Orbital point: ", dist, self.time_arr[idx], idx)
            if dist < eps:
                self.orbital_period = self.time_arr[idx]
                return (self.time_arr[idx], idx)
        
        print("could not calculate the orbital period.")
        exit()
    
    # the length of the astronomical/sidereal day based on the orbital period and num days/year
    # astronomical day is slightly shorter than the solar/synodic day if the two rotational
    # motions are in the same directions. e.g. for earch.
    def calculate_rotational_period(self, days_per_year=365.256, rotations_aligned=True):
        if self.orbital_period == 0: self.calculate_orbital_period()

        # calculate number of rotations:
        rotations_per_year = \
            days_per_year + 1 if rotations_aligned else days_per_year - 1

        self.rotational_period = self.orbital_period / rotations_per_year

        return self.rotational_period






