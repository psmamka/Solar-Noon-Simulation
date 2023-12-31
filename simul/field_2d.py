# Calculation of the trajectory in a given acceleration field using increments

import numpy as np

class Field2D:
    def __init__(self, r_0=[0, 0], v_0=[0, 0], acc_field=lambda x,y: [0, 0], 
                 dt=1, t_0=0, t_max=100, days_per_year=365.256, rotations_aligned=True,
                 axial_theta=0, axial_phi=0, psi_0=0,
                 observer_latitude=0, observer_longitude=0):
        self.r_0 = r_0
        self.x_0, self.y_0 = self.r_0

        self.v_0 = v_0
        self.v_x_0, self.v_y_0 = self.v_0

        self.acc_field = acc_field

        self.dt = dt
        self.t_max = t_max

        num_points = int(np.ceil((t_max - t_0) / dt) + 1)    # first and last points
        self.idx_ar = np.arange(start = 0, stop = num_points, step = 1)
        self.time_ar = self.idx_ar * dt

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

        # planetary rotation
        self.days_per_year = days_per_year
        self.rotations_aligned = rotations_aligned   # planetary and orbital rotations in the same direction
        self.rotational_period = 0.
        self.rotational_velocity = 0.
        self.psi_ar = np.zeros(num_points)
        self.solar_angle_ar = np.zeros(num_points)

        # planetary clock in units of simulation time
        self.avg_solar_day_dur = 0.
        self.hour_dur = 0.
        self.minute_dur = 0.
        self.second_dur = 0.

        # set up the axial tilt: angles of the rotation axis w.r.t ecliptic
        # axial theta and axial phi are in spherical coordinates w.r.t
        # the x and y axes on the obital plane (the "ecliptic") and the 
        # z-axis perpendicular to the ecliptic
        self.ax_theta = np.deg2rad(axial_theta)
        self.ax_phi  = np.deg2rad(axial_phi)

        # Set up the observer lattitude and longitude
        # reference point w.r.t. the x-y axis
        self.obs_lat = np.deg2rad(observer_latitude)
        self.obs_lon = np.deg2rad(observer_longitude)
        # initial rotation phase at time 0
        self.psi_0 = np.deg2rad(psi_0)

        # observer noraml vector (n_hat) w.r.t. the planet surface
        # and east (e_hat) vector w.r.t. the north pole
        # the two vectors will be orthogonal
        self.obs_n_ar = np.zeros((num_points, 3))
        self.obs_east_ar = np.zeros((num_points, 3))

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
        self.calculate_angles()
    
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

    # angle-based approach (more robust)
    # calculate the orbital period based on the time it takes for the planet to
    # make a complete revolution (2 * pi angle) around the sun
    def calculate_orbital_period(self, verbose=False):
        if self.th_ar[1] == 0: self.calculate_angles()

        if self.th_ar[1] > 0 :  # counter-clockwise orbit
            for idx in range(1, self.num_points):
                if self.th_ar[idx] >= 2 * np.pi and self.th_ar[idx - 1] < 2 * np.pi:    # calculate crossing of 2pi
                    self.orbital_period = self.get_angle_crossing(
                        self.time_ar[idx - 1], self.time_ar[idx], self.th_ar[idx - 1], self.th_ar[idx], 2 * np.pi
                    )
                    return self.orbital_period, idx
        else:                   # clockwise orbit
            for idx in range(1, self.num_points):
                if self.th_ar[idx] <= -2 * np.pi and self.th_ar[idx - 1] > -2 * np.pi:    # calculate crossing of -2pi
                    self.orbital_period = self.get_angle_crossing(
                        self.time_ar[idx - 1], self.time_ar[idx], self.th_ar[idx - 1], self.th_ar[idx], -2 * np.pi
                    )
                    return self.orbital_period, idx

        print("could not calculate the angle-based orbital period.")
        exit()

    # calculate the orbital period (primitive)
    # epsilon is the criterion for "sufficiently close"
    # primitive approach based on the closest distance to the start point
    # angle based calculation above should be more reliable
    def calculate_orbital_period_distbased(self, eps=0, skip_points=20, verbose=False):
        if not self.orbit_calculated:
            self.calculate_motion()

        if eps == 0: # set to the half of distanced covered in dt with initial velocity
            eps = (self.v_0[0]**2 + self.v_0[1]**2)**0.5 * self.dt / 2.
        
        # calculate distance to the initial point 
        for idx in range(skip_points, self.num_points):   # skip first few points point
            dist = abs(self.x_ar[idx] - self.x_0) + abs(self.y_ar[idx] - self.y_0)
            if verbose and dist < 10 * eps: 
                print("Orbital point: ", dist, self.time_ar[idx], idx)
            if dist < eps:
                self.orbital_period = self.time_ar[idx]
                return (self.time_ar[idx], idx)
        
        print("could not calculate the distance-based orbital period.")
        exit()
    
    # the length of the astronomical/sidereal day based on the orbital period and num days/year
    # astronomical day is slightly shorter than the solar/synodic day if the two rotational
    # motions are in the same directions. e.g. for earch.
    def calculate_rotational_period(self):
        if self.orbital_period == 0: self.calculate_orbital_period()

        # calculate number of rotations:
        rotations_per_year = \
            self.days_per_year + 1 if self.rotations_aligned else self.days_per_year - 1

        self.rotational_period = self.orbital_period / rotations_per_year

        self.rotational_velocity = 2 * np.pi / self.rotational_period

        return self.rotational_period
    
    # When we have the initial ref (utc) time instead of the phi_0
    # Calculate the initial rotation phase of the planet
    # based on the utc (clock) time point at simul time of 0
    # as well as the rotation phase. 
    def update_psi_0(self, utc_hr=0, utc_min=0, utc_sec=0, utc_psi=0):
        if self.rotational_period == 0: self.calculate_rotational_period()
        self.setup_planetary_clock()
        # self.psi_0 = (utc_hr * 3600 + utc_min * 60 + utc_sec) / (24 * 3600) * 2 * np.pi
        self.psi_0 = self.rotational_velocity * (utc_hr * self.hour_dur + utc_min * self.minute_dur + utc_sec * self.second_dur) + utc_psi
        
        return self.psi_0
    
    # observer n-hat is the perpendicular to the planet surface at
    # the location of the observer. 
    # It depends on observer coordinates (latt-long) as well as the 
    # axial tilt of the planet and the planet rotation (psi)
    # for the north pole it is always aligned with the rotation axis
    # The 8 contants in the final rotated vector are only calculated once
    # See figures/n_hat_calculation
    # psi_0 is the rotation angle of the planet at time zero
    def calculate_observer_n(self):
        if self.rotational_period == 0: self.calculate_rotational_period()
        self.psi_ar = self.rotational_velocity * self.time_ar + self.psi_0   # roration phase
        # observer spherical coordinates
        theta = np.pi/2 - self.obs_lat
        phi_ar = self.obs_lon + self.psi_ar # total phase is longitude + rotation phase
        cos_phi_ar = np.cos(phi_ar)
        sin_phi_ar = np.sin(phi_ar)
        # matrix multiplacation results:
        # x component consts:
        c1 = np.cos(self.ax_phi) * np.cos(self.ax_theta) * np.sin(theta)
        c2 = np.cos(self.ax_phi) * np.sin(self.ax_theta) * np.cos(theta)
        c3 = -1 * np.sin(self.ax_phi) * np.sin(theta)
        # y:
        c4 = np.sin(self.ax_phi) * np.cos(self.ax_theta) * np.sin(theta)
        c5 = np.sin(self.ax_phi) * np.sin(self.ax_theta) * np.cos(theta)
        c6 = np.cos(self.ax_phi) * np.sin(theta) 
        # z:    (sanity check: ax_phi independent)
        c7 = -1 * np.sin(self.ax_theta) * np.sin(theta)
        c8 = np.cos(self.ax_theta) * np.cos(theta)

        self.obs_n_ar[:,0] = c1 * cos_phi_ar + c2 + c3 * sin_phi_ar         # x component of n_hat
        self.obs_n_ar[:, 1] = c4 * cos_phi_ar + + c5 + c6 * sin_phi_ar      # y
        self.obs_n_ar[:, 2] = c7 * cos_phi_ar + c8                          # z
        return
    
    # Observer east vector. 
    # It will be used for calculation of solar noons,
    # The technique is to do a vector product of axis and n, then normalize:
    # 
    #       obserber_east = rotation_axis x observer_n / sin(observer_theta)
    # 
    # Vector product refresher:
    #
    #       (x1, y1, z1) x (x2, y2, z2) = ( y1 z2 - z1 y2 , z1 x2 - x1 z2 , x1 y2 - y1 x2 )
    # 
    def calculate_observer_east(self):
        if all(self.obs_n_ar[0, :] == [0, 0, 0]): self.calculate_observer_n()
        # special cases: observer on north or south poles
        if self.obs_lat == np.pi/2 or self.obs_lat == -np.pi/2:
            # TODO: define east on poles in relation to the sun
            return

        axial_vector = [np.sin(self.ax_theta) * np.cos(self.ax_phi),
                        np.sin(self.ax_theta) * np.sin(self.ax_phi),
                        np.cos(self.ax_theta)]
        
        self.obs_east_ar = np.cross(axial_vector, self.obs_n_ar) / np.sin(np.pi/2 - self.obs_lat)
    
    # calculate the times (simul) and indices (approx) of solar noons
    # solar noons and midnights occurs when the observer's east unit vector is 
    # perpendicular to the vector from the sun to the planet r = (x, y, 0)
    # for noons (obs_east dot r) changes sign from negative to positive
    # i.e. the sun moves to west for positive rototation velocity
    # vice versa for negative rototation velocity
    def calculate_solar_noons(self):
        self.noon_times_ar = []
        e_dot_r = np.zeros(self.num_points)
        for idx in self.idx_ar:                             # get e dot r
            e_dot_r[idx] = np.dot(self.obs_east_ar[idx], 
                                  [self.x_ar[idx]/self.r_ar[idx], self.y_ar[idx]/self.r_ar[idx], 0])

        if self.rotations_aligned:
            for idx in range(1, self.num_points):
                if e_dot_r[idx - 1] < 0 and e_dot_r[idx] >= 0:  # find the zero crossing
                    self.noon_times_ar.append(
                            self.get_angle_crossing(
                                self.time_ar[idx - 1], self.time_ar[idx], e_dot_r[idx - 1], e_dot_r[idx], 0
                            )
                    )
        else:
            # TODO: finish the case for negative rotation velocity
            pass

        return self.noon_times_ar

    # without axial tilt (deprecated)
    # angle of the sun in the sky. psi minus theta. For now psi at t0 is 0
    # for now the initial point is at psi - theta = 0, i.e. midnight.
    def calculate_solar_angles_flat(self):
        self.rotational_velocity = 2 * np.pi / self.rotational_period
        if not self.rotations_aligned: self.rotational_velocity *= -1

        # psi - theta
        self.solar_angle_ar = (self.time_ar * self.rotational_velocity - self.th_ar) % (2 * np.pi)
        return self.solar_angle_ar
    
    # without axial tilt (deorecated)
    # calculate the exact time points where the observer sees the sun at
    # the highsest point in the sky
    def calculate_solar_noons_flat(self, noon_angle=np.pi):
        self.noon_times_ar = []
        if self.rotations_aligned:
            for idx in range(1, self.num_points):
                # solar_ang_1, solar_ang_2 = (self.solar_angle_ar[idx - 1], self.solar_angle_ar[idx])
                if (self.solar_angle_ar[idx - 1] < noon_angle and self.solar_angle_ar[idx] > noon_angle):
                    self.noon_times_ar.append(
                        self.get_angle_crossing(
                            self.time_ar[idx - 1], self.time_ar[idx], self.solar_angle_ar[idx - 1], self.solar_angle_ar[idx], noon_angle
                        )
                    )
        
        return self.noon_times_ar

    # basic computation of a linear segment from (x0, y0) to (x1, y1) crossing a y value
    def get_angle_crossing(self, time_1, time_2, angle_1, angle_2, noon_angle=np.pi):
        if not ((angle_1 - noon_angle) * (angle_2 - noon_angle) <= 0):
            print("no crossing in the interval")
            exit()
        slope = (angle_2 - angle_1) / (time_2 - time_1)
        # angle_1 + slope * Δt = noon_angle  ⇒  Δt = (noon_angle - angle_1) / slope
        delta_t = (noon_angle - angle_1) / slope

        return time_1 + delta_t
    
    # solar days are on avg. 24h, hours are 60 minutes, etc. 
    def setup_planetary_clock(self, hours_in_day=24.0):
        if self.orbital_period == 0: self.calculate_orbital_period()

        # calculate durations of average solar day, hour, minute and second in simulation time units
        self.avg_solar_day_dur = self.orbital_period / self.days_per_year
        self.hour_dur = self.orbital_period / (24 * self.days_per_year)
        self.minute_dur = self.orbital_period / (24 * 60 * self.days_per_year)
        self.second_dur = self.orbital_period / (24 * 3600 * self.days_per_year)

        return(self.avg_solar_day_dur, self.hour_dur, self.minute_dur, self.second_dur)

    # convert simulation time to observer clock
    # the observer datetime after applying the time offset
    # in format: (day, hour, minute second)
    # offset_in_hours is the local time offset w.r.t. the planetary clock (e.g. UTC -7 for Phoenix)
    # hr_0_utc etc. are the UTC time at simul time of 0. 
    # e.g. 16:17 for perihelion of 2023
    def calculate_observer_datetime(self, time, offset_in_hours=0, hr_0_utc=0, mn_0_utc=0, sc_0_utc=0):
        if self.avg_solar_day_dur == 0: self.setup_planetary_clock()

        date_day, day_rem = np.divmod(time + 
                                      (hr_0_utc + offset_in_hours) * self.hour_dur + mn_0_utc * self.minute_dur + sc_0_utc * self.second_dur, 
                                      self.avg_solar_day_dur
                                      )
        date_hour, hour_rem = np.divmod(day_rem, self.hour_dur)
        date_minute, minute_rem = np.divmod(hour_rem, self.minute_dur)
        date_second = minute_rem / self.second_dur 

        return  (int(date_day), int(date_hour), int(date_minute), date_second)
        
