import numpy as np
from field_2d import Field2D
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Simple circular orbit
def circular_orbit_test():
    acc_field = lambda x, y: [-x / (x*x + y*y)**1.5, -y / (x*x + y*y)**1.5]
    print(acc_field(1, 0),      #   [-1.0, 0.0]
          acc_field(0, 1),      #   [0.0, -1.0]
          acc_field(1, 1))      #   [-0.3535, -0.3535]

    circular_2d = Field2D(
        r_0 = [1, 0],
        v_0 = [0, 1],
        acc_field = acc_field,
        dt = 1e-5,
        t_max = 4.1 * np.pi,
        days_per_year=365.256, 
        rotations_aligned=True
    )

    circular_2d.calculate_motion()
    print("final x, y, vx, vy:")
    print(circular_2d.x_ar[-1], circular_2d.y_ar[-1], circular_2d.vx_ar[-1], circular_2d.vy_ar[-1])

    # points = list(zip(circular_2d.x_ar, circular_2d.y_ar))

    fig, ax = plt.subplots()
    # ax.scatter(circular_2d.x_ar, circular_2d.y_ar, s=np.ones(circular_2d.num_points))
    ax.plot(circular_2d.x_ar, circular_2d.y_ar)
    ax.scatter(x=[0], y=[0], s=1, c=['black'])
    fig.set_size_inches((5, 5))
    ax.set_title(f"Circular Orbit around A\n Central 1/r^2 Force")
    # plt.show()

    (orb_time, orb_idx) = circular_2d.calculate_orbital_period_distbased(eps=1e-4, skip_points=1000, verbose=True)
    print("Distance Based Orbital Period: %f, | Orbital Index: %i" % (orb_time, orb_idx))    
    # Distance Based Orbital Period: 6.283450, | Orbital Index: 628345

    (orb_time, orb_idx) = circular_2d.calculate_orbital_period()
    print("Angle Based Orbital Period: %f, | Orbital Index: %i" % (orb_time, orb_idx))
    # Angle Based Orbital Period: 6.283481, | Orbital Index: 628349

    ecc, idx_min, idx_max = circular_2d.calculate_elliptical_eccentricity(idx_first=0, idx_last=orb_idx)

    print(f"eccentricity: %f | min index: %i | max index: %i" % (ecc, idx_min, idx_max))
    # eccentricity: 0.000031 | min index: 6 | max index: 649999

    circular_2d.calculate_angles(verbose=True)

    circular_2d.calculate_rotational_period()

    print(f"rotational period: {circular_2d.rotational_period}")
    # rotational period: 0.017155896422174657

    circular_2d.calculate_solar_angles()
    # fig, ax = plt.subplots()
    # ax.plot(circular_2d.time_ar, circular_2d.solar_angle_ar)
    # plt.show()

    time_1, time_2 = (-0.1, 0.1)
    anlge_1, angle_2 = (np.pi - 0.001 , np.pi + 0.001)

    # basic noon crossing test:
    noon_time = circular_2d.get_angle_crossing(time_1, time_2, anlge_1, angle_2, noon_angle=np.pi)
    print(f"noon crossing for {time_1, time_2, anlge_1, angle_2} should be at 0: {noon_time}")

    circular_2d.calculate_solar_noons(noon_angle = np.pi)
    print(f"first few solar noon orbital times: {circular_2d.noon_times_ar[0:20]}")
    print(f"number of noon crossings: {len(circular_2d.noon_times_ar)}")

    # fig, ax = plt.subplots()
    # ax.plot(circular_2d.noon_times_ar)
    # plt.show()   

    circular_2d.setup_planetary_clock(hours_in_day=24.0)
    dy, hr, mn, sc = circular_2d.calculate_planetary_datetime(circular_2d.orbital_period / 2.0)
    print(f"planetary datetime for half-orbital period: day: {dy}, hour: {hr}, minute: {mn}, second: {sc}")
    # planetary datetime for half-orbital period: day: 182.0, hour: 15.0, minute: 4.0, second: 19.199999999668258 
    
    noon_datetimes = list(map(circular_2d.calculate_planetary_datetime, circular_2d.noon_times_ar))
    noon_str = "\n".join(str(t) for t in noon_datetimes)
    print(f"noon times in planetary date-time: {noon_str}")

    plt.show()


# Elliptical orbit
def elliptical_orbit_test():
    acc_field = lambda x, y: [-x / (x*x + y*y)**1.5, -y / (x*x + y*y)**1.5]

    elliptical_2d = Field2D(
        r_0 = [1, 0],
        v_0 = [0, 1.00835],      # <=== increase the initial velocity (earth orbit ecc. 1.0167)
        acc_field = acc_field,
        dt = 1e-5,
        t_max = 3.0 * np.pi,
        days_per_year=365.256, 
        rotations_aligned=True
    )

    elliptical_2d.calculate_motion()

    (orb_time, orb_idx) = elliptical_2d.calculate_orbital_period()  
    print("Angle Based Orbital Period: %f, | Orbital Index: %i" % (orb_time, orb_idx))
    # Angle-based Orbital Period: 6.444912, | Orbital Index: 644492
    # Distance based Orbital Period: 6.444880, | Orbital Index: 644488

    ecc, idx_min, idx_max = elliptical_2d.calculate_elliptical_eccentricity(idx_first=0, idx_last=orb_idx)

    print(f"eccentricity: %f | min index: %i | max index: %i" % (ecc, idx_min, idx_max))
    # eccentricity: 0.016785 | min index: 0 | max index: 322358

    fig, ax = plt.subplots()
    ax.plot(elliptical_2d.x_ar, elliptical_2d.y_ar)
    ax.scatter(x=[0], y=[0], s=1, c=['black'])  # gravity center 
    fig.set_size_inches((5, 5))
    ax.set_title(f"Elliptical Orbit around A\n Central 1/r^2 Force")
    # plt.show()

    elliptical_2d.calculate_rotational_period()
    print(f"rotational period: {elliptical_2d.rotational_period}")

    elliptical_2d.calculate_solar_angles()
    
    elliptical_2d.calculate_solar_noons(noon_angle = np.pi)

    elliptical_2d.setup_planetary_clock(hours_in_day=24.0)

    noon_datetimes = list(map(elliptical_2d.calculate_planetary_datetime, elliptical_2d.noon_times_ar))
    # print(noon_datetimes.)
    noon_str = "\n".join(str(t) for t in noon_datetimes)
    print(f"noon times in planetary date-time: {noon_str}")

    datetime_to_offset = lambda tpl: (tpl[1] - 12)*3600 + tpl[2]*60 + tpl[3]    # tpl = (day, hour, min, sec)
    print(f"datetime offset of 11:59:59 from Noon: {datetime_to_offset((1, 11, 59, 59))}")

    noon_offsets = list(map(datetime_to_offset, noon_datetimes))

    fig, ax = plt.subplots()
    ax.plot(noon_offsets[0:365])
    ax.set_title("Solar noon offsets from 12pm in seconds")
    ax.set_xlabel("days")
    ax.set_ylabel("offset (sec)")
    plt.show()

def circular_axial_tilt_test():
    acc_field = lambda x, y: [-x / (x*x + y*y)**1.5, -y / (x*x + y*y)**1.5]

    circular_tilt = Field2D(
        r_0 = [1, 0],
        v_0 = [0, 1],
        acc_field = acc_field,
        dt = 1e-5,
        t_max = 4.1 * np.pi,
        days_per_year=365.256, 
        rotations_aligned=True,
        axial_theta = 23.4,             # earth axial tilt (theta_a): 23.4
        axial_phi = -108.1,              # earh axial azimuth (phi_a) at perihelion: -108.1
        observer_latitude = 33.50,      # Phoenix, AZ Latitude 
        observer_longitude = -112.0,    # Phoenix, AZ Longitude
        # observer_longitude = 0.0,         # North Pole Longitude
        # observer_latitude = 90.0,         # North Pole Latitude 
    )

    circular_tilt.calculate_motion()
    circular_tilt.calculate_orbital_period()
    circular_tilt.calculate_rotational_period()
    circular_tilt.calculate_observer_n()
    circular_tilt.calculate_observer_east()

    circular_tilt.calculate_solar_noons()
    noon_datetimes = list(map(circular_tilt.calculate_observer_datetime, 
                              circular_tilt.noon_times_ar,
                              np.ones(len(circular_tilt.noon_times_ar)) * 10.0   # 10 hr offset
                              ))

    # print(circular_tilt.obs_n_ar[0:100_000:1000])
    # print(circular_tilt.obs_east_ar[0:100_000:1000])
    #   sanity check:
    # print("should be near zero: ", np.dot(circular_tilt.obs_n_ar[1000], circular_tilt.obs_east_ar[1000]))
    # should be near zero:  4.672239546380282e-17

    # print(circular_tilt.noon_times_ar[0:365:30])
    # print(noon_datetimes[0:365:30])

    datetime_to_offset = lambda tpl: (tpl[1] - 12)*3600 + tpl[2]*60 + tpl[3]    # tpl = (day, hour, min, sec)
    # print(f"datetime offset of 11:59:59 from Noon: {datetime_to_offset((1, 11, 59, 59))}")

    noon_offsets = list(map(datetime_to_offset, noon_datetimes))

    fig, ax = plt.subplots()
    ax.plot(noon_offsets[0:365])
    ax.set_title("Solar noon offsets from 12pm in seconds")
    ax.set_xlabel("days")
    ax.set_ylabel("offset (sec)")
    plt.show()


def elliptical_axial_tilt_test():
    acc_field = lambda x, y: [-x / (x*x + y*y)**1.5, -y / (x*x + y*y)**1.5]

    elliptic_tilt = Field2D(
        r_0 = [1, 0],
        v_0 = [0, 1.00835],      # <=== increase the initial velocity (earth orbit ecc. 1.0167)
        acc_field = acc_field,
        dt = 1e-5,
        t_max = 4.1 * np.pi,
        days_per_year=365.256, 
        rotations_aligned=True,
        axial_theta = 23.4,             # earth axial tilt (theta_a): 23.4
        axial_phi = -108.1,              # earh axial azimuth (phi_a) at perihelion: -108.1
        observer_latitude = 33.50,      # Phoenix, AZ Latitude 
        observer_longitude = -112.0,    # Phoenix, AZ Longitude
        # observer_longitude = 0.0,         # North Pole Longitude
        # observer_latitude = 90.0,         # North Pole Latitude 
    )

    elliptic_tilt.calculate_motion()
    elliptic_tilt.calculate_orbital_period()
    elliptic_tilt.calculate_rotational_period()
    elliptic_tilt.calculate_observer_n()
    elliptic_tilt.calculate_observer_east()

    elliptic_tilt.calculate_solar_noons()
    noon_datetimes = list(map(elliptic_tilt.calculate_observer_datetime, 
                              elliptic_tilt.noon_times_ar,
                              np.ones(len(elliptic_tilt.noon_times_ar)) * 10.0   # 10 hr offset
                              ))

    datetime_to_offset = lambda tpl: (tpl[1] - 12)*3600 + tpl[2]*60 + tpl[3]    # tpl = (day, hour, min, sec)
    # print(f"datetime offset of 11:59:59 from Noon: {datetime_to_offset((1, 11, 59, 59))}")

    noon_offsets = list(map(datetime_to_offset, noon_datetimes))

    fig, ax = plt.subplots()
    ax.plot(noon_offsets[0:365])
    ax.set_title("Solar noon offsets from 12pm in seconds\nEccentricity: 1.67%, Axial Tilt: 23.4°")
    ax.set_xlabel("days from perihelion")
    ax.set_ylabel("offset (sec)")
    plt.show()


def run_all_tests():
    # circular_orbit_test()
    # elliptical_orbit_test()
    # circular_axial_tilt_test()
    elliptical_axial_tilt_test()

if __name__ == "__main__":
    run_all_tests()
