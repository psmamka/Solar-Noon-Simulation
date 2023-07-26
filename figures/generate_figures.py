import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from util.read_data import read_noaa_data
from simul.field_2d import Field2D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def get_noaa_noon_offsets():
    cur_path = os.path.dirname(__file__)
    fname = os.path.join(cur_path, '..', 'data', 'noaa-phoenix-solar-noon.txt')

    noon_data_noaa = read_noaa_data(fname)
    # print(noon_data_noaa)
    # print(noon_data_noaa.shape)

    data_tpl_noaa = list(map(lambda elm: tuple([int(x) for x in elm.split(':')]), noon_data_noaa))
    # print(data_tpl_noaa)

    datetime_to_offset = lambda tpl: (tpl[0] - 12)*3600 + tpl[1]*60 + tpl[2]

    noaa_noon_offsets = list(map(datetime_to_offset, data_tpl_noaa))

    return noaa_noon_offsets

def get_simul_noon_offsets(days_offset=0, hours_offset=-7, a_phi=-18.1):
    acc_field = lambda x, y: [-x / (x*x + y*y)**1.5, -y / (x*x + y*y)**1.5]

    elliptic_tilt = Field2D(
        r_0 = [1, 0],
        v_0 = [0, 1.00835],      # <=== increase the initial velocity (earth orbit ecc. 1.0167)
        acc_field = acc_field,
        dt = 1e-5,                  # 1e-5
        t_max = 4.1 * np.pi,
        days_per_year=365.256, 
        rotations_aligned=True,
        axial_theta = 23.4,             # earth axial tilt (theta_a): 23.4
        axial_phi = a_phi,              # earh axial azimuth (phi_a) at perihelion: -18.1
        psi_0 = 262.35,                 # rotation at time 0 (perihelion) from 16:17 UTC on Jan 4 2023 244.25 + 18.1
        observer_latitude = 33.46,      # Phoenix, AZ Latitude 
        observer_longitude = -112.06,    # Phoenix, AZ Longitude
    )

    elliptic_tilt.calculate_motion()
    elliptic_tilt.calculate_orbital_period()
    elliptic_tilt.calculate_rotational_period()
    elliptic_tilt.calculate_observer_n()
    elliptic_tilt.calculate_observer_east()

    elliptic_tilt.calculate_solar_noons()
    simul_noon_datetimes = list(map(elliptic_tilt.calculate_observer_datetime, 
                                    elliptic_tilt.noon_times_ar,
                                    np.ones(len(elliptic_tilt.noon_times_ar)) * hours_offset,   # hr offset from UTC
                                    np.ones(len(elliptic_tilt.noon_times_ar)) * 16.0,   # hr_0_utc: 16:17 for perihelion of 2023
                                    np.ones(len(elliptic_tilt.noon_times_ar)) * 17.0    # mn_0_utc
                                    ))

    datetime_to_offset = lambda tpl: (tpl[1] - 12)*3600 + tpl[2]*60 + tpl[3]    # tpl = (day, hour, min, sec)
    # print(f"datetime offset of 11:59:59 from Noon: {datetime_to_offset((1, 11, 59, 59))}")

    simul_noon_offsets = list(map(datetime_to_offset, simul_noon_datetimes))

    return simul_noon_offsets[days_offset: days_offset + 365]


if __name__ == "__main__":

    noaa_noon_offsets = get_noaa_noon_offsets()
    simul_noon_offsets = get_simul_noon_offsets(days_offset=0, hours_offset=-7.0, a_phi=-18.1)    # offset: 9.8, a_phi: -18.1

    fig, ax = plt.subplots()
    simul_curve = ax.plot(simul_noon_offsets[0:362], label='simul')     # [0:361]   [0:357]
    noaa_curve = ax.plot(noaa_noon_offsets[3:365], label='noaa')        # [4:365]   [8:365]
    ax.set_title("Earth Solar noon offsets from 12pm in seconds\nLocation: Phoenix AZ, Lat: 33.64 Lon: -112.06")
    ax.set_xlabel("days since perihelion (Jan 4th 2023)")
    ax.set_ylabel("noon offset (sec)")
    ax.legend()
    plt.show()