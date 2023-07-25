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

def get_simul_noon_offsets(days_offset=0, hours_offset=10.0):
    acc_field = lambda x, y: [-x / (x*x + y*y)**1.5, -y / (x*x + y*y)**1.5]

    elliptic_tilt = Field2D(
        r_0 = [1, 0],
        v_0 = [0, 1.00835],      # <=== increase the initial velocity (earth orbit ecc. 1.0167)
        acc_field = acc_field,
        dt = 5e-5,                  # 1e-5
        t_max = 6.3 * np.pi,
        days_per_year=365.256, 
        rotations_aligned=True,
        axial_theta = 23.4,             # earth axial tilt (theta_a): 23.4
        axial_phi = -108.1,              # earh axial azimuth (phi_a) at perihelion: -108.1
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
                                    np.ones(len(elliptic_tilt.noon_times_ar)) * hours_offset   # 10 hr offset
                                    ))

    datetime_to_offset = lambda tpl: (tpl[1] - 12)*3600 + tpl[2]*60 + tpl[3]    # tpl = (day, hour, min, sec)
    # print(f"datetime offset of 11:59:59 from Noon: {datetime_to_offset((1, 11, 59, 59))}")

    simul_noon_offsets = list(map(datetime_to_offset, simul_noon_datetimes))

    return simul_noon_offsets[days_offset: days_offset + 365]


if __name__ == "__main__":

    noaa_noon_offsets = get_noaa_noon_offsets()
    simul_noon_offsets = get_simul_noon_offsets(days_offset=360, hours_offset=9.8)

    fig, ax = plt.subplots()
    simul_curve = ax.plot(simul_noon_offsets, label='simul')
    noaa_curve = ax.plot(noaa_noon_offsets, label='noaa')
    ax.set_title("Earth Solar noon offsets from 12pm in seconds")
    ax.set_xlabel("days")
    ax.set_ylabel("noon offset (sec)")
    ax.legend()
    plt.show()