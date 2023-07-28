import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from util.read_data import read_noaa_data, read_navy_data
from simul.field_2d import Field2D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def get_data_noon_offsets(data_fname, file_reader):
    cur_path = os.path.dirname(__file__)
    fname = os.path.join(cur_path, '..', 'data', data_fname)

    noon_data = file_reader(fname) # read_noaa_data(fname)
    data_tpl = list(map(lambda elm: tuple([int(x) for x in elm.split(':')]), noon_data))

    datetime_to_offset = lambda tpl: (tpl[0] - 12)*3600 + tpl[1]*60 + tpl[2]
    noon_offsets = list(map(datetime_to_offset, data_tpl))
    return noon_offsets

# def get_navy_noon_offsets():
#     cur_path = os.path.dirname(__file__)
#     fname = os.path.join(cur_path, '..', 'data', 'navy-phoenix-sunrise-sunset.txt')

#     noon_data_navy = read_navy_data(fname)
#     data_tpl_navy = list(map(lambda elm: tuple([int(x) for x in elm.split(':')]), noon_data_navy))
#     datetime_to_offset = lambda tpl: (tpl[0] - 12)*3600 + tpl[1]*60 + tpl[2]
#     noon_offsets = list(map(datetime_to_offset, data_tpl_navy))
#     return noon_offsets

def get_simul_noon_offsets(days_offset=0, hours_offset=-7, ecc=0.0167, a_theta=23.4, a_phi=-18.1):
    acc_field = lambda x, y: [-x / (x*x + y*y)**1.5, -y / (x*x + y*y)**1.5]

    elliptic_tilt = Field2D(
        r_0 = [1, 0],
        v_0 = [0, 1.0 + ecc/2.0],      # <=== increase the initial velocity (earth orbit ecc. 0.0167)
        acc_field = acc_field,
        dt = 1e-5,                  # 1e-5
        t_max = 2.1 * np.pi,
        days_per_year=365.256, 
        rotations_aligned=True,
        axial_theta = a_theta,             # earth axial tilt (theta_a): 23.4
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

    fig_size = (6, 6)   # fig size in inches

    noaa_noon_offsets = get_data_noon_offsets(data_fname='noaa-phoenix-solar-noon.txt', file_reader=read_noaa_data)
    navy_noon_offsets = get_data_noon_offsets(data_fname='navy-phoenix-sunrise-sunset.txt', file_reader=read_navy_data)
    # print(navy_noon_offsets)

    simul_noon_offsets = get_simul_noon_offsets(days_offset=0, hours_offset=-7.0, ecc=0.0167, a_theta=23.4, a_phi=-18.1)  

    fig, ax = plt.subplots()
    simul_curve = ax.plot(simul_noon_offsets[0:362], 'green', label='simul')
    noaa_curve = ax.plot(noaa_noon_offsets[3:365], 'orange', label='noaa')
    navy_curve = ax.scatter(np.arange(3, 365), navy_noon_offsets[3:365], s=1, c='blue', label='navy')
    ax.set_title("Earth Solar noon offsets from 12pm in seconds\nLocation: Phoenix AZ (UTC-7), Lat: 33.64° Lon: -112.06°")
    ax.set_xlabel("days since perihelion (Jan 4th 2023)")
    ax.set_ylabel("noon offset (sec)")
    ax.legend()


    # no eccentricity
    simul_no_ecc = get_simul_noon_offsets(days_offset=0, hours_offset=-7.0, ecc=0.0, a_theta=23.4, a_phi=-18.1)    

    fig2, ax2 = plt.subplots()
    fig2.set_size_inches(fig_size)
    simul_no_ecc = ax2.plot(simul_no_ecc[0:365], label='simul')
    ax2.set_title("Circular Orbit (ecc = 0.0)\nAxial Tilt 23.4°\nLocation: Phoenix AZ, Lat: 33.64° Lon: -112.06°")
    ax2.set_xlabel("days since perihelion (Jan 4th 2023)")
    ax2.set_ylabel("noon offset (sec)")


    # no axial tilt
    simul_no_tilt = get_simul_noon_offsets(days_offset=0, hours_offset=-7.0, ecc=0.0167, a_theta=0.0, a_phi=-18.1)    

    fig3, ax3 = plt.subplots()
    fig3.set_size_inches(fig_size)
    simul_no_tilt = ax3.plot(simul_no_tilt[0:365], label='simul')
    ax3.set_title("Elliptic Orbit (ecc = 0.0167)\nNo Axial Tilt (0.0°)\nLocation: Phoenix AZ, Lat: 33.64° Lon: -112.06°")
    ax3.set_xlabel("days since perihelion (Jan 4th 2023)")
    ax3.set_ylabel("noon offset (sec)")

    plt.show()

