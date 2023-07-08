import numpy as np
from field_2d import Field2D
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Simple circular orbit
def inverse_square_test():
    acc_field = lambda x, y: [-x / (x*x + y*y)**1.5, -y / (x*x + y*y)**1.5]
    print(acc_field(1, 0),      #   [-1.0, 0.0]
          acc_field(0, 1),      #   [0.0, -1.0]
          acc_field(1, 1))      #   [-0.3535, -0.3535]

    motion_2d = Field2D(
        r_0 = [1, 0],
        v_0 = [0, 1],
        acc_field = acc_field,
        dt = 2e-4,
        t_max = 2 * np.pi
    )

    motion_2d.calculate_motion()
    print("final x, y, vx, vy:")
    print(motion_2d.x_ar[-1], motion_2d.y_ar[-1], motion_2d.vx_ar[-1], motion_2d.vy_ar[-1])

    # points = list(zip(motion_2d.x_ar, motion_2d.y_ar))

    fig, ax = plt.subplots()
    # ax.scatter(motion_2d.x_ar, motion_2d.y_ar, s=np.ones(motion_2d.num_points))
    ax.plot(motion_2d.x_ar, motion_2d.y_ar)
    ax.scatter(x=[0], y=[0], s=1, c=['black'])
    fig.set_size_inches((5, 5))
    ax.set_title(f"Circular Orbit around A\n Central 1/r^2 Force")
    plt.show()



def run_all_tests():
    inverse_square_test()

if __name__ == "__main__":
    run_all_tests()
