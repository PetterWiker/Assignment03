import numpy as np

from Laminate import Laminate
from Ply import Ply
from OptimizedCompositePipe import OptimizedCompositePipe
import matlib
import matplotlib.pyplot as plt
from matplotlib import cm


class PlotPoint2D:
    def __init__(self, x: float, y: float, failure_mechanism: str) -> None:
        self.x = x
        self.y = y
        self.pos = (self.x, self.y)
        self.failure_mechanism = failure_mechanism
        pass


class PlotPoint3D(PlotPoint2D):
    def __init__(self, x: float, y: float, z: float, failure_mechanism: str) -> None:
        super().__init__(x, y, failure_mechanism)
        self.z = z
        self.pos = (self.x, self.y, self.z)
        pass


def type_A_optimized_pipe():
    # Choose a material
    UDGRP = matlib.get("UDGRP")

    # Define the cases
    service_loads = {"case_1": {"P": 8, "F_x": 0, "T": 0},
                     "case_2": {"P": 0, "F_x": 100E3, "T": 0},
                     "case_3": {"P": 0, "F_x": 0, "T": 20E6},
                     "case_4": {"P": 8, "F_x": 100E3, "T": 0},
                     "case_5": {"P": 0, "F_x": 100E3, "T": 20E6}}

    # Initialize a dictionary for storing the thickness data.
    optimized_thicknesses = {}
    for key in service_loads.keys():
        optimized_thicknesses[key] = []

    n_theta = 300
    for theta_1 in np.linspace(10, 88, n_theta):
        # Make a type A layup with an arbitrarily large starting thickness
        layup = [Ply(material=UDGRP, orientation=theta_1, thickness=100),
                 Ply(material=UDGRP, orientation=-theta_1, thickness=100)]

        # Initialize a laminate with this layup
        laminate = Laminate(layup=layup, name="type_A_lam_theta={}".format(theta_1))

        # Runt through the cases, calculating the optimized thickness of the pipe for each case.
        for key, value in service_loads.items():
            # Initialize an optimized composite pipe
            pipe = OptimizedCompositePipe(material=laminate,
                                          internal_radius=75,
                                          thickness=laminate.thickness,
                                          service_load=value,
                                          progressive_failure=True)
            # Initialize an object storing information relating the angle to the optimal thickness
            plot_point = PlotPoint2D(x=theta_1, y=pipe.optimized_thickness, failure_mechanism=pipe.failure_mechanism)
            optimized_thicknesses[key].append(plot_point)
    for case in [1, 2, 3, 4, 5]:
        plot_optimized_case_mpl_2D(optimized_thicknesses, case="case_{}".format(case))
    pass


def type_B_optimized_pipe():
    # Choose a material
    UDGRP = matlib.get("UDGRP")

    # Define the cases
    service_loads = {"case_1": {"P": 8, "F_x": 0, "T": 0},
                     "case_2": {"P": 0, "F_x": 100E3, "T": 0},
                     "case_3": {"P": 0, "F_x": 0, "T": 20E6},
                     "case_4": {"P": 8, "F_x": 100E3, "T": 0},
                     "case_5": {"P": 0, "F_x": 100E3, "T": 20E6}}

    # Initalize a dictionary for storing the thickness data.
    optimized_thicknesses_3D = {}
    for key in service_loads.keys():
        optimized_thicknesses_3D[key] = []

    n_theta = 200

    for theta_1 in np.linspace(10, 88, n_theta):
        for theta_2 in np.linspace(10, 88, n_theta):
            print("Computing for theta1_theta_2: {}_{}".format(round(theta_1, 1), round(theta_2, 1)))
            # Make a type B layup with an arbitrarily large starting thickness
            layup = [Ply(material=UDGRP, orientation=theta_1, thickness=100),
                     Ply(material=UDGRP, orientation=-theta_1, thickness=100),
                     Ply(material=UDGRP, orientation=theta_2, thickness=100),
                     Ply(material=UDGRP, orientation=-theta_2, thickness=100)]

            # Initialize a laminate with this layup
            laminate = Laminate(layup=layup, name="type_B_lam_thetas=({}_{})".format(theta_1, theta_2))

            # Runt through the cases, calculating the optimized thickness of the pipe for each case.
            for key, value in service_loads.items():
                # Initialize an optimized composite pipe
                pipe = OptimizedCompositePipe(material=laminate,
                                              internal_radius=75,
                                              thickness=laminate.thickness,
                                              service_load=value,
                                              progressive_failure=True)
                # Initialize an object storing information relating the angle to the optimal thickness
                plot_point = PlotPoint3D(x=theta_1, y=theta_2, z=pipe.optimized_thickness,
                                         failure_mechanism=pipe.failure_mechanism)
                optimized_thicknesses_3D[key].append(plot_point)
    for case in [1, 2, 3, 4, 5]:
        plot_optimized_case_mpl_3D(optimized_thicknesses_3D, case="case_{}".format(case), n=n_theta, plt_heat_map=True,
                                   plt_fracture_modes=False)
        plot_optimized_case_mpl_3D(optimized_thicknesses_3D, case="case_{}".format(case), n=n_theta, plt_heat_map=False,
                                   plt_fracture_modes=True)
    pass


def plot_optimized_case_mpl_2D(optimized_thicknesses: dict[str], case: str) -> None:
    points = optimized_thicknesses[case]
    possible_failure_mechanisms = ["FF", "FF+IFF", "max_deform", "max_deform+FF", "max_deform+FF+IFF", "max_deform+IFF"]

    data = {}
    for mechanism in possible_failure_mechanisms:
        data[mechanism] = []
    for point in points:
        data[point.failure_mechanism].append(point.pos)

    # Make plot objects for the total case and the individual failure mechanisms
    legend = ["All_thicknesses"]
    plt.plot(*zip(*[point.pos for point in points]), "k-", alpha=0.8)
    for mechanism in possible_failure_mechanisms:
        if len(data[mechanism]):
            legend.append(mechanism)
            plt.plot(*zip(*data[mechanism]), "o", alpha=0.5)
    plt.xlabel("Theta 1 (°)")
    plt.ylabel("Thickness (mm)")
    plt.legend(legend)
    plt.show()
    pass


def plot_optimized_case_mpl_3D(optimized_thicknesses: dict[str], case: str, n: int, plt_heat_map: bool = True,
                               plt_fracture_modes: bool = False) -> None:
    points = optimized_thicknesses[case]
    possible_failure_mechanisms = ["FF", "FF+IFF", "max_deform", "max_deform+FF", "max_deform+FF+IFF", "max_deform+IFF"]

    data = {}
    for mechanism in possible_failure_mechanisms:
        data[mechanism] = []
    for point in points:
        data[point.failure_mechanism].append(point.pos)

    fig = plt.figure(figsize=plt.figaspect(1))
    ax = fig.add_subplot(projection='3d')

    if plt_heat_map:
        x, y, z = zip(*[point.pos for point in points])
        z = list(z)
        z = np.array([z[(i*n):(i*n+n)] for i in range(n)])
        x, y = np.meshgrid([10+i*78/(n-1) for i in range(0, z.shape[0])], [10+i*78/(n-1) for i in range(0, z.shape[1])])
        surf = ax.plot_surface(x, y, z, cmap=cm.RdYlGn_r, alpha=1)
        ax.set_xlabel("Theta 1 (°)")
        ax.set_ylabel("Theta 2 (°)")
        ax.set_zticklabels([])

        cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
        cbar.set_label("Thickness")

    if plt_fracture_modes:
        legend = []
        for mechanism in possible_failure_mechanisms:
            if len(data[mechanism]):
                legend.append(mechanism)
                ax.scatter(*zip(*data[mechanism]), "o", s=1, alpha=1)
        ax.set_xlabel("Theta 1 (°)")
        ax.set_ylabel("Theta 2 (°)")
        ax.set_zticklabels([])
        ax.legend(legend)

    plt.show()
    pass

def main():
    type_A_optimized_pipe()
    type_B_optimized_pipe()
    pass


if __name__ == "__main__":
    main()