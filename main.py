import numpy as np

from Laminate import Laminate, Ply
from OptimizedCompositePipe import OptimizedCompositePipe
import matlib
import matplotlib.pyplot as plt


class PlotPoint:
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y
        pass



def main():
    #service_load = {"P": 8, "F_x": 100E3, "T": 20E6}
    service_load = {"P": 2, "F_x": 0, "T": 0}

    UDGRP = matlib.get("UDGRP")
    xs = np.linspace(10, 89, 100)

    ys = []
    for x in xs:
        layup = [Ply(material=UDGRP, orientation=x, thickness=100),
                 Ply(material=UDGRP, orientation=-x, thickness=100)]

        lam = Laminate(layup=layup, name="my_lam")
        pipe = OptimizedCompositePipe(material=lam, internal_radius=75, thickness=lam.thickness, service_load=service_load)
        print("thickness at theta={}:".format(x), pipe.optimized_thickness)
        ys.append(pipe.optimized_thickness)
    plt.plot(xs, ys)
    plt.show()
    pass


if __name__ == "__main__":
    main()