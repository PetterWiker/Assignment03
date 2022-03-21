import numpy as np
from Laminate import Laminate


class OptimizedCompositePipe:

    def __init__(self, material: Laminate, internal_radius: float, thickness: float,
                       service_load: dict, cases: list) -> None:
        self.material = material
        self.internal_radius = internal_radius
        self.external_radius = internal_radius + thickness
        self.thickness = thickness
        self.service_load = service_load

        self.cross_sectional_area = np.pi*(self.external_radius**2 - self.internal_radius**2)
        self.optimized_mass_density, self.optimized_thickness = 0, 0
        self.N = [0, 0, 0]
        self.compute_stress_state(cases)
        pass

    def compute_stress_state(self, cases: list[int]) -> None:
        """

        :param: cases: List of the cases that should be included in the calculation.
        :return:
        """
        for case in cases:
            # Include the relevant cases
            getattr(self, "case_{}".format(case))()

        self.material.update_state(Nx=self.N[0], Ny=self.N[1], Nxy=self.N[2], Mx=0, My=0, Mxy=0)

        if not self.material.is_broken:
            self.__init__(self.material, self.internal_radius, self.thickness*0.95, self.service_load, cases)
        else:
            self.optimized_mass_density = self.compute_mass_density()
            self.optimized_thickness = self.thickness
        pass

    def compute_mass_density(self):
        """
        :return: The mass density in the longitudinal direction (mass per unit length of the pipe) [g/mm].
        """
        mass_density = self.cross_sectional_area*self.material.mass
        return mass_density

    def case_1(self):
        # Based on relations from https://folk.ntnu.no/nilspv/TMM4175/cs-thin-walled-pipes.html
        self.N[0] += self.service_load["P"]*self.internal_radius/2
        self.N[1] += self.service_load["P"]*self.internal_radius
        pass

    def case_2(self):
        self.N[0] += self.service_load["F_x"]*self.thickness/self.cross_sectional_area
        pass

    def case_3(self):
        # Based on relation from https://folk.ntnu.no/nilspv/TMM4175/cs-light-weight-drive-shaft.html
        self.N[2] += self.service_load["T"]/(2*np.pi*self.internal_radius**2)
        pass

