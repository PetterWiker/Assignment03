import numpy as np
from Laminate import Laminate


class OptimizedCompositePipe:

    def __init__(self, material: Laminate, internal_radius: float, thickness: float, service_load: dict[str]) -> None:
        """
        Recursive class for optimizing a composite pipe with regards to material used, i.e. mass.
        :param material: The "Laminate" object making out the composite pipe.
        :param internal_radius: The internal radius of the pipe.
        :param thickness: The thickness of the pipe.
        :param service_load: The service load case that the pipe should be optimized for.
        """
        self.material = material
        self.internal_radius = internal_radius
        self.thickness = thickness
        self.thickness_history = [self.thickness]

        self.N = self.compute_N(service_load)
        self.material.update_state(load=self.N)

        self.optimized_thickness, self.failure_mechanism = self.find_optimized_thickness()
        self.optimized_mass_density = self.material.compute_mass()
        pass

    def compute_N(self, service_load):
        P, F_x, T = service_load["P"], service_load["F_x"], service_load["T"]
        # By using assumption self.thickness/self.cross_sectional_area = (r2-r1)/((r2+r1)(r2-r1)*pi)= approx 1/(2*pi*r1)
        N = [P*self.internal_radius/2 + F_x/(np.pi*2*self.internal_radius),
             P*self.internal_radius,
             T/(2*np.pi*self.internal_radius**2)]
        return N

    def find_optimized_thickness(self) -> tuple[float, str]:
        fEs = []
        for ply in self.material.layup:
            ply.compute_exposure_factors()
            fEs.append(ply.fE_FF)
            if not ply.is_IFF_damaged:
                fEs.append(ply.fE_IFF)

        # Calculate a "thickness adjustment factor"
        fE_deform = max(self.material.defor)/0.05
        fE_dam_max = max(fEs)
        taf = max([fE_deform, fE_dam_max])

        for ply in self.material.layup:
            ply.fE_FF *= 1/taf
            ply.fE_IFF *= 1/taf
            fE_deform *= 1/taf

            if round(ply.fE_FF, 5) >= 1: #> ply.fE_IFF:
                if ply.is_FF_damaged:
                    self.update_thickness(adjustment_factor=taf)
                    return min(self.thickness_history), "fiber_failure"
                else:  # if ply is not already broken
                    ply.is_FF_damaged = True
                    ply.degrade_stiffness()

            elif (round(ply.fE_IFF, 5) >= 1 > ply.fE_FF) and not ply.is_IFF_damaged:
                ply.is_IFF_damaged = True
                ply.degrade_compressive_strength()

            elif round(fE_deform, 5) >= 0.05:
                self.update_thickness(adjustment_factor=taf)
                return min(self.thickness_history), "max_deformation"

        self.update_thickness(adjustment_factor=taf)
        self.material.update_state(load=self.N)

        try:
            return self.find_optimized_thickness()
        except RecursionError:
            return min(self.thickness_history), "max_recursion_depth"

    def update_thickness(self, adjustment_factor):
        self.thickness *= adjustment_factor
        self.thickness_history.append(self.thickness)
        self.material.update_thickness(adjustment_factor)
        pass


