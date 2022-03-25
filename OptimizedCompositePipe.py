import numpy as np
import copy
from Laminate import Laminate


class OptimizedCompositePipe:

    def __init__(self,
                 material: Laminate,
                 internal_radius: float,
                 thickness: float,
                 service_load: dict[str],
                 progressive_failure: bool) -> None:
        """
        Recursive class for optimizing a composite pipe with regards to material used, i.e. mass.
        :param material: The "Laminate" object making out the composite pipe.
        :param internal_radius: The internal radius of the pipe.
        :param thickness: The thickness of the pipe.
        :param service_load: The service load case that the pipe should be optimized for.
        """
        self.material = copy.deepcopy(material)
        self.internal_radius = internal_radius
        self.thickness = thickness
        self.thickness_history = [self.thickness]

        self.N = self.compute_N(service_load)
        self.material.update_state(load=self.N)

        self.progressive_failure = progressive_failure

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
            if not ply.is_IFF_damaged and self.progressive_failure:
                fEs.append(ply.fE_IFF)

        # Calculate a "thickness adjustment factor"
        fE_dam_max = max(fEs)
        taf = fE_dam_max

        fE_deform = max(self.material.defor)/0.05
        min_deform_limited_thickness = self.thickness*fE_deform

        if self.thickness*taf < min_deform_limited_thickness:
            self.update_thickness(adjustment_factor=fE_deform)
            failure_mechanism = "max_deform"
            if self.material.is_FF_damaged:
                failure_mechanism += "+FF"
            if self.material.is_IFF_damaged:
                failure_mechanism += "+IFF"
            return min(self.thickness_history), failure_mechanism

        for ply in self.material.layup:
            ply.fE_FF = ply.fE_FF/taf
            ply.fE_IFF = ply.fE_IFF/taf

            if round(ply.fE_FF, 5) >= 1:
                if ply.is_FF_damaged or not self.progressive_failure:
                    self.update_thickness(adjustment_factor=taf)
                    failure_mechanism = "FF"
                    if self.material.is_IFF_damaged:
                        failure_mechanism += "+IFF"
                    return min(self.thickness_history), failure_mechanism
                else:
                    ply.is_FF_damaged = True
                    self.material.is_FF_damaged = True
                    ply.degrade_stiffness()

            elif (round(ply.fE_IFF, 5) >= 1) and self.progressive_failure and not ply.is_IFF_damaged:
                ply.is_IFF_damaged = True
                self.material.is_IFF_damaged = True
                ply.degrade_compressive_strength()

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


