import laminatelib
import copy
import numpy as np
from numpy import ndarray


class Ply:

    def __init__(self, material: dict, orientation: int, thickness: float) -> None:
        """
        Simple class for storing ply information in an object.
        :param material: The material comprising the ply. On the format used in "matlib.py".
        :param orientation: The ply orientation in degrees relative to the xyz-coordinate system.
        :param thickness: The thickness of the ply.
        """
        self.material = copy.deepcopy(material)  # Deepcopy to remove pointer to original material, allowing for change
        self.orientation = orientation
        self.thickness = thickness
        self.is_FF_damaged = False
        self.is_IFF_damaged = False
        self.stress_state, self.strain_state = [], []
        self.Te = laminatelib.T2De(self.orientation)
        self.Q = laminatelib.Q2D(self.material)
        self.fE_FF = 0
        self.fE_IFF = 0
        pass

    def degrade_stiffness(self, degradation_factor: float = 0.9) -> None:
        self.material["E2"] *= degradation_factor
        self.material["G12"] *= degradation_factor
        self.Q = laminatelib.Q2D(self.material)
        pass

    def degrade_compressive_strength(self, reduction_factor: float = 0.1) -> None:
        self.material["XC"] *= reduction_factor
        self.Q = laminatelib.Q2D(self.material)
        pass

    def compute_exposure_factors(self) -> None:
        self.fE_FF = self.fE_MS2D()
        self.fE_IFF = self.fE_Hashin2D()
        pass

    def update_state(self, stresses: ndarray, strains: ndarray) -> None:
        self.stress_state = stresses
        self.strain_state = strains
        pass

    def fE_MS2D(self, load_resistance_factor: float = 1) -> float:
        """
        Calculates the FF exposure factor using the Maximum Stress criterion.
        :param load_resistance_factor:
        :return: True if the load case trips the criterion.
        """
        s_1 = load_resistance_factor*self.stress_state[0]
        fE_MS = max([s_1/self.material["XT"], -s_1/self.material["XC"]])
        return fE_MS

    def fE_Hashin2D(self, load_resistance_factor: float = 1.2) -> float:
        """
        Calculates the IFF exposure factor using the 2D Hashin failure criterion.
        :param load_resistance_factor: The exposure factor for IFF shall be less than 1.0 when the applied load is the
                                       service load multiplied with the load resistance factor.
        :return: True if the load case trips the criterion.
        """
        stress_state = [load_resistance_factor*s for s in self.stress_state]

        if stress_state[1] >= 0:
            temp = (1/self.material["YT"]**2)*stress_state[1]**2 + \
                   (1/self.material["S12"]**2)*(stress_state[2]**2)
            fE_IFF = np.sqrt(temp)
        else:
            temp_b = (1/self.material["YC"])*((self.material["YC"]/(2*self.material["S23"]))**2 - 1)*(stress_state[1])
            temp_a = (1/(4*self.material["S23"]**2))*(stress_state[1])**2 + \
                     (1/self.material["S12"]**2)*(stress_state[2]**2)
            fE_IFF = 0 if temp_a == 0 else (2*temp_a)/(-temp_b + (temp_b**2 + 4*temp_a)**0.5)

        return fE_IFF


class Laminate:

    def __init__(self, layup: list[Ply], name: str) -> None:
        """
        Class for modeling the behaviour of a laminate comprised of a stack of plies (using the "Ply"-class).
        Active boundary conditions: The laminate cannot curve, moments are disregarded.
        :param layup: A list of "Ply"-objects.
        :param name: A chosen name for the laminate.
        """
        self.layup = layup
        self.name = name
        self.thickness = self.compute_thickness()
        self.mass = self.compute_mass()
        self.A = self.compute_A()
        self.load, self.defor = [], []
        pass

    def compute_thickness(self) -> float:
        return sum([ply.thickness for ply in self.layup])

    def update_thickness(self, adjustment_factor: float) -> None:
        for ply in self.layup:
            ply.thickness *= adjustment_factor
        self.thickness = self.compute_thickness()
        pass

    def compute_mass(self) -> float:
        """
        :return: The mass density through thickness (mass per unit width and height) [g/mm^2].
        """
        return sum([ply.material["rho"]*ply.thickness for ply in self.layup])

    def compute_A(self) -> ndarray:
        A = np.zeros((3, 3), float)
        h_bot = -self.thickness/2
        for ply in self.layup:
            Q = laminatelib.Q2D(ply.material)
            Qt = laminatelib.Q2Dtransform(Q, ply.orientation)
            h_top = h_bot + ply.thickness
            A += Qt*(h_top-h_bot)
            h_bot = h_top
        return A

    def update_state(self, load: list[float]) -> None:
        self.load = load
        self.A = self.compute_A()
        self.defor = self.solve_A_load_case(self.load)
        self.update_ply_states()
        pass

    def update_ply_states(self) -> None:
        # As the laminate cannot have curvature, the strain remains the same throughout the ply thickness
        for ply in self.layup:
            strn123 = np.dot(ply.Te, self.defor)
            strs123 = np.dot(ply.Q, strn123)
            ply.update_state(stresses=strs123, strains=strn123)
        pass

    def solve_A_load_case(self, load: list[float]) -> ndarray:
        return np.dot(np.linalg.inv(self.A), load)
