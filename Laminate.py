import laminatelib
import numpy as np
from numpy import ndarray
from Ply import Ply


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
        self.is_IFF_damaged = False
        self.is_FF_damaged = False
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
