import laminatelib
import copy
import numpy as np


class Ply:

    def __init__(self, material: dict, orientation: int, thickness: float) -> None:
        """
        Simple class for storing ply information in an object.
        :param material: The material comprising the ply. On the format used in "matlib.py".
        :param orientation: The ply orientation in degrees relative to the xyz-coordinate system.
        :param thickness: The thickness of the ply.
        """
        self.material = material
        self.orientation = orientation
        self.thickness = thickness
        self.state = "active"
        pass


class Laminate:

    def __init__(self, layup: list[Ply], name: str, cross_section: float) -> None:
        """
        Class for modeling the behaviour of a laminate comprised of a stack of plies (using the "Ply"-class).
        :param layup: A list of "Ply"-objects.
        :param name: A chosen name for the laminate.
        :param cross_section: The cross section of the laminate.
        """
        self.is_broken = False
        self.layup = layup
        self.name = name
        self.cross_section = cross_section
        self.thickness = self.compute_thickness()
        self.mass = self.compute_mass()
        self.A = self.compute_A()
        self.B = self.compute_B()
        self.D = self.compute_D()
        self.ABD = self.build_stiffness_matrix()
        self.loads, self.defor = [], []

        pass

    def compute_thickness(self) -> float:
        return sum([ply.thickness for ply in self.layup])

    def compute_mass(self) -> float:
        return sum([ply.material["rho"]*ply.thickness for ply in self.layup])

    def compute_A(self) -> np.array:
        A = np.zeros((3, 3), float)
        h_bot = -self.thickness/2
        for ply in self.layup:
            Q = laminatelib.Q2D(ply.material)
            Qt = laminatelib.Q2Dtransform(Q, ply.orientation)
            h_top = h_bot + ply.thickness
            A += Qt*(h_top-h_bot)
            h_bot = h_top
        return A

    def compute_B(self) -> np.array:
        B = np.zeros((3, 3), float)
        h_bot = -self.thickness/2
        for ply in self.layup:
            Q = laminatelib.Q2D(ply.material)
            Qt = laminatelib.Q2Dtransform(Q, ply.orientation)
            h_top = h_bot + ply.thickness
            B += (1/2)*Qt*(h_top**2-h_bot**2)
            h_bot = h_top
        return B

    def compute_D(self) -> np.array:
        D = np.zeros((3, 3), float)
        h_bot = -self.thickness/2
        for ply in self.layup:
            Q = laminatelib.Q2D(ply.material)
            Qt = laminatelib.Q2Dtransform(Q, ply.orientation)
            h_top = h_bot + ply.thickness
            D += (1/3)*Qt*(h_top**3-h_bot**3)
            h_bot = h_top
        return D

    def build_stiffness_matrix(self) -> np.array:
        """
        Build the stiffness matrix from the precomputed(!) 3 x 3 A, B, and D matrices.
        :return: The ABD stiffness matrix.
        """
        ABD = np.zeros((6, 6), float)
        ABD[0:3, 0:3] = self.A
        ABD[0:3, 3:6] = self.B
        ABD[3:6, 0:3] = self.B
        ABD[3:6, 3:6] = self.D
        return ABD

    def update_state(self, **kwargs) -> None:
        self.loads, self.defor = laminatelib.solveLaminateLoadCase(self.ABD, **kwargs)
        pass

