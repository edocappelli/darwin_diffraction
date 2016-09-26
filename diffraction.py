"""
The following code follows the treatment in "Theory of X-Ray Diffraction in Crystals" by William H. Zachariasen
"""

import numpy as np
from scipy import constants


class Crystal(object):
    """
    This class stores crystal data read from a data file (e.g. Si.txt, Ge.txt, etc.)

    """
    # TODO: asymmetry angle not needed here, used in  Diffraction
    def __init__(self, crystal_type, asymmetry_angle):

        """
        Initializer.
        :param crystal_type: <class 'str'> string specifying which crystal is being used
        :param asymmetry_angle: <class 'float'> angle in degrees
        """
        self._crystal_type = crystal_type
        self._asymmetry_angle = asymmetry_angle * np.pi / 180  # converting into radians

        """
        The files are structured in the following way (each row has 11 elements):
        (0, 0)-(0, 4) -> 5x{a_i} | (0, 5) -> c | (0, 6)-(0, 10) -> 5x{b_i}     goal:f0 calculation
        (1, 0)-(1, 5) -> a b c alpha beta gamma | (1, 6)-(1, 10) -> 5x0     goal:unit cell volume calculation
        (2-9, 0-10) -> atomic_number fraction x y z 6x0     goal:structure factor calculation
        """

        self._data_array = np.loadtxt("{}.txt".format(self._crystal_type), dtype='float')
        self._a = []  # parameters a_i (x5)
        self._b = []  # parameters b_i (x5)
        self._unit_cell = []  # unit cell sides (x3) and angles (x3)
        self._atom_coordinates = np.zeros((8, 3))  # each row holds the x, y, z coordinates of one atom in the unit cell
        # there are 8 atoms in the silicon unit cell

        for i in np.arange(self._data_array.shape[0]):

            for j in np.arange(self._data_array.shape[1]):

                if i == 0:

                    if j < 5:
                        self._a.append(self._data_array[(i, j)])

                    elif 5 < j < self._data_array.shape[1]:
                        self._b.append(self._data_array[(i, j)])

                    elif j == 5:
                        self._c = self._data_array[(i, j)]

                elif i == 1 and j < 6:

                    self._unit_cell.append(self._data_array[(i, j)])

                else:

                    if j in [2, 3, 4]:
                        self._atom_coordinates[(i - 2, j - 2)] = self._data_array[(i, j)]

    def unit_cell_volume(self):

        """
        calculates the unit cell volume.
        Conventions: alpha is the angle between b and c, beta between a and c, gamma between a and b,
        with a, b being on the horizontal plane.
        Volume = a * b * c * sin(alpha) * sin(beta) * sin(gamma)
        :return: unit cell volume in Angstroms^3
        """
        a = self._unit_cell[0]
        b = self._unit_cell[1]
        c = self._unit_cell[2]
        alpha = self._unit_cell[3] * np.pi / 180  # converting into radians
        beta = self._unit_cell[4] * np.pi / 180  # converting into radians
        gamma = self._unit_cell[5] * np.pi / 180  # converting into radians

        return a * b * c * np.sin(alpha) * np.sin(beta) * np.sin(gamma)  # Angstroms ^ 3


class Diffraction(Crystal):
    """
    This class defines a diffraction geometry for the crystal and the incoming x-ray beam.
    In particular, this class is initialized by specifying the wavelength of the beam and then this is used to calculate
    the Bragg angle, the b parameter and the structure factors.
    """

    # todo: not here: global variables?? use "import scipy.constants as constants" and use "constants.c" etc
    e = constants.elementary_charge  # charge of the proton
    c = constants.speed_of_light  # speed of light in vacuum
    m = constants.electron_mass  # mass of the electron
    h = constants.Planck  # Planck's constant

    def __init__(self, crystal_type, asymmetry_angle, miller_h, miller_k, miller_l, wavelength, polarization):

        """
        :param miller_h: < class 'int'> Miller index H
        :param miller_k: < class 'int'> Miller index K
        :param miller_l: < class 'int'> Miller index L
        :param wavelength: < class 'float'> wavelength of the x-ray beam in keV
        :param polarization: < class 'bool'> True = normal polarization (s), False = parallel polarization (p)
        """

        super(Diffraction, self).__init__(crystal_type, asymmetry_angle)

        self._wavelength = self.h * self.c * 1e+7 / (self.e * wavelength)  # conversion keV -> Angstrom
        self._polarization = polarization
        self._miller_h = miller_h
        self._miller_k = miller_k
        self._miller_l = miller_l

    def d_spacing(self):

        """
        calculates the spacing between the Miller planes in the special case of a cubic unit cell
        :return: d spacing in Angstroms
        """

        return self._unit_cell[0] / np.sqrt(self._miller_h ** 2 + self._miller_k ** 2 + self._miller_l ** 2)

    def bragg_angle(self):

        """
        calculates the Bragg angle using Bragg's law: lambda = 2 * d_spacing * sin(theta_Bragg) (Zachariasen [3.5])
        :return: Bragg angle in radians
        """

        d_spacing = self.d_spacing()  # in Angstroms
        return np.arcsin(0.5 * self._wavelength / d_spacing)

    def parameter_f_0(self):

        """
        calculates f0 following the 5-gaussian interpolation by D. Waasmaier & A. Kirfel
        :return: atomic scattering power, f0 @ specified wavelength
        """

        f = self._c
        s = np.sin(self.bragg_angle()) / self._wavelength  # in Angstrom ^ -1
        for i in np.arange(5):

            f += self._a[i] * np.exp(-self._b[i] * s ** 2)

        return f

    def parameter_b(self):

        """
        calculates b = gamma_0 / gamma_h (Zachariasen [3.115])
        :return: b
        """

        gamma_h = - np.sin(self.bragg_angle())
        gamma_0 = np.sin(self.bragg_angle() + self._asymmetry_angle)

        return gamma_0 / gamma_h

    def parameter_k(self):

        """
        k = 1 for normal polarization; k = |cos(2*theta_Bragg)| for parallel polarization (Zachariasen [3.141])
        :return: k
        """

        if self._polarization:
            return 1

        if not self._polarization:
            return np.absolute(np.cos(2*self.bragg_angle()))

    def parameter_psi_0(self):

        """
        psi_0 = num_atoms_in_unit_cell * f_0 (Zachariasen [3.95])
        :return: psi_0
        """
        f_0 = self.parameter_f_0()
        v = self.unit_cell_volume() * 1e-30  # Angstrom ^3 -> m ^3
        wavelength = self._wavelength * 1e-10  # Angstrom -> metre
        e_radius = 2.8179403267e-15  # cgs units
        psi_0 = - e_radius * wavelength ** 2 / (v * np.pi) * f_0

        return (self._data_array.shape[0] - 2) * psi_0  # this only holds for crystals with atoms all of the same kind

    def parameter_f_h(self):

        """
        F_H = f_0 * SUM(i) { exp[ (i*2*pi) * (h*a_i + k*b_i + l*c_i) ] }
        where a_i, b_i, c_i are the coordinates of the i-th atom in the unit cell (Zachariasen [3.52])
        and h, k, l are the Miller indices.
        :return: F_H
        """

        f_h = 0 + 0j

        for i in np.arange(self._data_array.shape[0] - 2):

            j = 0

            f_h += np.exp(1j * 2 * np.pi * (self._miller_h * self._atom_coordinates[(i, j)]
                                            + self._miller_k * self._atom_coordinates[(i, j + 1)]
                                            + self._miller_l * self._atom_coordinates[(i, j + 2)]))

        return f_h * self.parameter_f_0()

    def parameter_psi_h(self):

        """
        calculates psi_H starting from F_H, the wavelength and the unit cell volume (Zachariasen [3.95])
        :return: psi_H
        """

        v = self.unit_cell_volume() * 1e-30  # Angstrom ^3 -> m ^3
        wavelength = self._wavelength * 1e-10  # Angstrom -> metre
        f_h = self.parameter_f_h()
        e_radius = 2.8179403267e-15  # cgs units

        return - e_radius * (wavelength ** 2 / (v * np.pi)) * f_h
