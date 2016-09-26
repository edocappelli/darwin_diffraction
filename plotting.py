"""
The following code follows the treatment in "Theory of X-Ray Diffraction in Crystals" by William H. Zachariasen
"""

import numpy as np
import matplotlib.pyplot as plt
import diffraction


#def create_instance():
def create_diffraction_setup(use_defaults=False):

    """
    asks user to insert the diffraction setup parameters and creates an instance of the DiffractionSetup subclass
    :return instance
    """

    if use_defaults:
        crystal_type = "Si"
        asymmetry_angle = 0.0
        miller_h = 1
        miller_k = 1
        miller_l = 1
        wavelength = 8.0
        polarization = 's'
        polarization_bool = True
    else:
        crystal_type = input("Please insert crystal type (e.g. Si): ")
        asymmetry_angle = float(input("Please insert the asymmetry angle in degrees: "))
        miller_h = int(input("Please insert Miller H: "))
        miller_k = int(input("Please insert Miller K: "))
        miller_l = int(input("Please insert Miller L: "))
        # todo: do not call wavelength to photon_energy, very confusing!!
        wavelength = float(input("Please insert the x-ray energy in keV: "))
        polarization = input("Please insert the x-ray polarization ('s' or 'p'): ")
        polarization_bool = True

    if polarization == "s":
        polarization_bool = True

    elif polarization == "p":
        polarization_bool = False

    diffraction_experiment = diffraction.Diffraction(crystal_type, asymmetry_angle,
                                       miller_h, miller_k, miller_l, wavelength, polarization_bool)

    return diffraction_experiment


def plotter(diffraction_experiment,start=-200, stop=200, num=1000):

    # diffraction_experiment = create_instance()


    theta_diff = np.linspace(start, stop, num)  # in micro radians


    # TODO: this part should not be here, as this is only for plotting, not for calculation

    # everything here is in Diffraction instance, except the inputs theta_diff
    # thus should be a new method in Diffraction:
    # def darwin_reflectivity(self,theta_diff):
    #   ...
    #   return reflectivity
    theta_bragg = diffraction_experiment.bragg_angle()
    alpha = 2 * theta_diff * 1e-6 * np.sin(2 * theta_bragg)  # Zachariasen [3.116]
    b = diffraction_experiment.parameter_b()
    k = diffraction_experiment.parameter_k()
    psi_0 = np.ones(num) * diffraction_experiment.parameter_psi_0()
    psi_h = diffraction_experiment.parameter_psi_h()

    y = ((0.5 * (1-b) * psi_0) + (0.5 * b * alpha)) / (np.sqrt(np.absolute(b)) * k * np.absolute(psi_h))

    reflectivity = np.ones(num)

    for i in np.arange(num):  # Zachariasen [3.192b]

        if y[i] < -1:
            reflectivity[i] = (-y[i] - np.sqrt(y[i] ** 2 - 1)) ** 2

        elif y[i] > 1:
            reflectivity[i] = (y[i] - np.sqrt(y[i] ** 2 - 1)) ** 2

        else:
            y[i] = 1

    ##################

    plt.figure()
    plt.plot(theta_diff, reflectivity)
    plt.title("Darwin curve")
    plt.xlabel("angular displacement from Bragg angle [urad]")
    plt.ylabel("power reflectivity")
    plt.show()


