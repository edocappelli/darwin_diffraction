

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

    # crystal_data = Crystal()...
    diffraction_experiment = diffraction.Diffraction(crystal_type, asymmetry_angle,
                                       miller_h, miller_k, miller_l, wavelength, polarization_bool)

    return diffraction_experiment


def plotter(theta_diff,reflectivity):





    plt.figure()
    plt.plot(theta_diff, reflectivity)
    plt.title("Darwin curve")
    plt.xlabel("angular displacement from Bragg angle [urad]")
    plt.ylabel("power reflectivity")
    plt.show()





if __name__ == "__main__":

    diffraction_experiment = create_diffraction_setup(use_defaults=True)


    theta_diff = np.linspace(-200,  200, 100)  # in micro radians


    reflectivity = diffraction_experiment.darwin_profile_theta(theta_diff)


    plotter(theta_diff,reflectivity)


