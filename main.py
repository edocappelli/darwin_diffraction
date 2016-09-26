import plotting

#todo: in fact, all plotting.py can be copied here



# then

#
# input Crystal
#
diffraction_experiment = plotting.create_diffraction_setup(use_defaults=False)


# todo: rename plotter to something more readable like "scan_theta_microrads"
plotting.plotter(diffraction_experiment,start=-200, stop=200, num=1000)


