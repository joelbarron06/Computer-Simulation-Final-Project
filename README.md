# Computer Simulation Final Project

This repository contains code to simulate a 2D N-body gravitational system (/src) and experiments on it (/test).

System is simulated accoring to a parameters JSON file. See files in /parameters for examples. Must include "timestep" and "simulation_period". Must also include "bodies", specifing a collection of bodies in the system, each with "name", "body_type", "mas", "orbital_radius", and "color". There must be exactly one body of type "star".

See /test/basic_simulation.py for simple use case. General structure is as follows:
* Initialise 'SolarSystem' class
* Use 'read_parameters' method to read parameters file (above)
* Use 'run_simulation' method to run
* (optional) use 'animate' method to create simulation animation

'SolarSystemEulerCromer' and 'SolarSystemDirectEuler' are used in the same way, with alternative integration methods.

To search for launch conditions for martian probes, use `SolarSystemSatelliteGrid' and use the 'add_satellite_grid' method after reading parameters (before simulating) to add an $m \times n$ grid of velocities in specified range. After simulating, specific satellites are found based on velovity, using 'close_to_mars_satellites' to extract velocities that allow close approch, and 'stats_from_velocity' to return summary statistics of a specfic satellite (identified by velocity). Use 'plot_heatmap' to create a heatmap showing how close to Mars each satellite gets. Note that to animate a satellite, you must use the root class, and add a custom 'Body' with desired mass and trajectory using the 'add_body' method of 'SolarSystem' before simulating and animating. 
