# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 16:53:18 2026

@author: JoelM
"""
import numpy as np

class Body:
    """
    Class to create any body in solar system simulation.
    """
    
    def __init__(self, name, body_type, mass, orbital_radius, color, initial_velocity=None):
        """
         Initialise class.

        Parameters
        ----------
        name : str
            Unique identifing name for the body.
        body_type : str
            What type of body instance is. Generally one of "star", "planet", or "satellite".
        mass : float
            Mass of body in earth masses.
        orbital_radius : float
            Orbital radius of body in AU.
        color : str
            Colour of body for visualisation.
        initial_velocity : ndarray of shape (2, 1), dtype float64, optional
            Initial velocity of body. The default is None.

        Returns
        -------
        None.

        """
        
        
        # makes unique to each instance, otherwise an edit would edit for all carrying default value
        if initial_velocity is None:
            initial_velocity = np.zeros(2)
        
        self.name = name
        self.body_type = body_type
        self.mass = mass
        self.color = color
        self.initial_position = np.array([orbital_radius, 0]) # bodies initalised on x-axis
        self.initial_velocity = initial_velocity
        
        # current vectors
        self.current_position = self.initial_position
        self.current_velocity = self.initial_velocity
        self.current_acceleration = np.zeros(2)
        
        self.previous_acceleration = np.zeros(2) # for beeman integration
        
        self.position_history = self.initial_position.reshape(1, 2)
        
        # values for orbital periods (used in simulation class)
        self.orbital_periods = [] 
        self.last_crossing_time = 0 # point when body last crossed helio-centric y-axis
        
    def calculate_kinetic_energy(self):
        """
        Calculate the current kinetic energy of body.

        Returns
        -------
        float
            Current kinetic energy in (EarthMass)AU/yr.

        """
        
        v = self.current_velocity
        return 0.5 * self.mass * (v @ v)
 
    
def main():
    print("Root File for Body Class")
    
if __name__ == "__main__":
    main()