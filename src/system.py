# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 16:54:44 2026

@author: JoelM
"""

from src.body import Body
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import json
import csv

class SolarSystem:
    """
    Parent class for solar system simulation. Uses Beeman Integration.
    """
    def __init__(self):
        """
        Initialise class.

        Returns
        -------
        None.

        """
        
        self.bodies = []
 
    def read_parameters(self, file_path):
        """
        Reads parameters from json file.
        Unpacks and adds bodies present to instance.

        Parameters
        ----------
        file_path : str
            Path to parameters file.

        Raises
        ------
        ValueError
            When values in file are incorrect for simulation.

        Returns
        -------
        None.

        """        
        with open(file_path) as f:
            parameters_solar = json.load(f)
        
        # read timestep with validation
        try:
            timestep = float(parameters_solar["timestep"])
            if timestep <= 0:
                raise ValueError
        except ValueError:
            raise ValueError("Timestep must be a positive real number.")
            
        self.timestep = timestep
        
        # read period with validation
        try:
            period = float(parameters_solar["simulation_period"])
            if timestep <= 0:
                raise ValueError
        except ValueError:
            raise ValueError("Simulation period must be a positive real number.")
        
        # set iterations in accordance to period and timestep 
        self.num_iterations = round(period / timestep)
        
        # check exactly one star in bodies of file
        file_stars = list(filter(lambda body: body["body_type"] == "star", parameters_solar["bodies"]))
        if len(file_stars) != 1:
            raise ValueError("System must include exactly one star.")
        
        # store mass to initalise velocity
        star_mass = file_stars[0]["mass"]
        # custom graviatational constant by unit conversion
        self.G = (4 * np.pi**2) / (star_mass + 1)
        
        # add bodies in file into system
        for file_body in parameters_solar["bodies"]:
            name = file_body["name"]
            body_type = file_body["body_type"]
            mass = file_body["mass"]
            orbital_radius = file_body["orbital_radius"]
            color = file_body["color"]
            
            # check for unique name
            if any(b.name == name for b in self.bodies):
                raise ValueError("Multiple Bodies cannot have the same name.")
            
            # add initial velocity if body is a planet
            if body_type == "planet":
                v = np.sqrt(self.G * star_mass / orbital_radius) # initial velocity estimate
                initial_v = np.array([0, v])
                
                system_body = Body(name, body_type, mass, orbital_radius, color, initial_velocity=initial_v)
                self.add_body(system_body)
            else:
                system_body = Body(name, body_type, mass, orbital_radius, color)
                self.add_body(system_body)

            
    def add_body(self, body):
        """
        Add a body object to instance.

        Parameters
        ----------
        body : Body
            Body object.

        Returns
        -------
        None.

        """
        self.bodies.append(body)
    
    def get_body(self, name):
        """
        Create a pointer to body in instance.

        Parameters
        ----------
        name : str
            Unique name of body in instance.

        Raises
        ------
        KeyError
            Raises if body is not found in instance.

        Returns
        -------
        Body
            Pointer to specified body.

        """
        
        # excract body from system by name
        body = list(filter(lambda b: b.name == name, self.bodies))
        if len(body) == 0:
            raise KeyError("Body does not exist in system.") 
        else:
            return body[0]
    
    
    def run_simulation(self, print_periods=True, write_energy=True, energy_file_name='energy'): 
        """
        Run full simulation of system.
        Update vectors at each timestep for whole time period.
        Detects orbital periods of each body.
        Optionally print average orbital period and/or write energy periodically to csv file.

        Parameters
        ----------
        print_periods : bool, optional
            Switch for printing average orbital period of bodies. The default is True.
        write_energy : bool, optional
            Switch for writing system energy to file periodically. The default is True.
        energy_file_name : str, optional
            Name for .csv file of system energy. The default is 'energy'.

        Returns
        -------
        None.

        """
        
        # initialise acceleration for simulations
        for body in self.bodies:
            a = self.calculate_acceleration(body.name)
            body.current_acceleration = a
            body.previous_acceleration = a # assume preious acceleration is current acceleration (introduces small error but simple)
        
        # pointer to system star
        star = list(filter(lambda b: b.body_type == "star", self.bodies))[0]
        
        if write_energy:
            BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # project root directory
            # create csv file to store timeseries of energy; overwrites if same name (re-run simulation)
            with open(os.path.join(BASE_DIR, "data", f"{energy_file_name}.csv"), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['time', 'energy']) # header
        
        # run simulation for number iterations
        for i in range(self.num_iterations):
                
            self.update_vectors()
                
            # check orbital period completed
            y0 = star.current_position[1] # set y axis to be helio-centric; star shifts over time
            for body in self.bodies:
                if body.name == star.name: # skip star as star should be held constant
                    continue
                
                # obrit completetion condition
                if body.position_history[-2, 1] < y0 and body.position_history[-1, 1] >= y0:
                    t_now = i * self.timestep
                    # keep time periodic (since last orbit achieved)
                    if body.last_crossing_time is not None:
                        orbit_time = t_now - body.last_crossing_time
                        body.orbital_periods.append(orbit_time)
                    body.last_crossing_time = t_now
               
            if write_energy:
                if i % 30 == 0: # every 30 iterations
                    t_now = i * self.timestep
                    energy = self.calculate_total_system_energy()
                    # reopen file to append
                    with open(os.path.join(BASE_DIR, "data", f"{energy_file_name}.csv"), 'a', newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow([t_now, energy]) # write new row
        
        if print_periods: 
            for body in self.bodies:
                if body.name == star.name: # skip star as star should be held constant
                    continue
                # print mean of body's orbital period
                print(f"Orbital period of {body.name}: {np.mean(body.orbital_periods):.4f}")
    
    def update_vectors(self):
        """
        Update vectors of every massive body in instance.
        Update done by Beeman Algorithm.

        Returns
        -------
        None.

        """
        # update position of all bodies
        for body in self.bodies:
            new_position = body.current_position + body.current_velocity*self.timestep + (1/6)*(4*body.current_acceleration - body.previous_acceleration)*(self.timestep**2)
            body.position_history = np.vstack((body.position_history, new_position)) # add new position to history for visualisation
            body.current_position = new_position
            
        # calculate new accelerations and velocity
        for body in self.bodies:
            new_acceleration = self.calculate_acceleration(body.name)
            new_velocity = body.current_velocity + (1/6)*(2*new_acceleration + 5*body.current_acceleration - body.previous_acceleration)*self.timestep
            
            # update for next timestep
            body.current_velocity = new_velocity
            body.previous_acceleration = body.current_acceleration
            body.current_acceleration = new_acceleration
    
    def calculate_acceleration(self, body_name):
        """
        Calculate the acceleration on specified body in N-body system.

        Parameters
        ----------
        body_name : str
            Unique identifier of body in instance.

        Returns
        -------
        acceleration : ndarray, shape (2, 1)
            Calculated acceleration of body.

        """
        # pointer to body we are calculating for
        body_j = list(filter(lambda body: body.name == body_name, self.bodies))[0]
        
        # summation as defined
        summation = 0
        for body_i in self.bodies:
            
            if body_i.name == body_name: # skip if i=j
                continue
            
            if body_i.body_type == "satellite": # skip if body is a satellite (negligible mass)
                continue
            
            mass_i = body_i.mass
            
            r_ij = body_i.current_position - body_j.current_position
            
            length_r_ij = np.linalg.norm(r_ij)
            
            if length_r_ij == 0:
                length_r_ij = 1e-16 # minimum distance to resolve (softening)
            
            term = ((mass_i) / (length_r_ij ** 3)) * r_ij 
            summation += term
        
        
        acceleration = summation * self.G  
        
        return acceleration
    
    def calculate_total_potential_energy(self):
        """
        Calculate total potential enrgy of system at current time.

        Returns
        -------
        U : float
            Potential energy in (EarthMass)AU/yr.

        """
        U = 0
        
        for body_i in self.bodies:
            for body_j in self.bodies:
                # i not equal j
                if body_i.name == body_j.name:
                    continue
                
                r_ij = body_i.current_position - body_j.current_position
                dist = np.linalg.norm(r_ij)
                
                if dist == 0:
                    dist = 1e-16 # minimum distance to resolve (softening)
                
                U += self.G * body_i.mass * body_j.mass / dist
                
        U *= -0.5 # overcounting correction
        
        return U
    
    def calculate_total_kinetic_energy(self):
        """
        Calculate kintetic energy across all bodies.

        Returns
        -------
        K : float
            Kinetic energy in (EarthMass)AU/yr.

        """       
        K = 0
        
        for body in self.bodies:
            K += body.calculate_kinetic_energy()
            
        return K
        
    def calculate_total_system_energy(self):
        """
        Calculate total system energy at current time.
        Sum of potenital and kinetic energies.

        Returns
        -------
        E : float
            Total system energy in (EarthMass)AU/yr.

        """
        
        K = self.calculate_total_kinetic_energy() # kinetic energy
        U = self.calculate_total_potential_energy() # potential energy
        E = K + U
        
        return E

    def animation(self, title="Solar System Simulation"):
        """
        Create animation of system based on bodies stored position histories.
        Bodies are sized on a nromalised log-scale of their masses.
        Animation is helio-centric as sun drifts slowly over time.
        5ms between frames.

        Parameters
        ----------
        title : str, optional
            Figure title. The default is "Solar System Simulation".

        Returns
        -------
        fig : matplotlib.figure
            Figure of animation.
        ax : matplotlib.axes
            Axes of animation.
        anim : matplotlib.animation.FuncAnimation
            Animation object.

        """
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_facecolor('black')
        ax.set_aspect('equal')
        
        # star pointer to make animation heliocentric
        star = list(filter(lambda b: b.body_type == "star", self.bodies))[0]
        star_positions = star.position_history
        
        # work out axis limits from position history
        all_positions = np.vstack([body.position_history - star_positions for body in self.bodies])
        max_r = np.max(np.abs(all_positions)) * 1.2 # maximum radius with buffer for plot
        ax.set_xlim(-max_r, max_r)
        ax.set_ylim(-max_r, max_r)
        
        # all masses for scaling patch sizes
        masses = np.array([body.mass for body in self.bodies])
        

        # create a patch for each body 
        self.patches = {}
        for body in self.bodies:
            # normalise mass on log scale to get relative patch sizes
            norm = (np.log(body.mass) - np.log(masses.min())) / (np.log(masses.max()) - np.log(masses.min()))
            patch_radius = max_r * 0.5 * (0.005 + norm * 0.025)
            # initial patch position
            x0 = body.position_history[0, 0] - star_positions[0, 0]
            y0 = body.position_history[0, 1] - star_positions[0, 1]
            # create the patch on the axes
            patch = plt.Circle((x0, y0), patch_radius, color=body.color, animated=True, label=body.name)
            ax.add_patch(patch)
            self.patches[body.name] = patch
        
        # customise plot
        ax.set_title(title, fontsize=20)
        ax.legend(title='Bodies', bbox_to_anchor=(1, 0.5), loc="center left", title_fontsize=13, fontsize=10)
        
        # remove axes ticks
        ax.set_xticks([])
        ax.set_yticks([])
        
        # create animator
        self.anim = FuncAnimation(fig, self.animate_step, frames=len(self.bodies[0].position_history), interval=5, blit=True)
        
        plt.tight_layout()
        
        return fig, ax, self.anim
    
    def animate_step(self, i):
        """
        Animation update at each step.

        Parameters
        ----------
        i : int
            Current step of animation.

        Returns
        -------
        list
            Values for next frame.

        """
        
        # star pointer to make animation heliocentric
        star = list(filter(lambda b: b.body_type == "star", self.bodies))[0]
        star_positions = star.position_history
        
        for body in self.bodies:
            x = body.position_history[i, 0] - star_positions[i, 0]
            y = body.position_history[i, 1] - star_positions[i, 1]
            self.patches[body.name].center = (x, y) # update patch position
        
        return list(self.patches.values()) # return values to plot in animation
        

class SolarSystemEulerCromer(SolarSystem):
    """
    Subclass of SolarSystem.
    Overwrites update_vectors method to use Euler-Cromer algorithm.
    """
    def update_vectors(self):
        """
        Update vectors of every massive body in instance.
        Update done by Euler-Cromer Algorithm.

        Returns
        -------
        None.

        """
        for body in self.bodies:
            # calcule velocity and position
            new_velocity = body.current_velocity + body.current_acceleration*self.timestep
            new_position = body.current_position + new_velocity*self.timestep
            # update velocity and position
            body.current_velocity = new_velocity
            body.position_history = np.vstack((body.position_history, new_position))
            body.current_position = new_position
            
        # update acceleration for next timestep
        for body in self.bodies:
            body.current_acceleration = self.calculate_acceleration(body.name)


class SolarSystemDirectEuler(SolarSystem):
    """
    Subclass of SolarSystem.
    Overwrites update_vectors method to use Direct Euler algorithm.
    """
    def update_vectors(self):
        """
        Update vectors of every massive body in instance.
        Update done by Direct Euler Algorithm.

        Returns
        -------
        None.

        """
        for body in self.bodies:
            # calcule velocity and position
            new_position = body.current_position + body.current_velocity*self.timestep
            new_velocity = body.current_velocity + body.current_acceleration*self.timestep
            # update velocity and position
            body.position_history = np.vstack((body.position_history, new_position))
            body.current_position = new_position
            body.current_velocity = new_velocity
            
        # update acceleration for next timestep
        for body in self.bodies:
            body.current_acceleration = self.calculate_acceleration(body.name)
            

class SolarSystemSatelliteGrid(SolarSystem):
    """
    Subclass of SolarSystem to simulate multiple satellite paths to mars.
    """
    
    def add_satellite_grid(self, vx_range, vy_range, initial_separation=0.001):
        """
        Add (n x m) grid of satellites based on intial velocities to system.
        Initialises required variables.

        Parameters
        ----------
        vx_range : tuple of (float, float, float)
            Elements are (minimum x speed (AU/yr), maximum x speed (AU/yr), number of points).
        vy_range : tuple of (float, float, float)
            Elements are (minimum y speed (AU/yr), maximum y speed (AU/yr), number of points)
        initial_separation : float, optional
            How far off earth on positive x-axis that paths begin in AU. The default is 0.001.

        Returns
        -------
        None.

        """
        # unpack tuples to create vector of velocities
        vx = np.linspace(vx_range[0], vx_range[1], vx_range[2])
        vy = np.linspace(vy_range[0], vy_range[1], vy_range[2])
        
        # store to label heatmap
        self.vx_vals = vx
        self.vy_vals = vy
        
        # create grid
        VX, VY = np.meshgrid(vx, vy)
        
        self.grid_shape = VX.shape # to reshape later
        N = VX.size # number of satellites
        
        # repeat initial position for all satellites; creates (N, 2) array - each position is (1, 2)
        self.sat_positions = np.tile([1 + initial_separation, 0.0], (N, 1)) 
        # form (N, 2) array of initial velocities
        self.sat_velocities = np.column_stack([VX.flatten(), VY.flatten()])
        
        # Beeman requires previous and current acceleration
        # initialise both to the same value (same small error as parent class)
        a = self.calculate_acceleration_grid()
        self.sat_accelerations = a.copy()  # current
        self.sat_previous_accelerations = a.copy()  # previous

        
        # track only min distance and time it takes
        self.sat_min_distance_mars = np.full(N, np.inf)
        self.sat_min_distance_earth = np.full(N, np.inf)
        self.sat_min_time_mars = np.zeros(N)
        self.sat_min_time_earth = np.zeros(N)
        
    def calculate_acceleration_grid(self):
        """
        Vectorised acceleration calculation for all satellites.

        Returns
        -------
        a : ndarray, shape (N, 2)
            acceleration vectors for each satellite in grid.

        """
        a = np.zeros_like(self.sat_positions)
        
        # calculate acceleration casued by massive bodies
        for body in self.bodies:
            
            if body.body_type == "satellite": # skip if body is a satellite (negligible mass); redundancy
                continue
            
            r = body.current_position - self.sat_positions
            
            length_r = np.linalg.norm(r, axis=1, keepdims=True)
            
            a += self.G * body.mass / length_r**3 * r
            
        return a
    
    def update_step_satellites(self):
        """
        Update satellite vectors (vectorised) by the Beeman Algorithm

        Returns
        -------
        None.

        """
        self.sat_positions += (self.sat_velocities * self.timestep + (1/6) * (4*self.sat_accelerations - self.sat_previous_accelerations) * self.timestep**2)
        
        a_new = self.calculate_acceleration_grid()

        self.sat_velocities += (1/6) * (2*a_new + 5 * self.sat_accelerations - self.sat_previous_accelerations) * self.timestep
        
        # update accelerations for next iteration
        self.sat_previous_accelerations = self.sat_accelerations
        self.sat_accelerations = a_new
    
    def update_min_distance_to(self, target_name, step_index):
        """
        Update the minimum distance to target if closer than previous closest distance.

        Parameters
        ----------
        target_name : str
            Unique name identifier of target body.
        step_index : int
            Current step of iteration.

        Returns
        -------
        None.

        """
        
        target = self.get_body(target_name)
        r = self.sat_positions - target.current_position # position vector from target to satellites
        distances = np.linalg.norm(r, axis=1)
        
        if target_name == "Mars":
            # vectorised updates of cloesrr to target and minimum time
            closer = distances < self.sat_min_distance_mars # mask for closer position
            self.sat_min_distance_mars = np.where(closer, distances, self.sat_min_distance_mars) 
            self.sat_min_time_mars = np.where(closer, step_index * self.timestep, self.sat_min_time_mars)
            
        if target_name == "Earth":
            # vectorised updates of cloesrr to target and minimum time
            condition = (distances < self.sat_min_distance_earth) & (step_index*self.timestep > self.sat_min_time_mars) # mask for closer position
            self.sat_min_distance_earth = np.where(condition, distances, self.sat_min_distance_earth) 
            self.sat_min_time_earth = np.where(condition, step_index * self.timestep, self.sat_min_time_earth)
    
    def run_simulation(self):
        """
        Overwrites parent method to include updating satellites.
        Removes writing energy and updating orbital periods.

        Returns
        -------
        None.

        """
        # initialise acceleration for simulations
        for body in self.bodies:
            a = self.calculate_acceleration(body.name)
            body.current_acceleration = a
            body.previous_acceleration = a # assume preious acceleration is current acceleration (introduces small error but simple)
       
        for i in range(self.num_iterations):
            # update massive bodies
            self.update_vectors()
            # update all satellites in one pass
            self.update_step_satellites()
            self.update_min_distance_to("Mars", i)
            self.update_min_distance_to("Earth", i)
    
    def closest_approach_mars_conditions(self):
        """
        Find the closests approach to mars and the conditions of it.

        Returns
        -------
        vx_star : float
            Optimal x speed for cloests approach in AU/yr.
        vy_star : float
            Optimal y speed for cloests approach in AU/yr.
        min_distance : float
            Minimum separation in AU.
        min_distance_time : float
            Time min_distance occured in yr

        """
        # find satellite with closest approach
        flat_best_idx = np.argmin(self.sat_min_distance_mars)
        min_distance = self.sat_min_distance_mars[flat_best_idx]
        min_distance_time = self.sat_min_time_mars[flat_best_idx]
        
        # find corresponding vx and vy values
        j, i = np.unravel_index(flat_best_idx, self.grid_shape)
        vx_star = self.vx_vals[i]
        vy_star = self.vy_vals[j]
        
        return vx_star, vy_star, min_distance, min_distance_time
    
    def close_to_mars_satellites(self, tol=0.006):
        """
        Find all initial velocities which allow satellite to get within 'tol' of Mars.

        Parameters
        ----------
        tol : float, optional
            Value minimum distance must be within in AU. The default is 0.006.

        Returns
        -------
        velocities : ndarray, shape (n, 2)
            Array of all (n) velocities which allow closest appraoch within tolerance.

        """ 
        close = self.sat_min_distance_mars < tol # mask for close values
        close_idx = np.where(close)[0] 
        
        # find corresponding vx and vy
        j, i = np.unravel_index(close_idx, self.grid_shape)
        vx_star = self.vx_vals[i]
        vy_star = self.vy_vals[j]
        velocities = np.column_stack((vx_star, vy_star)) # put back into pairs
        
        return velocities
    
    def stats_from_velocity(self, vx, vy):
        """
        Finds statistics from satellite with initial velocity closest to input.

        Parameters
        ----------
        vx : float
            Initial x-speed to search for.
        vy : float
            Initial x-speed to search for.

        Returns
        -------
        dist_mars : float
            Minimum distance from Mars in AU.
        time_mars : TYPE
            Time which minimum distance from Mars occurs in yr.
        dist_earth : TYPE
            Minimum distance from Earth in AU.
        time_earth : TYPE
            Time which minimum distance from Earth occurs in yr.

        """
        # find satellite with closest to input velocity
        i = np.argmin(np.abs(self.vx_vals - vx))
        j = np.argmin(np.abs(self.vy_vals - vy))
        flat_idx = np.ravel_multi_index((j, i), self.grid_shape)
        
        # stats to return
        dist_mars = self.sat_min_distance_mars[flat_idx]
        time_mars = self.sat_min_time_mars[flat_idx]
        
        dist_earth = self.sat_min_distance_earth[flat_idx]
        time_earth = self.sat_min_time_earth[flat_idx]
        
        return dist_mars, time_mars, dist_earth, time_earth
        
    def plot_heatmap(self):
        """
        Plots heatmap of all initial velocities minimum distance from mars.

        Returns
        -------
        fig : matplotlib.figure
            Figure of plot.
        ax : matplotlig.axes
            Axes of plot.

        """
        # grid of minium distances for plot
        grid = self.sat_min_distance_mars.reshape(self.grid_shape)
    
        fig, ax = plt.subplots()
        # plot heatmap with velocity axes
        im = ax.imshow(grid, cmap="plasma_r", origin='lower', aspect='auto', extent=[self.vx_vals[0], self.vx_vals[-1], self.vy_vals[0], self.vy_vals[-1]])
        plt.colorbar(im, ax=ax, label='Minimum Distance from Mars (AU)')
        ax.set_xlabel(r'$v_x$ (AU/yr)')
        ax.set_ylabel(r'$v_y$ (AU/yr)')
    
        return fig, ax
    
    

def main():
    print("Root File for SolarSystem Class and Subclasses")
    
if __name__ == "__main__":
    main()