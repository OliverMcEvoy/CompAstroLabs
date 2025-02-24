import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
from astroquery.jplhorizons import Horizons
from tqdm import tqdm
import argparse
from MassAndObjectInfo import get_masses_and_object_info

def calculate_orbit(objects, initial_time_step=3600, max_steps=8766):
      '''
      Calculate the orbit of multiple objects in 3D space.
      @param objects: A list of dictionaries, each containing the
            name, mass, initial x, y, z positions, initial
            x, y, z velocities of an object.
      @param initial_time_step: The initial time step for the simulation.
      @param min_time_step: The minimum time step for the simulation.
      @param max_steps: The maximum number of steps to simulate.
      @param tolerance: The tolerance for the simulation.
      @return: A dictionary containing the positions of each object
            at each time step, the time at each step, and the time    
            at which each object crosses the y-axis.
      '''
      # Definine G
      g = 6.67430e-11  

      # Initialise arrays for each object
      positions = {obj['name']: {'x': np.zeros(max_steps), 'y': np.zeros(max_steps), 'z': np.zeros(max_steps)} for obj in objects}
      velocities = {obj['name']: {'x': np.zeros(max_steps), 'y': np.zeros(max_steps), 'z': np.zeros(max_steps)} for obj in objects}
      time = np.zeros(max_steps)
      
      # Object initial conditions
      for obj in objects:
            positions[obj['name']]['x'][0] = obj['initial_x']
            positions[obj['name']]['y'][0] = obj['initial_y']
            positions[obj['name']]['z'][0] = obj['initial_z']
            velocities[obj['name']]['x'][0] = obj['initial_x_velocity']
            velocities[obj['name']]['y'][0] = obj['initial_y_velocity']
            velocities[obj['name']]['z'][0] = obj['initial_z_velocity']
      
      time[0] = 0

      time_step = initial_time_step

      # TODO Track radial velocity changes as a better indication of an orbit
      previous_radial_velocity = {obj['name']: None for obj in objects}
      radial_velocity_sign_changes = {obj['name']: [] for obj in objects}
      y_crossings = {obj['name']: [] for obj in objects}
      
      for i in tqdm(range(0, max_steps-1), desc="Calculating orbits", mininterval=1.0):
            for obj in objects:
                  name = obj['name']
                  # Radius not needed anymore
                  # radius = np.sqrt(positions[name]['x'][i]**2 + positions[name]['y'][i]**2 + positions[name]['z'][i]**2)
                  
                  x_acceleration = 0
                  y_acceleration = 0
                  z_acceleration = 0

                  # A bit much indentation but it'll do for now, but what this chunk does is calculate the acceleration of the object due to all other objects.
                  for other_obj in objects:
                        # Cursed check but itll do.
                        if other_obj['name'] != name:
                              other_radius = ((positions[other_obj['name']]['x'][i] - positions[name]['x'][i])**2 + 
                                                      (positions[other_obj['name']]['y'][i] - positions[name]['y'][i])**2 + 
                                                      (positions[other_obj['name']]['z'][i] - positions[name]['z'][i])**2) ** (0.5)
                              x_acceleration += -g * other_obj['mass'] * (positions[name]['x'][i] - positions[other_obj['name']]['x'][i]) / other_radius**3
                              y_acceleration += -g * other_obj['mass'] * (positions[name]['y'][i] - positions[other_obj['name']]['y'][i]) / other_radius**3
                              z_acceleration += -g * other_obj['mass'] * (positions[name]['z'][i] - positions[other_obj['name']]['z'][i]) / other_radius**3
                             
                  # Update velocities and positions
                  velocities[name]['x'][i+1] = velocities[name]['x'][i] + x_acceleration * time_step
                  velocities[name]['y'][i+1] = velocities[name]['y'][i] + y_acceleration * time_step
                  velocities[name]['z'][i+1] = velocities[name]['z'][i] + z_acceleration * time_step
                  positions[name]['x'][i+1] = positions[name]['x'][i] + velocities[name]['x'][i+1] * time_step
                  positions[name]['y'][i+1] = positions[name]['y'][i] + velocities[name]['y'][i+1] * time_step
                  positions[name]['z'][i+1] = positions[name]['z'][i] + velocities[name]['z'][i+1] * time_step

                  # Skip the Sun itself
                  if name == 'Sun':
                        continue

                  # Calculate relative position and velocity with respect to the Sun
                  sun_position = np.array([positions['Sun']['x'][i], positions['Sun']['y'][i], positions['Sun']['z'][i]])
                  sun_velocity = np.array([velocities['Sun']['x'][i], velocities['Sun']['y'][i], velocities['Sun']['z'][i]])
                  relative_position = np.array([positions[name]['x'][i], positions[name]['y'][i], positions[name]['z'][i]]) - sun_position
                  relative_velocity = np.array([velocities[name]['x'][i], velocities[name]['y'][i], velocities[name]['z'][i]]) - sun_velocity

                  # Calculate radial velocity relative to the Sun. Note the magnitude will be wrong as this is not normalised by radius but this will not change the sign which is what we are after!
                  radial_velocity = (
                        relative_position[0] * relative_velocity[0] +
                        relative_position[1] * relative_velocity[1] +
                        relative_position[2] * relative_velocity[2]
                        )


                  # Check for radial velocity sign change
                  # I have tried to think of better ways of doing the is not None. None I can think off do not introduce more errors/save processing

                  if previous_radial_velocity[name] is not None and radial_velocity * previous_radial_velocity[name] < 0:
                        radial_velocity_sign_changes[name].append((time[i]))
                  
                  previous_radial_velocity[name] = radial_velocity

                  # Check for y-axis crossing
                  if positions[name]['y'][i] * positions[name]['y'][i+1] < 0:
                        y_crossings[name].append(time[i])
            
            time[i+1] = time[i] + time_step

      # Only return used values
      return {name: {key: val[:i+1] for key, val in positions[name].items()} for name in positions}, time[:i+1], y_crossings,radial_velocity_sign_changes


def plot_orbit(positions, formatting, img_name, resolution=1000,just_3d = False):
      
      # 3D plot
      if just_3d:
            fig = plt.figure(figsize=(300,90))
            ax1 = fig.add_subplot(111,projection = '3d')
      else:
            fig = plt.figure(figsize=(20, 10))
            ax1 = fig.add_subplot(121, projection='3d')

      
      for name, pos in positions.items():
            # Get 1000 evenly spaced indices
            if len(pos['x']) > resolution:
                  indices = np.linspace(0, len(pos['x']) - 1, resolution, dtype=int)
                  pos_slice = {
                  'x': np.array(pos['x'])[indices],
                  'y': np.array(pos['y'])[indices],
                  'z': np.array(pos['z'])[indices]
                  }
            else:
                  pos_slice = pos  # If there are fewer than 1000 points, use all available points
            
            fmt = formatting.get(name, {'linewidth': 0.5, 'alpha': 0.75, 'color': 'black'})  # Default formatting
            ax1.plot(pos_slice['x'], pos_slice['y'], pos_slice['z'], label=name, **fmt)
            # Add marker at the last point
            if name == 'Sun':
                  ax1.scatter(pos['x'][-1], pos['y'][-1], pos['z'][-1], color='yellow', s=100, edgecolor='k', zorder=5 ,alpha=1)

            else :
                  ax1.scatter(pos['x'][-1], pos['y'][-1], pos['z'][-1], color=fmt['color'], s=25, edgecolor='k', zorder=5 ,alpha=0.8)
      
      ax1.set_title("Orbit Simulation (3D)")
      ax1.legend(loc='upper right')
      ax1.view_init(elev=10, azim=18)
      #Ensure the z doesnt look too weird
      ax1.set_aspect('equal')
      cool_photo = True
      if cool_photo:
            ax1.grid(False)
            ax1.set_xticks([])  # Remove x ticks
            ax1.set_yticks([])  # Remove y ticks
            ax1.set_zticks([])  # Remove z ticks
            ax1.xaxis.pane.fill = False
            ax1.yaxis.pane.fill = False
            ax1.zaxis.pane.fill = False
            ax1.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))  # Hide x axis line
            ax1.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))  # Hide y axis line
            ax1.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))  # Hide z axis line
      else:
            ax1.set_xlabel('x (m)')
            ax1.set_ylabel('y (m)')
            ax1.set_zlabel('z (m)')

      if just_3d:
            plt.savefig(f'{img_name}.png') 
            plt.show()
            return
      
      # 2D plot
      ax2 = fig.add_subplot(122)
      for name, pos in positions.items():
            if len(pos['x']) > resolution:
                  indices = np.linspace(0, len(pos['x']) - 1, resolution, dtype=int)
                  pos_slice = {
                  'x': np.array(pos['x'])[indices],
                  'y': np.array(pos['y'])[indices]
                  }
            else:
                  pos_slice = pos  
            
            fmt = formatting.get(name, {'linewidth': 0.5, 'alpha': 0.75, 'color': 'black'})  # Default formatting
            ax2.plot(pos_slice['x'], pos_slice['y'], label=name, **fmt)
            # Add marker at the last point
            ax2.scatter(pos['x'][-1], pos['y'][-1], color=fmt['color'], s=50, edgecolor='k', zorder=5)
      
      ax2.set_title("Orbit Simulation (2D)")
      ax2.set_xlabel('x (m)')
      ax2.set_ylabel('y (m)')
      ax2.set_aspect('equal', 'box') 

      ax2.legend(loc='upper right')
      
      plt.show()
      plt.savefig(f'{img_name}.png') 

def calculate_orbital_periods(y_crossings):
    for name, crossings in y_crossings.items():
        if len(crossings) > 1:
            # Get the half-periods by taking the difference between consecutive crossings.
            half_periods = [crossings[i+1] - crossings[i] for i in range(len(crossings) - 1)]

            # Convert half-periods from seconds to years.
            half_periods_in_years = [period / (365.25 * 24 * 3600) for period in half_periods]
            
            average_half_period = np.mean(half_periods_in_years)
            std_dev_half_period = np.std(half_periods_in_years)
            
            # Half to full.
            average_period = average_half_period * 2
            std_dev_period = std_dev_half_period * 2
            
            print(f"{name} - Average Period: {average_period:.2f} years, Standard Deviation: {std_dev_period:.2f} years")
        else:
            print(f"{name} - Not enough crossings to calculate period")

def calculate_cloest_and_furthest_approach(positions):
      for name, pos in positions.items():
            distances = np.sqrt(pos['x']**2 + pos['y']**2 + pos['z']**2)
            closest_approach = np.min(distances)
            furthest_approach = np.max(distances)

            #Convert to AU
            closest_approach = closest_approach / 1.496e11
            furthest_approach = furthest_approach / 1.496e11

            print(f"{name} - Closest Approach: {closest_approach:.2e} AU, Furthest Approach: {furthest_approach:.2e} AU")


#TODO I should probably move this to a seperate file or a class but cba rn.


def api_call(example,time_step,total_steps,plot_count,three_d_graph_only):


      masses, objects_info, formatting = get_masses_and_object_info(example)

      au = 1.496e11  # 1 Astronomical Unit in meters.
      day_in_seconds = 24 * 3600  # Number of seconds in a day.



      # Check if the CSV file exists
      csv_file = 'objects_data.csv'
      if os.path.exists(csv_file):
            # Load data from CSV
            print('CSV file already exists, loading data from CSV. NOTE if you have changed the example use the --refresh flag to rewrite the csv')
            objects = pd.read_csv(csv_file).to_dict(orient='records')
      else:
            # If the CSV does not exist, fetch data from the api
            objects = []
            for obj in objects_info:
                  horizons_obj = Horizons(
                  id=obj['id'],
                  location='500@0',  # Centered on Sun
                  epochs={'start': '2024-04-01', 'stop': '2025-04-02', 'step': '1d'}
                  )
                  eph = horizons_obj.vectors() 

                  print(masses[obj['name']])
                  obj_data = {
                  'name': obj['name'],
                  'mass': masses[obj['name']], 
                  'initial_x': eph['x'][0] * au,
                  'initial_y': eph['y'][0] * au,
                  'initial_z': eph['z'][0] * au,

                  'initial_x_velocity': eph['vx'][0] * (au / day_in_seconds),
                  'initial_y_velocity': eph['vy'][0] * (au / day_in_seconds),
                  'initial_z_velocity': eph['vz'][0] * (au / day_in_seconds)
                  }
                  objects.append(obj_data)
            
            # Save data to CSV as I dont wanna make more API calls than I need to 
            pd.DataFrame(objects).to_csv(csv_file, index=False)

      positions, time, y_crossings , radial_velocity_changes = calculate_orbit(objects,time_step,total_steps)
      plot_orbit(positions, formatting, 'plot_actual_data',plot_count,three_d_graph_only)

      print('Periods using y crossings')
      calculate_orbital_periods(y_crossings) 
      print()
      print('Periods using radiual velcoity changes')
      calculate_orbital_periods(radial_velocity_changes)
      print()
      print('cloests approaches')
      calculate_cloest_and_furthest_approach(positions)
     

if __name__ == "__main__":
      # Check if anything is passsed in.
      parser = argparse.ArgumentParser(
            description="Simulate solar system"
      )
      parser.add_argument("--refresh", action="store_true", help="Delete current CSV ")
      parser.add_argument("--single_graph", action="store_true", help="Delete current CSV ")

      parser.add_argument(
            "--example",
            type=str,
            default="solar_system",
            help="Show an exmaple system, currently implemented 'solar_system' and 'earth_and_moon'",
      )

      parser.add_argument(
            "--time_step",
            type=int,
            default=3600,
            help="Specify in seconds how often the time interval is to be taken"
      )

      parser.add_argument(
            "--time_total",
            type=int,
            default=3,
            help="specify in earth years how much time to take"                    
      )
      parser.add_argument(
            "--plot_count",
            type=int,
            default=1000,
            help="How many points should be plotted in the diagram"
      )


      args = parser.parse_args()

      # convert years to max time interval.
      # Takes how many seconds in a year and divies by the time step to get how many time steps in the given year. the plus one is purely for the tqdm loading bar to look nicer.
      total_steps = int(args.time_total*3600*24*365/args.time_step) + 1

      if args.refresh:
            # Delete existing CSV if it exists
            if os.path.exists('objects_data.csv'):
                  os.remove('objects_data.csv')

      api_call(args.example,args.time_step,total_steps,args.plot_count,args.single_graph)
 
