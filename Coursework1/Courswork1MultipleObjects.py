import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm

def calculate_orbit(objects, initial_time_step=100, min_time_step=2000, max_steps=10000001, tolerance=1e-6):
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
      # Defining constants
      G = 6.67430e-11  

      # Initialise arrays for each object
      positions = {obj['name']: {'x': np.zeros(max_steps), 'y': np.zeros(max_steps), 'z': np.zeros(max_steps)} for obj in objects}
      velocities = {obj['name']: {'x': np.zeros(max_steps), 'y': np.zeros(max_steps), 'z': np.zeros(max_steps)} for obj in objects}
      time = np.zeros(max_steps)
      
      # Object initial conditions
      for obj in objects:
            positions[obj['name']]['x'][0] = obj['initial_x']
            positions[obj['name']]['y'][0] = obj['initial_y']
            positions[obj['name']]['z'][0] = obj['initial_z']
            velocities[obj['name']]['x'][0] = 0
            velocities[obj['name']]['y'][0] = obj['initial_velocity']
            velocities[obj['name']]['z'][0] = obj['initial_z_velocity']
      
      time[0] = 0

      time_step = initial_time_step
      y_crossings = {obj['name']: [] for obj in objects}
      
      for i in tqdm(range(0, max_steps-1), desc="Calculating orbits"):
            for obj in objects:
                  name = obj['name']
                  #radius = np.sqrt(positions[name]['x'][i]**2 + positions[name]['y'][i]**2 + positions[name]['z'][i]**2)
                  
                  # Calculate acceleration due to all other objects
                  x_acceleration = 0
                  y_acceleration = 0
                  z_acceleration = 0

                  # A bit much indentation but it'll do for now
                  for other_obj in objects:
                        if other_obj['name'] != name:
                              other_radius = np.sqrt((positions[other_obj['name']]['x'][i] - positions[name]['x'][i])**2 + 
                                                      (positions[other_obj['name']]['y'][i] - positions[name]['y'][i])**2 + 
                                                      (positions[other_obj['name']]['z'][i] - positions[name]['z'][i])**2)
                              x_acceleration += -G * other_obj['mass'] * (positions[name]['x'][i] - positions[other_obj['name']]['x'][i]) / other_radius**3
                              y_acceleration += -G * other_obj['mass'] * (positions[name]['y'][i] - positions[other_obj['name']]['y'][i]) / other_radius**3
                              z_acceleration += -G * other_obj['mass'] * (positions[name]['z'][i] - positions[other_obj['name']]['z'][i]) / other_radius**3
                        
                  velocities[name]['x'][i+1] = velocities[name]['x'][i] + x_acceleration * time_step
                  velocities[name]['y'][i+1] = velocities[name]['y'][i] + y_acceleration * time_step
                  velocities[name]['z'][i+1] = velocities[name]['z'][i] + z_acceleration * time_step
                  positions[name]['x'][i+1] = positions[name]['x'][i] + velocities[name]['x'][i+1] * time_step
                  positions[name]['y'][i+1] = positions[name]['y'][i] + velocities[name]['y'][i+1] * time_step
                  positions[name]['z'][i+1] = positions[name]['z'][i] + velocities[name]['z'][i+1] * time_step
                  
                  # Check for y-axis crossing
                  if positions[name]['y'][i] * positions[name]['y'][i+1] < 0:
                        y_crossings[name].append(time[i])
            
            time[i+1] = time[i] + time_step

      # Only return used values
      return {name: {key: val[:i+1] for key, val in positions[name].items()} for name in positions}, time[:i+1], y_crossings


def plot_orbit(positions, formatting):
    fig = plt.figure(figsize=(20, 10))
    
    # 3D plot
    ax1 = fig.add_subplot(121, projection='3d')
    for name, pos in positions.items():
        fmt = formatting.get(name, {'linewidth': 0.5, 'alpha': 0.75, 'color': 'black'})  # Default formatting
        ax1.plot(pos['x'], pos['y'], pos['z'], label=name, **fmt)
    
    ax1.set_title("Orbit Simulation (3D)")
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    ax1.set_zlabel('z (m)')
    ax1.legend(loc='upper right')
    
    # 2D plot
    ax2 = fig.add_subplot(122)
    for name, pos in positions.items():
        fmt = formatting.get(name, {'linewidth': 0.5, 'alpha': 0.75, 'color': 'black'})  # Default formatting
        ax2.plot(pos['x'], pos['y'], label=name, **fmt)
    
    ax2.set_title("Orbit Simulation (2D)")
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel('y (m)')
    ax2.set_aspect('equal', 'box')  # Ensure the x and y axes have the same scale

    ax2.legend(loc='upper right')
    
    plt.savefig('orbit_plot_combined.png')  # Save the plot as an image file
    plt.show()

def main():
      # Initial conditions for objects
      objects = [
      {'name': 'Halleys_Comet', 'mass': 2.2e14, 'initial_x': 5.2e12, 'initial_y': 0, 'initial_z': 0, 'initial_velocity': 880, 'initial_z_velocity': 5},
      {'name': 'Sun', 'mass': 1.989e30, 'initial_x': 0, 'initial_y': 0, 'initial_z': 0, 'initial_velocity': 0, 'initial_z_velocity': 5},
      {'name': 'Mercury', 'mass': 3.301e23, 'initial_x': 5.791e10, 'initial_y': 0, 'initial_z': 0, 'initial_velocity': 47400, 'initial_z_velocity': 5},
      {'name': 'Venus', 'mass': 4.867e24, 'initial_x': 1.082e11, 'initial_y': 0, 'initial_z': 0, 'initial_velocity': 35020, 'initial_z_velocity': 5},
      {'name': 'Earth', 'mass': 5.972e24, 'initial_x': 1.471e11, 'initial_y': 0, 'initial_z': 0, 'initial_velocity': 29780, 'initial_z_velocity': 5},
      {'name': 'Mars', 'mass': 6.417e23, 'initial_x': 2.279e11, 'initial_y': 0, 'initial_z': 0, 'initial_velocity': 24070, 'initial_z_velocity': 5},
      {'name': 'Jupiter', 'mass': 1.898e27, 'initial_x': 7.785e11, 'initial_y': 0, 'initial_z': 0, 'initial_velocity': 13070, 'initial_z_velocity': 5},
      {'name': 'Saturn', 'mass': 5.683e26, 'initial_x': 1.433e12, 'initial_y': 0, 'initial_z': 0, 'initial_velocity': 9680, 'initial_z_velocity': 5},
      {'name': 'Uranus', 'mass': 8.681e25, 'initial_x': 2.872e12, 'initial_y': 0, 'initial_z': 0, 'initial_velocity': 6800, 'initial_z_velocity': 5},
      {'name': 'Neptune', 'mass': 1.024e26, 'initial_x': 4.495e12, 'initial_y': 0, 'initial_z': 0, 'initial_velocity': 5430, 'initial_z_velocity': 5}
      ]

      # Define custom formatting for each object
      formatting = {
            'Sun': {'linewidth': 4, 'alpha': 1.0, 'color': 'yellow'},
            'Halleys_Comet': {'linewidth': 0.5, 'alpha': 0.75, 'color': 'blue'},
            'Mercury': {'linewidth': 1, 'alpha': 0.75, 'color': 'gray'},
            'Venus': {'linewidth': 1, 'alpha': 0.75, 'color': 'orange'},
            'Earth': {'linewidth': 1, 'alpha': 0.75, 'color': 'green'},
            'Mars': {'linewidth': 1, 'alpha': 0.75, 'color': 'red'},
            'Jupiter': {'linewidth': 1.5, 'alpha': 0.75, 'color': 'orange'},
            'Saturn': {'linewidth': 1.5, 'alpha': 0.75, 'color': 'gold'},
            'Uranus': {'linewidth': 1.5, 'alpha': 0.75, 'color': 'lightblue'},
            'Neptune': {'linewidth': 1.5, 'alpha': 0.75, 'color': 'blue'}
      }

      # Run simulation
      positions, time, y_crossings = calculate_orbit(objects)
      plot_orbit(positions, formatting)

      # Calculate orbital period for each object
      for name, crossings in y_crossings.items():
            if len(crossings) > 1:
                  # Calculate half-periods by taking the difference between consecutive crossings
                  half_periods = [crossings[i+1] - crossings[i] for i in range(len(crossings) - 1)]
                  
                  # Convert half-periods from seconds to years
                  half_periods_in_years = [period / (365.25 * 24 * 3600) for period in half_periods]
                  
                  # Calculate average half-period and standard deviation
                  average_half_period = np.mean(half_periods_in_years)
                  std_dev_half_period = np.std(half_periods_in_years)
                  
                  # Adjust to get full orbital period
                  average_period = average_half_period * 2
                  std_dev_period = std_dev_half_period * 2
                  
                  # Print results
                  print(f"{name} - Average Period: {average_period:.2f} years, Standard Deviation: {std_dev_period:.2f} years")
            else:
                  print(f"{name} - Not enough crossings to calculate period")



if __name__ == "__main__":
    main()