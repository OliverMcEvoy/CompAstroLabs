import numpy as np 
from tqdm import tqdm 

def calculate_orbit_com(objects, initial_time_step=3600, max_steps=8766):
      """
      Calculate the orbit of multiple objects in 3D space using the Velocity Verlet method,
      using the dynamically shifting center of mass as the reference frame.
      """
      G = 6.67430e-11  # Gravitational constant
      epsilon = 1
      
      # Initialize arrays
      positions, velocities, accelerations = {}, {}, {}
      for obj in objects:
            name = obj['name']
            positions[name] = {'x': np.zeros(max_steps), 'y': np.zeros(max_steps), 'z': np.zeros(max_steps)}
            velocities[name] = {'x': np.zeros(max_steps), 'y': np.zeros(max_steps), 'z': np.zeros(max_steps)}
            accelerations[name] = {'x': np.zeros(max_steps), 'y': np.zeros(max_steps), 'z': np.zeros(max_steps)}
      
      time = np.zeros(max_steps)
      
      def compute_center_of_mass(i):
            total_mass = sum(obj['mass'] for obj in objects)
            com_position = np.zeros(3)
            com_velocity = np.zeros(3)
            for obj in objects:
                  name = obj['name']
                  mass = obj['mass']
                  com_position += mass * np.array([positions[name]['x'][i], positions[name]['y'][i], positions[name]['z'][i]])
                  com_velocity += mass * np.array([velocities[name]['x'][i], velocities[name]['y'][i], velocities[name]['z'][i]])
            return com_position/total_mass, com_velocity/total_mass
      
      # Initialize positions and velocities
      for obj in objects:
            name = obj['name']
            positions[name]['x'][0] = obj['initial_x']
            positions[name]['y'][0] = obj['initial_y']
            positions[name]['z'][0] = obj['initial_z']
            velocities[name]['x'][0] = obj['initial_x_velocity']
            velocities[name]['y'][0] = obj['initial_y_velocity']
            velocities[name]['z'][0] = obj['initial_z_velocity']
      
      # Compute initial center of mass
      com_position, com_velocity = compute_center_of_mass(0)
      
      # Adjust positions and velocities to COM frame
      for obj in objects:
            name = obj['name']
            positions[name]['x'][0] -= com_position[0]
            positions[name]['y'][0] -= com_position[1]
            positions[name]['z'][0] -= com_position[2]
            velocities[name]['x'][0] -= com_velocity[0]
            velocities[name]['y'][0] -= com_velocity[1]
            velocities[name]['z'][0] -= com_velocity[2]
      
      time[0] = 0
      time_step = initial_time_step
      
      def compute_acceleration(name, i):
            ax, ay, az = 0, 0, 0
            for other in objects:
                  if other['name'] != name:
                        dx = positions[other['name']]['x'][i] - positions[name]['x'][i]
                        dy = positions[other['name']]['y'][i] - positions[name]['y'][i]
                        dz = positions[other['name']]['z'][i] - positions[name]['z'][i]
                        radius_squared = dx**2 + dy**2 + dz**2
                        factor = G * other['mass'] / ((radius_squared+epsilon)**1.5)
                        ax += factor * dx
                        ay += factor * dy
                        az += factor * dz
            return ax, ay, az
      
      # Compute initial accelerations
      for obj in objects:
            name = obj['name']
            ax, ay, az = compute_acceleration(name, 0)
            accelerations[name]['x'][0] = ax
            accelerations[name]['y'][0] = ay
            accelerations[name]['z'][0] = az
      
      previous_radial_velocity = {obj['name']: None for obj in objects}
      radial_velocity_sign_changes = {obj['name']: [] for obj in objects}
      y_crossings = {obj['name']: [] for obj in objects}
      
      for i in tqdm(range(max_steps - 1), desc="Calculating orbits", mininterval=1.0):
            time[i+1] = time[i] + time_step  # Update time at the start of the step

            com_position, com_velocity = compute_center_of_mass(i)
            # First update positions using current acceleration
            for obj in objects:
                  name = obj['name']
                  positions[name]['x'][i+1] = positions[name]['x'][i] + velocities[name]['x'][i] * time_step + 0.5 * accelerations[name]['x'][i] * time_step**2 - com_position[0]
                  positions[name]['y'][i+1] = positions[name]['y'][i] + velocities[name]['y'][i] * time_step + 0.5 * accelerations[name]['y'][i] * time_step**2 - com_position[1]
                  positions[name]['z'][i+1] = positions[name]['z'][i] + velocities[name]['z'][i] * time_step + 0.5 * accelerations[name]['z'][i] * time_step**2 - com_position[2]

            # Then compute new accelerations for all objects
            new_accelerations = {}
            for obj in objects:
                  name = obj['name']
                  new_accelerations[name] = compute_acceleration(name, i+1)

            # Now update velocities
            for obj in objects:
                  name = obj['name']
                  ax_new, ay_new, az_new = new_accelerations[name]
                  velocities[name]['x'][i+1] = velocities[name]['x'][i] + 0.5 * (accelerations[name]['x'][i] + ax_new) * time_step - com_velocity[0]
                  velocities[name]['y'][i+1] = velocities[name]['y'][i] + 0.5 * (accelerations[name]['y'][i] + ay_new) * time_step - com_velocity[1]
                  velocities[name]['z'][i+1] = velocities[name]['z'][i] + 0.5 * (accelerations[name]['z'][i] + az_new) * time_step - com_velocity[2]

                  # Store new acceleration
                  accelerations[name]['x'][i+1] = ax_new
                  accelerations[name]['y'][i+1] = ay_new
                  accelerations[name]['z'][i+1] = az_new

            
      return {name: {key: val[:i+1] for key, val in positions[name].items()} for name in positions}, time[:i+1], y_crossings, radial_velocity_sign_changes
