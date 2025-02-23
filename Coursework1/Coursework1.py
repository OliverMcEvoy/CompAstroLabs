import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11  
M = 1.989e30     
m = 2.2e14       

# Initial conditions
initial_radius = 5.2e12      
initial_velocity = 880       

def calculate_orbit(initial_time_step=86400,min_time_step = 2000, max_steps=100000, tolerance=1e-6):
      # Initialise arrays
      x = np.zeros(max_steps)
      y = np.zeros(max_steps)
      x_velocity = np.zeros(max_steps)
      y_velocity = np.zeros(max_steps)
      time = np.zeros(max_steps)
      
      # Set initial conditions
      x[0] = initial_radius
      y[0] = 0
      x_velocity[0] = 0
      y_velocity[0] = initial_velocity
      time[0] = 0
      
      # Orbit completion parameters
      orbit_completion_tolerance = initial_radius * 0.9999
      min_points = 10000
      time_step = initial_time_step
      
      for i in range(0,max_steps-1):
            radius = np.sqrt(x[i]**2 + y[i]**2)
                  
            if i > min_points and x[i] > orbit_completion_tolerance:
                        return x[:i+1], y[:i+1], time[:i+1]
            
            # Calculate acceleration
            x_acceleration = -G * M * x[i] / radius**3
            y_acceleration = -G * M * y[i] / radius**3
            
            # Trial step
            x_velocity[i+1] = x_velocity[i] + x_acceleration * time_step
            y_velocity [i+1 ] = y_velocity[i] + y_acceleration * time_step
            x[i+1] = x[i] + x_velocity[i+1] * time_step
            y[i+1] = y[i] + y_velocity[i+1] * time_step
            time[i+1] = time[i] + time_step

            i += 1

      # Only return used values
      return x[:i+1], y[:i+1], time[:i+1]

def plot_orbit(x, y):
    plt.figure(figsize=(10, 10))
    plt.plot(x, y)
    plt.plot(0, 0, 'yo', label='Sun')
    plt.axis('equal')
    plt.title("Halley's Comet Orbit")
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.legend()
    plt.show()

# Run simulation
x, y, time = calculate_orbit()
plot_orbit(x, y)

# Calculate orbital period
period = time[-1] / (365.25 * 24 * 3600)  # Convert to years
print(f"Calculated orbital period: {period:.2f} years")