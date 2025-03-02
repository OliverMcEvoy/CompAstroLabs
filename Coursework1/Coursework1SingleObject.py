import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def calculate_orbit(
    initial_radius,
    initial_velocity,
    initial_z_velocity,
    max_orbit_count,
    initial_time_step=1000,
    min_time_step=2000,
    max_steps=10000000,
    tolerance=1e-6,
):
    # Constants
    G = 6.67430e-11
    star_mass = 1.989e30
    comet_mass = 2.2e14

    # Initialise arrays
    x = np.zeros(max_steps)
    y = np.zeros(max_steps)
    z = np.zeros(max_steps)
    x_velocity = np.zeros(max_steps)
    y_velocity = np.zeros(max_steps)
    time = np.zeros(max_steps)

    # Initialise Sun's z-position and z-velocity
    sun_z = np.zeros(max_steps)
    sun_z_velocity = np.zeros(max_steps)

    # Set initial conditions
    x[0] = initial_radius
    y[0] = 0
    z[0] = 0
    x_velocity[0] = 0
    y_velocity[0] = initial_velocity
    z_velocity = initial_z_velocity
    time[0] = 0

    # Orbit completion parameters
    orbit_completion_tolerance = initial_radius * 0.9999
    min_points = 10000
    time_step = initial_time_step
    orbit_completions = 0

    for i in range(0, max_steps - 1):
        radius = np.sqrt(x[i] ** 2 + y[i] ** 2)

        if i > min_points and x[i] > orbit_completion_tolerance:
            orbit_completions += 1
            min_points += i

            if orbit_completions == max_orbit_count:
                return (
                    x[: i + 1],
                    y[: i + 1],
                    z[: i + 1],
                    sun_z[: i + 1],
                    time[: i + 1],
                    orbit_completions,
                )

        # Calculate acceleration
        x_acceleration = -G * star_mass * x[i] / radius**3
        y_acceleration = -G * star_mass * y[i] / radius**3

        # Trial step
        x_velocity[i + 1] = x_velocity[i] + x_acceleration * time_step
        y_velocity[i + 1] = y_velocity[i] + y_acceleration * time_step
        x[i + 1] = x[i] + x_velocity[i + 1] * time_step
        y[i + 1] = y[i] + y_velocity[i + 1] * time_step
        z[i + 1] = z[i] + z_velocity * time_step
        sun_z[i + 1] = sun_z[i] + z_velocity * time_step
        time[i + 1] = time[i] + time_step

        i += 1

    # Only return used values
    return (
        x[: i + 1],
        y[: i + 1],
        z[: i + 1],
        sun_z[: i + 1],
        time[: i + 1],
        orbit_completions,
    )


def plot_orbit(x, y, z, sun_z):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(x, y, z, label="Comet")
    ax.plot([0], [0], sun_z[-1], "yo", label="Sun")
    ax.plot([0], [0], sun_z, "r--", label="Sun z-position")
    ax.set_title("Halley's Comet Orbit")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    ax.legend()
    plt.savefig("orbit_plot_3d.png")
    plt.show()


def main():
    # Initial conditions
    initial_radius = 5.2e12
    initial_velocity = 880
    initial_z_velocity = 5  # Initial z-velocity
    max_orbit_count = 10  # Number of orbits to simulate

    # Run simulation
    x, y, z, sun_z, time, orbit_completions = calculate_orbit(
        initial_radius, initial_velocity, initial_z_velocity, max_orbit_count
    )
    plot_orbit(x, y, z, sun_z)

    # Calculate orbital period
    period = time[-1] / (
        365.25 * 24 * 3600 * orbit_completions
    )  # Convert to years and divide by number of orbits
    print(f"Calculated orbital period: {period:.2f} years")


if __name__ == "__main__":
    main()
