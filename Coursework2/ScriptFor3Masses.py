import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt


# Ill preface that this file is not the finest code I have ever written, it will do for now so I can spend more time on the cool stuff


def acceleration(name, pos, all_positions):
    """
    Calculate acceleration due to gravity from all other objects.
    """
    ax, ay = 0, 0
    for other_obj in objects:
        if other_obj["name"] == name:
            continue  # Skip self

        # Calculate change in position
        dx = all_positions[other_obj["name"]]["x"] - pos[0]
        dy = all_positions[other_obj["name"]]["y"] - pos[1]

        r = (dx**2 + dy**2) ** 0.5
        factor = G * other_obj["mass"] / (r**3)
        ax += factor * dx
        ay += factor * dy

    return np.array([ax, ay])


SOLAR_MASS = 1.988e30
G = 6.67430e-11  # Gravitational constant


def intial_velocity(sma):
    # Assuming mass of planet negligable to mass of the sun
    period_squared = (4 * np.pi**2 * sma**3) / (G * SOLAR_MASS)
    return 2 * np.pi * sma / (period_squared**0.5)


au = 1.496e11  # 1 Astronomical Unit in meters.
day_in_seconds = 24 * 3600  # Number of seconds in a day.


objects = []

objects.append(
    {
        "name": "Star",
        "mass": SOLAR_MASS,
        "initial_x": 0,
        "initial_y": 0,
        "initial_x_velocity": 0,
        "initial_y_velocity": 0,
    }
)

objects.append(
    {
        "name": "Object_1",
        "mass": 10e-3 * SOLAR_MASS,
        "initial_x": -2.52 * au,
        "initial_y": 0,
        "initial_x_velocity": 0,
        "initial_y_velocity": intial_velocity(
            2.52 * au,
        ),
    }
)

objects.append(
    {
        "name": "Object_2",
        "mass": 4e-2 * SOLAR_MASS,
        "initial_x": 5.24 * au,
        "initial_y": 0,
        "initial_x_velocity": 0,
        "initial_y_velocity": -intial_velocity(
            5.24 * au,
        ),
    }
)


# Initialise position and velocity arrays.
max_steps = 20000000 + 1
positions = {
    obj["name"]: {"x": np.zeros(max_steps), "y": np.zeros(max_steps)} for obj in objects
}
velocities = {
    obj["name"]: {"x": np.zeros(max_steps), "y": np.zeros(max_steps)} for obj in objects
}
time = np.zeros(max_steps)

# Set initial conditions
for obj in objects:
    positions[obj["name"]]["x"][0] = obj["initial_x"]
    positions[obj["name"]]["y"][0] = obj["initial_y"]
    velocities[obj["name"]]["x"][0] = obj["initial_x_velocity"]
    velocities[obj["name"]]["y"][0] = obj["initial_y_velocity"]

time[0] = 0
time_step = 3600

# Tracking radial velocity sign changes and y-crossings
previous_radial_velocity = {obj["name"]: None for obj in objects}
radial_velocity_sign_changes = {obj["name"]: [] for obj in objects}
y_crossings = {obj["name"]: [] for obj in objects}


for i in tqdm(range(0, max_steps - 1), desc="Calculating orbits", mininterval=1.0):
    all_positions = {
        obj["name"]: {
            "x": positions[obj["name"]]["x"][i],
            "y": positions[obj["name"]]["y"][i],
        }
        for obj in objects
    }

    all_velocities = {
        obj["name"]: {
            "x": velocities[obj["name"]]["x"][i],
            "y": velocities[obj["name"]]["y"][i],
        }
        for obj in objects
    }

    # RK4 integration step
    for obj in objects:
        name = obj["name"]
        pos = np.array(
            [
                positions[name]["x"][i],
                positions[name]["y"][i],
            ]
        )
        vel = np.array(
            [
                velocities[name]["x"][i],
                velocities[name]["y"][i],
            ]
        )

        # Skip calculations for the star so it remains alwasyy stationary
        if name == "Star":
            positions[name]["x"][i + 1] = positions[name]["x"][i]
            positions[name]["y"][i + 1] = positions[name]["y"][i]

            velocities[name]["x"][i + 1] = velocities[name]["x"][i]
            velocities[name]["y"][i + 1] = velocities[name]["x"][i]
            continue
        # Compute RK4 coefficients
        k1v = acceleration(name, pos, all_positions) * time_step
        k1x = vel * time_step

        k2v = acceleration(name, pos + k1x / 2, all_positions) * time_step
        k2x = (vel + k1v / 2) * time_step

        k3v = acceleration(name, pos + k2x / 2, all_positions) * time_step
        k3x = (vel + k2v / 2) * time_step

        k4v = acceleration(name, pos + k3x, all_positions) * time_step
        k4x = (vel + k3v) * time_step

        # Update position and velocity using RK4 weighted sum
        new_vel = vel + (k1v + 2 * k2v + 2 * k3v + k4v) / 6
        new_pos = pos + (k1x + 2 * k2x + 2 * k3x + k4x) / 6

        positions[name]["x"][i + 1] = new_pos[0]
        positions[name]["y"][i + 1] = new_pos[1]

        velocities[name]["x"][i + 1] = new_vel[0]
        velocities[name]["y"][i + 1] = new_vel[1]

        # Calculate relative position and velocity with respect to the Sun
        sun_pos = np.array([positions["Star"]["x"][i], positions["Star"]["y"][i]])
        sun_vel = np.array([velocities["Star"]["x"][i], velocities["Star"]["y"][i]])
        relative_position = pos - sun_pos
        relative_velocity = vel - sun_vel

        # Calculate radial velocity relative to the Sun
        radial_velocity = (
            relative_position[0] * relative_velocity[0]
            + relative_position[1] * relative_velocity[1]
        )

        # Track radial velocity sign changes
        if (
            previous_radial_velocity[name] is not None
            and radial_velocity * previous_radial_velocity[name] < 0
        ):
            radial_velocity_sign_changes[name].append(time[i])

        previous_radial_velocity[name] = radial_velocity

        # Check for y-axis crossings
        if positions[name]["y"][i] * positions[name]["y"][i + 1] < 0:
            y_crossings[name].append(time[i])

    time[i + 1] = time[i] + time_step

# Plot a graph of the results.

plt.figure(figsize=(10, 6))

# Plotting orbits or trajectories
for obj in objects:
    name = obj["name"]
    # This method will have a couple of thousand points just for the sun, oh well the compute time is still absolutely minimal compare to the bigger programmes so I will spend more time on them.
    # This code is just to get a graph produced for the report,
    # I'll use the more important and significant pieces of code as demonstration I can code?
    # Hopefully this isnt the first piece of code looked at that would be unlucky.
    x = positions[name]["x"]
    y = positions[name]["y"]

    plt.plot(x, y, label=f"{name}", linewidth=1)
    plt.scatter(x[-1], y[-1], s=20)
    plt.xlabel("X Position")
    plt.ylabel("Y Position")
    plt.title(f"Starting on opposite sides, both clockwise, long time")

plt.legend(loc="upper right")
# The forward slash here might break on windows as they have a file system like \ Theres a libary I could use and implement but I would not be able to test as I have got rid of Windows off every device I own
plt.savefig("2a/clockwise_opposite_long_time.png")
plt.show()
