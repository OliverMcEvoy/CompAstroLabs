import numpy as np
from tqdm import tqdm


def calculate_orbit(objects, initial_time_step=3600, max_steps=8766):
    """
    Calculate the orbit of multiple objects in 3D space using the Velocity Verlet method.
    """
    G = 6.67430e-11  # Gravitational constant
    # Initialise arrays
    positions = {
        obj["name"]: {
            "x": np.zeros(max_steps),
            "y": np.zeros(max_steps),
            "z": np.zeros(max_steps),
        }
        for obj in objects
    }
    velocities = {
        obj["name"]: {
            "x": np.zeros(max_steps),
            "y": np.zeros(max_steps),
            "z": np.zeros(max_steps),
        }
        for obj in objects
    }
    accelerations = {
        obj["name"]: {
            "x": np.zeros(max_steps),
            "y": np.zeros(max_steps),
            "z": np.zeros(max_steps),
        }
        for obj in objects
    }
    time = np.zeros(max_steps)

    # Initialise positions and velocities
    for obj in objects:
        name = obj["name"]
        positions[name]["x"][0] = obj["initial_x"]
        positions[name]["y"][0] = obj["initial_y"]
        positions[name]["z"][0] = obj["initial_z"]
        velocities[name]["x"][0] = obj["initial_x_velocity"]
        velocities[name]["y"][0] = obj["initial_y_velocity"]
        velocities[name]["z"][0] = obj["initial_z_velocity"]

    time[0] = 0
    time_step = initial_time_step

    def compute_acceleration(name, i):
        ax, ay, az = 0, 0, 0
        for other in objects:
            if other["name"] != name:
                dx = positions[other["name"]]["x"][i] - positions[name]["x"][i]
                dy = positions[other["name"]]["y"][i] - positions[name]["y"][i]
                dz = positions[other["name"]]["z"][i] - positions[name]["z"][i]
                radius_squared = dx**2 + dy**2 + dz**2
                factor = G * other["mass"] / (radius_squared**1.5)
                ax += factor * dx
                ay += factor * dy
                az += factor * dz
        return ax, ay, az

    # Compute initial accelerations
    for obj in objects:
        name = obj["name"]
        ax, ay, az = compute_acceleration(name, 0)
        accelerations[name]["x"][0] = ax
        accelerations[name]["y"][0] = ay
        accelerations[name]["z"][0] = az

    previous_radial_velocity = {obj["name"]: None for obj in objects}
    radial_velocity_sign_changes = {obj["name"]: [] for obj in objects}
    y_crossings = {obj["name"]: [] for obj in objects}

    # TODO pass in bootstap iteration number
    for i in tqdm(range(max_steps - 1), desc="Calculating orbits", mininterval=1.0):
        time[i + 1] = time[i] + time_step

        # First update positions using current acceleration
        for obj in objects:
            name = obj["name"]
            positions[name]["x"][i + 1] = (
                positions[name]["x"][i]
                + velocities[name]["x"][i] * time_step
                + 0.5 * accelerations[name]["x"][i] * time_step**2
            )
            positions[name]["y"][i + 1] = (
                positions[name]["y"][i]
                + velocities[name]["y"][i] * time_step
                + 0.5 * accelerations[name]["y"][i] * time_step**2
            )
            positions[name]["z"][i + 1] = (
                positions[name]["z"][i]
                + velocities[name]["z"][i] * time_step
                + 0.5 * accelerations[name]["z"][i] * time_step**2
            )

        # Then compute new accelerations for all objects
        new_accelerations = {}
        for obj in objects:
            name = obj["name"]
            new_accelerations[name] = compute_acceleration(name, i + 1)

        # Now update velocities
        for obj in objects:
            name = obj["name"]
            ax_new, ay_new, az_new = new_accelerations[name]
            velocities[name]["x"][i + 1] = (
                velocities[name]["x"][i]
                + 0.5 * (accelerations[name]["x"][i] + ax_new) * time_step
            )
            velocities[name]["y"][i + 1] = (
                velocities[name]["y"][i]
                + 0.5 * (accelerations[name]["y"][i] + ay_new) * time_step
            )
            velocities[name]["z"][i + 1] = (
                velocities[name]["z"][i]
                + 0.5 * (accelerations[name]["z"][i] + az_new) * time_step
            )

            # Store new acceleration
            accelerations[name]["x"][i + 1] = ax_new
            accelerations[name]["y"][i + 1] = ay_new
            accelerations[name]["z"][i + 1] = az_new

            # Skip the Sun itself for calculating the radial velocity
            if name == "Sun":
                continue

            # Calculate relative position and velocity with respect to the Sun
            sun_position = np.array(
                [
                    positions["Sun"]["x"][i],
                    positions["Sun"]["y"][i],
                    positions["Sun"]["z"][i],
                ]
            )
            sun_velocity = np.array(
                [
                    velocities["Sun"]["x"][i],
                    velocities["Sun"]["y"][i],
                    velocities["Sun"]["z"][i],
                ]
            )
            relative_position = (
                np.array(
                    [
                        positions[name]["x"][i],
                        positions[name]["y"][i],
                        positions[name]["z"][i],
                    ]
                )
                - sun_position
            )
            relative_velocity = (
                np.array(
                    [
                        velocities[name]["x"][i],
                        velocities[name]["y"][i],
                        velocities[name]["z"][i],
                    ]
                )
                - sun_velocity
            )

            # Calculate radial velocity relative to the Sun. Note the magnitude will be wrong as this is not normalised by radius but this will not change the sign which is what we are after!
            radial_velocity = (
                relative_position[0] * relative_velocity[0]
                + relative_position[1] * relative_velocity[1]
                + relative_position[2] * relative_velocity[2]
            )

            # Check for radial velocity sign change
            # I have tried to think of better ways of doing the is not None. None I can think off do not introduce more errors/save processing

            if (
                previous_radial_velocity[name] is not None
                and radial_velocity * previous_radial_velocity[name] < 0
            ):
                radial_velocity_sign_changes[name].append((time[i]))

            previous_radial_velocity[name] = radial_velocity

            # Check for y-axis crossing
            if positions[name]["y"][i] * positions[name]["y"][i + 1] < 0:
                y_crossings[name].append(time[i])

    # Only return used values
    return (
        {
            name: {key: val[: i + 1] for key, val in positions[name].items()}
            for name in positions
        },
        time[: i + 1],
        y_crossings,
        radial_velocity_sign_changes,
    )
