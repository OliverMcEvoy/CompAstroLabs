import numpy as np
from tqdm import tqdm


def calculate_orbit_rk4(objects, initial_time_step=3600, max_steps=8766):
    g = 6.67430e-11  # Gravitational constant
    """
      Simulates the orbits of multiple objects in 3D space using the Runge-Kutta 4th order method.
      """

    def acceleration(name, pos, all_positions):
        """
        Calculate acceleration due to gravity from all other objects.
        Its minor but the scope of this function will only encompass the calcuate orbit function so I see no reason to define this outside the claculate orbit function.
        Defining this in the calculate orbit function also allows the refrencing of objects without explicitly passing it in which is a minor efficency save
        With how much this function is ran it could be done in C and inported but alas thats beyond the scope of my plans
        It might seem off to have this defined in this function but it was done with intent. Typically I would not do this.
        """
        ax, ay, az = 0, 0, 0
        for other_obj in objects:
            if other_obj["name"] == name:
                continue  # Skip self

            # Calculate change in position
            dx = all_positions[other_obj["name"]]["x"] - pos[0]
            dy = all_positions[other_obj["name"]]["y"] - pos[1]
            dz = all_positions[other_obj["name"]]["z"] - pos[2]

            r = (dx**2 + dy**2 + dz**2) ** 0.5
            factor = g * other_obj["mass"] / (r**3)
            ax += factor * dx
            ay += factor * dy
            az += factor * dz

        return np.array([ax, ay, az])

    # Initialise position and velocity arrays
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
    time = np.zeros(max_steps)

    # Set initial conditions
    for obj in objects:
        positions[obj["name"]]["x"][0] = obj["initial_x"]
        positions[obj["name"]]["y"][0] = obj["initial_y"]
        positions[obj["name"]]["z"][0] = obj["initial_z"]
        velocities[obj["name"]]["x"][0] = obj["initial_x_velocity"]
        velocities[obj["name"]]["y"][0] = obj["initial_y_velocity"]
        velocities[obj["name"]]["z"][0] = obj["initial_z_velocity"]

    time[0] = 0
    time_step = initial_time_step

    # Tracking radial velocity sign changes and y-crossings
    previous_radial_velocity = {obj["name"]: None for obj in objects}
    radial_velocity_sign_changes = {obj["name"]: [] for obj in objects}
    y_crossings = {obj["name"]: [] for obj in objects}

    for i in tqdm(range(0, max_steps - 1), desc="Calculating orbits", mininterval=1.0):
        all_positions = {
            obj["name"]: {
                "x": positions[obj["name"]]["x"][i],
                "y": positions[obj["name"]]["y"][i],
                "z": positions[obj["name"]]["z"][i],
            }
            for obj in objects
        }

        all_velocities = {
            obj["name"]: {
                "x": velocities[obj["name"]]["x"][i],
                "y": velocities[obj["name"]]["y"][i],
                "z": velocities[obj["name"]]["z"][i],
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
                    positions[name]["z"][i],
                ]
            )
            vel = np.array(
                [
                    velocities[name]["x"][i],
                    velocities[name]["y"][i],
                    velocities[name]["z"][i],
                ]
            )

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
            positions[name]["z"][i + 1] = new_pos[2]

            velocities[name]["x"][i + 1] = new_vel[0]
            velocities[name]["y"][i + 1] = new_vel[1]
            velocities[name]["z"][i + 1] = new_vel[2]

            # Skip further calculations for the Sun
            if name == "Sun":
                continue

            # Calculate relative position and velocity with respect to the Sun
            sun_pos = np.array(
                [
                    positions["Sun"]["x"][i],
                    positions["Sun"]["y"][i],
                    positions["Sun"]["z"][i],
                ]
            )
            sun_vel = np.array(
                [
                    velocities["Sun"]["x"][i],
                    velocities["Sun"]["y"][i],
                    velocities["Sun"]["z"][i],
                ]
            )
            relative_position = pos - sun_pos
            relative_velocity = vel - sun_vel

            # Calculate radial velocity relative to the Sun
            radial_velocity = (
                relative_position[0] * relative_velocity[0]
                + relative_position[1] * relative_velocity[1]
                + relative_position[2] * relative_velocity[2]
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

    # Return only used values
    return (
        {
            name: {key: val[: i + 1] for key, val in positions[name].items()}
            for name in positions
        },
        time[: i + 1],
        y_crossings,
        radial_velocity_sign_changes,
    )
