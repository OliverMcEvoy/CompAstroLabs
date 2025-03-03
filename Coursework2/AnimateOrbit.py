import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D


def animate_orbit(
    positions,
    formatting,
    img_name,
    resolution=500,
    fps=30,
    duration=10,
    step=50,
    num_frames=1000,
):
    """
    Creates an animated 3D plot of orbital positions with aggressive optimizations.

    :param positions: Dictionary of positions for each body with keys "x", "y", "z".
    :param formatting: Dictionary of formatting options for each body.
    :param img_name: Filename to save the animation.
    :param resolution: Number of points to use for the orbit line.
    :param fps: Frames per second for animation.
    :param duration: Duration of animation in seconds.
    :param step: Step size for downsampling the data (e.g., take every 50th point).
    :param num_frames: Total number of frames in the animation.
    """
    num_frames = fps * duration
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Determine the longest dataset among the positions
    max_points = max(len(pos["x"]) for pos in positions.values())

    # Downsample the total number of frames
    frame_indices = np.linspace(0, max_points - 1, num_frames, dtype=int)

    lines = {}
    markers = {}

    # Initialise plot elements
    for name, pos in positions.items():
        fmt = formatting.get(name, {"linewidth": 1, "alpha": 0.75, "color": "black"})

        # Precompute downsampled marker positions
        marker_pos = {
            "x": np.array(pos["x"])[frame_indices],
            "y": np.array(pos["y"])[frame_indices],
            "z": np.array(pos["z"])[frame_indices],
        }

        # Initialise an empty line for the orbit
        (line,) = ax.plot([], [], [], label=name, **fmt)  # Empty orbit line
        marker = ax.scatter(
            [], [], [], color=fmt["color"], s=50, edgecolor="k"
        )  # Moving marker

        lines[name] = (line, marker_pos)  # Store line and marker positions
        markers[name] = (marker, marker_pos)

    # Auto-scaling the axis limits for the first 10000 points
    all_x = np.concatenate([pos["x"][:10000] for pos in positions.values()])
    all_y = np.concatenate([pos["y"][:10000] for pos in positions.values()])
    all_z = np.concatenate([pos["z"][:10000] for pos in positions.values()])
    ax.set_xlim(np.min(all_x), np.max(all_x) * 1.5)
    ax.set_ylim(np.min(all_y), np.max(all_y) * 1.5)
    ax.set_zlim(np.min(all_z), np.max(all_z) * 1.5)

    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

    def update(frame):
        # Update orbit lines and markers
        for name, (line, marker_pos) in lines.items():
            # Get the current marker position
            current_x = marker_pos["x"][: frame + 1]
            current_y = marker_pos["y"][: frame + 1]
            current_z = marker_pos["z"][: frame + 1]

            # Update the line to include points up to the current frame
            line.set_data(current_x, current_y)
            line.set_3d_properties(current_z)

        # Update markers
        for name, (marker, marker_pos) in markers.items():
            marker._offsets3d = (
                np.array([marker_pos["x"][frame]]),
                np.array([marker_pos["y"][frame]]),
                np.array([marker_pos["z"][frame]]),
            )

        ax.view_init(elev=10, azim=30)
        return list(lines.values()) + list(markers.values())

    # Create animation
    ani = animation.FuncAnimation(
        fig,
        update,
        frames=num_frames,
        interval=1000 / fps,  # Interval in milliseconds
        blit=False,
    )

    # Save animation
    ani.save(f"{img_name}.mp4", writer="ffmpeg", fps=fps)
    plt.show()
