# List of examples to choose from:
# 'solar_system' - The Sun, Moon, Halley's Comet, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune
# 'earth_and_moon' - The Sun, Earth, Moon
# 'rocky_planets' - The Sun, Mercury, Venus, Earth, Moon, Mars


def get_masses_and_object_info(example="solar_system"):
    # Define the objects and their IDs in the JPL Horizons system along with the mass seperately as thats not stored for some reason... unless it is and I cannot find it.
    masses = {
        "Sun": 1.989e30,
        "Halleys_Comet": 2.2e14,
        "Mercury": 3.301e23,
        "Moon": 7.34767309e22,
        "Venus": 4.867e24,
        "Earth": 5.972e24,
        "Mars": 6.417e23,
        "Jupiter": 1.898e27,
        "Saturn": 5.683e26,
        "Uranus": 8.681e25,
        "Neptune": 1.024e26,
    }
    # Custom formatting no harm in setting this as a default and if I want a seperate one for an example then I can just override it.
    formatting = {
        "Sun": {"linewidth": 1, "alpha": 1.0, "color": "yellow"},
        "Halleys_Comet": {"linewidth": 1.5, "alpha": 0.75, "color": "purple"},
        "Mercury": {"linewidth": 1, "alpha": 0.75, "color": "gray"},
        "Venus": {"linewidth": 1, "alpha": 0.75, "color": "orange"},
        "Earth": {"linewidth": 1, "alpha": 0.75, "color": "green"},
        "Moon": {"linewidth": 0.5, "alpha": 0.75, "color": "black"},
        "Mars": {"linewidth": 1, "alpha": 0.75, "color": "red"},
        "Jupiter": {"linewidth": 1.5, "alpha": 0.75, "color": "orange"},
        "Saturn": {"linewidth": 1.5, "alpha": 0.75, "color": "gold"},
        "Uranus": {"linewidth": 1.5, "alpha": 0.75, "color": "lightblue"},
        "Neptune": {"linewidth": 1.5, "alpha": 0.75, "color": "blue"},
    }

    # As Python does not have switch statements this will do for now.
    # TODO Make this a disctionary or something this code makes me ashamed.
    if example == "solar_system":
        objects_info = [
            {"name": "Mercury", "id": "199"},
            {"name": "Sun", "id": "10"},
            {"name": "Moon", "id": "301"},
            {"name": "Halleys_Comet", "id": "90000001"},
            {"name": "Venus", "id": "299"},
            {"name": "Earth", "id": "399"},
            {"name": "Mars", "id": "499"},
            {"name": "Jupiter", "id": "599"},
            {"name": "Saturn", "id": "699"},
            {"name": "Uranus", "id": "799"},
            {"name": "Neptune", "id": "899"},
        ]

    elif example == "earth_and_moon":
        objects_info = [
            {"name": "Sun", "id": "10"},
            {"name": "Earth", "id": "399"},
            {"name": "Moon", "id": "301"},
        ]

        formatting = {
            "Sun": {"linewidth": 4, "alpha": 1.0, "color": "yellow"},
            "Moon": {
                "linewidth": 0.5,
                "alpha": 0.5,
                "color": "black",
                "linestyle": "-",
            },
            "Earth": {"linewidth": 2, "alpha": 0.5, "color": "green", "linestyle": ":"},
        }
    elif example == "rocky_planets":
        objects_info = [
            {"name": "Sun", "id": "10"},
            {"name": "Mercury", "id": "199"},
            {"name": "Venus", "id": "299"},
            {"name": "Moon", "id": "301"},
            {"name": "Earth", "id": "399"},
            {"name": "Mars", "id": "499"},
        ]
        formatting = {
            "Sun": {"linewidth": 4, "alpha": 1.0, "color": "yellow"},
            "Mercury": {"linewidth": 1, "alpha": 1, "color": "gray"},
            "Venus": {"linewidth": 1, "alpha": 0.75, "color": "orange"},
            "Moon": {"linewidth": 0.5, "alpha": 0.8, "color": "black"},
            "Earth": {"linewidth": 1, "alpha": 0.5, "color": "green"},
            "Mars": {"linewidth": 1, "alpha": 0.75, "color": "red"},
        }
    elif example == "comet_coursework":
        objects_info = [
            {"name": "Sun", "id": "10"},
            {"name": "Halleys_Comet", "id": "90000001"},
            {"name": "Earth", "id": "399"},
            {"name": "Jupiter", "id": "599"},
            {"name": "Saturn", "id": "699"},
            {"name": "Uranus", "id": "799"},
            {"name": "Neptune", "id": "899"},
        ]
    elif example == "just_comet":
        objects_info = [
            {"name": "Sun", "id": "10"},
            {"name": "Halleys_Comet", "id": "90000001"},
        ]
        formatting = {
            "Sun": {"linewidth": 2, "alpha": 1, "color": "orange"},
            "Halleys_Comet": {"linewidth": 2, "alpha": 1, "color": "lightblue"},
        }

    else:
        # Throw exception
        print(
            "Nothing selected currently existing, Check MassAndObjectInfo.py for existing examples"
        )

    return masses, objects_info, formatting
