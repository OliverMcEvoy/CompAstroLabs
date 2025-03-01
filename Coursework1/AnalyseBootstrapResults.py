import numpy as np
from datetime import datetime

# Compute statistics helper function
def compute_stats(data_dict):
      return {key: (np.mean(vals), np.std(vals)) for key, vals in data_dict.items()}

# Print results with formatting
def print_stats(stat_dict, label, unit=""):
      for key, (mean, std) in stat_dict.items():
            print(f"{key} & {mean:.8f} {unit} & ±{std:.8f} {unit}")

def analyse_bootstrap_results(bootstrap_results):
      # Initialize storage dictionaries
      named_values, halley_values = {}, {}
      approach_values = {'closest_approach': [], 'furthest_approach': []}
      time_values = {'closest_time': [], 'furthest_time': []}
      rk_values, vv_values = {}, {}
      additional_approach_values = {'closest_approach': [], 'furthest_approach': []}
      additional_time_values = {'closest_time': [], 'furthest_time': []}
      
      def update_values(data_dict, new_data):
            """Helper function to update dictionary values."""
            for key, value in new_data.items():
                  data_dict.setdefault(key, []).append(value)

      # Refrence data
      ref_date = datetime(datetime.now().year, 1, 1)

      for entry in bootstrap_results:
            # Unpack with defaults for missing values
            start_date_str, obj_dict, halley_period_list, earth_and_halley_approaches, obj_dict_if_both, halley_period_list_if_both, earth_and_halley_period_list_if_both = (
                  entry + (None,) * (7 - len(entry))
            )

            # Convert start_date to datetime and compute year difference, due to the variance in start dates
            start_date = datetime.strptime(start_date_str, "%Y-%m-%d")
            year_difference = (start_date - ref_date).days / 365.25


            # Process main dictionary values
            update_values(named_values, obj_dict)
            update_values(halley_values, halley_period_list)

            # Process approach and time values with adjusted times
            if earth_and_halley_approaches:
                  approach_values['closest_approach'].append(earth_and_halley_approaches[0])
                  approach_values['furthest_approach'].append(earth_and_halley_approaches[2])
                  time_values['closest_time'].append(earth_and_halley_approaches[1] + year_difference)
                  time_values['furthest_time'].append(earth_and_halley_approaches[3] + year_difference)

            # Process RK and VV values if they exist
            if obj_dict_if_both:
                  update_values(rk_values, obj_dict_if_both)
            if halley_period_list_if_both:
                  update_values(vv_values, halley_period_list_if_both)

            # Process additional approaches like earth_and_halley_approaches if it exists
            if earth_and_halley_period_list_if_both:
                  additional_approach_values['closest_approach'].append(earth_and_halley_period_list_if_both[0])
                  additional_approach_values['furthest_approach'].append(earth_and_halley_period_list_if_both[2])
                  additional_time_values['closest_time'].append(earth_and_halley_period_list_if_both[1] + year_difference)
                  additional_time_values['furthest_time'].append(earth_and_halley_period_list_if_both[3] + year_difference)

      # Compute statistics for each dataset
      stats = {
            "final_stats": compute_stats(named_values),
            "halley_stats": compute_stats(halley_values),
            "approach_stats": compute_stats(approach_values),
            "time_stats": compute_stats(time_values),
            "rk_stats": compute_stats(rk_values) if rk_values else {},
            "vv_stats": compute_stats(vv_values) if vv_values else {},
            "additional_approach_stats": compute_stats(additional_approach_values) if additional_approach_values["closest_approach"] else {},
            "additional_time_stats": compute_stats(additional_time_values) if additional_time_values["closest_time"] else {}
      }


      print_stats(stats["final_stats"], "Final Values")
      print_stats(stats["halley_stats"], "Halley Values")
      print_stats(stats["approach_stats"], "Approaches", "AU")
      print_stats(stats["time_stats"], "Times", "years")
      
      if stats["rk_stats"]:
            print("\nRK results above, VV below\n")
            print_stats(stats["rk_stats"], "RK Values")
      if stats["vv_stats"]:
            print_stats(stats["vv_stats"], "VV Values")
      if stats["additional_approach_stats"]:
            print_stats(stats["additional_approach_stats"], "Additional Approaches", "AU")
      if stats["additional_time_stats"]:
            print_stats(stats["additional_time_stats"], "Additional Times", "years")