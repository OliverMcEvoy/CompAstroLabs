import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd # Just used to save csv for cacheing api call.
from astroquery.jplhorizons import Horizons
import argparse
from MassAndObjectInfo import get_masses_and_object_info
from CalculateOrbitRK import calculate_orbit_rk4
from CalculateOrbitVelocityVerlet import calculate_orbit
from AnalyseBootstrapResults import analyse_bootstrap_results
from multiprocessing import Pool
import datetime
import random


def plot_orbit(positions, formatting, img_name, resolution=1000,just_3d = False):
      
      if just_3d:
            fig = plt.figure(figsize=(300,90))
            ax1 = fig.add_subplot(111,projection = '3d')
      else:
            fig = plt.figure(figsize=(20, 10))
            ax1 = fig.add_subplot(121, projection='3d')

      
      for name, pos in positions.items():
            # Get 1000 evenly spaced indices
            if len(pos['x']) > resolution:
                  indices = np.linspace(0, len(pos['x']) - 1, resolution, dtype=int)
                  pos_slice = {
                  'x': np.array(pos['x'])[indices],
                  'y': np.array(pos['y'])[indices],
                  'z': np.array(pos['z'])[indices]
                  }
            else:
                  pos_slice = pos  # If there are fewer than 1000 points, use all available points
            
            fmt = formatting.get(name, {'linewidth': 0.5, 'alpha': 0.75, 'color': 'black'})  # Default formatting if for whatever reason none is passed in.
            ax1.plot(pos_slice['x'], pos_slice['y'], pos_slice['z'], label=name, **fmt)
            # Add marker at the last point
            if name == 'Sun':
                  ax1.scatter(pos['x'][-1], pos['y'][-1], pos['z'][-1], color='yellow', s=25, edgecolor='k', zorder=5 ,alpha=1)

            else :
                  ax1.scatter(pos['x'][-1], pos['y'][-1], pos['z'][-1], color=fmt['color'], s=25, edgecolor='k', zorder=5 ,alpha=0.8)
      
      ax1.set_title("Orbit Simulation (3D)")
      ax1.view_init(elev=10, azim=18)
      #Ensure the z doesnt look too weird
      ax1.set_aspect('equal')
      if just_3d:
            ax1.grid(False)
            ax1.set_xticks([])
            ax1.set_yticks([]) 
            ax1.set_zticks([]) 
            ax1.xaxis.pane.fill = False
            ax1.yaxis.pane.fill = False
            ax1.zaxis.pane.fill = False
            #Cant figure out how to hide so just setting colour to white
            ax1.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0)) 
            ax1.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
            ax1.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))  
      else:
            ax1.legend(loc='upper right')
            ax1.set_xlabel('x (m)')
            ax1.set_ylabel('y (m)')
            ax1.set_zlabel('z (m)')

      if just_3d:
            plt.savefig(f'{img_name}.png') 
            plt.show()
            return
      
      # 2D plot
      ax2 = fig.add_subplot(122)
      for name, pos in positions.items():
            if len(pos['x']) > resolution:
                  indices = np.linspace(0, len(pos['x']) - 1, resolution, dtype=int)
                  pos_slice = {
                  'x': np.array(pos['x'])[indices],
                  'y': np.array(pos['y'])[indices]
                  }
            else:
                  pos_slice = pos  
            
            fmt = formatting.get(name, {'linewidth': 0.5, 'alpha': 0.75, 'color': 'black'})  # Default formatting
            ax2.plot(pos_slice['x'], pos_slice['y'], label=name, **fmt)
            # Add marker at the last point
            ax2.scatter(pos['x'][-1], pos['y'][-1], color=fmt['color'], s=50, edgecolor='k', zorder=5)
      
      ax2.set_title("Orbit Simulation (2D)")
      ax2.set_xlabel('x (m)')
      ax2.set_ylabel('y (m)')
      ax2.set_aspect('equal', 'box') 

      ax2.legend(loc='upper right')
      
      plt.savefig(f'{img_name}.png') 
      plt.show()

def random_past_date(days_ago):
    # Get the current date
    current_date = datetime.datetime.now()
    
    # Generate a random number of days between 0 and days_ago
    random_days = random.randint(0, days_ago)
    
    past_date = current_date - datetime.timedelta(days=random_days)
    past_date_minus_one_day = past_date - datetime.timedelta(days=1)
    
    # Format the date as "YYYY-mm-dd"
    formatted_date = past_date.strftime("%Y-%m-%d")
    formatted_date_minus_one_day = past_date_minus_one_day.strftime("%Y-%m-%d")
    return  formatted_date_minus_one_day , formatted_date

def calculate_orbital_periods(crossings, bootstrapping):
      # TODO Not sure if needed but the days in a year is not exactly 365.25
      # But will it make a differnce taking account of it as a leap year is only skipped each 100 years (but not when divisible by 400)...

      average_periods = {}
      halleys_periods = {}
            

      for name, crossings in crossings.items():
            if len(crossings) > 1:
                  # Get the half-periods by taking the difference between consecutive crossings.
                  half_periods = [crossings[i+1] - crossings[i] for i in range(len(crossings) - 1)]

                  # Convert half-periods from seconds to full periods in years.
                  period_in_years = [period * 2 / (365.25 * 24 * 3600) for period in half_periods]

                  if (name == 'Halleys_Comet'):
                        i = 1
                        for period in period_in_years:
                              halleys_periods[f"{i}th half crossing"] =period
                              if bootstrapping == False:
                                    #Im lazy and this alows me to copy and paste into latex directly 
                                    print(f"{i}, & {period} \\\\")
                              i += 1
                  
                  average_period = np.mean(period_in_years)
                  average_periods[name] = average_period
                  std_dev_period = np.std(period_in_years)
                  
                  #print(f"{name} - Average Period: {average_period:.2f} years, Standard Deviation: {std_dev_period:.2f} years")
                  if bootstrapping == False:
                        print(f"{name} % {average_period:.2f} % {std_dev_period:.2f} \\\\")

            else:
                  if bootstrapping == False:
                        print(f"{name} - Not enough crossings to calculate period")
      
      return average_periods, halleys_periods

# normal function name length
def calculate_closest_and_furthest_approach_between_earth_and_halleys_comet(positions, times, bootstrapping):
      if 'Earth' in positions and 'Halleys_Comet' in positions:
            earth_pos = positions['Earth']
            comet_pos = positions['Halleys_Comet']
            
            distances = ((earth_pos['x'] - comet_pos['x'])**2 +
                              (earth_pos['y'] - comet_pos['y'])**2 +
                              (earth_pos['z'] - comet_pos['z'])**2) **0.5
            
            closest_idx = np.argmin(distances)
            furthest_idx = np.argmax(distances)
            
            closest_approach = distances[closest_idx] / 1.496e11  # Convert to AU
            furthest_approach = distances[furthest_idx] / 1.496e11
            closest_time = times[closest_idx] / (60 * 60 * 24 * 365.25)  # Convert seconds to years
            furthest_time = times[furthest_idx] / (60 * 60 * 24 * 365.25)
            
            if bootstrapping:
                  return [closest_approach,closest_time,furthest_approach,furthest_time]
            print(f"Closest Approach: {closest_approach:.2e} AU at {closest_time:.2f} years, "
                  f"Furthest Approach: {furthest_approach:.2e} AU at {furthest_time:.2f} years")

      else:
            print('Position values for ethier Earth or Halleys comet missing')      

def main(example,time_step,total_steps,plot_count,three_d_graph_only,method,output_img, bootstrapping = False):

      masses, objects_info, formatting = get_masses_and_object_info(example)

      au = 1.496e11  # 1 Astronomical Unit in meters.
      day_in_seconds = 24 * 3600  # Number of seconds in a day.



      # Check if the CSV file exists
      csv_file = 'objects_data.csv'
      if os.path.exists(csv_file) and bootstrapping == False:
            # Load data from CSV
            print('CSV file already exists, loading data from CSV. NOTE if you have changed the example use the --refresh flag to rewrite the csv')
            objects = pd.read_csv(csv_file).to_dict(orient='records')
      else:
            if bootstrapping is not False:
                  # Get a random date in the past 5 years to allow for different initial starting conditions of the model to see how the physics change! This is how we work out the uncertainit from bootstrapping.
                  start_date,end_date = random_past_date(5*365)
            else:
                  start_date,end_date = '2024-03-17','2024-03-18'

            # If the CSV does not exist, fetch data from the api
            objects = []
            for obj in objects_info:
                  horizons_obj = Horizons(
                  id=obj['id'],
                  location='500@0',  # Centered on Sun @399 to be centered on earth
                  # St patricks day is as good as any to get data from 
                  epochs={'start': start_date, 'stop': end_date, 'step': '1d'}
                  )
                  # Extract the stuff I want 
                  eph = horizons_obj.vectors()

                  obj_data = {
                  'name': obj['name'],
                  'mass': masses[obj['name']], 
                  'initial_x': eph['x'][0] * au,
                  'initial_y': eph['y'][0] * au,
                  'initial_z': eph['z'][0] * au,

                  'initial_x_velocity': eph['vx'][0] * (au / day_in_seconds),
                  'initial_y_velocity': eph['vy'][0] * (au / day_in_seconds),
                  'initial_z_velocity': eph['vz'][0] * (au / day_in_seconds)
                  }
                  objects.append(obj_data)
            
            # Save data to CSV as I dont wanna make more API calls than I need to 
            pd.DataFrame(objects).to_csv(csv_file, index=False)

      if bootstrapping == False:
            print(f"Method is {method}" )

      if method == "both" or method =="rk":
            #Do the calculations to process the system
            positions_rk, time, y_crossings_rk , radial_velocity_changes_rk = calculate_orbit_rk4(objects,time_step,total_steps)
            calculate_orbital_periods(y_crossings_rk,bootstrapping) 
            rk_avgerage_periods , rk_halley_periods = calculate_orbital_periods(radial_velocity_changes_rk,bootstrapping)
            rk_earth_halley_appraoches  = calculate_closest_and_furthest_approach_between_earth_and_halleys_comet(positions_rk,time,bootstrapping)

      if method == "both" or method == "vv":
            #Do the calculations to process the system
            positions_euler, time, y_crossings_euler , radial_velocity_changes_euler = calculate_orbit(objects,time_step,total_steps)

            calculate_orbital_periods(y_crossings_euler,bootstrapping) 
            vv_average_periods , vv_halley_periods = calculate_orbital_periods(radial_velocity_changes_euler,bootstrapping)
            vv_earth_halley_appraoches  = calculate_closest_and_furthest_approach_between_earth_and_halleys_comet(positions_euler,time,bootstrapping)

      # Becuase of the nature of Matplotlib with vs code the code will pause while an image is able to be viewed, I want to be able to run my code and do something once its ran and be back once ONLY its all ran. 
      # Having the plotting seperatly allows this while also making sure 

      # If the bootstrapping method is used these are the results that need returned for further processing, False is used for if the result should not be processes
      if bootstrapping:
            if method == "rk":
                  return start_date , rk_avgerage_periods ,rk_halley_periods,rk_earth_halley_appraoches,False ,  False, False
            if method == "vv":
                  return start_date , vv_average_periods, vv_halley_periods , vv_earth_halley_appraoches, False , False , False
            if method == "both":
                  return start_date ,rk_avgerage_periods,rk_halley_periods ,rk_earth_halley_appraoches,vv_average_periods ,vv_halley_periods, vv_earth_halley_appraoches

      # If a user is boot strapping to find values it is unlikely theywant a large amount of graphs produced
      if method == "both" or method == "rk":
            plot_orbit(positions_rk, formatting,output_img + '_rk' ,plot_count,three_d_graph_only)
      if method == "both" or method =="vv":
            plot_orbit(positions_euler, formatting,output_img + '_vv',plot_count,three_d_graph_only)
     

if __name__ == "__main__":
      # Check if anything is passsed in.
      parser = argparse.ArgumentParser(
            description="Simulate solar system"
      )
      parser.add_argument("--refresh", action="store_true", help="Delete current CSV ")
      parser.add_argument("--single_graph", action="store_true", help="Just a single pretty 3D grpah ")
      parser.add_argument("--vv" ,  action="store_true", help="Euler method")
      parser.add_argument("--both" , action="store_true", help="both RK and Euler")

      parser.add_argument(
            "--example",
            type=str,
            default="solar_system",
            help="Show an exmaple system, currently implemented 'solar_system' and 'earth_and_moon'",
      )

      parser.add_argument(
            "--bootstrap_count",
            type= int,
            default=0,
            help = "work out uncertainity via bootstraping, This technique is taking a random value between the given uncertainities of any of the provided values and seeing how the simulation differs"
      )

      parser.add_argument(
            "--thread_count",
            type=int,
            default=0,
            help="Select how many threads you want the bootstrapping to run on"
      )

      parser.add_argument(
            "--time_step",
            type=int,
            default=3600,
            help="Specify in seconds how often the time interval is to be taken"
      )

      parser.add_argument(
            "--time_total",
            type=int,
            default=3,
            help="specify in earth years how much time to take"                    
      )
      parser.add_argument(
            "--plot_count",
            type=int,
            default=1000,
            help="How many points should be plotted in the diagram"
      )
      parser.add_argument(
            "--output_file_name",
            type= str,
            default="output_img",
            help="specify name of the image, note: dont specify file type as that is handelled by the code itself"
      )

      args = parser.parse_args()

      if args.bootstrap_count >  0 and args.thread_count <  0:
            raise Exception( " When using bootstrapping please ensure the amount of threads the process will run on is specified. do this with the --thread_count flag e.g --thread_count 4 (if your unsure what this means just go for 8 and if your computer chugs reduce this as needed)")

      # Convert years to the amount of intervals needed based on the time_step.
      # Takes how many seconds in a year and divies by the time step to get how many time steps in the given year. the plus one is purely for the tqdm loading bar to look nicer.
      total_steps = int(args.time_total*3600*24*365/args.time_step) + 1

      if (args.vv):
            method = "vv"
      elif (args.both):
            method = "both"
      else:
            method = "rk"

      if args.refresh:
            # Delete existing CSV if it exists
            if os.path.exists('objects_data.csv'):
                  os.remove('objects_data.csv')
      if args.bootstrap_count == 0:
            main(args.example,args.time_step,total_steps,args.plot_count,args.single_graph,method,args.output_file_name)
      else:
            print("Using a bootstrap method to get uncertainties, Might take a while. Note if the API call has a 503 error you have made too many requests, lower the thread count or bootstrap_count( or even decrease the time step to reduce time between api calls ). In an ideal world Id set up a database to store this information and refrence it to eliminate the reliance on third parties completely but ah well it do for now")
            with Pool(args.thread_count) as pool:
                  bootstrap_results = pool.starmap(main,[
                        (
                        args.example,args.time_step,total_steps,args.plot_count,args.single_graph,method,args.output_file_name,True             
                        ) for i in range (0,args.bootstrap_count)
                  ])
            
            analyse_bootstrap_results(bootstrap_results)

 
