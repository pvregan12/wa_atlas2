import pandas as pd
from pathlib import Path
import numpy as np
from scipy.interpolate import griddata
import os
import tkinter as tk
from tkinter import ttk
import traceback



class AtlasApp:
    def __init__(self, root):
        self.root = root
        self.root.title("NOAA Atlas 2 Rainfall Tool")
        self.root.geometry("600x475")
        self.root.resizable(True, True)
        
        self.root.eval('tk::PlaceWindow . center')

        self.setup_gui()

    def setup_gui(self):

        # create main frame
        main_frame = ttk.Frame(root, padding="20")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # input fields
        self.state_var = tk.StringVar(value="WA")
        self.lat_var = tk.StringVar(value="47.115078")
        self.long_var = tk.StringVar(value="-123.754755")
        self.region_var = tk.StringVar(value="3")
        self.elev_var = tk.StringVar(value="341")
        
        state_label = ttk.Label(main_frame, text="State:")
        state_label.grid(row=1, column=0, sticky=tk.W, pady=5, padx=(0,10))
        self.state_entry = ttk.Entry(main_frame, textvariable=self.state_var, width=25)
        self.state_entry.grid(row=1, column=1, sticky=(tk.W, tk.E))

        lat_label = ttk.Label(main_frame, text="Latitude (decimal):")
        lat_label.grid(row=2, column=0, sticky=tk.W, pady=5, padx=(0,10))
        self.lat_entry = ttk.Entry(main_frame, textvariable=self.lat_var, width=25)
        self.lat_entry.grid(row=2, column=1, sticky=(tk.W, tk.E))

        long_label = ttk.Label(main_frame, text="Longitude (decimal):")
        long_label.grid(row=3, column=0, sticky=tk.W, pady=5, padx=(0,10))
        self.long_entry = ttk.Entry(main_frame, textvariable=self.long_var, width=25)
        self.long_entry.grid(row=3, column=1, sticky=(tk.W, tk.E))

        region_label = ttk.Label(main_frame, text="Region (1-4):")
        region_label.grid(row=4, column=0, sticky=tk.W, pady=5, padx=(0,10))
        self.region_entry = ttk.Entry(main_frame, textvariable=self.region_var, width=25)
        self.region_entry.grid(row=4, column=1, sticky=(tk.W, tk.E))

        elev_label = ttk.Label(main_frame, text="Elevation (ft NGVD29):")
        elev_label.grid(row=5, column=0, sticky=tk.W, pady=5, padx=(0,10))
        self.elev_entry = ttk.Entry(main_frame, textvariable=self.elev_var, width=25)
        self.elev_entry.grid(row=5, column=1, sticky=(tk.W, tk.E))

        # run button 
        run_button = ttk.Button(main_frame, text="Run", command=self.run_action)
        run_button.grid(row=6, column=0, columnspan=2, pady=2)

        # status
        status_label = ttk.Label(main_frame, text="Status:")
        status_label.grid(row=7, column=0, sticky=tk.W, pady=(10,5))

        self.status_text = tk.Text(main_frame, height=10, width=40)
        self.status_text.grid(row=8, column=0, columnspan=2, pady=5, sticky=(tk.W, tk.E, tk.N, tk.S))

        # scrollbar
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=self.status_text.yview)
        scrollbar.grid(row=8, column=2, sticky=(tk.N, tk.S))
        self.status_text.configure(yscrollcommand=scrollbar.set)

        # grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)

    def run_action(self):
        try:
            state = self.state_entry.get()
            lat = float(self.lat_entry.get())
            long = float(self.long_entry.get())
            region = int(self.region_entry.get())
            elev = float(self.elev_entry.get())

            self.status_text.delete(1.0, tk.END)
            result_df = self.generate_full_data(state, region, lat, long, elev)
            self.update_status_with_dataframe(result_df)
            
        except ValueError as e:
            self.update_status(f"Error: Please enter valid numbers - {e}")
        except Exception as e:
            self.update_status(f"Error: {e}")

    def update_status(self, message:str):
        self.status_text.insert(tk.END, message + "\n")
        self.status_text.see(tk.END)

    def update_status_with_dataframe(self, df):
        # clear
        #self.status_text.delete(1,0, tk.END)

        # convert dataframe to string for nice format
        df_string = df.to_string(
            index=True,
            justify="left",
            float_format=lambda x: f"{x:.2f}",
            max_cols=10,
            max_rows=50
        )
        self.status_text.insert(tk.END, df_string)

    def generate_ascii_name(self, state:str, recurrence:int, duration:int) -> str:
        state = state.lower()
        asc_nm = "na2_"+state+"_"+str(recurrence)+"yr"+str(duration)+"hr.asc"        
        print(f"filename: {asc_nm}")
        return asc_nm

    def clean_ascii_line(self, ascii_line:str) -> str:
        return ascii_line.split(" ")[1].split("\n")[0]

    def read_ascii_grid(self, state:str, recurrence:int, duration:int, lat:float, long:float) -> dict:
        output_dict = {
            'state': state,
            'recurrence': recurrence,
            'duration': duration,
            'latitude': lat,
            'longitude': long
        }
        
        file_name = self.generate_ascii_name(state, recurrence, duration)
        directory = os.getcwd()

        try:
            file_path = Path(directory, file_name)
            with open(file_path, 'r') as data_dump:
                lines_dump = data_dump.readlines()

            # Get numeric values from header
            num_cols = int(self.clean_ascii_line(lines_dump[0]))
            num_rows = int(self.clean_ascii_line(lines_dump[1]))
            xllc = float(self.clean_ascii_line(lines_dump[2]))
            yllc = float(self.clean_ascii_line(lines_dump[3]))
            cell_size = float(self.clean_ascii_line(lines_dump[4]))
            nodata_val = float(self.clean_ascii_line(lines_dump[5]))  # Convert to float

        except Exception as e:
            print(f"Error getting ascii grid: {e}")
            return None  # Fixed: was missing None

        try:
            # Find the rows and cols nearest to the given lat/long
            row_from_bottom = int((lat - yllc) // cell_size)
            row_below = num_rows - row_from_bottom - 1 + 6 # flip + header offset
            row_above = row_below + 1
            
            col_left = int((long - xllc) // cell_size)
            col_right = col_left + 1

            # Get the data lines
            row_lines = [lines_dump[row_below], lines_dump[row_above]]

            # Calculate coordinates
            latitudes = [yllc + (row_below - 6) * cell_size,
                        yllc + (row_below - 6) * cell_size,
                        yllc + (row_above - 6) * cell_size,
                        yllc + (row_above - 6) * cell_size]
            
            longitudes = [xllc + col_left * cell_size,
                        xllc + col_right * cell_size,
                        xllc + col_right * cell_size,
                        xllc + col_left * cell_size]
            
            # Extract and convert values to float immediately
            try:
                values = [float(row_lines[0].split()[col_left]),
                        float(row_lines[0].split()[col_right]),
                        float(row_lines[1].split()[col_right]),
                        float(row_lines[1].split()[col_left])]
            except (IndexError, ValueError) as e:
                print(f"Error extracting values: {e}")
                return None

            # Remove nodata values
            valid_indices = [i for i, val in enumerate(values) if val != nodata_val]
            values = [values[i] for i in valid_indices]
            latitudes = [latitudes[i] for i in valid_indices]
            longitudes = [longitudes[i] for i in valid_indices]
            
            if len(values) == 0:
                print(f"Error: all points are nodata value")
                return None
            
            # Interpolation logic
            if len(values) == 1:
                val_interpd = float(values[0])
            
            elif len(values) == 2:
                p0 = np.array([longitudes[0], latitudes[0]])
                p1 = np.array([longitudes[1], latitudes[1]])
                query_point = np.array([long, lat])

                line_vec = p1 - p0
                line_length_sq = np.dot(line_vec, line_vec)
                
                if line_length_sq == 0:
                    val_interpd = float(values[0])
                else:
                    t = np.dot(query_point - p0, line_vec) / line_length_sq
                    if t < 0 or t > 1:
                        val_interpd = float(np.mean(values))
                    else:
                        val_interpd = float(values[0] + t * (values[1] - values[0]))

            elif len(values) >= 3:
                try:
                    points = np.array([[longitudes[i], latitudes[i]] for i in range(len(values))])
                    z_vals = np.array(values, dtype=float)
                    new_point = np.array([[long, lat]])

                    result = griddata(points, z_vals, new_point, method='linear')
                    
                    if len(result) > 0 and not np.isnan(result[0]):
                        val_interpd = float(result[0])
                    else:
                        print(f"Griddata failed, using mean of {len(values)} points")
                        val_interpd = float(np.mean(values))
                        
                except Exception as e:
                    print(f"Griddata exception: {e}, using mean of {len(values)} points")
                    val_interpd = float(np.mean(values))

            # Final validation
            if np.isnan(val_interpd):
                print(f"Warning: final interpolated value is NaN, using mean")
                val_interpd = float(np.mean(values))

            output_dict['depth'] = val_interpd / 100000 # values are given in 100,000 inches
            output_dict['intensity'] = val_interpd / duration
            
            return output_dict
            
        except Exception as e:
            print(f"Error in coordinate processing: {e}")
            import traceback
            traceback.print_exc()
            return None

    def general_recurrence_nomogram(self, two_year:float, hundred_year:float) -> dict:
        print('a')
        mapping_array = [1, 0.75073, 0.58722, 0.37898, 0.18518, 0]
        print('b')
        mapping_array_comp = [1 - i for i in mapping_array]
        print('c')
        print(10)
        output_dict = {
                        'two_year': two_year,
                        'five_year': (two_year*mapping_array[1]) + (hundred_year*mapping_array_comp[1]),
                        'ten_year': (two_year*mapping_array[2]) + (hundred_year*mapping_array_comp[2]),
                        'twentyfive_year': (two_year*mapping_array[3]) + (hundred_year*mapping_array_comp[3]),
                        'fifty_year': (two_year*mapping_array[4]) + (hundred_year*mapping_array_comp[4]),
                        'hundred_year': hundred_year
        }

        return output_dict

    def compute_one_hour_data(self, six_hour_data:dict, twentyfour_hour_data:dict, region:int, lat:float, long:float, elev:float) -> dict:
            
        x_one = six_hour_data['two_year']
        x_two = twentyfour_hour_data['two_year']

        x_three = six_hour_data['hundred_year']
        x_four = twentyfour_hour_data['hundred_year']

        x_five = lat - 40
        x_six = abs(long) - 100
        zee = elev / 100
            
        if region == 1:
            two_one = 0.019 + 0.711*(x_one*(x_one/x_two)) + 0.001*zee
            hundred_one = 0.338 + 0.670*(x_three*(x_three/x_four)) + 0.001*zee
        elif region == 2:
            two_one = 0.077 + 0.715*(x_one*(x_two/x_three)) - 0.0004*(x_five/x_six)
            hundred_one = 0.187 + 0.833*(x_three*(x_three/x_four))
        elif region == 3:
            two_one = 0.157 + 0.513*(x_one*(x_one/x_two))
            hundred_one = 0.324 + 0.752*(two_one*(x_three/x_one))
        elif region == 4:
            two_one = 0.160 + 0.520*(x_one*(x_one/x_two))
            hundred_one = 0.177 + 0.965*(two_one*(x_three/x_one))
        else:
            print(f"Invalide region number: {region}")
        one_hour = self.general_recurrence_nomogram(two_one, hundred_one)
        return one_hour

    def compute_twelve_hour_depths(self, six_hours:dict, twentyfour_hours:dict) -> dict:
        output_dict = {
                        'two_year': (six_hours['two_year'] + twentyfour_hours['two_year'])/2,
                        'five_year': (six_hours['five_year'] + twentyfour_hours['five_year'])/2,
                        'ten_year': (six_hours['ten_year'] + twentyfour_hours['ten_year'])/2,
                        'twentyfive_year': (six_hours['twentyfive_year'] + twentyfour_hours['twentyfive_year'])/2,
                        'fifty_year': (six_hours['fifty_year'] + twentyfour_hours['fifty_year'])/2,
                        'hundred_year': (six_hours['hundred_year'] + twentyfour_hours['hundred_year'])/2
        }
        return output_dict

    def compute_two_hour_depths(self, one_hours:dict, six_hours:dict, region:int) -> dict:
        if region == 1:
            return {
                "two_year": 0.250*six_hours['two_year']+ 0.750*one_hours['two_year'],
                "five_year": 0.250*six_hours['five_year'] + 0.750*one_hours['five_year'],
                "ten_year": 0.250*six_hours['ten_year'] + 0.750*one_hours['ten_year'],
                "twentyfive_year": 0.250*six_hours['twentyfive_year'] + 0.750*one_hours['twentyfive_year'],
                "fifty_year": 0.250*six_hours['fifty_year'] + 0.750*one_hours['fifty_year'],
                "hundred_year": 0.250*six_hours['hundred_year'] + 0.750*one_hours['hundred_year']
            }
        elif region == 2:
            return {
                    "two_year": 0.278*six_hours['two_year'] + 0.722*one_hours['two_year'],
                    "five_year": 0.278*six_hours['five_year'] + 0.722*one_hours['five_year'],
                    "ten_year": 0.278*six_hours['ten_year'] + 0.722*one_hours['ten_year'],
                    "twentyfive_year": 0.278*six_hours['twentyfive_year'] + 0.722*one_hours['twentyfive_year'],
                    "fifty_year": 0.278*six_hours['fifty_year'] + 0.722*one_hours['fifty_year'],
                    "hundred_year": 0.278*six_hours['hundred_year'] + 0.722*one_hours['hundred_year']
            }
        elif region == 3 or region == 4:
            return {
                    "two_year": 0.240*six_hours['two_year'] + 0.760*one_hours['two_year'],
                    "five_year": 0.240*six_hours['five_year'] + 0.760*one_hours['five_year'],
                    "ten_year": 0.240*six_hours['ten_year'] + 0.760*one_hours['ten_year'],
                    "twentyfive_year": 0.240*six_hours['twentyfive_year'] + 0.760*one_hours['twentyfive_year'],
                    "fifty_year": 0.240*six_hours['fifty_year'] + 0.760*one_hours['fifty_year'],
                    "hundred_year": 0.240*six_hours['hundred_year'] + 0.760*one_hours['hundred_year']
            }
        else:
            print(f"Error: invalid region: {region}")
        
    def compute_three_hour_depths(self, one_hours:dict, six_hours:dict, region:int) -> dict:
        if region == 1:
            return {
                "two_year": 0.467*six_hours['two_year'] + 0.533*one_hours['two_year'],
                "five_year": 0.467*six_hours['five_year'] + 0.533*one_hours['five_year'],
                "ten_year": 0.467*six_hours['ten_year'] + 0.533*one_hours['ten_year'],
                "twentyfive_year": 0.467*six_hours['twentyfive_year'] + 0.533*one_hours['twentyfive_year'],
                "fifty_year": 0.467*six_hours['fifty_year'] + 0.533*one_hours['fifty_year'],
                "hundred_year": 0.467*six_hours['hundred_year'] + 0.533*one_hours['hundred_year']
            }
        elif region == 2:
            return {
                    "two_year": 0.503*six_hours['two_year'] + 0.497*one_hours['two_year'],
                    "five_year": 0.503*six_hours['five_year'] + 0.497*one_hours['five_year'],
                    "ten_year": 0.503*six_hours['ten_year'] + 0.497*one_hours['ten_year'],
                    "twentyfive_year": 0.503*six_hours['twentyfive_year'] + 0.497*one_hours['twentyfive_year'],
                    "fifty_year": 0.503*six_hours['fifty_year'] + 0.497*one_hours['fifty_year'],
                    "hundred_year": 0.503*six_hours['hundred_year'] + 0.497*one_hours['hundred_year']
            }
        elif region == 3 or region == 4:
            return {
                    "two_year": 0.468*six_hours['two_year'] + 0.532*one_hours['two_year'],
                    "five_year": 0.468*six_hours['five_year'] + 0.532*one_hours['five_year'],
                    "ten_year": 0.468*six_hours['ten_year'] + 0.532*one_hours['ten_year'],
                    "twentyfive_year": 0.468*six_hours['twentyfive_year'] + 0.532*one_hours['twentyfive_year'],
                    "fifty_year": 0.468*six_hours['fifty_year'] + 0.532*one_hours['fifty_year'],
                    "hundred_year": 0.468*six_hours['hundred_year'] + 0.532*one_hours['hundred_year']
            }
        else:
            print(f"Error: invalid region: {region}")

    def generate_full_data(self, state:str, region:int, lat:float, long:float, elev:float):

        # first, get 2-year and 100-year 6-hour and 24-hour depths
        
        twosix = self.read_ascii_grid(state, 2, 6, lat, long)
        twotwentyfour = self.read_ascii_grid(state, 2, 24, lat, long)
        hundredsix = self.read_ascii_grid(state, 100, 6, lat, long)
        hundredtwentyfour = self.read_ascii_grid(state, 100, 24, lat, long)
        self.update_status(f"2yr6hr: {twosix['depth']}")
        self.update_status(f"2yr24hr: {twotwentyfour['depth']}")
        self.update_status(f"100yr6hr: {hundredsix['depth']}")
        self.update_status(f"100yr24hr: {hundredtwentyfour['depth']}")

        # get all 6-hour depths, including storing 2 and 24 hour
        six_hours = self.general_recurrence_nomogram(twosix['depth'], hundredsix['depth'])

        # get all 24-hour depths
        twentyfour_hours = self.general_recurrence_nomogram(twotwentyfour['depth'], hundredtwentyfour['depth'])

        # get 1-hour depths
        one_hours = self.compute_one_hour_data(six_hours, twentyfour_hours, region, lat, long, elev)

        # compute the 2 and 3 hour depths
        two_hours = self.compute_two_hour_depths(one_hours, six_hours, region)
        three_hours = self.compute_three_hour_depths(one_hours, six_hours, region)

        # get all 12-hour depths
        twelve_hours = self.compute_twelve_hour_depths(six_hours, twentyfour_hours)
        
        out_df = pd.DataFrame({
            '1-hr': one_hours,
            '2-hr': two_hours,
            '3-hr': three_hours,
            '6-hr': six_hours,
            '12-hr': twelve_hours,
            '24-hr': twentyfour_hours
        })

        #self.update_status_with_dataframe(out_df)
        return out_df
    
# testing
#test = generate_full_data("wa", 3, 47.115047, -123.754755, 341.0)

if __name__ == "__main__":
    root = tk.Tk()
    app = AtlasApp(root)
    root.mainloop()

