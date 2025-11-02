# implement the get_data_from() function

from datetime import datetime
from qcodes_loop.data.data_set import load_data
import os 

def get_data_from(start_time, end_time, fig_dir, date_folder):
    data_path = fig_dir + '\\' + date_folder
    # --- Convert string times to datetime objects for comparison ---
    time_format = "%Y-%m-%d\\%H-%M-%S"
    start_time_obj = datetime.strptime(start_time, time_format)
    end_time_obj = datetime.strptime(end_time, time_format)
    print('Start time: ', start_time_obj, '; End time: ', end_time_obj)


    # --- Main Logic ---
    datfiles = [] # Initialize an empty list to store the loaded data

    # Check if the data path actually exists
    if os.path.exists(data_path):
        # List all the experiment folders in the date directory
        experiment_folders = os.listdir(data_path)
        # print('Experiment Folders', os.listdir(data_path))

        for folder_name in experiment_folders:
            try:
                # Extract the time part from the folder name (e.g., "11-09-42" from "11-09-42_...")
                time_part = folder_name.split('_')[0]
                # print('time_part: ', time_part)
                
                # Create a full timestamp string for the current folder
                folder_time_str = f"{date_folder}\\{time_part}"
                
                # Convert the folder's timestamp string to a datetime object
                folder_time_obj = datetime.strptime(folder_time_str, "%Y-%m-%d\\%H-%M-%S")
                # print('folder_time_obj: ', folder_time_obj)

                # Check if the folder's time is within the desired range
                if start_time_obj <= folder_time_obj <= end_time_obj:
                    # Construct the full path to the data folder
                    full_path = os.path.join(data_path, folder_name)
                    
                    print(f"Loading data from: {full_path}")
                    
                    # Load the data and append it to the list
                    # Make sure you have the 'load_data' function available in your environment
                    dataset = load_data(full_path)
                    datfiles.append(dataset)

            except (ValueError, IndexError):
                # This will catch any folders that don't match the expected naming format
                print(f"Could not parse time from folder: {folder_name}")
                continue

    else:
        print(f"Error: Directory not found at {data_path}")


    # Now 'datafiles' contains all the datasets from the specified time range
    print(f"\nSuccessfully loaded {len(datfiles)} datasets.")

    return datfiles