# download_data.py
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import os
from datetime import datetime, timedelta
import time
import pandas as pd

# Create a dictionary of known launches during Oct-Dec 2022
# Source: Space Launch Schedule and Wikipedia
launches = [
    # October 2022
    {"mission": "SpaceX Falcon 9 - Crew-5", "date": datetime(2022, 10, 5, 16, 0, 0), "site": "Kennedy Space Center"},
    {"mission": "Rocket Lab Electron - Argos-4", "date": datetime(2022, 10, 7, 17, 9, 0), "site": "Wallops Island"},
    {"mission": "SpaceX Falcon 9 - Starlink 4-36", "date": datetime(2022, 10, 8, 23, 5, 0), "site": "Cape Canaveral"},
    {"mission": "SpaceX Falcon 9 - Hotbird 13F", "date": datetime(2022, 10, 15, 5, 22, 0), "site": "Cape Canaveral"},
    {"mission": "ULA Atlas V - JPSS-2", "date": datetime(2022, 10, 27, 8, 49, 0), "site": "Vandenberg SFB"},
    {"mission": "SpaceX Falcon Heavy - USSF-44", "date": datetime(2022, 10, 31, 13, 41, 0), "site": "Kennedy Space Center"},
    
    # November 2022
    {"mission": "SpaceX Falcon 9 - Hotbird 13G", "date": datetime(2022, 11, 3, 1, 22, 0), "site": "Cape Canaveral"},
    {"mission": "Rocket Lab Electron - 'Catch Me If You Can'", "date": datetime(2022, 11, 4, 17, 27, 0), "site": "Mahia Peninsula, NZ"},
    {"mission": "SpaceX Falcon 9 - Galaxy 31 & 32", "date": datetime(2022, 11, 12, 16, 6, 0), "site": "Cape Canaveral"},
    {"mission": "NASA Artemis I - SLS", "date": datetime(2022, 11, 16, 6, 47, 0), "site": "Kennedy Space Center"},
    {"mission": "SpaceX Falcon 9 - Eutelsat 10B", "date": datetime(2022, 11, 23, 2, 57, 0), "site": "Cape Canaveral"},
    {"mission": "SpaceX Falcon 9 - Dragon CRS-26", "date": datetime(2022, 11, 26, 19, 20, 0), "site": "Kennedy Space Center"},
    
    # December 2022
    {"mission": "SpaceX Falcon 9 - OneWeb #15", "date": datetime(2022, 12, 8, 22, 27, 0), "site": "Cape Canaveral"},
    {"mission": "ULA Delta IV Heavy - NROL-91", "date": datetime(2022, 12, 10, 23, 46, 0), "site": "Vandenberg SFB"},
    {"mission": "SpaceX Falcon 9 - SWOT", "date": datetime(2022, 12, 16, 11, 46, 0), "site": "Vandenberg SFB"},
    {"mission": "SpaceX Falcon 9 - Starlink 4-37", "date": datetime(2022, 12, 28, 9, 34, 0), "site": "Cape Canaveral"}
]

def date_to_year_doy(date):
    """Convert datetime to year and day of year."""
    year = date.strftime('%Y')
    doy = date.strftime('%j')
    return year, doy

def list_goes_files_for_hour(product, year, doy, hour):
    """List files in the S3 bucket for specific hour."""
    s3_client = boto3.client('s3', 
                          region_name='us-east-1', 
                          config=Config(signature_version=UNSIGNED))  # Use UNSIGNED from botocore
    
    prefix = f"{product}/{year}/{doy}/{hour:02d}/"
    
    response = s3_client.list_objects_v2(
        Bucket='noaa-goes16',
        Prefix=prefix
    )
    
    if 'Contents' in response:
        return [item['Key'] for item in response['Contents']]
    return []

def download_goes_file(file_key, local_dir='data'):
    """Download a file from S3 bucket."""
    s3_client = boto3.client('s3', 
                          region_name='us-east-1', 
                          config=Config(signature_version=UNSIGNED))
    
    # Create directory if it doesn't exist
    os.makedirs(local_dir, exist_ok=True)
    
    # Extract filename from key
    filename = os.path.basename(file_key)
    local_path = os.path.join(local_dir, filename)
    
    # Don't re-download if file exists
    if os.path.exists(local_path):
        print(f"File {filename} already exists, skipping")
        return local_path
    
    # Download file
    print(f"Downloading {filename}...")
    s3_client.download_file('noaa-goes16', file_key, local_path)
    return local_path

def download_launch_data(launches, window_hours=2, product='ABI-L1b-RadC'):
    """Download data for all launches within specified time window."""
    all_files = []
    
    for launch in launches:
        print(f"\nProcessing launch: {launch['mission']}")
        print(f"Launch date/time: {launch['date']}")
        print(f"Launch site: {launch['site']}")
        
        # Create folder name from mission (sanitize for filesystem)
        mission_folder = launch['mission'].replace('/', '_').replace(' ', '_')
        
        # Calculate window start and end times
        start_time = launch['date'] - timedelta(hours=window_hours)
        end_time = launch['date'] + timedelta(hours=window_hours)
        
        print(f"Downloading data from {start_time} to {end_time}")
        
        # Download data for each hour in the window
        current_time = start_time
        while current_time <= end_time:
            year, doy = date_to_year_doy(current_time)
            hour = current_time.hour
            
            # Create folder for this launch
            folder_path = os.path.join('data', mission_folder)
            
            # List files for this hour
            files = list_goes_files_for_hour(product, year, doy, hour)
            
            if files:
                # Count total files and filter files
                total_files = len(files)
                filtered_files = [f for f in files if "C07" in f or "C14" in f]
                
                print(f"  Found {len(filtered_files)} relevant files out of {total_files} total for {current_time.strftime('%Y-%m-%d %H:%M')}")
                
                # Download each filtered file
                for file_key in filtered_files:
                    local_path = download_goes_file(file_key, folder_path)
                    all_files.append(local_path)
                    # Sleep briefly to avoid overwhelming connection
                    time.sleep(0.5)
            else:
                print(f"  No files found for {current_time.strftime('%Y-%m-%d %H:%M')}")
            
            # Move to next hour
            current_time += timedelta(hours=1)
    
    return all_files

if __name__ == "__main__":
    # Create data directory
    os.makedirs('data', exist_ok=True)
    
    # Download data for all launches
    files = download_launch_data(launches, window_hours=2, product='ABI-L1b-RadC')
    
    print(f"\nDownloaded {len(files)} files")
    print(f"Estimated storage used: {len(files) * 20 / 1024:.2f} GB")