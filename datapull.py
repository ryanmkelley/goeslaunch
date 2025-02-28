import boto3
import os
from datetime import datetime, timedelta
import time

# Function to convert date to day of year
def date_to_day_of_year(date):
    return date.strftime('%Y/%j')

# Function to list files for a specific date
def list_goes_files_for_date(product, year, doy, hour=None):
    s3_client = boto3.client('s3', region_name='us-east-1', 
                            config=boto3.session.Config(signature_version=boto3.session.UNSIGNED))
    
    prefix = f"{product}/{year}/{doy}/"
    if hour is not None:
        prefix += f"{hour:02d}/"
    
    response = s3_client.list_objects_v2(
        Bucket='noaa-goes16',
        Prefix=prefix
    )
    
    if 'Contents' in response:
        return [item['Key'] for item in response['Contents']]
    return []

# Function to download a file
def download_goes_file(file_key, local_dir='data'):
    s3_client = boto3.client('s3', region_name='us-east-1', 
                            config=boto3.session.Config(signature_version=boto3.session.UNSIGNED))
    
    # Create directory if it doesn't exist
    os.makedirs(local_dir, exist_ok=True)
    
    # Extract filename from key
    filename = os.path.basename(file_key)
    local_path = os.path.join(local_dir, filename)
    
    # Download file
    print(f"Downloading {filename}...")
    s3_client.download_file('noaa-goes16', file_key, local_path)
    return local_path

# Example: Download FDCC data for a specific date range 
def download_data_for_date_range(start_date, end_date, product='ABI-L2-FDCC'):
    current_date = start_date
    while current_date <= end_date:
        year = current_date.strftime('%Y')
        doy = current_date.strftime('%j')
        
        # Process each hour of the day
        for hour in range(24):
            print(f"Processing {current_date.strftime('%Y-%m-%d')} hour {hour:02d}...")
            
            # List files for this hour
            files = list_goes_files_for_date(product, year, doy, hour)
            
            # Download each file
            for file_key in files:
                download_goes_file(file_key)
                # Sleep briefly to avoid overwhelming the connection
                time.sleep(0.5)
        
        # Move to next day
        current_date += timedelta(days=1)

# Example usage: Download data for a known launch day
# SpaceX Falcon 9 launch on January 31, 2022 from Cape Canaveral
start_date = datetime(2022, 1, 31)
end_date = datetime(2022, 1, 31)
download_data_for_date_range(start_date, end_date)