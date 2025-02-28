import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import argparse

def parse_time_from_filename(filename):
    """Extract start time from GOES filename."""
    # OR_ABI-L2-FDCC-M6_G16_s20221051600000_e20221051609000_c20221051609157.nc
    parts = filename.split('_')
    if len(parts) >= 4:
        timestr = parts[3][1:14]  # Extract s20221051600000 -> 20221051600000
        try:
            year = int(timestr[0:4])
            doy = int(timestr[4:7])
            hour = int(timestr[7:9])
            minute = int(timestr[9:11])
            second = int(timestr[11:13])
            
            # Convert day of year to date
            date = datetime(year, 1, 1) + pd.Timedelta(days=doy-1)
            return datetime(date.year, date.month, date.day, hour, minute, second)
        except:
            return None
    return None

def process_fire_detection_file(file_path, output_dir):
    """Process a single GOES fire detection NetCDF file."""
    
    file_name = os.path.basename(file_path)
    file_time = parse_time_from_filename(file_name)
    time_str = file_time.strftime("%Y-%m-%d %H:%M:%S") if file_time else "Unknown Time"
    
    # Create output name base (without extension)
    output_base = os.path.splitext(file_name)[0]
    
    try:
        # Open the NetCDF file
        ds = xr.open_dataset(file_path)
        
        # Print basic information
        print(f"Processing: {file_name}")
        print(f"  Time: {time_str}")
        print(f"  Variables: {list(ds.data_vars)}")
        
        # Extract metadata about the scan
        if 'goes_imager_projection' in ds:
            print(f"  Projection: GOES imager projection")
        
        # 1. Create visualization image of fire detections
        plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ccrs.Mercator())
        
        # Add geographic features
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.add_feature(cfeature.STATES, linewidth=0.25)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        
        # Add grid lines
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        
        # Extract fire mask and coordinates
        fire_mask = ds.Mask == 2  # 2 = fire detected
        
        # Set map extent - focus on continental US
        ax.set_extent([-125, -66, 24, 50], crs=ccrs.PlateCarree())
        
        # Check if any fires are detected
        if fire_mask.sum() > 0:
            # Get latitude and longitude coordinates
            if 'x' in ds.dims and 'y' in ds.dims:
                # Convert x,y to lat,lon if needed
                if 'goes_imager_projection' in ds:
                    # Complex conversion from x,y to lat,lon using projection info
                    # This is simplified - proper conversion would use pyproj
                    pass
            
            # For simplicity, we'll use the built-in coordinates from the dataset
            # Plot fire radiative power with color indicating intensity
            fire_power = ds.Power.where(fire_mask)
            img = fire_power.plot(ax=ax, transform=ccrs.PlateCarree(), 
                                  cmap='hot_r', vmin=0, vmax=200,
                                  add_colorbar=True)
            plt.colorbar(img, label='Fire Radiative Power (MW)')
            
            plt.title(f"GOES-16 Fire Detection - {time_str}")
        else:
            plt.title(f"GOES-16 Fire Detection - {time_str}\nNo fires detected")
        
        # Save visualization
        viz_path = os.path.join(output_dir, f"{output_base}_visualization.png")
        plt.savefig(viz_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved visualization to: {viz_path}")
        
        # 2. Extract fire detection data to CSV
        if fire_mask.sum() > 0:
            # Extract all fire pixels
            fire_df = pd.DataFrame({
                'time': [file_time] * fire_mask.sum(),
                'power_MW': ds.Power.values[fire_mask.values],
                'temp_K': ds.Temp.values[fire_mask.values],
                'area_km2': ds.Area.values[fire_mask.values]
            })
            
            # Add coordinates if available
            if 'lon' in ds and 'lat' in ds:
                fire_df['longitude'] = ds.lon.values[fire_mask.values]
                fire_df['latitude'] = ds.lat.values[fire_mask.values]
            
            # Save to CSV
            csv_path = os.path.join(output_dir, f"{output_base}_fires.csv")
            fire_df.to_csv(csv_path, index=False)
            print(f"  Saved {len(fire_df)} fire detections to: {csv_path}")
            
            return {
                'file': file_name,
                'time': file_time,
                'fire_count': fire_mask.sum(),
                'max_power': float(fire_df['power_MW'].max()) if len(fire_df) > 0 else 0,
                'visualization': viz_path,
                'csv': csv_path if fire_mask.sum() > 0 else None
            }
        else:
            return {
                'file': file_name,
                'time': file_time,
                'fire_count': 0,
                'max_power': 0,
                'visualization': viz_path,
                'csv': None
            }
    
    except Exception as e:
        print(f"  Error processing {file_name}: {str(e)}")
        return {
            'file': file_name,
            'error': str(e)
        }

def process_radiance_file(file_path, output_dir):
    """Process a raw radiance file to look for rocket plumes."""
    ds = xr.open_dataset(file_path)
    
    # Print variables to see what's available
    print(f"Variables: {list(ds.data_vars)}")
    
    # For ABI-L1b-RadC files, we want to extract Band 7 (3.9μm)
    # The variable name will be 'Rad' and it will have a band dimension
    
    # Check if we have the right product
    if 'Rad' in ds:
        # For Band 7, we'll need to select it from the bands
        # This might need adjustment based on actual structure
        if 'band' in ds.Rad.dims:
            # Get band 7 if possible
            try:
                rad_band7 = ds.Rad.sel(band=7)
            except:
                # If direct selection doesn't work, try other approaches
                if ds.Rad.shape[0] >= 7:  # If we have at least 7 bands
                    rad_band7 = ds.Rad[6]  # 0-indexed, so band 7 is index 6
                else:
                    print("  Could not find Band 7")
                    return
        
        # Convert radiance to brightness temperature (simplified)
        # This is a complex conversion, but for detection we can use a simple threshold
        
        # Create visualization
        plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ccrs.Mercator())
        
        # Add geographic features
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.add_feature(cfeature.STATES, linewidth=0.25)
        
        # Focus on Cape Canaveral area
        ax.set_extent([-82, -80, 28, 29.5], crs=ccrs.PlateCarree())
        
        # Plot the radiance (higher values = hotter)
        img = rad_band7.plot(ax=ax, transform=ccrs.PlateCarree(),
                           cmap='inferno', add_colorbar=True)
        
        # Save visualization
        file_time = parse_time_from_filename(os.path.basename(file_path))
        time_str = file_time.strftime("%Y-%m-%d %H:%M:%S") if file_time else "Unknown"
        plt.title(f"GOES-16 IR (3.9μm) - {time_str}")
        
        viz_path = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(file_path))[0]}_radiance.png")
        plt.savefig(viz_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        # Look for extremely high values that could be rocket plumes
        # This would need refinement for an actual detection algorithm
        
        # Return information
        return {
            'file': os.path.basename(file_path),
            'time': file_time,
            'max_radiance': float(rad_band7.max())
        }

def process_directory(input_dir, output_dir=None):
    """Process all NetCDF files in a directory."""
    
    # If no output directory specified, create one based on input
    if output_dir is None:
        output_dir = input_dir + "_processed"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Processing files from: {input_dir}")
    print(f"Saving results to: {output_dir}")
    
    # Find all .nc files in the directory
    nc_files = [f for f in os.listdir(input_dir) if f.endswith('.nc')]
    print(f"Found {len(nc_files)} NetCDF files")
    
    # Process each file
    results = []
    for filename in nc_files:
        file_path = os.path.join(input_dir, filename)
        result = process_radiance_file(file_path, output_dir)
        results.append(result)
    
    # Create a summary file
    summary_df = pd.DataFrame(results)
    summary_df.sort_values('time', inplace=True)
    summary_path = os.path.join(output_dir, "processing_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    
    print(f"\nProcessing complete. Summary saved to: {summary_path}")
    
    # Create a time series of fire counts
    if 'fire_count' in summary_df.columns and 'time' in summary_df.columns:
        plt.figure(figsize=(12, 6))
        plt.plot(summary_df['time'], summary_df['fire_count'], 'ro-')
        plt.title('Fire Detections Over Time')
        plt.xlabel('Time')
        plt.ylabel('Number of Fire Pixels Detected')
        plt.grid(True)
        plt.xticks(rotation=45)
        plt.tight_layout()
        timeseries_path = os.path.join(output_dir, "fire_timeseries.png")
        plt.savefig(timeseries_path)
        plt.close()
        print(f"Time series visualization saved to: {timeseries_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process GOES-16 NetCDF files')
    parser.add_argument('input_dir', help='Directory containing NetCDF files')
    parser.add_argument('--output_dir', help='Directory to save processed files (default: input_dir + "_processed")')
    
    args = parser.parse_args()
    
    process_directory(args.input_dir, args.output_dir)