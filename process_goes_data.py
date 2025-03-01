import os
import pyproj
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import math
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# ===============================================================
# Configuration and Helper Functions
# ===============================================================

def parse_time_from_filename(filename):
    """Extract start time from GOES filename."""
    # OR_ABI-L1b-RadC-M6C07_G16_s20221051600000_e20221051609000_c20221051609157.nc
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
            date = datetime(year, 1, 1)
            date = date.replace(day=1) + pd.Timedelta(days=doy-1)
            return datetime(date.year, date.month, date.day, hour, minute, second)
        except:
            return None
    return None

def haversine_distance(lat1, lon1, lat2, lon2):
    """Calculate the great circle distance between two points in kilometers."""
    # Convert degrees to radians
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
    
    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    r = 6371  # Radius of Earth in kilometers
    return c * r

def create_circular_mask(lon_grid, lat_grid, center_lon, center_lat, radius_km):
    """Create a mask for points within radius_km of the center coordinates."""
    mask = np.zeros_like(lon_grid, dtype=bool)
    
    # Calculate distances for each point
    for i in range(lon_grid.shape[0]):
        for j in range(lon_grid.shape[1]):
            distance = haversine_distance(lat_grid[i, j], lon_grid[i, j], center_lat, center_lon)
            if distance <= radius_km:
                mask[i, j] = True
    
    return mask

# ===============================================================
# Core Data Processing Functions
# ===============================================================

def process_radiance_file(file_path, output_dir, center_coords, radius_km, std_threshold=3.0):
    """
    Process a GOES radiance file to detect potential rocket launches.
    
    Args:
        file_path: Path to the NetCDF file
        output_dir: Directory to save outputs
        center_coords: (lat, lon) tuple of the area of interest center
        radius_km: Radius in kilometers around the center to analyze
        std_threshold: Number of standard deviations above mean to consider a hotspot
    
    Returns:
        Dictionary with processing results
    """
    # Extract file information
    file_name = os.path.basename(file_path)
    is_band7 = 'C07' in file_name
    is_band14 = 'C14' in file_name
    
    if not (is_band7 or is_band14):
        print(f"Unknown band in file: {file_name}, skipping.")
        return None
    
    band_name = "3.9μm (Band 7)" if is_band7 else "11.2μm (Band 14)"
    file_time = parse_time_from_filename(file_name)
    time_str = file_time.strftime("%Y-%m-%d %H:%M:%S") if file_time else "Unknown Time"
    
    print(f"Processing: {file_name}")
    print(f"  Time: {time_str}")
    print(f"  Band: {band_name}")
    
    # Open the file and extract the data
    try:
        import pyproj
        from pyproj import Transformer
        
        ds = xr.open_dataset(file_path)
        
        # Display available variables
        print(f"  Variables: {list(ds.data_vars)}")
        
        if 'geospatial_lat_lon_extent' in ds:
            geo_extent = ds.geospatial_lat_lon_extent
            print(f"  File covers: {geo_extent.geospatial_westbound_longitude}° to {geo_extent.geospatial_eastbound_longitude}° longitude")
            print(f"               {geo_extent.geospatial_southbound_latitude}° to {geo_extent.geospatial_northbound_latitude}° latitude")
            
            # Check if our area of interest is within this extent
            center_lat, center_lon = center_coords
            is_in_bounds = (
                geo_extent.geospatial_westbound_longitude <= center_lon <= geo_extent.geospatial_eastbound_longitude and
                geo_extent.geospatial_southbound_latitude <= center_lat <= geo_extent.geospatial_northbound_latitude
            )
            print(f"  Kennedy Space Center coordinates within file bounds: {is_in_bounds}")
        
        # Extract the radiance data
        if 'Rad' not in ds:
            print("  No radiance data found in file")
            return None
        
        rad_data = ds.Rad
        
        # Get the projection information from the file
        if 'goes_imager_projection' not in ds:
            print("  No projection information found in file")
            return None
            
        # Get the GOES projection parameters
        proj_info = ds.goes_imager_projection
        h = proj_info.perspective_point_height
        lambda_0 = proj_info.longitude_of_projection_origin
        
        # Create a proper projection transformer
        goes_proj = f"+proj=geos +h={h} +lon_0={lambda_0} +sweep=x"
        transformer = Transformer.from_crs(
            "+proj=latlong +datum=WGS84",
            goes_proj,
            always_xy=True
        )
        
        # Convert our center coordinates to satellite projection coordinates
        center_lat, center_lon = center_coords
        x_center, y_center = transformer.transform(center_lon, center_lat)
        
        # Get the x and y coordinates from the dataset
        x = ds.x.values
        y = ds.y.values
        
        # Create coordinate grids
        x_grid, y_grid = np.meshgrid(x, y)
        
        # Create a mask for pixels within our radius (in meters)
        # The satellite coordinates are in meters from the projection center
        sat_radius = radius_km * 1000  # Convert km to meters
        distance = np.sqrt((x_grid - x_center)**2 + (y_grid - y_center)**2)
        mask = distance <= sat_radius
        
        # Apply the mask to the radiance data
        rad_values = rad_data.values
        masked_rad = np.where(mask, rad_values, np.nan)
        
        # Create full image visualization to understand what we're looking at
        plt.figure(figsize=(15, 10))
        plt.pcolormesh(x_grid, y_grid, rad_values, cmap='inferno', shading='auto')
        plt.colorbar(label=f'Radiance (W/m²/sr/μm) - {band_name}')
        plt.title(f"GOES-16 Full Image - {band_name}\n{time_str}")
        
        # Mark where we think Kennedy Space Center is
        plt.scatter([x_center], [y_center], c='red', marker='x', s=100, label='Kennedy Space Center')
        
        # Draw a circle for our area of interest
        theta = np.linspace(0, 2*np.pi, 100)
        circle_x = x_center + sat_radius * np.cos(theta)
        circle_y = y_center + sat_radius * np.sin(theta)
        plt.plot(circle_x, circle_y, 'r--', linewidth=1.5, label=f'{radius_km}km Radius')
        
        plt.legend()
        
        # Save this visualization
        debug_viz_path = os.path.join(output_dir, f"{os.path.splitext(file_name)[0]}_full_image.png")
        plt.savefig(debug_viz_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved full image visualization to: {debug_viz_path}")
        
        # Check if we have any valid data points in our mask
        valid_data = masked_rad[~np.isnan(masked_rad)]
        if len(valid_data) == 0:
            print(f"  No valid data points in the specified radius. Check coordinates or increase radius.")
            
            # As a fallback, let's create a visualization of the whole image to debug
            plt.figure(figsize=(12, 8))
            plt.pcolormesh(x_grid, y_grid, rad_values, cmap='inferno', shading='auto')
            plt.colorbar(label=f'Radiance (W/m²/sr/μm) - {band_name}')
            plt.title(f"GOES-16 Full Image - {band_name}\n{time_str}")
            plt.scatter([x_center], [y_center], c='red', marker='x', s=100)  # Mark our target point
            
            viz_path = os.path.join(output_dir, f"{os.path.splitext(file_name)[0]}_debug.png")
            plt.savefig(viz_path, dpi=150, bbox_inches='tight')
            plt.close()
            
            return {
                'file': file_name,
                'time': file_time,
                'band': 'C07' if is_band7 else 'C14',
                'error': 'No valid data in radius',
                'debug_viz': viz_path
            }
        
        # Calculate statistics on the masked data
        data_mean = np.nanmean(valid_data)
        data_std = np.nanstd(valid_data)
        data_max = np.nanmax(valid_data)
        data_min = np.nanmin(valid_data)
        
        # Find hotspots (points above threshold)
        hotspot_threshold = data_mean + (std_threshold * data_std)
        hotspots = (masked_rad > hotspot_threshold) & (~np.isnan(masked_rad))
        hotspot_count = np.sum(hotspots)
        
        print(f"  Data summary - Mean: {data_mean:.2f}, Std: {data_std:.2f}, Min: {data_min:.2f}, Max: {data_max:.2f}")
        print(f"  Found {hotspot_count} potential hotspots (>{hotspot_threshold:.2f})")
        
        # Create visualization
        plt.figure(figsize=(12, 8))
        
        # First, plot the masked radiance data
        plt.pcolormesh(x_grid, y_grid, masked_rad, cmap='inferno', shading='auto')
        
        # Add a colorbar
        cbar = plt.colorbar()
        cbar.set_label(f'Radiance (W/m²/sr/μm) - {band_name}')
        
        # Draw a circle for the area of interest
        theta = np.linspace(0, 2*np.pi, 100)
        circle_x = x_center + sat_radius * np.cos(theta)
        circle_y = y_center + sat_radius * np.sin(theta)
        plt.plot(circle_x, circle_y, 'w--', linewidth=2)
        
        # Highlight hotspots if any
        if hotspot_count > 0:
            hotspot_y, hotspot_x = np.where(hotspots)
            hotspot_x_coords = x[hotspot_x]
            hotspot_y_coords = y[hotspot_y]
            plt.scatter(hotspot_x_coords, hotspot_y_coords, c='cyan', 
                       s=50, edgecolors='white', marker='o')
        
        plt.title(f"GOES-16 Radiance - {band_name}\n{time_str}")
        
        # Save the visualization
        output_base = os.path.splitext(file_name)[0]
        viz_path = os.path.join(output_dir, f"{output_base}_radiance.png")
        plt.savefig(viz_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        # Return the results
        return {
            'file': file_name,
            'time': file_time,
            'band': 'C07' if is_band7 else 'C14',
            'mean_radiance': float(data_mean),
            'std_radiance': float(data_std),
            'max_radiance': float(data_max),
            'min_radiance': float(data_min),
            'hotspot_threshold': float(hotspot_threshold),
            'hotspot_count': int(hotspot_count),
            'visualization': viz_path
        }
    
    except Exception as e:
        import traceback
        print(f"  Error processing file: {str(e)}")
        print(traceback.format_exc())
        return {
            'file': file_name,
            'time': file_time,
            'error': str(e)
        }

def visualize_hotspots(results_df, output_dir, center_coords, radius_km, selected_band=None):
    """
    Create visualizations of hotspot analysis for review.
    
    Args:
        results_df: DataFrame with processing results
        output_dir: Directory to save outputs
        center_coords: (lat, lon) tuple of the area of interest center
        radius_km: Radius in kilometers
        selected_band: Optional band to filter by ('C07' or 'C14')
    """
    if len(results_df) == 0:
        print("No results to visualize")
        return
    
    # Filter by band if specified
    if selected_band:
        band_df = results_df[results_df['band'] == selected_band]
    else:
        band_df = results_df
    
    if len(band_df) == 0:
        print(f"No results for band {selected_band}")
        return
    
    # Create time series of hotspot counts
    plt.figure(figsize=(14, 7))
    plt.plot(band_df['time'], band_df['hotspot_count'], 'ro-', linewidth=2, markersize=8)
    plt.title(f'Hotspot Detections Over Time - {selected_band if selected_band else "All Bands"}')
    plt.xlabel('Time')
    plt.ylabel('Number of Hotspots Detected')
    plt.grid(True)
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # Save the time series
    timeseries_path = os.path.join(output_dir, f"hotspot_timeseries{'_'+selected_band if selected_band else ''}.png")
    plt.savefig(timeseries_path, dpi=150)
    plt.close()
    print(f"Saved hotspot time series to: {timeseries_path}")
    
    # Create time series of maximum radiance
    plt.figure(figsize=(14, 7))
    plt.plot(band_df['time'], band_df['max_radiance'], 'bo-', linewidth=2, markersize=8)
    plt.title(f'Maximum Radiance Over Time - {selected_band if selected_band else "All Bands"}')
    plt.xlabel('Time')
    plt.ylabel('Maximum Radiance (W/m²/sr/μm)')
    plt.grid(True)
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    # Save the radiance time series
    radiance_path = os.path.join(output_dir, f"max_radiance_timeseries{'_'+selected_band if selected_band else ''}.png")
    plt.savefig(radiance_path, dpi=150)
    plt.close()
    print(f"Saved max radiance time series to: {radiance_path}")

def process_directory(input_dir, output_dir, center_coords, radius_km, band_choice=None, std_threshold=3.0):
    """
    Process all NetCDF files in a directory to detect rocket launches.
    
    Args:
        input_dir: Directory containing NetCDF files
        output_dir: Directory to save outputs
        center_coords: (lat, lon) tuple of the area of interest center
        radius_km: Radius in kilometers around the center to analyze
        band_choice: Which band to process ('C07', 'C14', or None for both)
        std_threshold: Number of standard deviations for hotspot detection
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Processing files from: {input_dir}")
    print(f"Saving results to: {output_dir}")
    print(f"Area of interest: {center_coords} with radius {radius_km}km")
    print(f"Hotspot threshold: Mean + {std_threshold} standard deviations")
    
    # Find all .nc files in the directory
    all_nc_files = [f for f in os.listdir(input_dir) if f.endswith('.nc')]
    
    # Filter files based on band choice
    if band_choice == 'C07':
        nc_files = [f for f in all_nc_files if 'C07' in f]
        print(f"Found {len(nc_files)} Band 7 (3.9μm) files")
    elif band_choice == 'C14':
        nc_files = [f for f in all_nc_files if 'C14' in f]
        print(f"Found {len(nc_files)} Band 14 (11.2μm) files")
    else:
        nc_files = all_nc_files
        print(f"Found {len(nc_files)} total NetCDF files")
    
    # Process each file
    results = []
    for filename in nc_files:
        file_path = os.path.join(input_dir, filename)
        result = process_radiance_file(file_path, output_dir, center_coords, radius_km, std_threshold)
        if result:
            results.append(result)
    
    # Create a summary file
    if results:
        summary_df = pd.DataFrame(results)
        # Filter out rows with errors
        valid_results = summary_df[~summary_df.get('error', pd.Series([False] * len(summary_df))).notna()]
        
        if len(valid_results) > 0:
            # Sort by time
            valid_results.sort_values('time', inplace=True)
            
            # Save to CSV
            summary_path = os.path.join(output_dir, "processing_summary.csv")
            valid_results.to_csv(summary_path, index=False)
            print(f"\nProcessing complete. Summary saved to: {summary_path}")
            
            # Count hotspots
            total_hotspots = valid_results['hotspot_count'].sum()
            print(f"Found a total of {total_hotspots} hotspots across all files")
            
            # Ask if user wants visualizations
            if total_hotspots > 0:
                resp = input(f"Found {total_hotspots} hotspots. Generate visualizations? (Y/N): ")
                if resp.lower() == 'y':
                    # Generate additional visualizations for Band 7
                    if 'C07' in valid_results['band'].unique():
                        visualize_hotspots(valid_results, output_dir, center_coords, radius_km, 'C07')
                    
                    # Generate additional visualizations for Band 14
                    if 'C14' in valid_results['band'].unique():
                        visualize_hotspots(valid_results, output_dir, center_coords, radius_km, 'C14')
                    
                    # Generate visualizations for all bands together
                    visualize_hotspots(valid_results, output_dir, center_coords, radius_km)
        else:
            print("No valid results to summarize")
    else:
        print("No results generated")

# ===============================================================
# Main script execution
# ===============================================================

if __name__ == "__main__":
    # ===== INPUT PARAMETERS =====
    
    # Input directory containing NetCDF files
    input_dir = "data/NASA_Artemis_I_-_SLS"
    
    # Output directory for results (will be created if it doesn't exist)
    output_dir = "processed_results/NASA_Artemis_I_-_SLS"
    
    # Band choice: 
    # 'C07' - Band 7 (3.9μm shortwave IR) - Best for detecting hot objects
    # 'C14' - Band 14 (11.2μm longwave IR) - Good for general thermal features
    # None - Process both bands
    band_choice = 'C07'
    
    # Area of interest:
    # Kennedy Space Center (Launch Complex 39) coordinates:
    center_coords = (28.608, -80.604)  # (latitude, longitude)
    
    # Radius in kilometers around the center point to analyze
    radius_km = 75
    
    # Statistical threshold for hotspot detection
    # Lower values (2-3) will detect more potential hotspots but may include false positives
    # Higher values (5-8) will be more selective but might miss subtle signatures
    std_threshold = 3.0
    
    # Process the data
    process_directory(input_dir, output_dir, center_coords, radius_km, band_choice, std_threshold)