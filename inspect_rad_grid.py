import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.colors import LogNorm
import argparse

def visualize_radiance_grid(file_path, output_dir=None):
    """
    Load a GOES NetCDF file and visualize its radiance grid.
    
    Args:
        file_path: Path to the NetCDF file
        output_dir: Directory to save output visualization (default: current directory)
    """
    # Create output directory if it doesn't exist
    if output_dir is None:
        output_dir = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)
    
    # Create output filename based on input filename
    base_name = os.path.basename(file_path)
    output_base = os.path.splitext(base_name)[0]
    
    print(f"Loading file: {file_path}")
    try:
        # Open the NetCDF file
        ds = xr.open_dataset(file_path)
        
        # Extract metadata
        file_time = None
        if 'time_bounds' in ds:
            # Convert from epoch time to datetime
            start_time = ds.time_bounds.values[0]
            file_time = pd.to_datetime(start_time)
            print(f"File time: {file_time}")
        
        # Check if this is a radiance file with the expected variable
        if 'Rad' not in ds:
            print("Error: This file does not contain a 'Rad' variable.")
            return False
        
        # Get the radiance data
        rad_data = ds.Rad
        
        # Print basic statistics about the data
        valid_data = rad_data.values[~np.isnan(rad_data.values)]
        print(f"Radiance data shape: {rad_data.shape}")
        print(f"Valid data points: {len(valid_data)} out of {np.prod(rad_data.shape)}")
        print(f"Min value: {np.nanmin(rad_data.values):.6f}")
        print(f"Max value: {np.nanmax(rad_data.values):.6f}")
        print(f"Mean value: {np.nanmean(rad_data.values):.6f}")
        print(f"Median value: {np.nanmedian(rad_data.values):.6f}")
        print(f"Standard deviation: {np.nanstd(rad_data.values):.6f}")
        
        # Check if we have geographic extent information
        if 'geospatial_lat_lon_extent' in ds:
            geo_extent = ds.geospatial_lat_lon_extent
            print(f"Geographic coverage: {geo_extent.geospatial_westbound_longitude}° to {geo_extent.geospatial_eastbound_longitude}° longitude")
            print(f"                     {geo_extent.geospatial_southbound_latitude}° to {geo_extent.geospatial_northbound_latitude}° latitude")
        
        # Create multiple visualizations with different scales to better understand the data
        
        # 1. Basic linear scale visualization
        plt.figure(figsize=(12, 8))
        plt.imshow(rad_data.values, cmap='inferno')
        plt.colorbar(label='Radiance (mW m⁻² sr⁻¹ (cm⁻¹)⁻¹)')
        plt.title(f"GOES-16 Band 7 (3.9μm) Radiance - Linear Scale\n{base_name}")
        plt.axis('off')  # Turn off axis for cleaner view
        
        # Save linear visualization
        linear_path = os.path.join(output_dir, f"{output_base}_linear.png")
        plt.savefig(linear_path, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Saved linear scale visualization to: {linear_path}")
        
        # 2. Log scale visualization to better see variations across orders of magnitude
        plt.figure(figsize=(12, 8))
        # Using LogNorm to enhance visibility of small variations
        plt.imshow(rad_data.values, cmap='inferno', norm=LogNorm(vmin=max(0.001, np.nanmin(rad_data.values)), 
                                                              vmax=np.nanmax(rad_data.values)))
        plt.colorbar(label='Radiance (mW m⁻² sr⁻¹ (cm⁻¹)⁻¹) - Log Scale')
        plt.title(f"GOES-16 Band 7 (3.9μm) Radiance - Log Scale\n{base_name}")
        plt.axis('off')
        
        # Save log visualization
        log_path = os.path.join(output_dir, f"{output_base}_log.png")
        plt.savefig(log_path, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Saved log scale visualization to: {log_path}")
        
        # 3. Histogram of radiance values to understand the distribution
        plt.figure(figsize=(12, 6))
        plt.hist(valid_data, bins=100, log=True)
        plt.xlabel('Radiance (mW m⁻² sr⁻¹ (cm⁻¹)⁻¹)')
        plt.ylabel('Frequency (log scale)')
        plt.title(f"Distribution of Radiance Values\n{base_name}")
        plt.grid(True, alpha=0.3)
        
        # Save histogram
        hist_path = os.path.join(output_dir, f"{output_base}_histogram.png")
        plt.savefig(hist_path, dpi=150)
        plt.close()
        print(f"Saved histogram to: {hist_path}")
        
        # 4. Enhanced visualization to highlight potential hotspots
        # Calculate threshold for highlighting (e.g., mean + 3*std)
        threshold = np.nanmean(rad_data.values) + 3 * np.nanstd(rad_data.values)
        print(f"Hotspot threshold (mean + 3*std): {threshold:.6f}")
        
        # Create a masked array for values above threshold
        plt.figure(figsize=(12, 8))
        plt.imshow(rad_data.values, cmap='gray', vmax=np.nanmean(rad_data.values) * 1.5)
        
        # Overlay potential hotspots in red
        hotspot_mask = rad_data.values > threshold
        if np.any(hotspot_mask):
            print(f"Found {np.sum(hotspot_mask)} potential hotspots above threshold")
            # Create a highlighted overlay for hotspots
            highlight = np.zeros_like(rad_data.values)
            highlight[hotspot_mask] = rad_data.values[hotspot_mask]
            plt.imshow(highlight, cmap='hot', vmin=threshold, vmax=np.nanmax(rad_data.values))
        else:
            print("No hotspots found above threshold")
        
        plt.colorbar(label='Radiance (mW m⁻² sr⁻¹ (cm⁻¹)⁻¹)')
        plt.title(f"GOES-16 Band 7 (3.9μm) Radiance - Hotspot Enhanced\n{base_name}")
        plt.axis('off')
        
        # Save hotspot visualization
        hotspot_path = os.path.join(output_dir, f"{output_base}_hotspots.png")
        plt.savefig(hotspot_path, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Saved hotspot visualization to: {hotspot_path}")
        
        return True
        
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize GOES-16 radiance grid from NetCDF file')
    parser.add_argument('file_path', help='Path to NetCDF file')
    parser.add_argument('--output_dir', help='Directory to save output visualizations')
    
    args = parser.parse_args()
    
    # Ensure pandas is imported for datetime conversion
    import pandas as pd
    
    visualize_radiance_grid(args.file_path, args.output_dir)