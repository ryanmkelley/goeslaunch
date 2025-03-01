import xarray as xr
import numpy as np
import json
import os
import sys
from datetime import datetime

def summarize_variable(var):
    """Create a summary of a variable that can be serialized to JSON."""
    
    # Basic information
    info = {
        "name": var.name,
        "dimensions": list(var.dims),
        "shape": list(var.shape),
        "dtype": str(var.dtype),
    }
    
    # Add attributes
    if hasattr(var, 'attrs'):
        info["attributes"] = {k: str(v) for k, v in var.attrs.items()}
    
    # For data arrays, add statistics instead of all values
    if hasattr(var, 'values'):
        # Check if it's a large array
        if np.prod(var.shape) > 100:  # Only summarize large arrays
            # Handle different data types
            if np.issubdtype(var.dtype, np.number):
                # For numeric data, provide statistics
                valid_data = var.values[~np.isnan(var.values)] if hasattr(var.values, 'mask') or np.isnan(var.values).any() else var.values
                if len(valid_data) > 0:
                    info["data_summary"] = {
                        "min": float(np.nanmin(var.values)),
                        "max": float(np.nanmax(var.values)),
                        "mean": float(np.nanmean(var.values)),
                        "std": float(np.nanstd(var.values)),
                        "non_nan_values": int(np.sum(~np.isnan(var.values))),
                        "total_elements": int(np.prod(var.shape))
                    }
                else:
                    info["data_summary"] = "No valid numeric data found"
            else:
                # For non-numeric data, just describe the size
                info["data_summary"] = f"Array of shape {var.shape} with dtype {var.dtype}"
        else:
            # For small arrays, include all values
            try:
                # Try to convert values to a serializable format
                values = var.values.tolist() if hasattr(var.values, 'tolist') else str(var.values)
                info["values"] = values
            except:
                info["values"] = "Data could not be serialized"
    
    return info

def inspect_netcdf(file_path, output_dir):
    """Inspect a NetCDF file and output its structure as JSON."""
    # Make sure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Create output filename based on input filename
    base_name = os.path.basename(file_path)
    output_file = os.path.join(output_dir, f"{os.path.splitext(base_name)[0]}_structure.json")
    
    # Load the NetCDF file
    print(f"Opening file: {file_path}")
    try:
        ds = xr.open_dataset(file_path)
        
        # Create a dictionary to store file information
        file_info = {
            "file_name": base_name,
            "dimensions": {dim: size for dim, size in ds.dims.items()},
            "coordinates": [],
            "data_variables": [],
            "global_attributes": {k: str(v) for k, v in ds.attrs.items()}
        }
        
        # Process coordinates
        for name, var in ds.coords.items():
            file_info["coordinates"].append(summarize_variable(var))
        
        # Process data variables
        for name, var in ds.data_vars.items():
            file_info["data_variables"].append(summarize_variable(var))
        
        # Write the JSON output
        with open(output_file, 'w') as f:
            json.dump(file_info, f, indent=2)
        
        print(f"File structure saved to: {output_file}")
        ds.close()
        return output_file
    
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        return None

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python inspect_netcdf.py <path_to_netcdf_file>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    output_dir = "processed_results/file_inspection"
    
    inspect_netcdf(file_path, output_dir)