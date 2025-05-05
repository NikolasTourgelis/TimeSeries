```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from collections import defaultdict

# File paths for the two NetCDF files
file1 = 'athensnew.nc'
file2 = 'thessalonikinew.nc'

# Coordinates for Athens and Thessaloniki
athens_coords = (37.9838, 23.7275)
thessaloniki_coords = (40.6401, 22.9444)

# Open the second NetCDF file
ds2 = nc.Dataset(file2, 'r')

# Extract latitude and longitude arrays from the file
latitudes = ds2.variables['latitude'][:]
longitudes = ds2.variables['longitude'][:]

# Function to find the nearest index for a given coordinate
def find_nearest_idx(array, value):
    return np.abs(array - value).argmin()

# Find the closest grid points for Athens and Thessaloniki
athens_lat_idx = find_nearest_idx(latitudes, athens_coords[0])
athens_lon_idx = find_nearest_idx(longitudes, athens_coords[1])

thessaloniki_lat_idx = find_nearest_idx(latitudes, thessaloniki_coords[0])
thessaloniki_lon_idx = find_nearest_idx(longitudes, thessaloniki_coords[1])

# Extract the 'u100' and 'v100' variables (3D data: time, lat, lon)
u100_2 = ds2.variables['u100'][:]
v100_2 = ds2.variables['v100'][:]

# Calculate wind speed from u100 and v100
wind_speed_2 = np.sqrt(u100_2**2 + v100_2**2)

# Extract time series for Athens and Thessaloniki
u100_athens = u100_2[:, athens_lat_idx, athens_lon_idx]
u100_thessaloniki = u100_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]
v100_athens = v100_2[:, athens_lat_idx, athens_lon_idx]
v100_thessaloniki = v100_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]
wind_speed_athens = wind_speed_2[:, athens_lat_idx, athens_lon_idx]
wind_speed_thessaloniki = wind_speed_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

# Extract the 'valid_time' variable and convert it to datetime objects
valid_time = ds2.variables['valid_time'][:]
time_units = ds2.variables['valid_time'].units

# Convert valid_time to datetime using num2date from netCDF4
dates = nc.num2date(valid_time, units=time_units)

# Convert cftime objects to standard datetime objects
dates = [datetime(dt.year, dt.month, dt.day) for dt in dates]

# Close the NetCDF file
ds2.close()

# Function to compute seasonal statistics across all years
def compute_seasonal_stats(dates, data):
    seasonal_data = defaultdict(list)
    for date, value in zip(dates, data):
        # Determine season based on month
        if date.month in [12, 1, 2]:  # Winter
            season = 'Winter'
        elif date.month in [3, 4, 5]:  # Spring
            season = 'Spring'
        elif date.month in [6, 7, 8]:  # Summer
            season = 'Summer'
        else:  # Autumn (9, 10, 11)
            season = 'Autumn'
        
        # Group data by season
        seasonal_data[season].append(value)
    
    # Compute statistics for each season
    seasonal_stats = {}
    for season, values in seasonal_data.items():
        values = np.array(values)
        seasonal_stats[season] = {
            'average': np.mean(values),
            'median': np.median(values),
            'max': np.max(values),
            'std_dev': np.std(values),
            'variance': np.var(values),
            'range': np.ptp(values)  # range = max - min
        }
    
    # Sort seasons for consistent plotting order
    sorted_seasons = ['Winter', 'Spring', 'Summer', 'Autumn']
    sorted_stats = {season: seasonal_stats[season] for season in sorted_seasons}
    
    return sorted_seasons, sorted_stats

# Calculate seasonal statistics for Athens and Thessaloniki
seasons_athens, u100_stats_athens = compute_seasonal_stats(dates, u100_athens)
seasons_thessaloniki, u100_stats_thessaloniki = compute_seasonal_stats(dates, u100_thessaloniki)
seasons_athens, v100_stats_athens = compute_seasonal_stats(dates, v100_athens)
seasons_thessaloniki, v100_stats_thessaloniki = compute_seasonal_stats(dates, v100_thessaloniki)
seasons_athens, wind_speed_stats_athens = compute_seasonal_stats(dates, wind_speed_athens)
seasons_thessaloniki, wind_speed_stats_thessaloniki = compute_seasonal_stats(dates, wind_speed_thessaloniki)

# Set up the subplots (1 row, 3 columns)
fig, axs = plt.subplots(1, 3, figsize=(18, 6))

# Set the width for the bars
bar_width = 0.3

# Set positions for Athens and Thessaloniki bars
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

# Get seasonal averages for u100, v100, and wind speed for plotting
u100_seasonal_athens = [u100_stats_athens[season]['average'] for season in seasons_athens]
u100_seasonal_thessaloniki = [u100_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

v100_seasonal_athens = [v100_stats_athens[season]['average'] for season in seasons_athens]
v100_seasonal_thessaloniki = [v100_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

wind_speed_seasonal_athens = [wind_speed_stats_athens[season]['average'] for season in seasons_athens]
wind_speed_seasonal_thessaloniki = [wind_speed_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Plot u100
axs[0].bar(positions_athens, u100_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')
axs[0].bar(positions_thessaloniki, u100_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')
axs[0].set_xticks(np.arange(len(seasons_athens)))
axs[0].set_xticklabels(seasons_athens)
axs[0].set_xlabel('Season')
axs[0].set_ylabel('Average U100 (m/s)')
axs[0].set_title('Seasonal Avg U100 for Athens and Thessaloniki')
axs[0].legend()

# Plot v100
axs[1].bar(positions_athens, v100_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')
axs[1].bar(positions_thessaloniki, v100_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')
axs[1].set_xticks(np.arange(len(seasons_athens)))
axs[1].set_xticklabels(seasons_athens)
axs[1].set_xlabel('Season')
axs[1].set_ylabel('Average V100 (m/s)')
axs[1].set_title('Seasonal Avg V100 for Athens and Thessaloniki')
axs[1].legend()

# Plot wind speed
axs[2].bar(positions_athens, wind_speed_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')
axs[2].bar(positions_thessaloniki, wind_speed_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')
axs[2].set_xticks(np.arange(len(seasons_athens)))
axs[2].set_xticklabels(seasons_athens)
axs[2].set_xlabel('Season')
axs[2].set_ylabel('Average Wind Speed (m/s)')
axs[2].set_title('Seasonal Avg Wind Speed for Athens and Thessaloniki')
axs[2].legend()

# Adjust layout for better spacing
plt.tight_layout()

# Display the plots
plt.show()

```


    
![png](output_0_0.png)
    



```python

```
