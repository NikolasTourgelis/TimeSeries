```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from collections import defaultdict

# File paths for the two NetCDF files
file1 = 'athensnewtotal.nc'
file2 = 'thessalonikinewtotal.nc'

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

# Extract the 'tp' variable (3D data: time, lat, lon)
tp_2 = ds2.variables['tp'][:]

# Extract time series for Athens and Thessaloniki
tp_athens = tp_2[:, athens_lat_idx, athens_lon_idx]
tp_thessaloniki = tp_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

# Extract the time variable and convert it to datetime objects
time = ds2.variables['valid_time'][:]
time_units = ds2.variables['valid_time'].units

# Convert time to datetime using num2date from netCDF4
dates = nc.num2date(time, units=time_units)

# Convert cftime objects to standard datetime objects
dates = [datetime(dt.year, dt.month, dt.day) for dt in dates]

# Close the NetCDF file
ds2.close()

# Function to compute monthly averages
def compute_monthly_average(dates, data):
    monthly_data = defaultdict(list)
    for date, value in zip(dates, data):
        month = date.replace(day=1)  # Set day to 1 to group by month
        monthly_data[month].append(value)
    
    monthly_avg = {month: np.mean(values) for month, values in monthly_data.items()}
    sorted_months = sorted(monthly_avg.keys())
    sorted_averages = [monthly_avg[month] for month in sorted_months]
    
    return sorted_months, sorted_averages

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

# Calculate monthly averages for Athens and Thessaloniki
months_athens, tp_monthly_athens = compute_monthly_average(dates, tp_athens)
months_thessaloniki, tp_monthly_thessaloniki = compute_monthly_average(dates, tp_thessaloniki)

# Calculate seasonal statistics for Athens and Thessaloniki
seasons_athens, tp_stats_athens = compute_seasonal_stats(dates, tp_athens)
seasons_thessaloniki, tp_stats_thessaloniki = compute_seasonal_stats(dates, tp_thessaloniki)

# Prepare seasonal averages for bar plot
tp_seasonal_athens = [tp_stats_athens[season]['average'] for season in seasons_athens]
tp_seasonal_thessaloniki = [tp_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Plot subplots
fig, axs = plt.subplots(1, 2, figsize=(15, 6))

# Monthly averages subplot
axs[0].plot(months_athens, tp_monthly_athens, label='Athens', color='blue')
axs[0].plot(months_thessaloniki, tp_monthly_thessaloniki, label='Thessaloniki', color='red')
axs[0].xaxis.set_major_locator(mdates.MonthLocator())
axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
axs[0].set_title('Monthly Avg Total Precipitation')
axs[0].set_xlabel('Month')
axs[0].set_ylabel('Monthly Avg TP (m)')
axs[0].legend()
axs[0].tick_params(axis='x', rotation=45)

# Seasonal averages subplot
bar_width = 0.3
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2
axs[1].bar(positions_athens, tp_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens')
axs[1].bar(positions_thessaloniki, tp_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')
axs[1].set_xticks(np.arange(len(seasons_athens)))
axs[1].set_xticklabels(seasons_athens)
axs[1].set_title('Seasonal Avg Total Precipitation')
axs[1].set_xlabel('Season')
axs[1].set_ylabel('Average TP (m)')
axs[1].legend()

# Adjust layout
plt.tight_layout()
plt.show()

```


    
![png](output_0_0.png)
    



```python

```
