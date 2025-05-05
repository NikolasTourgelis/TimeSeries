```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from collections import defaultdict

# File paths for the two NetCDF files
file1 = 'Athina.nc'
file2 = 'Thessaloniki.nc'

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

# Extract the 't2m' variable (3D data: time, lat, lon)
t2m_2 = ds2.variables['t2m'][:]

# Extract time series for Athens and Thessaloniki
t2m_athens = t2m_2[:, athens_lat_idx, athens_lon_idx]
t2m_thessaloniki = t2m_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

# Extract the time variable and convert it to datetime objects
time = ds2.variables['time'][:]
time_units = ds2.variables['time'].units

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
        seasonal_stats[season] = {'average': np.mean(values)}
    
    # Sort seasons for consistent plotting order
    sorted_seasons = ['Winter', 'Spring', 'Summer', 'Autumn']
    sorted_stats = {season: seasonal_stats[season] for season in sorted_seasons}
    
    return sorted_seasons, sorted_stats

# Convert Kelvin to Celsius for plotting
def kelvin_to_celsius(data):
    return data - 273.15

# Calculate monthly averages for Athens and Thessaloniki
months_athens, t2m_monthly_athens = compute_monthly_average(dates, kelvin_to_celsius(t2m_athens))
months_thessaloniki, t2m_monthly_thessaloniki = compute_monthly_average(dates, kelvin_to_celsius(t2m_thessaloniki))

# Calculate seasonal statistics for Athens and Thessaloniki
seasons_athens, t2m_stats_athens = compute_seasonal_stats(dates, kelvin_to_celsius(t2m_athens))
seasons_thessaloniki, t2m_stats_thessaloniki = compute_seasonal_stats(dates, kelvin_to_celsius(t2m_thessaloniki))

# Get seasonal averages for plotting
t2m_seasonal_athens = [t2m_stats_athens[season]['average'] for season in seasons_athens]
t2m_seasonal_thessaloniki = [t2m_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(24, 6))

# Subplot 1: Monthly Average Time Series
axes[0].plot(months_athens, t2m_monthly_athens, label='Athens (Athina)', color='blue')
axes[0].plot(months_thessaloniki, t2m_monthly_thessaloniki, label='Thessaloniki', color='red')
axes[0].xaxis.set_major_locator(mdates.MonthLocator())
axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
axes[0].set_xlabel('Month')
axes[0].set_ylabel('Monthly Avg T2M (2-meter Temperature) (°C)')
axes[0].set_title('Monthly Avg T2M Time Series')
axes[0].legend()
axes[0].tick_params(axis='x', rotation=45)

# Subplot 2: Seasonal Statistics (Zoomed In)
bar_width = 0.3
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2
axes[1].bar(positions_athens, t2m_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')
axes[1].bar(positions_thessaloniki, t2m_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')
axes[1].set_xticks(np.arange(len(seasons_athens)))
axes[1].set_xticklabels(seasons_athens)
axes[1].set_xlabel('Season')
axes[1].set_ylabel('Average T2M (2-meter Temperature) (°C)')
axes[1].set_title('Seasonal Avg T2M ')
axes[1].legend()

# Zoom y-axis: Adjust the y-limits to zoom in
min_temp = min(min(t2m_seasonal_athens), min(t2m_seasonal_thessaloniki))
max_temp = max(max(t2m_seasonal_athens), max(t2m_seasonal_thessaloniki))
axes[1].set_ylim(min_temp - 1, max_temp + 1)  # Zoom in by tightening the range

# Adjust layout and show the plot
plt.tight_layout()
plt.show()

```


    
![png](output_0_0.png)
    



```python

```
