```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from collections import defaultdict

# File paths for the two NetCDF files
file1 = 'athinacloudcover.nc'
file2 = 'thessalonikicloudcover.nc'

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

# Extract the 'tcc' variable (3D data: time, lat, lon)
tcc_2 = ds2.variables['tcc'][:]

# Extract time series for Athens and Thessaloniki
tcc_athens = tcc_2[:, athens_lat_idx, athens_lon_idx]
tcc_thessaloniki = tcc_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

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
            'average': np.mean(values)
        }
    
    # Sort seasons for consistent plotting order
    sorted_seasons = ['Winter', 'Spring', 'Summer', 'Autumn']
    sorted_stats = {season: seasonal_stats[season] for season in sorted_seasons}
    
    return sorted_seasons, sorted_stats

# Calculate monthly averages for Athens and Thessaloniki
months_athens, tcc_monthly_athens = compute_monthly_average(dates, tcc_athens)
months_thessaloniki, tcc_monthly_thessaloniki = compute_monthly_average(dates, tcc_thessaloniki)

# Calculate seasonal statistics for Athens and Thessaloniki
seasons_athens, tcc_stats_athens = compute_seasonal_stats(dates, tcc_athens)
seasons_thessaloniki, tcc_stats_thessaloniki = compute_seasonal_stats(dates, tcc_thessaloniki)

# Prepare seasonal averages for plotting
tcc_seasonal_athens = [tcc_stats_athens[season]['average'] * 100 for season in seasons_athens]  # Convert to %
tcc_seasonal_thessaloniki = [tcc_stats_thessaloniki[season]['average'] * 100 for season in seasons_thessaloniki]  # Convert to %

# Convert monthly averages to percentages
tcc_monthly_athens = [value * 100 for value in tcc_monthly_athens]
tcc_monthly_thessaloniki = [value * 100 for value in tcc_monthly_thessaloniki]

# Subplot to display both analyses
fig, axs = plt.subplots(1, 2, figsize=(16, 6))

# First subplot: Monthly averages time series
axs[0].plot(months_athens, tcc_monthly_athens, label='Athens (Athina)', color='blue')
axs[0].plot(months_thessaloniki, tcc_monthly_thessaloniki, label='Thessaloniki', color='red')
axs[0].xaxis.set_major_locator(mdates.MonthLocator())
axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
axs[0].set_xlabel('Month')
axs[0].set_ylabel('Monthly Avg TCC (%)')
axs[0].set_title('Monthly Avg Total Cloud Cover Time Series')
axs[0].legend()
axs[0].tick_params(axis='x', rotation=45)

# Second subplot: Seasonal averages bar plot
bar_width = 0.3
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

axs[1].bar(positions_athens, tcc_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')
axs[1].bar(positions_thessaloniki, tcc_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')
axs[1].set_xticks(np.arange(len(seasons_athens)))
axs[1].set_xticklabels(seasons_athens)
axs[1].set_xlabel('Season')
axs[1].set_ylabel('Seasonal Avg TCC (%)')
axs[1].set_title('Seasonal Avg Total Cloud Cover')
axs[1].legend()

# Adjust layout and show plot
plt.tight_layout()
plt.show()


```


    
![png](output_0_0.png)
    



```python

```
