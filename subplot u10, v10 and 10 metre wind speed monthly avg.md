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

# Extract the 'u10' and 'v10' variables (3D data: time, lat, lon)
u10_2 = ds2.variables['u10'][:]
v10_2 = ds2.variables['v10'][:]

# Extract time series for Athens and Thessaloniki for u10 and v10
u10_athens = u10_2[:, athens_lat_idx, athens_lon_idx]
v10_athens = v10_2[:, athens_lat_idx, athens_lon_idx]

u10_thessaloniki = u10_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]
v10_thessaloniki = v10_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

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

# Calculate wind speed (sqrt(u10^2 + v10^2)) for Athens and Thessaloniki
wind_speed_athens = np.sqrt(u10_athens**2 + v10_athens**2)
wind_speed_thessaloniki = np.sqrt(u10_thessaloniki**2 + v10_thessaloniki**2)

# Calculate monthly averages for Athens and Thessaloniki
months_athens, u10_monthly_athens = compute_monthly_average(dates, u10_athens)
months_thessaloniki, u10_monthly_thessaloniki = compute_monthly_average(dates, u10_thessaloniki)

months_athens, v10_monthly_athens = compute_monthly_average(dates, v10_athens)
months_thessaloniki, v10_monthly_thessaloniki = compute_monthly_average(dates, v10_thessaloniki)

months_athens, wind_speed_monthly_athens = compute_monthly_average(dates, wind_speed_athens)
months_thessaloniki, wind_speed_monthly_thessaloniki = compute_monthly_average(dates, wind_speed_thessaloniki)

# Create subplots for u10, v10, and wind speed
fig, axs = plt.subplots(3, 1, figsize=(12, 18), sharex=True)

# Plot u10 data
axs[0].plot(months_athens, u10_monthly_athens, label='Monthly Avg U10 - Athens (Athina)', color='blue')
axs[0].plot(months_thessaloniki, u10_monthly_thessaloniki, label='Monthly Avg U10 - Thessaloniki', color='red')
axs[0].set_ylabel('Monthly Avg U10 (10-meter U wind Component) (m/s)')
axs[0].set_title('Monthly Avg U10 (10-meter U Wind Component) Time Series for Athens and Thessaloniki')
axs[0].legend()
axs[0].tick_params(axis='x', rotation=45)

# Plot v10 data
axs[1].plot(months_athens, v10_monthly_athens, label='Monthly Avg V10 - Athens (Athina)', color='blue')
axs[1].plot(months_thessaloniki, v10_monthly_thessaloniki, label='Monthly Avg V10 - Thessaloniki', color='red')
axs[1].set_ylabel('Monthly Avg V10 (10-meter V wind Component) (m/s)')
axs[1].set_title('Monthly Avg V10 (10-meter V Wind Component) Time Series for Athens and Thessaloniki')
axs[1].legend()
axs[1].tick_params(axis='x', rotation=45)

# Plot wind speed data
axs[2].plot(months_athens, wind_speed_monthly_athens, label='Monthly Avg Wind Speed - Athens (Athina)', color='blue')
axs[2].plot(months_thessaloniki, wind_speed_monthly_thessaloniki, label='Monthly Avg Wind Speed - Thessaloniki', color='red')
axs[2].set_xlabel('Month')
axs[2].set_ylabel('Monthly Avg Wind Speed (m/s)')
axs[2].set_title('Monthly Avg Wind Speed (10-meter) Time Series for Athens and Thessaloniki')
axs[2].legend()
axs[2].tick_params(axis='x', rotation=45)

# Format the x-axis to show months
for ax in axs:
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))

plt.tight_layout()
plt.show()

```


    
![png](output_0_0.png)
    



```python

```
