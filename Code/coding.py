import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# --------------Coding---------------


# Open the NetCDF file
data = xr.open_dataset("data_stream-oper_stepType-accum.nc")
# View basic info
print("Dataset loaded successfully!")
print(data)

time = data["valid_time"].values

# Compute mean across all lat/lon for each time
mean_tp = data["tp"].mean(dim=["latitude", "longitude"])
mean_ssrd = data["ssrd"].mean(dim=["latitude", "longitude"])

# Plot
plt.figure(figsize=(12, 6))
plt.plot(time, mean_tp, label="Mean Total Precipitation", color="blue")
plt.plot(time, mean_ssrd, label="Mean Solar Radiation Downward", color="orange")
plt.xlabel("Time")
plt.ylabel("Mean Value")
plt.title("Average Precipitation & Solar Radiation Over Time (Pakistan)")
plt.legend()
plt.tight_layout()
plt.show()


print("Available variables:", list(data.data_vars))

# If u10, v10, and t2m exist
if {'u10', 'v10', 't2m'}.issubset(data.data_vars):
    # Compute wind speed
    wind_speed = np.sqrt(data['u10']**2 + data['v10']**2)
    mean_wind = wind_speed.mean(dim=["latitude", "longitude"])
    mean_temp = data['t2m'].mean(
        dim=["latitude", "longitude"]) - 273.15  # Convert K to °C

    # Plot
    plt.figure(figsize=(12, 6))
    plt.plot(data["valid_time"], mean_temp,
             label="Mean 2m Temperature (°C)", color="red")
    plt.plot(data["valid_time"], mean_wind,
             label="Mean Wind Speed (m/s)", color="green")
    plt.xlabel("Time")
    plt.ylabel("Mean Value")
    plt.title("Average Temperature & Wind Speed Over Time (Pakistan)")
    plt.legend()
    plt.tight_layout()
    plt.show()
else:
    print("The dataset does not contain u10, v10, or t2m. Please re-download with these variables.")


# Load datasets
accum = xr.open_dataset(
    "data_stream-oper_stepType-accum.nc").drop_vars("expver", errors="ignore")
avg = xr.open_dataset(
    "data_stream-oper_stepType-avg.nc").drop_vars("expver", errors="ignore")
instant = xr.open_dataset(
    "data_stream-oper_stepType-instant.nc").drop_vars("expver", errors="ignore")

# Merge safely
merged = xr.merge([accum, avg, instant])

print("Datasets merged successfully!")
print("Variables in merged dataset:", list(merged.data_vars.keys()))

# Save merged dataset
merged.to_netcdf("merged_era5_pakistan.nc")

# Compute mean variables
mean_tp = merged['tp'].mean(dim=["latitude", "longitude"])
mean_t2m = merged['t2m'].mean(dim=["latitude", "longitude"])
wind_speed = ((merged['u10']**2 + merged['v10']**2) **
              0.5).mean(dim=["latitude", "longitude"])

# Plot
plt.figure(figsize=(12, 6))
plt.plot(merged.valid_time, mean_tp,
         label='Mean Total Precipitation (m)', color='blue')
plt.plot(merged.valid_time, mean_t2m - 273.15,
         label='Mean 2m Temperature (°C)', color='red')
plt.plot(merged.valid_time, wind_speed,
         label='Mean Wind Speed (m/s)', color='green')
plt.xlabel('Time')
plt.ylabel('Value')
plt.title('Average ERA5 Climate Variables Over Time (Pakistan)')
plt.legend()
plt.show()

# Load merged dataset
merged = xr.open_dataset("merged_era5_pakistan.nc")

# Compute average values (spatial mean)
mean_tp = merged['tp'].mean(dim=["latitude", "longitude"])
mean_t2m = merged['t2m'].mean(
    dim=["latitude", "longitude"]) - 273.15  # Convert K → °C
mean_ssrd = merged['ssrd'].mean(dim=["latitude", "longitude"])
wind_speed = np.sqrt(merged['u10']**2 + merged['v10']
                     ** 2).mean(dim=["latitude", "longitude"])

# Create plot
fig, ax1 = plt.subplots(figsize=(12, 6))

# Primary Y-axis
ax1.plot(merged.valid_time, mean_t2m, color='r', label='2m Temperature (°C)')
ax1.plot(merged.valid_time, mean_tp, color='black', label='Precipitation (m)')
ax1.plot(merged.valid_time, wind_speed, color='blue', label='Wind Speed (m/s)')
ax1.set_xlabel('Time')
ax1.set_ylabel('Temperature / Precipitation / Wind Speed')
ax1.tick_params(axis='y')

# Secondary Y-axis (for solar radiation)
ax2 = ax1.twinx()
ax2.plot(merged.valid_time, mean_ssrd, color='gray',
         label='Solar Radiation (J/m²)', alpha=0.6)
ax2.set_ylabel('Solar Radiation (J/m²)')

# Combine legends
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='upper right')

plt.title('Weather Trends Over Time (ERA5 Pakistan)')
plt.tight_layout()
plt.show()
