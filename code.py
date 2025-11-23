import xarray as xr
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split


# Load merged file
ds = xr.open_dataset("merged_era5_pakistan.nc")

# Select key weather variables
df = ds[["t2m", "tp", "u10", "v10", "ssrd", "msl"]].mean(
    dim=["latitude", "longitude"]).to_dataframe().reset_index()

# Convert temperature from Kelvin to Celsius
df["t2m"] = df["t2m"] - 273.15

# Calculate wind speed
df["wind_speed"] = np.sqrt(df["u10"]**2 + df["v10"]**2)

# Drop missing values (if any)
df = df.dropna()

print(df.head())


# Add month and season (helps model learn periodic patterns)
df["month"] = pd.to_datetime(df["valid_time"]).dt.month
df["day"] = pd.to_datetime(df["valid_time"]).dt.day

# Seasonal sin/cos encoding (helps capture cyclic trends)
df["month_sin"] = np.sin(2 * np.pi * df["month"]/12)
df["month_cos"] = np.cos(2 * np.pi * df["month"]/12)

# Lag features – yesterday's values to predict today's
df["t2m_lag1"] = df["t2m"].shift(1)
df["tp_lag1"] = df["tp"].shift(1)

# Drop missing rows after shift
df = df.dropna()
print(df.head())

# Define features and target
X = df[["tp", "ssrd", "wind_speed", "msl",
        "month_sin", "month_cos", "t2m_lag1", "tp_lag1"]]
y = df["t2m"]

scaler = MinMaxScaler()
X_scaled = scaler.fit_transform(X)
X_train, X_test, y_train, y_test = train_test_split(
    X_scaled, y, test_size=0.2, shuffle=False)
print(df.head())