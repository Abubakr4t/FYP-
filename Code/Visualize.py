import xarray as xr

data = xr.open_dataset("data_stream-oper_stepType-accum.nc")
print("Dataset loaded successfully")
print(data)
