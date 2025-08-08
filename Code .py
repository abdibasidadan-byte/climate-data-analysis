import xarray as xr

file_path = r"C:\Users\hp\OneDrive\Desktop\T2m.nc"
ds = xr.open_dataset(file_path)

print(ds)
print(ds['t2m'].coords)
print(ds['t2m'].coords['valid_time'].values[:5])
print(ds['t2m'].coords['valid_time'].dtype)

jjas_2023 = ds.sel(valid_time=(
    (ds['valid_time'].dt.year == 2023) & 
    (ds['valid_time'].dt.month.isin([6, 7, 8, 9]))
))

t2m_mean_jjas_2023 = jjas_2023['t2m'].mean(dim='valid_time')
t2m_mean_jjas_2023_celsius = t2m_mean_jjas_2023 - 273.15
t2m_mean_jjas_2023_celsius.name = 't2m_mean_jjas_2023_celsius'


import os
output_path = r"C:\Temp\T2m_moyenne_JJAS_2023_celsius.nc"
os.makedirs(os.path.dirname(output_path), exist_ok=True)
t2m_mean_jjas_2023_celsius.to_netcdf(output_path)


import matplotlib.pyplot as plt
t2m_mean_jjas_2023_celsius.plot(cmap='coolwarm')
plt.title("Température moyenne 2m – JJAS 2023 (°C)")
plt.show()

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

t2m_mean_jjas_2023_celsius.plot(
    ax=ax,
    transform=ccrs.PlateCarree(),
    cmap='coolwarm',
    cbar_kwargs={'label': 'Temperature (°C)'},
    add_colorbar=True,
    add_labels=False
)

ax.coastlines()
ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='black')

ax.set_extent([
    float(t2m_mean_jjas_2023_celsius.longitude.min()),
    float(t2m_mean_jjas_2023_celsius.longitude.max()),
    float(t2m_mean_jjas_2023_celsius.latitude.min()),
    float(t2m_mean_jjas_2023_celsius.latitude.max())
], crs=ccrs.PlateCarree())

gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linewidth=0.5, color='gray', alpha=0.5)
gl.top_labels = False
gl.right_labels = False

plt.title("Mean 2m Air Temperature – JJAS 2023 (°C)")
plt.show()
