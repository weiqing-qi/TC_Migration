# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 15:13:09 2024
绘制风场图形
@author: 29585
"""
from Functions import *
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import seaborn as sns


import pandas as pd
import rasterio
from rasterio.transform import from_origin
from osgeo import gdal
from datetime import datetime

hPa = 500
nc_path_U = f'E:\\DATA\\1.Reanalysis\\EAR5\\ERA5-UCW\\{hPa}\\'
nc_file_name_U = f'ERA5_U_Component_Of_Wind_on_{hPa}hpa_levels_daily'
nc_path_V = f'E:\\DATA\\1.Reanalysis\\EAR5\\ERA5-VCW\\{hPa}\\'
nc_file_name_V = f'ERA5_V_Component_Of_Wind_on_{hPa}hpa_levels_daily'

year = 1982
mon = 1
day = 1
year_str = str(year)
mon_str = str(mon).zfill(2)
day_str = str(day).zfill(2)

# ----读取数据----
mon_UCW_DS = Dataset("".join([nc_path_U, nc_file_name_U, year_str, mon_str, '.nc']), 'r')
mon_VCW_DS = Dataset("".join([nc_path_V, nc_file_name_V, year_str, mon_str, '.nc']), 'r')
mon_UCW = mon_UCW_DS.variables['u'][:]
mon_VCW = mon_VCW_DS.variables['v'][:]

if 1980 <= year <= 2000:
    lats = mon_UCW_DS.variables['latitude'][:] # UV这lonlat是一样的
    lons = mon_VCW_DS.variables['longitude'][:]
    mon_UCW = np.squeeze(mon_UCW)
    mon_VCW = np.squeeze(mon_VCW)
    
elif 2001 <= year <= 2020:
    lats = mon_UCW_DS.variables['lat'][:] # UV这lonlat是一样的
    lons = mon_VCW_DS.variables['lon'][:]
    
mon_UCW_DS.close()
mon_VCW_DS.close()
# ----数据处理----
# 可以保留多的一列后面采样会处理
day_UCW_Ori = mon_UCW[day-1]
day_VCW_Ori = mon_VCW[day-1]
# 重采样到 0.5 同时自动处理extent到-180-180问题
day_UCW05 = resample_model_grid(lats, lons, day_UCW_Ori, 360, 720, np.nan, method = 'nearest')
day_VCW05 = resample_model_grid(lats, lons, day_VCW_Ori, 360, 720, np.nan, method = 'nearest')

lats05 = np.linspace(90 - 0.5/2, -90 + 0.5/2, 360)
lons05 = np.linspace(-180 + 0.5/2, 180 - 0.5/2, 720)
# 创建经纬度网格
lon_grid, lat_grid = np.meshgrid(lons05, lats05)

#%% 只画一张全球的日风场

sns.set_style('white') # 禁用 seaborn 的默认样式以避免网格线影响
plt.figure(figsize=(12, 6), dpi = 400)
ax = plt.axes(projection=ccrs.PlateCarree())  # 使用 PlateCarree 投影绘制地图
ax.coastlines(zorder=1)  # 添加海岸线
ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=1)  # 可选：添加国家边界线
ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=1)
ax.set_extent([-180, 180, -90, 90.], crs=ccrs.PlateCarree()) # 设置地图的显示范围
gridlines = ax.gridlines(draw_labels=True, linewidth=0)  # 仅显示标签，不显示网格线
# gridlines.right_labels = False
# gridlines.top_labels = False

# 绘制风场箭头
quiver_step = 10  # 控制箭头的间隔
q = plt.quiver(lon_grid[::quiver_step, ::quiver_step], lat_grid[::quiver_step, ::quiver_step],
               day_UCW05[::quiver_step, ::quiver_step], day_VCW05[::quiver_step, ::quiver_step],
               scale=1900, width=0.0015, pivot='middle', color='#1f4e78', transform=ccrs.PlateCarree())
# 图例
rect = plt.Rectangle((0.015, 0.91), 0.1, 0.06, transform=plt.gca().transAxes, color='white', alpha=0.90, zorder=1)
plt.gca().add_patch(rect)
plt.quiverkey(q, 0.04, 0.94, 30, '30 m/s', labelpos='E', coordinates='axes', color='black', zorder=2)

# 标题
plt.title(f'Wind Field at {hPa} hPa on {year_str}-{mon_str}-{day_str}')
plt.show()

#%% 指定区域
sns.set_style('white') # 禁用 seaborn 的默认样式以避免网格线影响
plt.figure(figsize=(12, 6), dpi = 400)
ax = plt.axes(projection=ccrs.PlateCarree())  # 使用 PlateCarree 投影绘制地图
ax.coastlines(zorder=1)  # 添加海岸线
ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=1)  # 可选：添加国家边界线
ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=1)
ax.set_extent([-180, 180, -90, 90.], crs=ccrs.PlateCarree()) # 设置地图的显示范围
gridlines = ax.gridlines(draw_labels=True, linewidth=0)  # 仅显示标签，不显示网格线

# 绘制指定经纬度的方框
regions = [
    [100, 170, 5, 60],     # WP: [lon_min, lon_max, lat_min, lat_max]
    [-160, -101, 5, 60],   # EP
    [-99, -40, 5, 60],     # NA
    [40, 76, 5, 40],       # NI L
    [77, 99, 5, 40],       # NI R
    [90, 134, -50, -5],     # AU L
    [136, 180, -50, -5],    # AU R
    [-180, -160, -50, -5],    # AU R
    [30, 70, -40, -5],       # SA
]

for region in regions:
    lon_min, lon_max, lat_min, lat_max = region
    rect = plt.Rectangle((lon_min, lat_min), lon_max - lon_min, lat_max - lat_min,
                         linewidth=1.5, edgecolor='black', facecolor='none', transform=ccrs.PlateCarree(), zorder=4)
    ax.add_patch(rect)

# 绘制风场箭头，仅显示在方框内的数据
quiver_step = 10  # 控制箭头的间隔
q = None  # 初始化 q
for region in regions:
    lon_min, lon_max, lat_min, lat_max = region
    lon_min += 1 # 控制溢出
    lon_max -= 1
    lat_min += 1
    lat_max -= 1
    mask = (lon_grid >= lon_min) & (lon_grid <= lon_max) & (lat_grid >= lat_min) & (lat_grid <= lat_max)
    lon_masked = np.ma.masked_where(~mask, lon_grid)
    lat_masked = np.ma.masked_where(~mask, lat_grid)
    u_masked = np.ma.masked_where(~mask, day_UCW05)
    v_masked = np.ma.masked_where(~mask, day_VCW05)
    q = ax.quiver(lon_masked[::quiver_step, ::quiver_step], lat_masked[::quiver_step, ::quiver_step],
                  u_masked[::quiver_step, ::quiver_step], v_masked[::quiver_step, ::quiver_step],
                  scale=1500, width=0.0015, pivot='middle', color='#1f4e78', transform=ccrs.PlateCarree(), zorder=3)

# 添加图例到图的右上角，并加上背景框
rect = plt.Rectangle((0.015, 0.91), 0.1, 0.06, transform=ax.transAxes, facecolor='white', alpha=0.90, zorder=1, edgecolor='black', linewidth=1.5)
ax.add_patch(rect)
q = ax.quiverkey(q, 0.04, 0.94, 30, '30 m/s', labelpos='E', coordinates='axes', color='black', zorder=2)

# 标题
plt.title(f'Wind Field at {hPa} hPa on {year_str}-{mon_str}-{day_str}')

# 显示图像
plt.show()

#%% 0-360
# 禁用 seaborn 的默认样式以避免网格线影响
sns.set_style('white')
# ----绘制风场图----
plt.figure(figsize=(12, 6), dpi=400)
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))  # 将中心经度设置为180度
ax.set_global()
ax.coastlines(zorder=1)  # 将海岸线的图层顺序设置为较低
ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=1)
ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=1)
ax.set_extent([0, 360, -90, 90], crs=ccrs.PlateCarree())  # 设置地图的显示范围为0-360度

# 添加经纬度数字标签并禁用网格线
gridlines = ax.gridlines(draw_labels=True, linewidth=0)  # 仅显示标签，不显示网格线

# 创建经纬度网格
lon_grid, lat_grid = np.meshgrid(lons05, lats05)
# 将经度从 -180 到 180 转换为 0 到 360
lon_grid = np.where(lon_grid < 0, lon_grid + 360, lon_grid)

# 绘制指定经纬度的方框
regions = [
    [101, 170, 5, 60],     # WP: [lon_min, lon_max, lat_min, lat_max]
    [200, 259, 5, 60],     # EP (将-160至-101度转换为200至259度)
    [261, 320, 5, 60],     # NA (将-99至-40度转换为261至320度)
    [40, 76, 5, 40],       # NI L
    [78, 99, 5, 40],       # NI R
    [90, 134, -50, -5],    # AU L
    [136, 200, -50, -5],   # AU R
    [30, 70, -40, -5],     # SA
]

for region in regions:
    lon_min, lon_max, lat_min, lat_max = region
    rect = plt.Rectangle((lon_min, lat_min), lon_max - lon_min, lat_max - lat_min,
                         linewidth=1.5, edgecolor='black', facecolor='none', transform=ccrs.PlateCarree(), zorder=4)
    ax.add_patch(rect)

# 绘制风场箭头，仅显示在方框内的数据
quiver_step = 10  # 控制箭头的间隔
q = None  # 初始化 q
for region in regions:
    lon_min, lon_max, lat_min, lat_max = region
    lon_min += 1  # 控制溢出
    lon_max -= 1
    lat_min += 1
    lat_max -= 1
    mask = (lon_grid >= lon_min) & (lon_grid <= lon_max) & (lat_grid >= lat_min) & (lat_grid <= lat_max)
    lon_masked = np.ma.masked_where(~mask, lon_grid)
    lat_masked = np.ma.masked_where(~mask, lat_grid)
    u_masked = np.ma.masked_where(~mask, day_UCW05)
    v_masked = np.ma.masked_where(~mask, day_VCW05)
    q = ax.quiver(lon_masked[::quiver_step, ::quiver_step], lat_masked[::quiver_step, ::quiver_step],
                  u_masked[::quiver_step, ::quiver_step], v_masked[::quiver_step, ::quiver_step],
                  scale=1500, width=0.0015, pivot='middle', color='#1f4e78', transform=ccrs.PlateCarree(), zorder=3)

# 添加图例到图的右上角，并加上背景框
rect = plt.Rectangle((0.015, 0.91), 0.1, 0.06, transform=ax.transAxes, facecolor='white', alpha=0.90, zorder=1, edgecolor='black', linewidth=1.5)
ax.add_patch(rect)
ax.quiverkey(q, 0.04, 0.94, 30, '30 m/s', labelpos='E', coordinates='axes', color='black', zorder=2)

# 标题
plt.title(f'Wind Field at {hPa} hPa on {year_str}-{mon_str}-{day_str}')

# 显示图像
plt.show()

#%% 全球路径图
from Functions import *
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.cm as cm
import matplotlib.colors as mcolors

SI, EI = 0, 4459                 #2001-2020 [) 0-BASED
record_len = EI - SI
IBTrACS_Path = r'D:\DATA\Segmentation map tutorial\IBTrACS.since1980.v04r00.240321.nc'
wind_dict, lonlat, times, otherdata = get_IBTrACK(IBTrACS_Path, SI, EI)

Years=[]
for i in range(record_len):
    year = datetime.strptime(''.join(times[i, 0, :]), '%Y-%m-%d %H:%M:%S').year
    Years.append(int(year))
    
# 假设数据在 lon_data 和 lat_data 数组中，形状为 (2000, 360)，填充值为 -9999
lon_data = lonlat['lons'] # 示例数据，随机生成
lat_data = lonlat['lats']     # 示例数据，随机生成

# 创建绘图
plt.figure(figsize=(12, 6), dpi=300)
ax = plt.axes(projection=ccrs.PlateCarree())

# 添加地图要素
ax.set_global()
ax.coastlines()
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, facecolor='lightgray')

# 设置颜色映射
norm = mcolors.Normalize(vmin=min(Years), vmax=max(Years))
cmap = cm.get_cmap('Blues')  # 使用颜色渐变 'viridis'，你也可以选择其他渐变

# 遍历每一行数据绘制路径
for i in range(lon_data.shape[0]):
    # 筛选有效数据（去除填充值 -9999）
    valid_mask = (lon_data[i] != -9999) & (lat_data[i] != -9999)
    valid_lons = lon_data[i][valid_mask]
    valid_lats = lat_data[i][valid_mask]
    
    # 获取对应年份的颜色
    year = Years[i]
    color = cmap(norm(year))
    
    # 绘制路径
    if len(valid_lons) > 1:  # 至少需要两个点才能绘制路径
        ax.plot(valid_lons, valid_lats, transform=ccrs.PlateCarree(), linewidth=0.5, alpha=0.6, color=color)

# 设置标题
plt.title('Global Paths of All Entries (Monochrome Color by Year)')

# 添加颜色条
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation='vertical', shrink=0.7, pad=0.02)
cbar.set_label('Year')

# 显示图像
plt.show()
