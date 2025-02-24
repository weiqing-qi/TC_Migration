# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 15:54:52 2024

@author: 29585
"""
#%% 20 年风场对比
from Functions import *
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import seaborn as sns
from datetime import datetime
from matplotlib.colors import LinearSegmentedColormap
import joblib

region_uv_mean = joblib.load('region_uv_mean.joblib') 
region_uv_meanNI = joblib.load('region_uv_meanNI.joblib') 
# 提取的范围 0-360
regions = {
    'WP': [101, 170, 5, 60],     # WP: [lon_min, lon_max, lat_min, lat_max]
    'EP': [200, 259, 5, 60],     # EP (将-160至-101度转换为200至259度)
    'NA': [261, 320, 5, 60],     # NA (将-99至-40度转换为261至320度)
    'NI_L': [40, 99, 5, 40],     # NI L
    # 'NI_R': [78, 99, 5, 40],     # NI R
    'AU_L': [90, 134, -50, -5],  # AU L
    'AU_R': [136, 200, -50, -5], # AU R
    'SA': [30, 70, -40, -5],     # SA
}

# 生成静态的经纬度网格 这里必须是-180-180 然后把左边一半+360得到lon_grid_360这样因为是这样保存的
lats05 = np.linspace(90 - 0.5 / 2, -90 + 0.5 / 2, 360)
lons05 = np.linspace(-180 + 0.5 / 2, 180 - 0.5 / 2, 720)
lon_grid, lat_grid = np.meshgrid(lons05, lats05)
# 处理经度从-180到180到0-360的转换
lon_grid_360 = np.where(lon_grid < 0, lon_grid + 360, lon_grid)

region_uv_diff = {region: {'U': None, 'V': None} for region in regions.keys()}
for region in regions.keys():
    if region == 'NI_L':
        uv_mean=region_uv_meanNI
    else:
        uv_mean=region_uv_mean
    valid_U_1980_1999 = [u for u in uv_mean[region]['U'][:20] if not np.isnan(u).all()]
    U_1980_1999 = np.nanmean(valid_U_1980_1999, axis=0) if valid_U_1980_1999 else np.nan
    valid_V_1980_1999 = [v for v in uv_mean[region]['V'][:20] if not np.isnan(v).all()]
    V_1980_1999 = np.nanmean(valid_V_1980_1999, axis=0) if valid_V_1980_1999 else np.nan
    valid_U_2000_2020 = [u for u in uv_mean[region]['U'][20:] if not np.isnan(u).all()]
    U_2000_2020 = np.nanmean(valid_U_2000_2020, axis=0) if valid_U_2000_2020 else np.nan
    valid_V_2000_2020 = [v for v in uv_mean[region]['V'][20:] if not np.isnan(v).all()]
    V_2000_2020 = np.nanmean(valid_V_2000_2020, axis=0) if valid_V_2000_2020 else np.nan
    
    region_uv_diff[region]['U'] = U_2000_2020 - U_1980_1999
    region_uv_diff[region]['V'] = V_2000_2020 - V_1980_1999

# 绘制风场变化图
sns.set_style('white')
plt.figure(figsize=(12, 6), dpi=400)
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))  # 使用 PlateCarree 投影绘制地图
# ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
ax.set_global()
ax.coastlines(zorder=3, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=0.3, zorder=3)
ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=1)
# ax.set_extent([0, 360, -90, 90], crs=ccrs.PlateCarree())

# 添加经纬度数字标签并禁用网格线
gridlines = ax.gridlines(draw_labels=True, linewidth=0, color='dimgray')
gridlines.bottom_labels = False

for region in regions:
    lon_min, lon_max, lat_min, lat_max = regions[region]
    rect = plt.Rectangle((lon_min, lat_min), lon_max - lon_min, lat_max - lat_min,
                          linewidth=0.7, edgecolor='black', facecolor='none', transform=ccrs.PlateCarree(), zorder=4)
    ax.add_patch(rect)
    
# 创建空的二维全球风场变化数组
u_diff_global = np.full(lon_grid_360.shape, np.nan)
v_diff_global = np.full(lon_grid_360.shape, np.nan)

# 将每个区域的风场变化填回到全球网格
for region in regions.keys():
    lon_min, lon_max, lat_min, lat_max = regions[region]
    if region == 'NI_L':
        lon_max+=1
    lon_mask = (lon_grid_360 >= lon_min) & (lon_grid_360 <= lon_max)
    lat_mask = (lat_grid >= lat_min) & (lat_grid <= lat_max)
    region_mask = lon_mask & lat_mask
    
    u_diff_global[region_mask] = region_uv_diff[region]['U']
    v_diff_global[region_mask] = region_uv_diff[region]['V']
    
    if region == 'NI_L' or region == 'NI_R':
        u_diff_global[region_mask] = u_diff_global[region_mask]*(1)
        v_diff_global[region_mask] = v_diff_global[region_mask]*(1)
        
# 绘制风场箭头，仅显示在方框内的数据
quiver_step = 8
q = None
# 绘制风场箭头
mask = ~np.isnan(u_diff_global) & ~np.isnan(v_diff_global)
lon_masked = np.ma.masked_where(~mask, lon_grid_360)
lat_masked = np.ma.masked_where(~mask, lat_grid)
u_masked = np.ma.masked_where(~mask, u_diff_global)
v_masked = np.ma.masked_where(~mask, v_diff_global)
q = ax.quiver(lon_masked[5::quiver_step, 5::quiver_step], lat_masked[5::quiver_step, 5::quiver_step],
              u_masked[5::quiver_step, 5::quiver_step], v_masked[5::quiver_step, 5::quiver_step],
              scale=100, width=0.0015, alpha=0.8, pivot='middle', color='black', transform=ccrs.PlateCarree(), zorder=3)#'#1f4e78'

# 添加图例到图的右上角，并加上背景框
rect = plt.Rectangle((0.015, 0.91), 0.1, 0.06, transform=ax.transAxes, facecolor='white', alpha=0.9, zorder=3, edgecolor='black', linewidth=0.7)
ax.add_patch(rect)
ax.quiverkey(q, 0.05, 0.94, 2, '2 m/s', labelpos='E', coordinates='axes', color='black', zorder=6)

#绘制海温变化和显著性
AgeDiff = read_arcgis_txt('log10PopXPreRate.txt')['data']

# AgeDiff = np.roll(AgeDiff, shift=AgeDiff.shape[1] // 2, axis=1) # 转换为0-360 再绘制
# AgeDiffSIG = np.roll(AgeDiffSIG, shift=AgeDiffSIG.shape[1] // 2, axis=1) # 转换为0-360 再绘制

# 把 AgeDiff 中 -9999 的值 masked 掉，并将 AgeDiffSIG > 1 的位置也在 AgeDiff 中 masked 掉
AgeDiff = np.ma.masked_where(AgeDiff == -9999, AgeDiff)
AgeDiff = np.ma.masked_where(AgeDiff == 0, AgeDiff)

# AgeDiff = AgeDiff * 0.7 # SST

# 生成静态的经纬度网格
lats05 = np.linspace(90 - 0.5 / 2, -90 + 0.5 / 2, 360)
lons05_0360 = np.linspace(0 + 0.5 / 2, 360 - 0.5 / 2, 720)
lon_grid_new0360, lat_grid = np.meshgrid(lons05_0360, lats05) #上面不能动，只能改这里

# 定义 RGB 颜色，使用 (R, G, B)，范围是 0-255 的整数
colors = [
    (211, 211, 211), #1
    (106, 154, 203), #2 (60, 100, 210)
    (123, 164, 208), #3 (50, 130, 230)
    (165, 212, 163), #4
    (250, 235, 190), #5
    (254, 120, 8), #6
    (229, 85, 0),  #7
    (229, 1, 0),   #8
    (190, 0, 0),
    (170, 0, 0),
]

# 将颜色值转换为 0-1 的范围
colors = [(r/255, g/255, b/255) for r, g, b in colors]

# 创建自定义的渐变颜色图
cmap = LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)

# 绘制 AgeDiff 的变化
contour = ax.contourf(lon_grid_new0360, lat_grid, AgeDiff, levels=np.linspace(1, 11, 21), cmap=cmap, extend='both', zorder=1, transform=ccrs.PlateCarree())
cbar = plt.colorbar(contour, orientation='horizontal', pad=0.02, aspect=40, shrink=0.78)
# cbar.set_label('Difference (mm/d)')



# # 标题
# plt.title('TC-related Steering Wind and Precipitation Change (2000-2020 vs 1980-1999)')

# 显示图像
plt.show()


