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

# 提取的范围 0-360
regions = {
    'WP': [101, 170, 5, 60],     # WP: [lon_min, lon_max, lat_min, lat_max]
    'EP': [200, 259, 5, 60],     # EP (将-160至-101度转换为200至259度)
    'NA': [261, 320, 5, 60],     # NA (将-99至-40度转换为261至320度)
    'NI': [40, 99, 5, 40],     # NI L
    # 'NI_R': [78, 99, 5, 40],     # NI R
    'AU_L': [90, 134, -50, -5],  # AU L
    'AU_R': [136, 200, -50, -5], # AU R
    'SA': [30, 70, -40, -5],     # SA
}

SI, EI = 0, 4459                 #2001-2020 [) 0-BASED
record_len = EI - SI
IBTrACS_Path = r'D:\DATA\Segmentation map tutorial\IBTrACS.since1980.v04r00.240321.nc'
wind_dict, lonlat, times, otherdata = get_IBTrACK(IBTrACS_Path, SI, EI)
    
# 筛选出lat位于44N-44S之内的位置 筛选陆地海洋格点
selected_lats = (lonlat["lats"] <= 44) & (lonlat["lats"] >= -44)

# 筛选出风速值大于34的位置的并集, 然后取44N-44S之内的值
# C5>=137 C4: 113~136 C3: 96~112 C2: 83~95 C1: 64~82 TS: 35~63 Knots=0.514 m/s
selected_winds = np.zeros((4459, 360), dtype=bool)
for key in wind_dict:
    selected_winds |= wind_dict[key] > 34 # 修改get_IBTrACK.wind_dict来选择需要合并的机构 
selected_winds &= selected_lats 

# 找到真值所在行：记录条目
selected_rows = np.where(np.any(selected_winds, axis=1))[0]

time_dict = {region: set() for region in regions.keys()}
for row in selected_rows:
    lat = lonlat["lats"][row]
    lon = lonlat["lons"][row]
    # 处理经度从-180-180到0-360，跳过-9999的值 小于-180的也可以
    lon = np.where(lon == -9999, lon, np.where(lon < 0, lon + 360, lon))
    
    # 对于每个非-9999的经度位置，检测当前行的经纬度是否属于各个分区范围内
    for idx in range(len(lon)):
        if lon[idx] == -9999:
            continue
        for region, bounds in regions.items():
            lon_min, lon_max, lat_min, lat_max = bounds
            if (lon_min <= lon[idx] <= lon_max) and (lat_min <= lat[idx] <= lat_max):
                # 如果在范围内，把对应的时间年月日提取放进字典
                date = datetime.strptime(''.join(times[row, idx, :]), '%Y-%m-%d %H:%M:%S').date()
                time_dict[region].add(date)

# 循环结束对time_dict，每个区域中的年月日信息进行去重
time_dict = {region: sorted(dates) for region, dates in time_dict.items()}

# 计算1980-2020年每个区域逐年的年平均UV风场
hPa = 500
nc_path_U = f'E:\\DATA\\1.Reanalysis\\EAR5\\ERA5-UCW\\{hPa}\\'
nc_file_name_U = f'ERA5_U_Component_Of_Wind_on_{hPa}hpa_levels_daily'
nc_path_V = f'E:\\DATA\\1.Reanalysis\\EAR5\\ERA5-VCW\\{hPa}\\'
nc_file_name_V = f'ERA5_V_Component_Of_Wind_on_{hPa}hpa_levels_daily'

# 生成静态的经纬度网格
lats05 = np.linspace(90 - 0.5 / 2, -90 + 0.5 / 2, 360)
lons05 = np.linspace(-180 + 0.5 / 2, 180 - 0.5 / 2, 720)
lon_grid, lat_grid = np.meshgrid(lons05, lats05)
# 处理经度从-180到180到0-360的转换
lon_grid_360 = np.where(lon_grid < 0, lon_grid + 360, lon_grid)

region_uv_mean = {region: {'U': [], 'V': []} for region in regions.keys()}

for year in range(1980, 2021):
    print(f'Processing year: {year}')
    yearly_uv_data = {region: {'U': [], 'V': []} for region in regions.keys()}
    for mon in range(1, 13):
        year_str = str(year)
        mon_str = str(mon).zfill(2)
        
        # 读取当月的U和V风场数据
        mon_UCW_DS = Dataset(f"{nc_path_U}{nc_file_name_U}{year_str}{mon_str}.nc", 'r')
        mon_VCW_DS = Dataset(f"{nc_path_V}{nc_file_name_V}{year_str}{mon_str}.nc", 'r')
        mon_UCW = mon_UCW_DS.variables['u'][:]
        mon_VCW = mon_VCW_DS.variables['v'][:]

        if 1980 <= year <= 2000:
            lats = mon_UCW_DS.variables['latitude'][:]
            lons = mon_VCW_DS.variables['longitude'][:]
            mon_UCW = np.squeeze(mon_UCW)
            mon_VCW = np.squeeze(mon_VCW)
        elif 2001 <= year <= 2020:
            lats = mon_UCW_DS.variables['lat'][:]
            lons = mon_VCW_DS.variables['lon'][:]

        mon_UCW_DS.close()
        mon_VCW_DS.close()

        # 重采样到0.5°并处理经度范围到-180到180
        for day in range(1, mon_UCW.shape[0] + 1):
            day_UCW_Ori = mon_UCW[day - 1]
            day_VCW_Ori = mon_VCW[day - 1]
            day_UCW05 = resample_model_grid(lats, lons, day_UCW_Ori, 360, 720, np.nan, method='nearest')
            day_VCW05 = resample_model_grid(lats, lons, day_VCW_Ori, 360, 720, np.nan, method='nearest')
            
            # 对于每个区域，检查是否有该天的数据，并提取相应的UV风场数据
            for region, dates in time_dict.items():
                current_date = datetime(year, mon, day).date()
                if current_date in dates:
                    lon_min, lon_max, lat_min, lat_max = regions[region]
                    lon_mask = (lon_grid_360 >= lon_min) & (lon_grid_360 <= lon_max)
                    lat_mask = (lat_grid >= lat_min) & (lat_grid <= lat_max)
                    region_mask = lon_mask & lat_mask
                    
                    region_U = day_UCW05[region_mask]
                    region_V = day_VCW05[region_mask]
                    
                    if len(region_U) > 0 and len(region_V) > 0:
                        yearly_uv_data[region]['U'].append(region_U)
                        yearly_uv_data[region]['V'].append(region_V)

    # 计算每个区域的年平均UV风场
    for region in regions.keys():
        if len(yearly_uv_data[region]['U']) > 0:
            region_uv_mean[region]['U'].append(np.nanmean(np.stack(yearly_uv_data[region]['U']), axis=0))
            region_uv_mean[region]['V'].append(np.nanmean(np.stack(yearly_uv_data[region]['V']), axis=0))
        else:
            region_uv_mean[region]['U'].append(np.nan)
            region_uv_mean[region]['V'].append(np.nan)

joblib.dump(region_uv_mean, 'region_uv_mean.joblib')  
  
#%% 计算风场变化 2000-2020 相比 1980-1999
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
plt.figure(figsize=(12, 6), dpi=500)
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))  # 使用 PlateCarree 投影绘制地图
# ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
ax.set_global()
ax.coastlines(zorder=3, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=0.3, zorder=3)
ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=1)
# ax.set_extent([0, 360, -90, 90], crs=ccrs.PlateCarree())
ax.spines['geo'].set_linewidth(1.2)# 设置地图外框线的宽度

# 添加经纬度数字标签并禁用网格线
gridlines = ax.gridlines(draw_labels=True, linewidth=0, color='dimgray')
gridlines.xlabel_style = {'size': 12}
gridlines.ylabel_style = {'size': 12}

gridlines.bottom_labels = False

# for region in regions:
#     lon_min, lon_max, lat_min, lat_max = regions[region]
#     rect = plt.Rectangle((lon_min, lat_min), lon_max - lon_min, lat_max - lat_min,
#                           linewidth=0.7, edgecolor='black', facecolor='none', transform=ccrs.PlateCarree(), zorder=4)
#     ax.add_patch(rect)
    
# # 创建空的二维全球风场变化数组
# u_diff_global = np.full(lon_grid_360.shape, np.nan)
# v_diff_global = np.full(lon_grid_360.shape, np.nan)

# # 将每个区域的风场变化填回到全球网格
# for region in regions.keys():
#     lon_min, lon_max, lat_min, lat_max = regions[region]
#     if region == 'NI_L':
#         lon_max+=1
#     lon_mask = (lon_grid_360 >= lon_min) & (lon_grid_360 <= lon_max)
#     lat_mask = (lat_grid >= lat_min) & (lat_grid <= lat_max)
#     region_mask = lon_mask & lat_mask
    
#     u_diff_global[region_mask] = region_uv_diff[region]['U']
#     v_diff_global[region_mask] = region_uv_diff[region]['V']
    
#     if region == 'NI_L' or region == 'NI_R':
#         u_diff_global[region_mask] = u_diff_global[region_mask]*(1)
#         v_diff_global[region_mask] = v_diff_global[region_mask]*(1)
        
# # 绘制风场箭头，仅显示在方框内的数据
# quiver_step = 8
# q = None
# # 绘制风场箭头
# mask = ~np.isnan(u_diff_global) & ~np.isnan(v_diff_global)
# lon_masked = np.ma.masked_where(~mask, lon_grid_360)
# lat_masked = np.ma.masked_where(~mask, lat_grid)
# u_masked = np.ma.masked_where(~mask, u_diff_global)
# v_masked = np.ma.masked_where(~mask, v_diff_global)
# q = ax.quiver(lon_masked[5::quiver_step, 5::quiver_step], lat_masked[5::quiver_step, 5::quiver_step],
#               u_masked[5::quiver_step, 5::quiver_step], v_masked[5::quiver_step, 5::quiver_step],
#               scale=100, width=0.0015, alpha=0.6, pivot='middle', color='black', transform=ccrs.PlateCarree(), zorder=3)#'#1f4e78'

# # 添加图例到图的右上角，并加上背景框'#1f4e78'
# rect = plt.Rectangle((0.015, 0.91), 0.1, 0.06, transform=ax.transAxes, facecolor='white', alpha=0.9, zorder=3, edgecolor='black', linewidth=0.7)
# ax.add_patch(rect)
# ax.quiverkey(q, 0.05, 0.94, 2, '2 m/s', labelpos='E', coordinates='axes', color='black', zorder=6)

#绘制海温变化和显著性
AgeDiff = read_arcgis_txt('TPWAgeDiffAve0020-8099_-180-180.txt')['data']
AgeDiffSIG = read_arcgis_txt('TPWAgeDiffSIG0020-8099_-180-180.txt')['data']
Landmask = read_arcgis_txt('countriesmasks.txt')['data']

AgeDiff = np.roll(AgeDiff, shift=AgeDiff.shape[1] // 2, axis=1) # 转换为0-360 再绘制
AgeDiffSIG = np.roll(AgeDiffSIG, shift=AgeDiffSIG.shape[1] // 2, axis=1) # 转换为0-360 再绘制


# 把 AgeDiff 中 -9999 的值 masked 掉，并将 AgeDiffSIG > 1 的位置也在 AgeDiff 中 masked 掉
AgeDiff = np.ma.masked_where(AgeDiff == -9999, AgeDiff)
AgeDiff = np.ma.masked_where(AgeDiff == 0, AgeDiff)
AgeDiff = np.ma.masked_where(AgeDiffSIG > 1, AgeDiff)
AgeDiffland = np.ma.masked_where(Landmask <0, AgeDiff)
AgeDiffsea = np.ma.masked_where(Landmask >= 0, AgeDiff)
AgeDiffland.mask[220:, 560:] = True # 非洲去掉
AgeDiffsea.mask[220:, 560:] = True # 非洲去掉
# AgeDiff = AgeDiff * 0.7 # SST

# 生成静态的经纬度网格
lats05 = np.linspace(90 - 0.5 / 2, -90 + 0.5 / 2, 360)
lons05_0360 = np.linspace(0 + 0.5 / 2, 360 - 0.5 / 2, 720)
lon_grid_new0360, lat_grid = np.meshgrid(lons05_0360, lats05) #上面不能动，只能改这里

# 定义 RGB 颜色，使用 (R, G, B)，范围是 0-255 的整数
# CL,CR = (96, 144, 193), (205, 115, 115) #低饱和
CL,CR = (50, 100, 230), (200, 70, 30) # cool warm

colors = [CL, (211, 211, 211),CR,]
colors = [(r/255, g/255, b/255) for r, g, b in colors]
cland = LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)

colors = [CL, (255, 255, 255),CR,]
colors = [(r/255, g/255, b/255) for r, g, b in colors]
csea = LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)

# 绘制 AgeDiff 的变化 coolwarm
contourland = ax.contourf(lon_grid_new0360, lat_grid, AgeDiffland, levels=np.linspace(-20, 20, 21), cmap=cland, extend='both', zorder=1, transform=ccrs.PlateCarree())
contoursea = ax.contourf(lon_grid_new0360, lat_grid, AgeDiffsea, levels=np.linspace(-20, 20, 21), cmap=csea, extend='both', zorder=1, transform=ccrs.PlateCarree())
cbar = plt.colorbar(contoursea, orientation='horizontal', pad=0.02, aspect=40, shrink=0.78)
cbar.set_label('Difference (mm/d)',fontsize = 12) #°C
cbar.ax.tick_params(labelsize=11.5)
cbar.outline.set_linewidth(1.2)
cbar.ax.tick_params(length=6, width=1.2)
# 绘制显著性区域的填充
# 在 AgeDiffSIG == 1 的位置绘制显著性标记
sig_locs = np.where(AgeDiffSIG == 1)
ax.plot(lon_grid_new0360[sig_locs], lat_grid[sig_locs], 'k.', markersize=0.15, transform=ccrs.PlateCarree(), zorder=2)




# 标题
# plt.title('TC-related Steering Wind and TPW Change (2000-2020 vs 1980-1999)')

# 显示图像
plt.show()


