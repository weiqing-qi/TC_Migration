# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 15:54:52 2024

@author: 29585
"""
from Functions import *
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import seaborn as sns
from datetime import datetime
import joblib

# 计算1980-2020年每个区域逐年的年平均UV风场
hPa = 500
nc_path_U = f'E:\\DATA\\1.Reanalysis\\EAR5\\ERA5-UCW\\{hPa}\\'
nc_file_name_U = f'ERA5_U_Component_Of_Wind_on_{hPa}hpa_levels_daily'
nc_path_V = f'E:\\DATA\\1.Reanalysis\\EAR5\\ERA5-VCW\\{hPa}\\'
nc_file_name_V = f'ERA5_V_Component_Of_Wind_on_{hPa}hpa_levels_daily'

yearly_uv_mean = {'U': [], 'V': []}

for year in range(1980, 2021):
    print(f'Processing year: {year}')
    monthly_uv_data = {'U': [], 'V': []}
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
        
        # 把月数据在第一个维度上平均得到月平均
        mon_UCW_mean = np.nanmean(mon_UCW, axis=0)
        mon_VCW_mean = np.nanmean(mon_VCW, axis=0)
        
        # 利用下面的语句重采样
        mon_UCW05_mean = resample_model_grid(lats, lons, mon_UCW_mean, 360, 720, np.nan, method='nearest')
        mon_VCW05_mean = resample_model_grid(lats, lons, mon_VCW_mean, 360, 720, np.nan, method='nearest')
        
        # 把月数据保存进monthly_uv_data
        monthly_uv_data['U'].append(mon_UCW05_mean)
        monthly_uv_data['V'].append(mon_VCW05_mean)
    
    # 计算每年的平均风场，把年平均风场存进yearly_uv_mean
    yearly_uv_mean['U'].append(np.nanmean(monthly_uv_data['U'], axis=0))
    yearly_uv_mean['V'].append(np.nanmean(monthly_uv_data['V'], axis=0))

# 使用 joblib 保存yearly_uv_mean 计算结果到文件
joblib.dump(yearly_uv_mean, 'yearly_global_uv_mean.joblib')

#%% 海温数据
from Functions import *
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import seaborn as sns
from datetime import datetime
import joblib

# 计算1980-2020年每个区域逐年的年平均UV风场
nc_path = 'E:\\DATA\\1.Reanalysis\\EAR5\\ERA5-Sea_surface_temperature\\'
nc_file_name = 'ERA5_Sea_Surface_Temperature_on_single_levels_daily'

yearly_sst_mean = []

for year in range(1980, 2021):
    print(f'Processing year: {year}')
    monthly_sst_data = []
    for mon in range(1, 13):
        year_str = str(year)
        mon_str = str(mon).zfill(2)
        
        # 读取当月的SST数据
        sst_DS = Dataset(f"{nc_path}{nc_file_name}{year_str}{mon_str}.nc", 'r')

        if 1950 <= year <= 1965:
            sst_data = sst_DS.variables['sst'][:]
            lats = sst_DS.variables['lat'][:]
            lons = sst_DS.variables['lon'][:]
            sst_data = np.squeeze(sst_data)
        elif 1966 <= year <= 2020:
            sst_data = sst_DS.variables['tos'][:]
            lats = sst_DS.variables['lat'][:]
            lons = sst_DS.variables['lon'][:]

        sst_DS.close()
        
        # 把SST数据从开尔文转换为摄氏度
        sst_data_celsius = sst_data - 273.15
        
        # 把月数据在第一个维度上平均得到月平均
        mon_sst_mean = np.nanmean(sst_data_celsius, axis=0)
        
        # 利用下面的语句重采样
        mon_sst05_mean = resample_model_grid(lats, lons, mon_sst_mean, 360, 720, np.nan, method='nearest')
        
        # 把月数据保存进monthly_sst_data
        monthly_sst_data.append(mon_sst05_mean)
    
    # 计算每年的平均海表温度，把年平均存进yearly_sst_mean
    yearly_sst_mean.append(np.nanmean(monthly_sst_data, axis=0))

# 使用 joblib 保存yearly_sst_mean 计算结果到文件
joblib.dump(yearly_sst_mean, 'yearly_global_sst_mean-180-180.joblib')
        
#%% TCW数据
from Functions import *
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import seaborn as sns
from datetime import datetime
import joblib

# 计算1980-2020年每个区域逐年的年平均UV风场
nc_path = 'E:\\DATA\\1.Reanalysis\\EAR5\\ERA5-Total-Column-Water\\'
nc_file_name = 'ERA5_Total_Column_Water_on_single_levels_daily'

yearly_tcw_mean = []

for year in range(1980, 2021):
    print(f'Processing year: {year}')
    monthly_tcw_data = []
    for mon in range(1, 13):
        year_str = str(year)
        mon_str = str(mon).zfill(2)
        
        # 读取当月的tcw数据
        tcw_DS = Dataset(f"{nc_path}{nc_file_name}{year_str}{mon_str}.nc", 'r')
        tcw_data = tcw_DS.variables['tcw'][:]
        lats = tcw_DS.variables['lat'][:]
        lons = tcw_DS.variables['lon'][:]
        tcw_data = np.squeeze(tcw_data)
        tcw_DS.close()
        
        # 把月数据在第一个维度上平均得到月平均
        mon_tcw_mean = np.nanmean(tcw_data, axis=0)
        
        # 利用下面的语句重采样
        mon_tcw05_mean = resample_model_grid(lats, lons, mon_tcw_mean, 360, 720, np.nan, method='nearest')
        
        # 把月数据保存进monthly_tcw_data
        monthly_tcw_data.append(mon_tcw05_mean)
    
    # 计算每年的平均海表温度，把年平均存进yearly_tcw_mean
    yearly_tcw_mean.append(np.nanmean(monthly_tcw_data, axis=0))

# 使用 joblib 保存yearly_tcw_mean 计算结果到文件
joblib.dump(yearly_tcw_mean, 'yearly_global_tcw_mean-180-180.joblib')

#%% 画图全球多年平均
from Functions import *
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import seaborn as sns
from datetime import datetime
import joblib

yearly_global_sst_mean = joblib.load('yearly_global_sst_mean-180-180.joblib')
yearly_global_uv_mean = joblib.load('yearly_global_uv_mean-180-180.joblib')

# 生成静态的经纬度网格
lats05 = np.linspace(90 - 0.5 / 2, -90 + 0.5 / 2, 360)
lons05 = np.linspace(0 + 0.5 / 2, 360 - 0.5 / 2, 720)
lon_grid_360, lat_grid = np.meshgrid(lons05, lats05)

# 计算1980-1999年和2000-2020年的SST平均和500hPa UV风场平均
yearly_sst_mean_1980_1999 = np.nanmean(yearly_global_sst_mean[:20], axis=0)
yearly_sst_mean_2000_2020 = np.nanmean(yearly_global_sst_mean[20:], axis=0)

yearly_u_mean_1980_1999 = np.nanmean(yearly_global_uv_mean['U'][:20], axis=0)
yearly_u_mean_2000_2020 = np.nanmean(yearly_global_uv_mean['U'][20:], axis=0)

yearly_v_mean_1980_1999 = np.nanmean(yearly_global_uv_mean['V'][:20], axis=0)
yearly_v_mean_2000_2020 = np.nanmean(yearly_global_uv_mean['V'][20:], axis=0)

# 将SST和风场数据转换为0-360度经度范围
sst_diff = yearly_sst_mean_2000_2020 - yearly_sst_mean_1980_1999
sst_diff = np.roll(sst_diff, shift=sst_diff.shape[1] // 2, axis=1) # 转换为0-360 再绘制

du = yearly_u_mean_2000_2020 - yearly_u_mean_1980_1999
dv = yearly_v_mean_2000_2020 - yearly_v_mean_1980_1999
du = np.roll(du, shift=du.shape[1] // 2, axis=1) # 转换为0-360 再绘制
dv = np.roll(dv, shift=dv.shape[1] // 2, axis=1) # 转换为0-360 再绘制

# 绘制1980-1999年和2000-2020年的SST和风场变化
plt.figure(figsize=(12, 6), dpi=600)
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))  # 使用 PlateCarree 投影绘制地图
# ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))  # 使用 PlateCarree 投影绘制地图
ax.coastlines(zorder=1)  # 添加海岸线
ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=0.3, zorder=2)  # 添加国家边界线
ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=1)
# ax.set_extent([0, 360, -90, 90], crs=ccrs.PlateCarree())  # 设置地图的显示范围为0-360 Robinson可以关闭
# 添加经纬度数字标签并禁用网格线
gridlines = ax.gridlines(draw_labels=True, linewidth=0, color='dimgray')
gridlines.bottom_labels = False

# 定义 RGB 颜色，使用 (R, G, B)，范围是 0-255 的整数
# CL,CR = (96, 144, 193), (205, 115, 115) #低饱和
CL,CR = (50, 100, 230), (200, 70, 30) # cool warm

colors = [CL, (211, 211, 211),CR,]
colors = [(r/255, g/255, b/255) for r, g, b in colors]
coolwarm = LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)

# 绘制SST变化
contour = ax.contourf(lon_grid_360, lat_grid, sst_diff, levels=np.linspace(-2, 2, 21), cmap=coolwarm, extend='both', zorder=0, transform=ccrs.PlateCarree())
cbar = plt.colorbar(contour, orientation='horizontal', pad=0.02, aspect=40, shrink=0.78)
cbar.set_label('SST Difference (°C)')

# 绘制风场箭头变化
quiver_step = 10  # 控制箭头的间隔
q = plt.quiver(lon_grid_360[5::quiver_step, 5::quiver_step], lat_grid[5::quiver_step, 5::quiver_step],
               du[5::quiver_step, 5::quiver_step], dv[5::quiver_step, 5::quiver_step],
               scale=40, width=0.0015, alpha=0.6, pivot='middle', color='black', transform=ccrs.PlateCarree())

# 图例
rect = plt.Rectangle((0.015, 0.91), 0.1, 0.06, transform=plt.gca().transAxes, facecolor='white', alpha=0.90, zorder=1, edgecolor='black', linewidth=0.7)
plt.gca().add_patch(rect)
plt.quiverkey(q, 0.045, 0.94, 1, '1 m/s', labelpos='E', coordinates='axes', color='black', zorder=1)

# 标题
plt.title('Steering Wind and SST Changes (2000-2020 vs 1980-1999)')
plt.show()

