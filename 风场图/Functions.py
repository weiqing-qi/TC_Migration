# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:47:15 2024

@author: 29585
"""
import os
import pandas as pd
import netCDF4 as nc
import numpy as np
import rasterio
from scipy.interpolate import griddata
from datetime import datetime


def col_letter_to_index(letter):
    """将Excel列字母转换为数字索引（0基）"""
    num = 0
    for c in letter:
        if 'A' <= c <= 'Z':
            num = num * 26 + (ord(c) - ord('A')) + 1
    return num - 1


def num2abc2(num):
    """
    将10进制列数转换为26进制字母excel表列名
    num: 1-based column number
    """
    string = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    
    strs = ""
    while num > 0:
        m = num % 26
        if m == 0:
            m = 26
        # Python中的字符串是不可变类型，因此每次实际上是创建了一个新的字符串
        strs = string[m - 1] + strs
        num = int((num - m) / 26)
        
    return strs



def get_Excel_DPRE(Excel_dir_Pre, oiorolos, skiprows = 1, nrows = 2273):
    # EXCEL之中对应的存贮列名
    Exc_oiorolos_ColName = [":".join([num2abc2(5  ),num2abc2(54 )]),
                            ":".join([num2abc2(55 ),num2abc2(104)]),
                            ":".join([num2abc2(105),num2abc2(154)]),
                            ":".join([num2abc2(155),num2abc2(204)]),
                            ":".join([num2abc2(205),num2abc2(254)])]
    
    Exc_oiorolos_ColName4Locat = [":".join([num2abc2(5  ),num2abc2(104 )]),
                            ":".join([num2abc2(105 ),num2abc2(204)]),
                            ":".join([num2abc2(205),num2abc2(304)]),
                            ":".join([num2abc2(305),num2abc2(404)]),
                            ":".join([num2abc2(405),num2abc2(504)])]
    # 读取
    PRate_df = pd.read_excel(Excel_dir_Pre, 
                          sheet_name = 0,                                 # 或者sheet_name='Sheet2'
                          usecols = Exc_oiorolos_ColName[oiorolos],                                      # 指定读取列
                          skiprows = skiprows,                                     # 跳过前两行
                          nrows = nrows)                                        # 读取行数 第一行作为key 1-based
    
    PAmount_df = pd.read_excel(Excel_dir_Pre, 
                          sheet_name=1,
                          usecols = Exc_oiorolos_ColName[oiorolos],
                          skiprows = skiprows,
                          nrows = nrows)
    
    PMax_df = pd.read_excel(Excel_dir_Pre, 
                          sheet_name=2,
                          usecols = Exc_oiorolos_ColName[oiorolos],
                          skiprows = skiprows,
                          nrows = nrows)
    
    PArea_df = pd.read_excel(Excel_dir_Pre, 
                          sheet_name=3,
                          usecols = Exc_oiorolos_ColName[oiorolos],
                          skiprows = skiprows,
                          nrows = nrows)
    
    MPCo_df = pd.read_excel(Excel_dir_Pre, 
                          sheet_name=4,
                          usecols = Exc_oiorolos_ColName4Locat[oiorolos],
                          skiprows = skiprows,
                          nrows = nrows)
    MPLon_df = MPCo_df.iloc[:,0:50]
    MPLat_df = MPCo_df.iloc[:,50:100]
    
    
    
    
    return {
        "PRate"   : PRate_df,
        "PAmount" : PAmount_df,
        "PMax"    : PMax_df,
        "PArea"   : PArea_df,
        "MPLon"   : MPLon_df,
        "MPLat"   : MPLat_df,
        }




def decode_element(element): 
    """
    定义一个函数，用于解码单个 bytes8 元素
    """
    return element.decode('utf-8')
vectorized_decode = np.vectorize(decode_element) # 使用 Numpy 的 vectorize 函数创建一个向量化的解码函数




def get_IBTrACK(IBTrACS_Path, SI, EI):
    """
    ---------------------------------------------------------------------------
    SI, EI  = start and end index of the TC in the IBTrACS IBTrACS_Path
    IBTrACS_Path = IBTrACS file
    
    USA_Agencies is using 1-min m.s.wind, TRANS into 10-min m.s.wind below
    One knot is 0.514 m s-1

    Returns a dictionary.
    ---------------------------------------------------------------------------
    """
    ibtracs = nc.Dataset(IBTrACS_Path,"r") 

    # Get time stamps, has the structure of YYYY-MM-DD HH:MM:SS,
    # access HH:MM with <dict>['times'][<index>][11:16]
    time_ori= ibtracs['iso_time'][SI:EI].data
    times = vectorized_decode(time_ori) # 应用向量化的解码函数到整个数组

    # Get other values
    dist2land = ibtracs['dist2land'][SI:EI].data
    usa_pres = ibtracs['usa_pres'][SI:EI].data
    landfall = ibtracs['landfall'][SI:EI].data
    storm_dir = ibtracs['storm_dir'][SI:EI].data
    storm_speed = ibtracs['storm_speed'][SI:EI].data
    
    lats = ibtracs['lat'][SI:EI].data
    lons = ibtracs['lon'][SI:EI].data
    # basin_raw = ibtracs['basin'][SI:EI].data
    # basin = vectorized_decode(basin_raw)    
    # sshs = ibtracs['usa_sshs'][SI:EI].data
    # usa_wind = ibtracs['usa_wind'][SI:EI].data / 1.12
    usa_wind = ibtracs['usa_wind'][SI:EI].data
    wmo_wind = ibtracs['wmo_wind'][SI:EI].data
    bom_wind = ibtracs['bom_wind'][SI:EI].data
    cma_wind = ibtracs['cma_wind'][SI:EI].data
    ds824_wind = ibtracs['ds824_wind'][SI:EI].data
    hko_wind = ibtracs['hko_wind'][SI:EI].data
    mlc_wind = ibtracs['mlc_wind'][SI:EI].data
    nadi_wind = ibtracs['nadi_wind'][SI:EI].data
    neumann_wind = ibtracs['neumann_wind'][SI:EI].data
    newdelhi_wind = ibtracs['newdelhi_wind'][SI:EI].data
    reunion_wind = ibtracs['reunion_wind'][SI:EI].data
    td9636_wind = ibtracs['td9636_wind'][SI:EI].data
    tokyo_wind = ibtracs['tokyo_wind'][SI:EI].data
    wellington_wind = ibtracs['wellington_wind'][SI:EI].data
   
    return {        
        "usa_wind": usa_wind,               # NOAA and JTWC 1 min -> 10 min？
        # "wmo_wind": wmo_wind,
        #"bom_wind": bom_wind,                    # 10 min
        # "cma_wind": cma_wind,             # 2 min
        # "ds824_wind": ds824_wind,         # 1877-1980 historical
        # "hko_wind": hko_wind,
        # "mlc_wind": mlc_wind,             # 1851-1898 historical
        #"nadi_wind": nadi_wind,                  # 10 min
        # "neumann_wind": neumann_wind,
        # "newdelhi_wind": newdelhi_wind,   # India IMD 3 min
        #"reunion_wind": reunion_wind,            # 10 min
        # "td9636_wind": td9636_wind,       # before 1980 historical
        #"tokyo_wind": tokyo_wind,                # 10 min
        #"wellington_wind": wellington_wind,      # 10 min
           }, {"lats": lats,"lons": lons, }, times,  {"dist2land": dist2land, "usa_wind": usa_wind,
                                                      "usa_pres": usa_pres, 'landfall':landfall,
                                                      "storm_dir":storm_dir, "storm_speed":storm_speed}




def get_NC_EC_data(NC_dir_EC, oiorolos = 0):
    """
    ---------------------------------------------------------------------------
    oiorolos  = 需要的气旋范围，有'O','I', 'OR', 'OL', 'OS'五个选项，默认为0

    Returns a dictionary containing the all EC on all 18 pressure levels, from
    950: 50: 100 hPa. Shape[18, 5, 50, 2272] 18 LEVELS, 5 SCALE, 50 DAYS, 2272
    RECORDS.
    ---------------------------------------------------------------------------
    """
    ALL_EC = nc.Dataset(NC_dir_EC,"r")
    
    HD_abs_mean = ALL_EC['HD_abs_mean'][:, oiorolos, :, :].data
    HD_mean = ALL_EC['HD_mean'][:, oiorolos, :, :].data
    RH_mean = ALL_EC['RH_mean'][:, oiorolos, :, :].data
    RV_abs_mean = ALL_EC['RV_abs_mean'][:, oiorolos, :, :].data
    RV_mean = ALL_EC['RV_mean'][:, oiorolos, :, :].data
    T_mean = ALL_EC['T_mean'][:, oiorolos, :, :].data
    UVW_abs_mean = ALL_EC['UVW_abs_mean'][:, oiorolos, :, :].data
    UVW_mean = ALL_EC['UVW_mean'][:, oiorolos, :, :].data
    UVW_mean_Azimuth = ALL_EC['UVW_mean_Azimuth'][:, oiorolos, :, :].data
    VV_abs_mean = ALL_EC['VV_abs_mean'][:, oiorolos, :, :].data
    VV_mean = ALL_EC['VV_mean'][:, oiorolos, :, :].data
    
    return{
        "HD_abs_mean": HD_abs_mean,
        "HD_mean": HD_mean,
        "RH_mean": RH_mean,
        "RV_abs_mean": RV_abs_mean,
        "RV_mean": RV_mean,
        "T_mean": T_mean,
        "UVW_abs_mean": UVW_abs_mean,
        "UVW_mean": UVW_mean,
        "UVW_mean_Azimuth": UVW_mean_Azimuth,
        "VV_abs_mean": VV_abs_mean,
        "VV_mean": VV_mean,
        }




def index_trans_3h_to_day(selected_winds, times, record_len = 2272, time_dim_3h = 360):
    """
    ---------------------------------------------------------------------------
    selected_winds = (:, 360) bool array; from get_IBTrACK.wind_dick; 
                     filtered by lats scale
    times = (:, 360, 19) str32 array; from get_IBTrACK.times
    
    Returns a pair array, which contain coordinates for (2272, 50) (records, day)
    daily TC PRE and EC data
    ---------------------------------------------------------------------------
    """
    day_coor = []  # 用于存储天数 0-based
    record_coor = [] # 用于存储记录编号
    unique_coords = set() # 初始化用于去重的集合
    
    for i in range(record_len):
        for j in range(time_dim_3h):
            
            if j == 0:  # 第一个时间点，初始化基准日期
                base_date = datetime.strptime(''.join(times[i, j, :]), '%Y-%m-%d %H:%M:%S').date()   
                
            if selected_winds[i, j]:
                date_str = ''.join(times[i, j, :])  # 将字符数组合并成字符串
                current_date = datetime.strptime(date_str, '%Y-%m-%d %H:%M:%S').date() #             
                
                day_diff = (current_date - base_date).days
                
                # 使用元组(i, day_diff)作为坐标对进行去重检查
                if (i, day_diff) not in unique_coords:
                    unique_coords.add((i, day_diff))  # 添加新坐标对到集合中进行去重
                    record_coor.append(i)
                    day_coor.append(day_diff)
    
    return day_coor, record_coor




def identify_basin(longitude, latitude):
    """
    Identify the basin for the given longitude and latitude.
    
    :param longitude: float, longitude of the location (-180 to 180 or 180 to 360)
    :param latitude: float, latitude of the location (-90 to 90)
    :return: list, [basin_number, basin_name]
    """
    if (longitude > 360) or (longitude < -180) or (latitude > 90) or (latitude < -90):
        raise ValueError("Check identify_basin's inputs")
        
    # Normalize longitude to -180 to 180 if necessary
    if longitude > 180:
        longitude -= 360
        
    basins = [
        (1, 'SI', 10, 135, 'Southern'),   # South Indian
        (2, 'SP', 135, -70, 'Southern'),  # South Pacific
        (3, 'SA', -70, 10, 'Southern'),   # South Atlantic
        (4, 'NI', 30, 100, 'Northern'),   # North Indian
        (5, 'WP', 100, 180, 'Northern'),  # Western Pacific
        (6, 'EP', -180, -105, 'Northern'),# Eastern Pacific
        (7, 'NA', -105, 30, 'Northern')   # North Atlantic
    ]

    for basin in basins:
        number, name, lon_min, lon_max, hemisphere = basin
        
        # Consider latitude 0 as Northern Hemisphere
        if (hemisphere == 'Southern') and (latitude >= 0):
            continue
        if (hemisphere == 'Northern') and (latitude < 0):
            continue
        
        if lon_min < lon_max:
            if lon_min < longitude <= lon_max:
                return int(number)
        else:  # Handle wrap-around for South Pacific (SP) and South Atlantic (SA)
            if (longitude > lon_min) or (longitude <= lon_max):
                return int(number)

    return None

# loop_results = []
# for lon, lat in zip(result_lons, result_lats):
#     if lon == -9999 or lat == -9999:
#         loop_results.append(np.nan)
#     else:
#         loop_results.append(identify_basin(lon, lat))
# 
# loop_results = np.array(loop_results)


def identify_basin_array(longitudes, latitudes):
    """
    Identify the basins for the given arrays of longitudes and latitudes.
    
    :param longitudes: np.ndarray, array of longitudes of the locations (-180 to 180 or 180 to 360)
    :param latitudes: np.ndarray, array of latitudes of the locations (-90 to 90)
    :return: np.ndarray, array of basin numbers, or np.nan if invalid coordinates (-9999)
    """
    # Ensure inputs are NumPy arrays
    longitudes = np.asarray(longitudes)
    latitudes = np.asarray(latitudes)

    # Initialize output array with np.nan
    basin_numbers = np.full(longitudes.shape, np.nan)
    
    # Check for invalid coordinates (-9999) and skip processing for those
    invalid_mask = (longitudes == -9999) | (latitudes == -9999)
    
    # Normalize longitudes to -180 to 180 if necessary
    longitudes = np.where(longitudes > 180, longitudes - 360, longitudes)
    
    # Check bounds and handle out of range values by retaining np.nan
    out_of_bounds_mask = (longitudes > 360) | (longitudes < -180) | (latitudes > 90) | (latitudes < -90)
    
    # Combined mask to identify positions to process
    valid_mask = ~(invalid_mask | out_of_bounds_mask)

    # Define basin boundaries
    basins = [
        (1, 'SI', 10, 135, 'Southern'),   # South Indian
        (2, 'SP', 135, -70, 'Southern'),  # South Pacific
        (3, 'SA', -70, 10, 'Southern'),   # South Atlantic
        (4, 'NI', 30, 100, 'Northern'),   # North Indian
        (5, 'WP', 100, 180, 'Northern'),  # Western Pacific
        (6, 'EP', -180, -105, 'Northern'),# Eastern Pacific
        (7, 'NA', -105, 30, 'Northern')   # North Atlantic
    ]
    
    # Process valid coordinates
    for basin in basins:
        number, name, lon_min, lon_max, hemisphere = basin
        
        # Create masks for latitude conditions
        if hemisphere == 'Southern':
            lat_mask = latitudes < 0
        else:
            lat_mask = latitudes >= 0
        
        # Handle wrap-around cases
        if lon_min < lon_max:
            lon_mask = (lon_min < longitudes) & (longitudes <= lon_max)
        else:
            lon_mask = (longitudes > lon_min) | (longitudes <= lon_max)
        
        # Combine valid, latitude, and longitude conditions
        mask = valid_mask & lat_mask & lon_mask
        
        # Assign basin number where mask is True
        basin_numbers[mask] = int(number)
    
    return basin_numbers

# identify_basin(111,50)
# identify_basin(-99,50)
# identify_basin(-80,-3)
# identify_basin_array(np.array([111, -99, -80]),np.array([50,50,-3]))

#%% 读取DEM    
def read_geotiff(file_path):
    """
    读取GeoTIFF文件并返回包含元数据和波段数据的字典。
    """
    with rasterio.open(file_path) as geotiff:
        metadata = {
            'width': geotiff.width,                   # 获取图像的宽度(int)
            'height': geotiff.height,                 # 获取图像的高度(int)
            'bands_num': geotiff.count,               # 获取波段数(int)
            'crs': geotiff.crs,                       # 获取投影信息(crs)
            'transform': geotiff.transform,           # 获取地理变换信息 (Affine)
            'bands': [geotiff.read(i) for i in range(1, geotiff.count + 1)]  # 读取所有波段数据(numpy.ndarray) 1-based
        }
    
    return metadata    
  
  
def read_arcgis_txt(file_path):
    """
    读取ArcGIS ASCII Grid格式的TXT文件并返回包含头文件字典和数据
    """
    header = {} 
    # 读取头文件信息
    with open(file_path, 'r') as file:
        for _ in range(6):
            line = file.readline().strip()
            key, value = line.split()
            header[key] = float(value)
        
    # 读取数据部分
    data = np.genfromtxt(file_path, skip_header=6)
    
    return {'header': header, 'data': data}    


def resample_model_grid(lats, lons, data, num_rows, num_cols, fill_value, method = 'nearest', ):
    """
    将非规则间隔的地理网格重采样到规则全球网格(超出部分为nan值)
    
    param:
        lats (np.array): 纬度数组 (n, 1)。
        lons (np.array): 经度数组 (n, 1)。
        data (np.array): (行,列)对应于 (lats, lons) 的数据网格。
        num_rows (int): 全球网格的行数。
        num_cols (int): 全球网格的列数。    
        method (str): nearest('默认，保正mask值不变')
    return: extend LL should be -90 -180 Masked array

    """
    # 检查输入形状
    if (lats.shape[0] != data.shape[0]) or (lons.shape[0] != data.shape[1]):
        raise ValueError("输入的 lats 和 lons 的长度必须分别对应 data 的行和列")
    #应对lons的不同范围0-360
    lons = np.where(lons > 180, lons - 360, lons)
    # 计算新网格的经纬度范围和分辨率
    lat_min, lat_max = -90, 90
    lon_min, lon_max = -180, 180
    lat_res = (lat_max - lat_min) / num_rows
    lon_res = (lon_max - lon_min) / num_cols

    # 创建规则的全球网格
    new_lats = np.linspace(lat_max - lat_res/2, lat_min + lat_res/2, num_rows)
    new_lons = np.linspace(lon_min + lon_res/2, lon_max - lon_res/2, num_cols)
    new_lats_grid, new_lons_grid = np.meshgrid(new_lats, new_lons, indexing='ij')

    # 将原始数据的网格化经纬度数据转换为一维数组
    lats_grid, lons_grid = np.meshgrid(lats, lons, indexing='ij')
    points = np.column_stack((lats_grid.ravel(), lons_grid.ravel()))
    values = data.ravel()

    # 使用最近邻插值填充新的数据网格
    resampled_grid = griddata(points, values, (new_lats_grid, new_lons_grid), method, fill_value)

    resampled_grid = np.ma.masked_where(resampled_grid == fill_value, resampled_grid) # 注意与np.ma.masked_values之间有不同，这个更合适
    
    return resampled_grid


def extent_handler(data_slice, first_lat, first_lon, ):
    """
    输入一张2D图像自动匹配到-180 180 90 -90
    
    :param: data_slice: 全球二维格网，要求经度比维度维度更长
    :param: first_lat: NC 中经常有的一列lat中的第一个值
    :param: first_lon: NC 中经常有的一列lon中的第一个值
    :return: -180 180 90 -90 的2d格网
    :example: extent_handler(data_slice, -89.14, 0)
    """
    num_rows, num_cols = data_slice.shape
    #是否需要转置，注意经度的维度要高于纬度纬度长度
    if (num_rows > num_cols) and (num_rows > 0) and (num_cols > 0):
        data_slice = data_slice.T
    elif (num_rows < num_cols) and (num_rows > 0) and (num_cols > 0):
        pass
    else:
        raise ValueError("please check the shape of 2d-grid")
        
    #是否需要上下翻转    
    if -90.5 < first_lat < -88:
        data_slice = data_slice[::-1, :]
    elif 90.5 > first_lat > 88:
        pass
    else: 
        raise ValueError("please check the lat of 2d-grid")
        
    #是否需要经度转换
    if -0.5 < first_lon < 2:
        num_rows, num_cols = data_slice.shape
        if num_cols % 2 != 0:
            raise ValueError("The number of longitude columns cannot be evenly divided by 2")
        mid_cols = num_cols // 2
        data_slice = np.hstack((data_slice[:, mid_cols:], data_slice[:, :mid_cols]))
    elif -178 > first_lon > -180.5:
        pass
    else: 
        raise ValueError("please check the lon of 2d-grid")
        
    return data_slice




#%% 经度平均处理方法
def mean_lon(longitudes):
    '''
    该脚本提供了一个函数，用于计算经度的平均值，考虑到经度跨越 180° 边界的问题。
    方法：该函数将经度值转换为单位圆上的复数表示，计算这些点的平均值，然后将平均值转换回经度。这种方法可以有效处理经度跨越 -180° 和 180° 的情况。
    输入：- 一个经度值的列表或 NumPy 数组，每个值在 -180° 到 180° 范围内。
    输出：- 一个表示平均经度的单一值，范围在 -180° 到 180° 之间。
    '''
    # 将经度转换为弧度
    radians = np.radians(longitudes)
    
    # 将经度转换为单位圆上的复数表示
    x = np.cos(radians)
    y = np.sin(radians)
    
    # 计算 x 和 y 的平均值
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    
    # 计算平均角度（弧度）
    mean_angle = np.arctan2(mean_y, mean_x)
    
    # 将平均角度转换回度数，并确保在 -180° 到 180° 范围内
    mean_long = np.degrees(mean_angle)
    mean_long = normalize_lon(mean_long)
    
    return mean_long

# 确保经度在 -180° 到 180° 范围内
def normalize_lon(longitude):

    longitude = (longitude + 180) % 360 - 180
    return longitude

# 将纬度和经度转换为格网的行列索引
def latlon_to_grid_indices(latitudes, longitudes, resolution=0.1): 
    row_indices = (90 - latitudes) / resolution # 行索引：纬度从 90 到 -90，每个格网间隔为 resolution 度
    col_indices = (longitudes + 180) / resolution # 列索引：经度从 -180 到 180，每个格网间隔为 resolution 度
    # 将 NaN 值替换为指定的 fill_value
    row_indices = np.where(np.isnan(row_indices), -9999, row_indices).astype(np.int32)
    col_indices = np.where(np.isnan(col_indices), -9999, col_indices).astype(np.int32)
    return row_indices, col_indices


def window_calculater(data, row_indices, col_indices, size=3, method='max', center=True):
    """
    # param: data 2d-masked.array 全球尺寸一帧格网
    # param: row_indices, col_indices 两列索引值
    # param: size window size
    # param: method 要计算的值
    # param: center 是否考虑中心点
    # return: 一维array
    """
    nrow, ncol = data.shape

    if size % 2 == 0 or size < 3:
        raise ValueError("Size must be an odd number and at least 3.")

    # 计算中间索引，生成偏移数组    
    mid_index = size // 2
    offsets = np.arange(-mid_index, mid_index + 1)
    row_offsets, col_offsets = np.meshgrid(offsets, offsets, indexing='ij')

    # 将2D数组转换为1D数组，且排除中心点
    row_offsets = row_offsets.flatten()
    col_offsets = col_offsets.flatten()

    if center == False:
        center_point_index = len(row_offsets) // 2
        row_offsets = np.delete(row_offsets, center_point_index)
        col_offsets = np.delete(col_offsets, center_point_index)
    else:
        pass

    # 计算邻域索引
    neighbor_rows = (row_indices.reshape(-1, 1) + row_offsets) % nrow
    neighbor_cols = (col_indices.reshape(-1, 1) + col_offsets) % ncol

    # Ensure the indices are integers
    neighbor_rows = neighbor_rows.astype(int)
    neighbor_cols = neighbor_cols.astype(int)

    neighbors = data[neighbor_rows, neighbor_cols]

    if method == 'max':
        window_result = np.nanmax(neighbors, axis=1)
    elif method == 'min':
        window_result = np.nanmin(neighbors, axis=1)        
    elif method == 'mean':
        window_result = np.nanmean(neighbors, axis=1)
    else:
        raise ValueError("Unsupported method. Choose from 'max', 'min', 'mean'.")

    return window_result




