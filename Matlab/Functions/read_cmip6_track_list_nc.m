function [lat, lon, observation_time, tos, psl_deficit, max_ta_deficit, max_wind850, rv850] = read_cmip6_track_list_nc(nc_file_path)
    % 读取 CMIP6 track_list NetCDF 文件中的所有变量
    %[lat, lon, observation_time] = read_cmip6_track_list_nc('D:\DATA\CMIP6\BCC-CSM2-MR\BCC-CSM2-MR_TC_tracks1980-2015.nc');
    % 读取变量数据
    lat = ncread(nc_file_path, 'lat');
    lon = ncread(nc_file_path, 'lon');
    tos = ncread(nc_file_path, 'tos');
    psl_deficit = ncread(nc_file_path, 'psl_deficit');
    max_ta_deficit = ncread(nc_file_path, 'max_ta_deficit');
    max_wind850 = ncread(nc_file_path, 'max_wind850');
    rv850 = ncread(nc_file_path, 'rv850');
    observation_time = ncread(nc_file_path, 'observation_time');
    
    % 转换 observation_time 为字符数组
    % observation_time = permute(observation_time, [1, 2, 3]);  % 调整维度顺序
    observation_time = char(observation_time);  % 转换为字符数组
end