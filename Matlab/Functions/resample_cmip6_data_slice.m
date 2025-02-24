function data_resampled = resample_cmip6_data_slice(data_slice, rows, cols, target_lat_reso, target_lon_reso, method)
    % 计算原始分辨率
    original_lat_reso = 180 / rows; % 计算纬度分辨率
    original_lon_reso = 360 / cols; % 计算经度分辨率

    % 创建原始的地理参考对象
    R_original = georefcells([-90, 90], [-180, 180], [rows, cols], 'ColumnsStartFrom', 'north', 'RowsStartFrom', 'west');

    % 计算重采样因子
    lat_scale = original_lat_reso / target_lat_reso;
    lon_scale = original_lon_reso / target_lon_reso;

    % 创建目标的地理参考对象
    [data_resampled, ~] = georesize(data_slice, R_original, lat_scale, lon_scale, method);
end