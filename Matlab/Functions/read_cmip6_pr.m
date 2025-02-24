function alldata = read_cmip6_pr(model_name, allday_S, search_path, target_reso, Rsmp_method, leap)
    % 搜索路径中的所有文件 注意还要处理单位
    files = dir(fullfile(search_path, '**', ['*pr*' model_name(end-6:end) '*.nc'])); %有的文件名不包含全名
    
    % 初始化变量
    alldata_initialized = false;
    days = size(allday_S,1);

    % 提取文件名和起止日期
    for f = 1:length(files)
        file = files(f);
        filename = file.name;
        filefolder =file.folder;
        file_path = fullfile(filefolder, filename);
        
        % 获取文件的起止日期
        start_date_str = filename(end-19:end-12);
        end_date_str = filename(end-10:end-3);
        start_date = datetime(start_date_str, 'InputFormat', 'yyyyMMdd');
        end_date = datetime(end_date_str, 'InputFormat', 'yyyyMMdd');

        if ~alldata_initialized %默认所有的pr文件这些信息是一样的
            % 读取纬度和经度数据
            lat = ncread(file_path, 'lat');
            lon = ncread(file_path, 'lon');

            % 获取 pr 的维度信息
            pr_info = ncinfo(file_path, 'pr');
            pr_dims = {pr_info.Dimensions.Name};
            pr_sizes = [pr_info.Dimensions.Length];
        
            % 找到时间、经度和纬度维度
            lon_dim_i = find(strcmp(pr_dims,'lon'));
            lat_dim_i = find(strcmp(pr_dims,'lat'));
            time_dim_i = find(strcmp(pr_dims,'time'));  

            % 获取第一个文件的 missing_value 和 _FillValue
            try
                missing_value = ncreadatt(file_path, 'pr', 'missing_value');
            catch
                missing_value = NaN; disp('不存在missing_value属性,已经设置为NaN')% 设置默认值
            end
            
            try
                fill_value = ncreadatt(file_path, 'pr', '_FillValue');
            catch
                error('The _FillValue attribute is missing in the NetCDF file');
            end

            % 初始化 alldata 变量
            rows = pr_sizes(lat_dim_i);
            cols = pr_sizes(lon_dim_i);
            alldata = NaN(round(180/target_reso), round(360/target_reso), days);
            alldata_initialized = true;
        end

        % 查找 allday_S 中的日期，并提取对应的 pr 数据
        for i = 1:days
            date_str = allday_S(i, :);
            date = datetime(date_str, 'InputFormat', 'yyyyMMdd');
            
            % 检查日期是否在文件的起止日期范围内
            if date >= start_date && date <= end_date

                % 计算该日期在文件中的索引
                day_index = daysact(start_date, date) + 1;

                % 如果不考虑闰年，调整 day_index
                if ~leap
                    % 计算 start_date 和 date 之间的闰日数（2月29日）
                    leap_days = 0;
                    for y = year(start_date):year(date)
                        if leapyear(y)
                            feb29 = datetime(y, 2, 29);
                            if feb29 >= start_date && feb29 <= date
                                leap_days = leap_days + 1;
                            end
                        end
                    end
                    day_index = day_index - leap_days;
                end

                % 读取特定日期的 pr 数据                
                if time_dim_i == 1
                    data_slice = ncread(file_path, 'pr', [day_index, 1, 1], [1, Inf, Inf]);
                    data_slice = permute(data_slice, [lat_dim_i lon_dim_i time_dim_i]);
                elseif time_dim_i == 2
                    data_slice = ncread(file_path, 'pr', [1, day_index, 1], [Inf, 1, Inf]);
                    data_slice = permute(data_slice, [lat_dim_i lon_dim_i time_dim_i]);
                elseif time_dim_i == 3
                    data_slice = ncread(file_path, 'pr', [1, 1, day_index], [Inf, Inf, 1]);
                end

                % 处理数据切片以确保其地理范围正确
                data_slice = extent_handler(data_slice, lat(1), lon(1));
                
                % 将缺失值替换为 NaN , 假设所有文件中的 missing_value 和 _FillValue 是相同的
                data_slice(data_slice == missing_value) = NaN;
                data_slice(data_slice == fill_value) = NaN;

                % 重采样到目标分辨率,默认cmip6数据经纬方向是等间距格网
                %if (rows ~= 180/target_reso) || (cols ~= 360/target_reso)
                resampled_data_slice = resample_cmip6_data_slice(data_slice, rows, cols, target_reso, target_reso, Rsmp_method);
                %end

                % 存储处理后的数据
                alldata(:, :, i) = resampled_data_slice;
            end
        end
    end
    

end

function data_slice = extent_handler(data_slice, first_lat, first_lon)
    % 函数确保输出的地理范围为 -180 到 180 和 90 到 -90
    
    [num_rows, num_cols] = size(data_slice);
    
    % 转置数据，使经度为列，纬度为行
    if num_rows > num_cols
        data_slice = data_slice';
    end
    
    % 检查并翻转纬度
    if -90.5 < first_lat && first_lat < -88
        data_slice = flipud(data_slice);
    elseif first_lat > 88 && first_lat < 90.5
        % already in correct order
    else
        error('请检查纬度的范围');
    end
    
    % 检查并重新排列经度
    if -0.5 < first_lon && first_lon < 2
        if mod(num_cols, 2) ~= 0
            error('经度列数不能被2整除');
        end
        mid_cols = num_cols / 2;
        data_slice = [data_slice(:, mid_cols+1:end), data_slice(:, 1:mid_cols)];
    elseif -178 > first_lon && first_lon > -180.5
        % already in correct order
    else
        error('请检查经度的范围');
    end
end

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