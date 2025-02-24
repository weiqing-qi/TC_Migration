function MaketxtLinks(start_date, end_date, output_folder)
% https://disc2.gesdisc.eosdis.nasa.gov/data/MERGED_IR/GPM_MERGIR.1/2001/001/merg_2001010123_4km-pixel.nc4
% https://data.chc.ucsb.edu/products/CHIRP/daily/1981/chirp.1981.01.02.tif.gz
% MaketxtLinks('2002.01.08', '2006.01.05', 'D:\Desktop2\Attri-project\Programming\Matlab\Functions')
% 类似链接制作
    % 转换起始和结束日期为日期向量
    start_date_vec = datevec(start_date, 'yyyy.mm.dd');
    end_date_vec = datevec(end_date, 'yyyy.mm.dd');
    
    % 计算起始和结束日期之间的天数
    num_days = datenum(end_date_vec) - datenum(start_date_vec);
   
    % 检查输出文件夹是否存在，如果不存在则创建
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    
    % 定义输出文件的路径
    output_file_path = fullfile(output_folder, ['links',start_date,'-',end_date,'.txt']);
    
    % 打开文件准备写入
    fileID = fopen(output_file_path, 'w');
    
    % 生成每一天的链接并写入文件
    for day = 0:num_days
        % 当前日期
        current_date = addtodate(datenum(start_date_vec), day, 'day');
        current_date_vec = datevec(current_date);
        %%
        link = sprintf('https://persiann.eng.uci.edu/CHRSdata/PERSIANN/daily/ms6s4_d%02d%03d.bin.gz', ...
            mod(current_date_vec(1), 100), ...
            YMD_num(current_date_vec(1),current_date_vec(2),current_date_vec(3)) ...
            );


        %%
        % link = sprintf('https://data.chc.ucsb.edu/products/CHIRP/daily/%04d/chirp.%04d.%02d.%02d.tif.gz', ...
        %     current_date_vec(1), ...
        %     current_date_vec(1), ...
        %     current_date_vec(2), ...
        %     current_date_vec(3) ...
        %     );

        %% for hour=1:24
        % % 修正为确保月份和日期在链接中正确显示为两位数imergir
        % link = sprintf(['https://disc2.gesdisc.eosdis.nasa.gov/data/MERGED_IR/GPM_MERGIR.1/' ...
        %     '%04d/%03d/merg_%04d%02d%02d%02d_4km-pixel.nc4'], ...
        %     current_date_vec(1), ...
        %     YMD_num(current_date_vec(1),current_date_vec(2),current_date_vec(3)), ...
        %     current_date_vec(1), ...
        %     current_date_vec(2), ...
        %     current_date_vec(3), ...
        %%     hour-1);
        % 写入文件

        fprintf(fileID, '%s\r\n', link);
        % end
    end
    
    % 关闭文件
    fclose(fileID);
end