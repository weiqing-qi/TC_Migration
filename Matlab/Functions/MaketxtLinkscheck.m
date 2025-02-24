function [link]=MaketxtLinkscheck(start_date, end_date, output_folder)
%link=MaketxtLinkscheck('2001.01.01', '2006.01.05', 'D:\DATA\MERGIR\20010101-20060105');
start_date_vec = datevec(start_date, 'yyyy.mm.dd');
end_date_vec = datevec(end_date, 'yyyy.mm.dd');
num_days = datenum(end_date_vec) - datenum(start_date_vec);
link=[];
for day = 0:num_days
    % 当前日期
    current_date = addtodate(datenum(start_date_vec), day, 'day');
    current_date_vec = datevec(current_date);
    if ~exist([output_folder,'\',sprintf('chirp.%04d.%02d.%02d.tif.gz',current_date_vec(1), current_date_vec(2), current_date_vec(3))])

        link = sprintf('https://data.chc.ucsb.edu/products/CHIRP/daily/%04d/chirp.%04d.%02d.%02d.tif.gz', ...
        current_date_vec(1), ...
        current_date_vec(1), ...
        current_date_vec(2), ...
        current_date_vec(3) ...
        );
    % for hour=1:24
    %     if ~exist([output_folder,'\',sprintf('merg_%04d%02d%02d%02d_4km-pixel.nc4', ...
    %             current_date_vec(1), current_date_vec(2), current_date_vec(3), hour-1)])
    % 
    %     link = [link;sprintf(['https://disc2.gesdisc.eosdis.nasa.gov/data/MERGED_IR/GPM_MERGIR.1/' ...
    %     '%04d/%03d/merg_%04d%02d%02d%02d_4km-pixel.nc4'], ...
    %     current_date_vec(1), ...
    %     YMD_num(current_date_vec(1),current_date_vec(2),current_date_vec(3)), ...
    %     current_date_vec(1), ...
    %     current_date_vec(2), ...
    %     current_date_vec(3), ...
    %     hour-1)];
    %     end
    end
end
disp(link)

% fileID = fopen('output.txt', 'w');
% 
% if fileID == -1
%     error('File could not be opened.');
% end
% 
% for idx = 1:size(link, 1)
%     fprintf(fileID, '%s\r\n', link(idx, :));
% end
% 
% fclose(fileID);
end