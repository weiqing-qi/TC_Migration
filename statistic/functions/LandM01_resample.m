function [re_mask] = LandM01_resample(mask, incs, ocs)
% AIGM01_resample 用最邻近像元重采样一个01AI掩膜
% mask = 90 90 -180 180  输出是 90 90 -180 180 二维数组
% rate>1 (downscale) 
% 目前最好ocs和incs能被180整除,不然采样会控制在之前范围内新精度最大能容纳的网格数
% 所以可能网格中心不会匹配

% 为输入栅格创建地理参考对象
latlim = [-90 90];  % 纬度限制
lonlim = [-180 180]; % 经度限制
rasterSize = round([180/incs 360/incs]); % 栅格尺寸
IN_R = georefcells(latlim, lonlim, rasterSize, 'ColumnsStartFrom', 'north', 'RowsStartFrom', 'west');

Preciser_rate = incs/ocs;

[re_mask,re_geoinfo] = georesize(mask,IN_R,Preciser_rate,'nearest');

if re_geoinfo.RasterSize(1) ~= rasterSize(1)*Preciser_rate || re_geoinfo.RasterSize(2) ~= rasterSize(2)*Preciser_rate
   re_mask="因分辨率不能整除而网格中心不匹配";
end

%不全如果要写的话需要分很多情况包括考虑到比之前更大的范围
% if re_geoinfo.RasterSize(1) < rasterSize(1)*Preciser_rate || re_geoinfo.RasterSize(2) < rasterSize(2)*Preciser_rate
%     H=[ IN_R.RasterSize(2);
%         IN_R.RasterSize(1);
%         IN_R.LongitudeLimits(1);
%         IN_R.LatitudeLimits(1);
%         ocs;
%         -9999];
%     [re_mask] = FillView2Global(H,re_mask);
% end

end

