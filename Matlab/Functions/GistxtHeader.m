function [header] = GistxtHeader(cs,ncols,nrows,xllcorner,yllcorner,NODATA_value)
%头文件 input number cs,ncols,nrows,xllcorner,yllcorner,NODATA_value
%最多5位小数
if nargin == 1
   ncols=360/cs;
   nrows=180/cs;
   xllcorner=-180;
   yllcorner=-90;
   NODATA_value=-9999;
end

header=['ncols         ',sprintf('% 10.0f',ncols);
        'nrows         ',sprintf('% 10.0f',nrows);
        'xllcorner     ',sprintf('% 10.5f',xllcorner);
        'yllcorner     ',sprintf('% 10.5f',yllcorner);
        'cellsize      ',sprintf('% 10.5f',cs);
        'NODATA_value  ',sprintf('% 10.0f',NODATA_value)];
end