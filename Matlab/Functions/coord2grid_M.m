function [WZ] = coord2grid_M(LON,LAT,cs)
%给定任意经纬度返回网格索引 INPUT -180-180 Nx1支持 
%上90下-90左-180右180
GC=ceil((LON+180)/cs);
GC(GC==0)=1;

GR=ceil((90-LAT)/cs);
GR(GR==0)=1;
%计算位置
WZ=(GC-1)*(180/cs)+GR;
end
