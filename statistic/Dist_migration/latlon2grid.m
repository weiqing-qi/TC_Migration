function [GR,GC] = latlon2grid(LON,LAT,size)
%给定任意经纬度返回网格索引LON-180-180 
%上90下-90左-180右180
GC=ceil((LON+180)/size);
if GC==0 
    GC=1; 
end

GR=ceil((90-LAT)/size);
if GR==0 
    GR=1; 
end
%如果变成随机就完美了
end
