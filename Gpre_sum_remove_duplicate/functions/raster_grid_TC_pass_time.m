function [weizhi,passyear,passmonth,passday,locationR,locationC] = raster_grid_TC_pass_time(header,raster,year,month,day,lon,lat)
%calculate time of each grid from TCRoutedata
%update20201020

[locationR,locationC]= find(raster==1);
weizhi=find(raster==1);
passyear=cell(length(weizhi),1);
passmonth=cell(length(weizhi),1);
passday=cell(length(weizhi),1);

for i= 1:length(weizhi) %每一个掩膜网格
    GCLON=locationC(i)*header(5)-0.5*header(5);%locationC是0-360可以这样算
    GCLAT=90-locationR(i)*header(5)+0.5*header(5);
    minDist=40076;%赤道周长
    
    for j=1:length(lon)   %每一个台风位点 找最近距离对应的时间
        Dist=SphereDist([GCLON,GCLAT],[lon(j),lat(j)]); %单位KM lon是-180-180但是里面的正余弦函数结果一致
        if Dist<minDist
            minDist=Dist;
            Syear=year(j);    %找最近距离对应的时间
            Smonth=month(j);
            Sday=day(j);
        end
    end
    
   passyear(i)= Syear;%应该还缺一个年月日到每年第几天的转化  
   passmonth(i)= Smonth;
   passday(i)= Sday;
end
