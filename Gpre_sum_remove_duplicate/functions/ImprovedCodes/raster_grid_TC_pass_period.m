function [weizhi,passyear,passmonth,passday,locationR,locationC] = raster_grid_TC_pass_period(header,raster,year,month,day,lon,lat,ran)
%This is an update version of calculate time of each grid from TCRoutedata 
%The passing time is calculated rigidly, including every period it appeared
%in 500 km. And we add the former & later 1 day.
%update2023.9.6
%Raster is 0-360 

[locationR,locationC]= find(raster==1);
weizhi=find(raster==1);

Ltrac = length(lon);
LG=length(weizhi);

%每一行一个网格，每列记录经过的时间，最多为Ltrac，取3Ltrac为了方便后面BF1day衔接
%拓展的2Ltrac保证附加BF1day时绝对不超出数组长度
passyear=cell(LG,3*Ltrac);                     
passmonth=cell(LG,3*Ltrac);
passday=cell(LG,3*Ltrac);
                         
for i= 1:LG                                                                %每一个掩膜网格
    %计算500km网格中心经纬度
    GCLON=locationC(i)*header(5)-0.5*header(5);
    GCLAT=90-locationR(i)*header(5)+0.5*header(5);
    MGLON=GCLON*ones(Ltrac,1);                                             %统一SphereDist_Matrix输入的维度
    MGLAT=GCLAT*ones(Ltrac,1);

    %计算每个轨迹点到当前网格中心距离(Spheredist这个函数0-360或-180-180都可以用结果不变)
    Dist=SphereDist_Matrix([MGLON,MGLAT],[lon,lat]);                       %单位KM 列向量
    
    
    Lpass=length(find(Dist<=ran));
    passyear (i,1:Lpass) = year (Dist<=ran)';                              %有可能依然找不到任何值，这行时间可能为空
    passmonth(i,1:Lpass) = month(Dist<=ran)';                              %Lpass=0时不会进行赋值
    passday  (i,1:Lpass) = day  (Dist<=ran)'; 
 
    %每行每个时间断点
    Index=find(Dist<=ran);                                                 %得到TC经过这个点位的日期的全部索引
    IndexB=unique(Index-(1:Lpass)');                                       %变成数值like:000333555，连续部分索引等差1,去重然后进行正反搜寻
    Ls=length(IndexB);
    
    BF1D = zeros(1,2*Ls);                                                  %每一段连续经历两端都要加一天，所以是两倍
    BFM  = zeros(1,2*Ls);
    BFY  = zeros(1,2*Ls);

    for j= 1:Ls
        %-------正序断点加Backday前一天
        Ord = find(index==IndexB(j),1,'first');                            
        [BF1D(2*j-1),BFM(2*j-1),BFY(2*j-1),~,~,~]=daybackdayforward( ...   
            str2double(char(day(Ord))), ...
            str2double(char(month(Ord))), ...
            str2double(char(year(Ord))) );
        %-------反序断点加Forwardday后一天
        inOrd= find(index==IndexB(j),1,'last' ); 
        [~,~,~,BF1D(2*j),BFM(2*j),BFY(2*j)]=daybackdayforward( ...         
            str2double(char(day(inOrd))), ...
            str2double(char(month(inOrd))), ...
            str2double(char(year(inOrd))) );
    end

    %将BF1day的日期连接到passdmy数组后，外部函数包含去重这里就不去重了
    passyear(i,Ltrac+1:Ltrac+2*Ls)=BFY;  %Lpass<Ltrac                   
    passmonth(i,Ltrac+1:Ltrac+2*Ls)=BFM;
    passday(i,Ltrac+1:Ltrac+2*Ls)=BF1D;
end
end
