%逐年全球台风500km影响区降水总量 前后五天相加的版本    
%update20201020
clear; clc; close all;

header=['ncols           720';
        'nrows           360';
        'xllcorner         0';
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value   -999'];
    
TC_idir='E:\DATA\IBTrACS\IBTrACS.since1980.v04r00.nc';
files_idir ='F:\GlobalTCMask1\TCraster500km\';
raster_files = dir([files_idir,'*.txt']);
MSWEP_idir='E:\DATA\MSWEP\MSWEP\MSWEP\';
odir='F:\GlobalTCMask1\sum\';

G_SUM_PRE =zeros(360,720); 
Y=0;  %Hist_val=zeros(1,360*720); %YEAR,MONTH,DAY,WEIZHI

[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

for NO=1:1 %length(raster_files)
    
time=ISO_TIME(:,:,NO);
lat=LAT(:,NO);       %从上到下是不是对应台风的时间从前到后还是反过来的？？？
lat(isnan(lat) )=[];
lon=LON(:,NO);
lon(isnan(lon) )=[]; 

%一场台风每点位(每三小时)的日尺度时间点数值形式get
year=cell(length(lat));%预分配内存
month=cell(length(lat));
day=cell(length(lat));

  for i=1:length(lat)
    year(i)=cellstr((time(1:4,i))');%将字符数组转换为cell型字符串数组:字符数组S中的每行分割成为cell细胞元组C的一个元素,并删除S的每行尾部空格
    month(i)=cellstr((time(6:7,i))');
    day(i)=cellstr((time(9:10,i))');
  end
  
%赋值一个Hist_val去重记录矩阵  每年重新开始去重，不然太多，在每年年分赋值之后开始      %把上一年的总降水输出
  if str2double(char(year(1))) > Y      %这个Y还留着上一场台风的年分,年分肯定是越来越大的
      
  % Hist_val=zeros(100,360*720); %YEAR,MONTH,DAY,WEIZHI  一年大概100场台风但是不是每一个台风都会在同一天经过同一个点
  
  OUTPUT(odir,['MSWEP_G_SUM_PRE_',num2str(Y)],header,G_SUM_PRE);  %把上一年的总降水输出
  end
  Y = str2double(char(year(1))); %本来是year(length(lat))但是这样如果最后一个台风跨年的话Hist_val会不更新，year(1)可能会有重复但是影响不大
          
%读acrtxt头文件和掩膜 %header   %1--ncols  %2--nrows  %3--xllcorner  %4--yllcorner  %5--cellsize  %6--NODATA_value
[~,inputheader,ras]= read_ARCascii ([files_idir,raster_files(NO).name]); 

%process 1 globalize the imcomplete ras(ter) 任意格网大小arctxt处理成全球的
raster = globalize_the_imcomplete_raster(inputheader,ras);
% TCpre =zeros(180/inputheader(5),360/inputheader(5)); %每一次都要重新赋值

%process 2 calculate time of each grid from TCRoutedata
[weizhi,passyear,passmonth,passday,locationR,locationC] = raster_grid_TC_pass_time(inputheader,raster,year,month,day,lon,lat);
%注意这里输出的是横向0-360而下面直接进数据集里面读取的是竖向-180-180

%process 3 calculate Precipitation of each grid 
p = parpool(15) %线程
% p.IdleTimeout = 100; %最大时间
%--------k开启进程池并设定
  parfor j=1:(length(weizhi))
      if j==1 
          Hist_val=zeros(100,360*720); 
          TCpre =zeros(180/inputheader(5),360/inputheader(5));
      end
          
      if locationC(j)>180/0.5    %0》360to-180》180回归到数据原始文档里排序 应该是inputheader(5)的为了parfor
        locationC(j)=locationC(j)-180/0.5;
      else
        locationC(j)=locationC(j)+180/0.5;
      end

      %后面已经行列互换了，所以这里不再行列交换(转置)
       
      %读全球pre1=ncread([MSWEP_idir,num2str(passyear(j)),num2str(passmonth(j)),'.nc'],'precipitation',[1,1,passday],[360/header(5),180/header(5),1]);
      [dayB,monthB,yearB,dayF,monthF,yearF] = daybackdayforward(str2double(char(passday(j))),str2double(char(passmonth(j))),str2double(char(passyear(j))));
      [dayBB,monthBB,yearBB,~,~,~] = daybackdayforward(str2double(dayB),str2double(monthB),str2double(yearB));
      [~,~,~,dayFF,monthFF,yearFF] = daybackdayforward(str2double(dayF),str2double(monthF),str2double(yearF));
   
      pre1=ncread([MSWEP_idir,yearBB,monthBB,'.nc'],'precipitation',[locationC(j),locationR(j),str2double(dayBB)],[1,1,1]);
      pre2=ncread([MSWEP_idir,yearB,monthB,'.nc'],'precipitation',[locationC(j),locationR(j),str2double(dayB)],[1,1,1]);%注意文件里面的行列维度顺序这里行列位置交换了因为转置之后行变为列列变为行
      pre3=ncread([MSWEP_idir,char(passyear(j)),char(passmonth(j)),'.nc'],'precipitation',[locationC(j),locationR(j),str2double(char(passday(j)))],[1,1,1]);%注意文件里面的行列维度顺序
      pre4=ncread([MSWEP_idir,yearF,monthF,'.nc'],'precipitation',[locationC(j),locationR(j),str2double(dayF)],[1,1,1]);%注意文件里面的行列维度顺序
      pre5=ncread([MSWEP_idir,yearFF,monthFF,'.nc'],'precipitation',[locationC(j),locationR(j),str2double(dayFF)],[1,1,1]);
      
      %去重
      [pre1,pre2,pre3,pre4,pre5,Hist_val] = historical_deduplicate(pre1,pre2,pre3,pre4,pre5,yearBB,monthBB,dayBB,yearB,monthB,dayB,str2double(char(passyear(j))),str2double(char(passmonth(j))),str2double(char(passday(j))),yearF,monthF,dayF,yearFF,monthFF,dayFF,Hist_val,weizhi(j));    
           
      TCpre(weizhi(j)) = pre1+pre2+pre3+pre4+pre5;
  end
  delete(p)

G_SUM_PRE =G_SUM_PRE + TCpre;   
end

%process 4 Final output
%OUTPUT(odir,'MSWEP_G_SUM_PRE_entir',header,G_SUM_PRE);