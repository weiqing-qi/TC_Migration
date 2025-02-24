%IMERG 3天里面最大降水强度所在经纬度
clear; clc; close all;
%------------Done--------
header=['ncols           720';
        'nrows           360';
        'xllcorner      -180';
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value  -9999'];
    
TC_idir='G:\DATA\IBTrACS\NEW\IBTrACS.ALL.v04r00.nc';
files_idir ='G:\DATA\GlobalTCMask1\TCraster500km\';
raster_files = dir([files_idir,'*.txt']);
IMERG_idir='G:\DATA\IMERG_FV6_OrigAW05d\';
odir='G:\DATA\GlobalTCMask1\single_pre\IMERG_Single_Precip3d_mmd\';

[lmheader1,lmheader2,Clandmask] = read_ARCascii('F:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 根据输入数据的extend进行修改
[Garea025,~] = Gridarea(0.25);
[Garea05,~] = Gridarea(0.5);
[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);


MaxPreLAT=zeros(8171,1);%每场台风的最大降水强度---都是位置 
MaxPre1DLAT=zeros(8171,1);%每场台风的最大日降水强度
MaxPreLON=zeros(8171,1);%每场台风的最大降水强度---都是位置 
MaxPre1DLON=zeros(8171,1);%每场台风的最大日降水强度
MP=zeros(8171,1);%最大降水

Max10PreLAT=zeros(8171,1);%每场台风的最大降水强度---都是位置 10个网格
Max10Pre1DLAT=zeros(8171,1);%每场台风的最大日降水强度
Max10PreLON=zeros(8171,1);%每场台风的最大降水强度---都是位置 
Max10Pre1DLON=zeros(8171,1);%每场台风的最大日降水强度

parfor NO=11416:13476 %注意年分2001开始
    
time=ISO_TIME(:,:,NO);
lat=LAT(:,NO);       %从上到下是不是对应台风的时间从前到后还是反过来的？？？
lat(isnan(lat) )=[];
lon=LON(:,NO);
lon(isnan(lon) )=[]; 

%一场台风每点位(每三小时)的日尺度时间点数值形式get
year=cell(length(lat));%预分配内存
month=cell(length(lat));
day=cell(length(lat));
numtime=zeros(length(lat),1);

  for i=1:length(lat)
    year(i)=cellstr((time(1:4,i))');%将字符数组转换为cell型字符串数组:字符数组S中的每行分割成为cell细胞元组C的一个元素,并删除S的每行尾部空格
    month(i)=cellstr((time(6:7,i))');
    day(i)=cellstr((time(9:10,i))');
    numtime(i,1)=str2double([time(1:4,i)',time(6:7,i)',time(9:10,i)']);%以数字形式储存时间方便去重
  end
  
            
%读acrtxt头文件和掩膜 %header   %1--ncols  %2--nrows  %3--xllcorner  %4--yllcorner  %5--cellsize  %6--NODATA_value
[~,inputheader,ras]= read_ARCascii ([files_idir,raster_files(NO-5305).name]); %注意文件名虽然不连续但是时间顺序是连续的就不会影响

%process 1 globalize the imcomplete ras(ter) 任意格网大小arctxt处理成全球的
raster = globalize_the_imcomplete_raster(inputheader,ras);
TCpre =ones(180/inputheader(5),360/inputheader(5))*-9999; %每一次都要重新赋值,要想加求年降水的话不能-999，全零就可以了
TCpreMAX1D =ones(180/inputheader(5),360/inputheader(5))*-9999; %最大的一天的降水强度

%process 2 calculate time of each grid from TCRoutedata
[weizhi,passyear,passmonth,passday,locationR,locationC] = raster_grid_TC_pass_time(inputheader,raster,year,month,day,lon,lat);
%注意这里输出的是横向0-360而下面直接进数据集里面读取的是竖向-180-180

%process 3 preread all txt--saving time
[DB,MB,YB,~,~,~] = daybackdayforward(str2double(char(day(1))),str2double(char(month(1))),str2double(char(year(1))));%第一天之前的两天
[DBB,MBB,YBB,~,~,~] = daybackdayforward(str2double(DB),str2double(MB),str2double(YB));
[~,~,~,DF,MF,YF] = daybackdayforward(str2double(char(day(length(numtime)))),str2double(char(month(length(numtime)))),str2double(char(year(length(numtime)))));%第一天之后的两天
[~,~,~,DFF,MFF,YFF] = daybackdayforward(str2double(DF),str2double(MF),str2double(YF));
allday=[str2double([YBB,MBB,DBB]);str2double([YB,MB,DB]);numtime;str2double([YF,MF,DF]);str2double([YFF,MFF,DFF])];%需要补充2021年1月1日的数据
allday=unique(allday);%unique默认去重之后升序排序
allpre=zeros(360,720,length(allday));%初始化
for k=1:length(allday)
    [~,~,allpre(:,:,k)] = read_ARCascii([IMERG_idir,'IMERG_FV6_AW05d_',num2str(allday(k)),'.txt']);%每一个allday对应一个allpre
end

%process 4 calculate Precipitation of each grid 
  for j=1:(length(weizhi))
          
      if locationC(j)>180/inputheader(5)    %0》360to-180》180回归到数据原始文档里排序
        locationC(j)=locationC(j)-180/inputheader(5);
      else
        locationC(j)=locationC(j)+180/inputheader(5);
      end

      %后面已经行列互换了，所以这里不再行列交换(转置)
       
      [dayB,monthB,yearB,dayF,monthF,yearF] = daybackdayforward(str2double(char(passday(j))),str2double(char(passmonth(j))),str2double(char(passyear(j))));
      [dayBB,monthBB,yearBB,~,~,~] = daybackdayforward(str2double(dayB),str2double(monthB),str2double(yearB));
      [~,~,~,dayFF,monthFF,yearFF] = daybackdayforward(str2double(dayF),str2double(monthF),str2double(yearF));
   
     NOdataexclute=zeros(3,1);
%       pre1 =allpre(:,:, allday==str2double([yearBB,monthBB,dayBB]) );
%       NOdataexclute(1,1)=pre1(locationR(j),locationC(j));
      
      pre2 =allpre(:,:, allday==str2double([yearB,monthB,dayB]) );
      NOdataexclute(1,1)=pre2(locationR(j),locationC(j));
      
      pre3 =allpre(:,:, allday==str2double([char(passyear(j)),char(passmonth(j)),char(passday(j))]) );
      NOdataexclute(2,1)=pre3(locationR(j),locationC(j));
      
      pre4 =allpre(:,:, allday==str2double([yearF,monthF,dayF]) );
      NOdataexclute(3,1)=pre4(locationR(j),locationC(j));
      
%       pre5 =allpre(:,:, allday==str2double([yearFF,monthFF,dayFF]) );
%       NOdataexclute(5,1)=pre5(locationR(j),locationC(j));
        
      NOdataexclute(NOdataexclute<0)=[];%降水坏点-9999去除
      if max(NOdataexclute)>=0
      TCpreMAX1D(locationR(j),locationC(j))=max(NOdataexclute);
      end
      if mean(NOdataexclute)>=0
      TCpre(locationR(j),locationC(j)) = mean(NOdataexclute); %-180180，位置不变，注意
      end 
  end
  
TCpre(isnan(TCpre))=-9999;TCpreMAX1D(isnan(TCpre))=-9999;
TCpre(TCpre<0)=-9999;TCpreMAX1D(TCpre<0)=-9999;
OUTPUT(odir,['IMERGV6_Ori_SingleTC_3dPrecip_Ummd','_SID',SID(:,NO)'],header,TCpre);%MM/DAY
MP(NO-5305)=max(max(TCpre));

[~,I1] = sort(reshape(TCpre,[],1),'descend');%变成一列之后顺序还是不变的
[~,I2] = sort(reshape(TCpreMAX1D,[],1),'descend');
[LonCenter,LatCenter] = GridCenterLocation(0.5);
MaxPreLON(NO-5305)=LonCenter(I1(1));
MaxPreLAT(NO-5305)=LatCenter(I1(1));
MaxPre1DLON(NO-5305)=LonCenter(I2(1));
MaxPre1DLAT(NO-5305)=LatCenter(I2(1));

Max10PreLON(NO-5305)=mean(LonCenter(I1(1:10)));
Max10PreLAT(NO-5305)=mean(LatCenter(I1(1:10)));
Max10Pre1DLON(NO-5305)=mean(LonCenter(I2(1:10)));
Max10Pre1DLAT(NO-5305)=mean(LatCenter(I2(1:10)));
end

output=[MaxPreLON,MaxPreLAT,MaxPre1DLON,MaxPre1DLAT,Max10PreLON,Max10PreLAT,Max10Pre1DLON,Max10Pre1DLAT];
%delete(p);