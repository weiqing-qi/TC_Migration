%逐年全球台风500km影响区降水总量 前后五天相加的版本  FOR IMERGAW05D
%update20201020 注意这个没有去除降水数据中的nodata value因为mswep没有坏点  跨年台风的降水会算进前一年  
clear; clc; close all;
header=['ncols           720';
        'nrows           360';
        'xllcorner      -180';
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value  -9999'];
    
TC_idir='E:\DATA\IBTrACS\NEW\IBTrACS.ALL.v04r00.nc';
files_idir ='F:\GlobalTCMask1\TCraster500km\';
raster_files = dir([files_idir,'*.txt']);
IMERG_idir='E:\IMERG_FV6_OrigAW05d\';
odir='F:\GlobalTCMask1\single_pre\IMERG\';

[lmheader1,lmheader2,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 根据输入数据的extend进行修改
[Garea025,~] = Gridarea(0.25);
[Garea05,~] = Gridarea(0.5);
[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

S_TC_Pre=zeros(8196,1);%每场台风的降水
S_Land_Pre=zeros(8196,1);%每场台风的陆地降水
MaximumPre=zeros(8196,1);%每场台风的最大降水
MaximumLPre=zeros(8196,1);%每场台风的最大陆地降水
MaximumPre10grid=zeros(8196,1);%每场台风的最大降水
MaximumLPre10grid=zeros(8196,1);%每场台风的最大陆地降水
PreRate=zeros(8196,1);%每场台风的全平均降水强度

parfor NO=11416:13476%2001年第一场开始 后面要改的
    
time=ISO_TIME(:,:,NO);
lat=LAT(:,NO);       %从上到下是不是对应台风的时间从前到后还是反过来的？？？
lat(isnan(lat) )=[];
lon=LON(:,NO);
lon(isnan(lon) )=[]; 

%一场台风每点位(每三小时)的日尺度时间点数值形式get
year=cell(length(lat),1);%预分配内存
month=cell(length(lat),1);
day=cell(length(lat),1);
numtime=zeros(length(lat),1);

  for i=1:length(lat)
    year(i)=cellstr(time(1:4,i)');%将字符数组转换为cell型字符串数组:字符数组S中的每行分割成为cell细胞元组C的一个元素,并删除S的每行尾部空格
    month(i)=cellstr(time(6:7,i)');
    day(i)=cellstr(time(9:10,i)');
    numtime(i,1)=str2double([time(1:4,i)',time(6:7,i)',time(9:10,i)']);%以数字形式储存时间方便去重
  end 
            
%读acrtxt头文件和掩膜 %header   %1--ncols  %2--nrows  %3--xllcorner  %4--yllcorner  %5--cellsize  %6--NODATA_value
[~,inputheader,ras]= read_ARCascii ([files_idir,raster_files(NO-5305).name]); %把2000年的放起来不要，从2001年开始

%process 1 globalize the imcomplete ras(ter) 任意格网大小arctxt处理成全球的
raster = globalize_the_imcomplete_raster(inputheader,ras);
TCpre =zeros(180/inputheader(5),360/inputheader(5)); %每一次都要重新赋值,要想加求年降水的话不能-999，全零就可以了

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
          
      if locationC(j)>180/inputheader(5)    %0》360to-180》180回归到数据原始文档里排序 这之后代表weizhi0360和location-180180
        locationC(j)=locationC(j)-180/inputheader(5);
      else
        locationC(j)=locationC(j)+180/inputheader(5);
      end
       
      %读全球pre1=ncread([MSWEP_idir,num2str(passyear(j)),num2str(passmonth(j)),'.nc'],'precipitation',[1,1,passday],[360/header(5),180/header(5),1]);貌似不对
      [dayB,monthB,yearB,dayF,monthF,yearF] = daybackdayforward(str2double(char(passday(j))),str2double(char(passmonth(j))),str2double(char(passyear(j))));
      [dayBB,monthBB,yearBB,~,~,~] = daybackdayforward(str2double(dayB),str2double(monthB),str2double(yearB));
      [~,~,~,dayFF,monthFF,yearFF] = daybackdayforward(str2double(dayF),str2double(monthF),str2double(yearF));
      
      NOdataexclute=zeros(5,1);
      pre1 =allpre(:,:, allday==str2double([yearBB,monthBB,dayBB]) );
      NOdataexclute(1,1)=pre1(locationR(j),locationC(j));
      
      pre2 =allpre(:,:, allday==str2double([yearB,monthB,dayB]) );
      NOdataexclute(2,1)=pre2(locationR(j),locationC(j));
      
      pre3 =allpre(:,:, allday==str2double([char(passyear(j)),char(passmonth(j)),char(passday(j))]) );
      NOdataexclute(3,1)=pre3(locationR(j),locationC(j));
      
      pre4 =allpre(:,:, allday==str2double([yearF,monthF,dayF]) );
      NOdataexclute(4,1)=pre4(locationR(j),locationC(j));
      
      pre5 =allpre(:,:, allday==str2double([yearFF,monthFF,dayFF]) );
      NOdataexclute(5,1)=pre5(locationR(j),locationC(j));
        
      NOdataexclute(NOdataexclute<0)=0;%降水坏点-9999去除
      TCpre(locationR(j),locationC(j)) = sum(NOdataexclute);%-180180位置
  end

  OUTPUT(odir,['IMERG_Ori_SingleTC_pre','_SID',SID(:,NO)'],header,TCpre);

Pcurrent=TCpre;
LPcurrent=TCpre;
LPcurrent(Clandmask<0)=0;%海洋去除
LPcurrent(LPcurrent<0)=0;%nonvalue变成0可以求和
Pcurrent(Pcurrent<0)=0;%nonvalue变成0可以求和
OnecolumnP=reshape(Pcurrent,[],1);
OnecolumnLP=reshape(LPcurrent,[],1);
OnecolumnP=sort(OnecolumnP,'descend');%因为网格面积不一样所以可能存在影响，最好是以降水量为依据
OnecolumnLP=sort(OnecolumnLP,'descend');

S_TC_Pre(NO-11415)=10^(-6)*sum(sum(Pcurrent.*Garea05));%每场台风的降水量 十亿m3
S_Land_Pre(NO-11415)=10^(-6)*sum(sum(LPcurrent.*Garea05));%每场台风的陆地降水量 十亿m3
MaximumPre(NO-11415)=max(max(Pcurrent));%每场台风的最大降水强度 0.5°网格mm/5d
MaximumLPre(NO-11415)=max(max(LPcurrent));%每场台风的陆地最大降水强度 0.5°网格mm/5d
MaximumPre10grid(NO-11415)=mean(OnecolumnP(1:10));%每场台风的最大10降水强度0.5°网格mm/5d
MaximumLPre10grid(NO-11415)=mean(OnecolumnLP(1:10));%每场台风的最大10陆地降水强度0.5°网格mm/5d
area=sum(sum(Garea05(Pcurrent>1)));
t=3*length(lat);
PreRate(NO-11415)=S_TC_Pre(NO-11415)*10^6/area/t;
end
output=[S_TC_Pre,S_Land_Pre,MaximumPre,MaximumLPre,MaximumPre10grid,MaximumLPre10grid];