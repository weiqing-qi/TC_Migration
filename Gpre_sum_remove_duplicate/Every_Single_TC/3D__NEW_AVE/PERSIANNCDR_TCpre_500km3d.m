%逐年全球台风500km影响区降水总量 前后五天相加的版本  FOR PERSIANN
%update20210602 注意去除降水数据中的nodata 跨年台风的降水会算进前一年  
clear; clc; close all;
header=['ncols           720';
        'nrows           360';
        'xllcorner      -180';
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value  -9999'];
    
TC_idir='E:\DATA\IBTrACS\NEW\IBTrACS.since1980.v04r00.nc';
files_idir ='F:\GlobalTCMask1\TCraster500km\';
raster_files = dir([files_idir,'*.txt']);
PERSIANN_idir='G:\PERSIANNCDR83_20\PERSIANNCDR05Daily_Ori_fullsize\';
odir='F:\GlobalTCMask1\NODEduplicate_sum\PERSIANNCDRorigin\';
Sodir='F:\GlobalTCMask1\NODEduplicate_sum\PERSIANNCDRorigin\Season\';%季节降水
h=waitbar(0,'please wait');%progress bar

G_SUM_PRE =zeros(360,720); 
UntilLastyearG =G_SUM_PRE;
sensonprebase=0;
S1=0;S2=0;S3=0;S4=0;

Y=1983;  

[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

for NO=340:length(raster_files)%1983年第一场开始
    
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

%把每一年的总降水和20年的整体降水输出
  if str2double(char(year(1))) > Y       %这个Y还留着上一场台风的年分,年分肯定是越来越大,这样最后场跨年台风的降水会算进前一年
  out = G_SUM_PRE - UntilLastyearG;
  outO=[out(:,(180/0.5)+1:360/0.5),out(:,1:180/0.5)];%输出-180-180
  OUTPUT(odir,['PERSIANNCDR_GSUMTCP_Origin',num2str(Y)],header,outO);
  UntilLastyearG = G_SUM_PRE;
  
      if str2double(char(year(1)))==2021 %20年总计降水输出 但是理论上到不了2021最后还是要输出一下20年的G_SUM_PRE
         G_O=[G_SUM_PRE(:,(180/0.5)+1:360/0.5),G_SUM_PRE(:,1:180/0.5)];%输出-180-180
         %OUTPUT('F:\GlobalTCMask1\NODEduplicate_sum\SUMandAVE\','PERSIANNCDR_AVE_GTCP_Origin_1983_2020_38YEARS',header,G_O/38); %平均38年!!
         break %跳出循环
      end
  Y = str2double(char(year(1)));
  
  S1=0;S2=0;S3=0;S4=0;%每年复位季节开关
  end
  
%每季度降水输出 按台风第一次出现的时间算
if str2double(char(month(1))) == 6 && S1==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%输出-180-180
   OUTPUT(Sodir,['PERSIANNCDR_GSUMTCP_Origin',char(year(1)),'_S1'],header,SoutO); 
   sensonprebase= G_SUM_PRE;
   S1=1;
end
if str2double(char(month(1))) == 9 && S2==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%输出-180-180
   OUTPUT(Sodir,['PERSIANNCDR_GSUMTCP_Origin',char(year(1)),'_S2'],header,SoutO);
   sensonprebase= G_SUM_PRE;
   S2=1; 
end
if str2double(char(month(1))) == 12 && S3==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%输出-180-180
   OUTPUT(Sodir,['PERSIANNCDR_GSUMTCP_Origin',char(year(1)),'_S3'],header,SoutO);
   sensonprebase= G_SUM_PRE;
   S3=1;    
end
if str2double(char(month(1))) == 3 && S4==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%输出-180-180
   OUTPUT(Sodir,['PERSIANNCDR_GSUMTCP_Origin',num2str(Y-1),'_S4'],header,SoutO);
   sensonprebase= G_SUM_PRE;
   S4=1;    
end
            
%读acrtxt头文件和掩膜 %header   %1--ncols  %2--nrows  %3--xllcorner  %4--yllcorner  %5--cellsize  %6--NODATA_value
[~,inputheader,ras]= read_ARCascii ([files_idir,raster_files(NO).name]); 

%process 1 globalize the imcomplete ras(ter) 任意格网大小arctxt处理成全球的
raster = globalize_the_imcomplete_raster(inputheader,ras);
TCpre =zeros(180/inputheader(5),360/inputheader(5)); %每一次都要重新赋值,要想加求年降水的话不能-999，全零就可以了

%process 2 calculate time of each grid from TCRoutedata
[weizhi,passyear,passmonth,passday,locationR,locationC] = raster_grid_TC_pass_time(inputheader,raster,year,month,day,lon,lat);
%注意这里输出的是横向0-360而下面直接进数据集里面读取的是-180-180

%process 3 preread all txt--saving time
[DB,MB,YB,~,~,~] = daybackdayforward(str2double(char(day(1))),str2double(char(month(1))),str2double(char(year(1))));%台风第一天之前的两天
[DBB,MBB,YBB,~,~,~] = daybackdayforward(str2double(DB),str2double(MB),str2double(YB));
[~,~,~,DF,MF,YF] = daybackdayforward(str2double(char(day(length(numtime)))),str2double(char(month(length(numtime)))),str2double(char(year(length(numtime)))));%台风最后一天之后的两天
[~,~,~,DFF,MFF,YFF] = daybackdayforward(str2double(DF),str2double(MF),str2double(YF));
allday=[str2double([YBB,MBB,DBB]);str2double([YB,MB,DB]);numtime;str2double([YF,MF,DF]);str2double([YFF,MFF,DFF])];%需要补充2021年1月1日的数据
allday=unique(allday);%unique默认去重之后升序排序
allpre=zeros(360,720,length(allday));%初始化
for k=1:length(allday)
    [~,~,allpre(:,:,k)] = read_ARCascii([PERSIANN_idir,'PERSIANN_CDR_AW05d_',num2str(floor(allday(k)/10000)),YMD_num(allday(k)),'.txt']);%每一个allday对应一个allpre
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
      TCpre(weizhi(j)) = sum(NOdataexclute);%-180180位置提出来放在360位置里了
  end

G_SUM_PRE =G_SUM_PRE + TCpre;   

    str=['运行中...',num2str(NO/length(raster_files)*100,'%0.2f'),'%'];%progress bar sprintf('%8.2f',a) or '%03d'
    waitbar(NO/length(raster_files),h,str)%progress bar
end

