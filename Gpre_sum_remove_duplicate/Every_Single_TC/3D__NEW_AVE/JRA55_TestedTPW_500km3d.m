%����ȫ��̨��500kmӰ�����ܿɽ�ˮ ǰ��������ӵİ汾  ����̨��Ľ�ˮ�����ǰһ��  ע����������ERA5-TCW
%update20201020 ע�����û��ȥ����ˮ�����е�nodata value JRA55�����TC issue����ERA5���ݽ������� 
clear; clc; close all;
header=['ncols           720';
        'nrows           360';
        'xllcorner      -180';
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value  -9999'];
ErrorDate=load('E:\DATA\doc\JRA-55\JRA55_TCerror_DATE_NA.txt');
    
TC_idir='E:\DATA\IBTrACS\NEW\IBTrACS.ALL.v04r00.nc';
files_idir ='F:\GlobalTCMask1\TCraster500km\';
raster_files = dir([files_idir,'*.txt']);
JRA5_idir='E:\JRA55TPW\JRA55TPW_Bilinear05\';
ERA5_idir='E:\DATA\ERA5\ERA5-Total-Column-Water\';
JRAfilename='JRA55_Daily_Bilinear_Ori05_TPW';
ERAfilename='ERA5_Total_Column_Water_on_single_levels_daily';

S_TC_Pre=zeros(8171,1);%ÿ��̨��Ľ�ˮ
NAmask = Basinmasks(0.5);  %EPmask&NAmask:6 
[Garea025,~] = Gridarea(0.25);
[Garea05,~] = Gridarea(0.5);
[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

parfor NO=8899:13476  %ע�����1979��ʼ
    
time=ISO_TIME(:,:,NO);
lat=LAT(:,NO);       %���ϵ����ǲ��Ƕ�Ӧ̨���ʱ���ǰ�����Ƿ������ģ�����
lat(isnan(lat) )=[];
lon=LON(:,NO);
lon(isnan(lon) )=[]; 

%һ��̨��ÿ��λ(ÿ��Сʱ)���ճ߶�ʱ�����ֵ��ʽget
year=cell(length(lat));%Ԥ�����ڴ�
month=cell(length(lat));
day=cell(length(lat));
numtime=zeros(length(lat),1);

  for i=1:length(lat)
    year(i)=cellstr((time(1:4,i))');%���ַ�����ת��Ϊcell���ַ�������:�ַ�����S�е�ÿ�зָ��Ϊcellϸ��Ԫ��C��һ��Ԫ��,��ɾ��S��ÿ��β���ո�
    month(i)=cellstr((time(6:7,i))');
    day(i)=cellstr((time(9:10,i))');
    numtime(i,1)=str2double([time(1:4,i)',time(6:7,i)',time(9:10,i)']);%��������ʽ����ʱ�䷽��ȥ��
  end
           
%��acrtxtͷ�ļ�����Ĥ %header   %1--ncols  %2--nrows  %3--xllcorner  %4--yllcorner  %5--cellsize  %6--NODATA_value
[~,inputheader,ras]= read_ARCascii ([files_idir,raster_files(NO-5305).name]); %ע���ļ�����Ȼ����������ʱ��˳���������ľͲ���Ӱ��  �����50�꿪ʼ����5305

%process 1 globalize the imcomplete ras(ter) ���������Сarctxt�����ȫ���
raster = globalize_the_imcomplete_raster(inputheader,ras);
TCpre =zeros(180/inputheader(5),360/inputheader(5)); %ÿһ�ζ�Ҫ���¸�ֵ,Ҫ������꽵ˮ�Ļ�����-999��ȫ��Ϳ�����

%process 2 calculate time of each grid from TCRoutedata
[weizhi,passyear,passmonth,passday,locationR,locationC] = raster_grid_TC_pass_time(inputheader,raster,year,month,day,lon,lat);
%ע������������Ǻ���0-360������ֱ�ӽ����ݼ������ȡ��������-180-180

%process 3 preread all txt--saving time
[DB,MB,YB,~,~,~] = daybackdayforward(str2double(char(day(1))),str2double(char(month(1))),str2double(char(year(1))));%��һ��֮ǰ������
[DBB,MBB,YBB,~,~,~] = daybackdayforward(str2double(DB),str2double(MB),str2double(YB));
[~,~,~,DF,MF,YF] = daybackdayforward(str2double(char(day(length(numtime)))),str2double(char(month(length(numtime)))),str2double(char(year(length(numtime)))));%��һ��֮�������
[~,~,~,DFF,MFF,YFF] = daybackdayforward(str2double(DF),str2double(MF),str2double(YF));
allday=[str2double([YBB,MBB,DBB]);str2double([YB,MB,DB]);numtime;str2double([YF,MF,DF]);str2double([YFF,MFF,DFF])];%��Ҫ����2021��1��1�յ�����
allday=unique(allday);%uniqueĬ��ȥ��֮����������
allpre=zeros(360,720,length(allday));%��ʼ��
for k=1:length(allday)
    iserror = ismember(allday(k),ErrorDate);
    STRdate = num2str(allday(k));
    
    if iserror == 1
        Amonthdata = ncread([ERA5_idir,ERAfilename,STRdate(1:6),'.nc'],'tcw');%ÿһ��allday��Ӧһ��allpre
        patch=ERA5dailyresample( flipud(Amonthdata(:,:,str2double(STRdate(7:8)))') ,Garea025); %��ȡ������ERA5ԭ���ݾ����ز�������׼����
        errorpre=imread([JRA5_idir,JRAfilename,STRdate,'.tif']);%��λmm
        errorpre(NAmask==6)=patch(NAmask==6); %����ERA5����JRA55
        allpre(:,:,k)=errorpre;
    else 
        allpre(:,:,k)=imread([JRA5_idir,JRAfilename,STRdate,'.tif']);%��λmm
    end

end

%process 4 calculate Precipitation of each grid 
  for j=1:(length(weizhi))
          
      if locationC(j)>180/inputheader(5)    %0��360to-180��180�ع鵽����ԭʼ�ĵ�������, ��֮�����weizhi0360��location-180180
        locationC(j)=locationC(j)-180/inputheader(5);
      else
        locationC(j)=locationC(j)+180/inputheader(5);
      end

      %�����Ѿ����л����ˣ��������ﲻ�����н���(ת��)
       
      %��ȫ��
      [dayB,monthB,yearB,dayF,monthF,yearF] = daybackdayforward(str2double(char(passday(j))),str2double(char(passmonth(j))),str2double(char(passyear(j))));
%       [dayBB,monthBB,yearBB,~,~,~] = daybackdayforward(str2double(dayB),str2double(monthB),str2double(yearB));
%       [~,~,~,dayFF,monthFF,yearFF] = daybackdayforward(str2double(dayF),str2double(monthF),str2double(yearF));
   
      NOdataexclute=zeros(3,1);
      
      pre2 =allpre(:,:, allday==str2double([yearB,monthB,dayB]) );
      NOdataexclute(1,1)=pre2(locationR(j),locationC(j));
      
      pre3 =allpre(:,:, allday==str2double([char(passyear(j)),char(passmonth(j)),char(passday(j))]) );
      NOdataexclute(2,1)=pre3(locationR(j),locationC(j));
      
      pre4 =allpre(:,:, allday==str2double([yearF,monthF,dayF]) );
      NOdataexclute(3,1)=pre4(locationR(j),locationC(j));
        
      NOdataexclute(NOdataexclute<=0)=[];%��ˮ����-9999ȥ��
      TCpre(locationR(j),locationC(j)) = mean(NOdataexclute);%-180180λ�����������360λ������
  end
Pcurrent=TCpre;%180-180
Pcurrent(Pcurrent<0)=0;%nonvalue���0�������
Pcurrent(isnan(Pcurrent))=0;%nonvalue���0�������
S_TC_Pre(NO-5305)=sum(sum(Pcurrent.*Garea05))/sum(sum(Garea05(Pcurrent>0)));%ÿ��̨��Ľ�ˮ�� ʮ��m3
% G_SUM_PRE =G_SUM_PRE + TCpre;   
end

%process 4 Final output
%OUTPUT(odir,'MSWEP_G_SUM_PRE_entir',header,G_SUM_PRE);
