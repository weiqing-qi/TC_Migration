%Jra55 3���������ˮǿ�����ھ�γ��
%update20221129
clear; clc; close all;
%------------Done--------
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
JRA5_idir='E:\JRA55Daily05\';
ERA5_idir='E:\DATA\ERA5\ERA5-Totalprecip\';
odir='F:\GlobalTCMask1\single_pre\JRA55_Single_Precip3d_mmd\';
JRAfilename='JRA55_Daily_Bilinear_Ori05_TotalPre';
ERAfilename='ERA5_Total_precipitation_on_single_levels_daily';

NAmask = Basinmasks(0.5);  %EPmask&NAmask:6 
[Garea,~] = Gridarea(0.25);
[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

S_TC_MPre=zeros(8171,1);%ÿ��̨������ˮ
MaxPreLAT=zeros(8171,1);%ÿ��̨������ˮǿ��---����λ�� 
MaxPre1DLAT=zeros(8171,1);%ÿ��̨�������ս�ˮǿ��
MaxPreLON=zeros(8171,1);%ÿ��̨������ˮǿ��---����λ�� 
MaxPre1DLON=zeros(8171,1);%ÿ��̨�������ս�ˮǿ��

Max10PreLAT=zeros(8171,1);%ÿ��̨������ˮǿ��---����λ�� 10������
Max10Pre1DLAT=zeros(8171,1);%ÿ��̨�������ս�ˮǿ��
Max10PreLON=zeros(8171,1);%ÿ��̨������ˮǿ��---����λ�� 
Max10Pre1DLON=zeros(8171,1);%ÿ��̨�������ս�ˮǿ��

parfor NO=6121:13476 %ע�����1950��ʼ
% for NO=6121:6121    
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
[~,inputheader,ras]= read_ARCascii ([files_idir,raster_files(NO-5305).name]); %ע���ļ�����Ȼ����������ʱ��˳���������ľͲ���Ӱ��

%process 1 globalize the imcomplete ras(ter) ���������Сarctxt�����ȫ���
raster = globalize_the_imcomplete_raster(inputheader,ras);
TCpre =ones(180/inputheader(5),360/inputheader(5))*-9999; %ÿһ�ζ�Ҫ���¸�ֵ,Ҫ������꽵ˮ�Ļ�����-999��ȫ��Ϳ�����
TCpreMAX1D =ones(180/inputheader(5),360/inputheader(5))*-9999; %����һ��Ľ�ˮǿ��
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
        Amonthdata = ncread([ERA5_idir,ERAfilename,STRdate(1:6),'.nc'],'tprate');%ÿһ��allday��Ӧһ��allpre
        patch=ERA5dailyresample( flipud(Amonthdata(:,:,str2double(STRdate(7:8)))')*1000*86400 ,Garea); %��ȡ������ERA5ԭ���ݾ����ز�������׼����
        errorpre=imread([JRA5_idir,JRAfilename,STRdate,'.tif'])/8;%��λmm/8d
        errorpre(NAmask==6)=patch(NAmask==6); %����ERA5����JRA55
        allpre(:,:,k)=errorpre;
    else 
        allpre(:,:,k)=imread([JRA5_idir,JRAfilename,STRdate,'.tif'])/8;%��λmm/8d
    end
end

%process 4 calculate Precipitation of each grid 
  for j=1:(length(weizhi))
          
      if locationC(j)>180/inputheader(5)    %0��360to-180��180�ع鵽����ԭʼ�ĵ�������
        locationC(j)=locationC(j)-180/inputheader(5);
      else
        locationC(j)=locationC(j)+180/inputheader(5);
      end

      %�����Ѿ����л����ˣ��������ﲻ�����н���(ת��)
       
      [dayB,monthB,yearB,dayF,monthF,yearF] = daybackdayforward(str2double(char(passday(j))),str2double(char(passmonth(j))),str2double(char(passyear(j))));
%       [dayBB,monthBB,yearBB,~,~,~] = daybackdayforward(str2double(dayB),str2double(monthB),str2double(yearB));
%       [~,~,~,dayFF,monthFF,yearFF] = daybackdayforward(str2double(dayF),str2double(monthF),str2double(yearF));
   
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
        
      NOdataexclute(NOdataexclute<0)=[];%��ˮ����-9999ȥ��
      TCpreMAX1D(locationR(j),locationC(j))=max(NOdataexclute);
      TCpre(locationR(j),locationC(j)) = mean(NOdataexclute);%-180180��λ�ò��䣬ע��
  end
  
TCpre(isnan(TCpre))=-9999;TCpreMAX1D(isnan(TCpre))=-9999;
TCpre(TCpre<0)=-9999;TCpreMAX1D(TCpre<0)=-9999;
%OUTPUT(odir,['JRA555_Ori_SingleTC_3dPrecip_Ummd','_SID',SID(:,NO)'],header,TCpre);%MM/DAY

[~,I1] = sort(reshape(TCpre,[],1),'descend');%���һ��֮��˳���ǲ����
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

S_TC_MPre(NO-5305)=max(max(TCpre));
end

output=[S_TC_MPre,MaxPreLON,MaxPreLAT,MaxPre1DLON,MaxPre1DLAT,Max10PreLON,Max10PreLAT,Max10Pre1DLON,Max10Pre1DLAT];
%delete(p);