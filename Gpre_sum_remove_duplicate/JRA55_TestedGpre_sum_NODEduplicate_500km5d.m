%����ȫ��̨��500kmӰ������ˮ���� ǰ��������ӵİ汾  ����̨��Ľ�ˮ�����ǰһ��  
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
JRA5_idir='E:\JRA55Daily05\';
ERA5_idir='E:\DATA\ERA5-Totalprecip\';
odir='F:\GlobalTCMask1\NODEduplicate_sum\JRA55origin\';
Sodir='F:\GlobalTCMask1\NODEduplicate_sum\JRA55origin\JRA55seasonOri\';
JRAfilename='JRA55_Daily_Bilinear_Ori05_TotalPre';
ERAfilename='ERA5_Total_precipitation_on_single_levels_daily';

G_SUM_PRE =zeros(360,720); 
UntilLastyearG =G_SUM_PRE;
Y=2020   ;  
sensonprebase=0;
S1=0;S2=0;S3=0;S4=0;
NAmask = Basinmasks(0.5);  %EPmask&NAmask:6 
[Garea,~] = Gridarea(0.25);

[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

for NO=6121:13501  %ע�����1958��ʼ
    
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
  
% %��ÿһ����ܽ�ˮ��40������彵ˮ���
%   if str2double(char(year(1))) > Y       %���Y��������һ��̨������,��ֿ϶���Խ��Խ��,������󳡿���̨��Ľ�ˮ�����ǰһ��
%   out = G_SUM_PRE - UntilLastyearG;
%   outO=[out(:,(180/0.5)+1:360/0.5),out(:,1:180/0.5)];
%   OUTPUT(odir,['JRA55_GSUMTCP_Origin',num2str(Y)],header,outO);
%   UntilLastyearG = G_SUM_PRE;
%       if str2double(char(year(1)))==2021 %40���ܼƽ�ˮ��� ����ֻ��37�������ϵ�����2020�����Ҫ���һ��37���G_SUM_PRE
%          break
%           %OUTPUT(odir,'MSWEP_G_SUM_PRE_1980_2019_40YEARS',header,G_SUM_PRE); 
%       end
%   Y = str2double(char(year(1)));
%   S1=0;S2=0;S3=0;S4=0;%ÿ�긴λ���ڿ���
%   end
%   
% %SeasonOri Calculate
% if str2double(char(month(1))) == 6 && S1==0
%    Sout = G_SUM_PRE - sensonprebase;
%    SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
%    %OUTPUT(Sodir,['JRA55_GSUMTCP_Origin',char(year(1)),'_S1'],header,SoutO); 
%    sensonprebase= G_SUM_PRE;
%    S1=1;
% end
% if str2double(char(month(1))) == 9 && S2==0
%    Sout = G_SUM_PRE - sensonprebase;
%    SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
%    %OUTPUT(Sodir,['JRA55_GSUMTCP_Origin',char(year(1)),'_S2'],header,SoutO);
%    sensonprebase= G_SUM_PRE;
%    S2=1; 
% end
% if str2double(char(month(1))) == 12 && S3==0
%    Sout = G_SUM_PRE - sensonprebase;
%    SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
%    %OUTPUT(Sodir,['JRA55_GSUMTCP_Origin',char(year(1)),'_S3'],header,SoutO);
%    sensonprebase= G_SUM_PRE;
%    S3=1;    
% end
% if str2double(char(month(1))) == 3 && S4==0
%    Sout = G_SUM_PRE - sensonprebase;
%    SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
%    %OUTPUT(Sodir,['JRA55_GSUMTCP_Origin',num2str(Y-1),'_S4'],header,SoutO);
%    sensonprebase= G_SUM_PRE;
%    S4=1;    
% end
            
%��acrtxtͷ�ļ�����Ĥ %header   %1--ncols  %2--nrows  %3--xllcorner  %4--yllcorner  %5--cellsize  %6--NODATA_value
[~,inputheader,ras]= read_ARCascii ([files_idir,raster_files(NO-5305).name]); %ע���ļ�����Ȼ����������ʱ��˳���������ľͲ���Ӱ��

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
          
      if locationC(j)>180/inputheader(5)    %0��360to-180��180�ع鵽����ԭʼ�ĵ�������, ��֮�����weizhi0360��location-180180
        locationC(j)=locationC(j)-180/inputheader(5);
      else
        locationC(j)=locationC(j)+180/inputheader(5);
      end

      %�����Ѿ����л����ˣ��������ﲻ�����н���(ת��)
       
      %��ȫ��pre1=ncread([MSWEP_idir,num2str(passyear(j)),num2str(passmonth(j)),'.nc'],'precipitation',[1,1,passday],[360/header(5),180/header(5),1]);ò�Ʋ���
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
        
      NOdataexclute(NOdataexclute<0)=0;%��ˮ����-9999ȥ��
      TCpre(weizhi(j)) = sum(NOdataexclute);%-180180λ�����������360λ������
  end

% G_SUM_PRE =G_SUM_PRE + TCpre;   
end

%process 4 Final output
%OUTPUT(odir,'MSWEP_G_SUM_PRE_entir',header,G_SUM_PRE);
