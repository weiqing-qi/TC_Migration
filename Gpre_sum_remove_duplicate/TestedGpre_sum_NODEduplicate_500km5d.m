%����ȫ��̨��500kmӰ������ˮ���� ǰ��������ӵİ汾  ����̨��Ľ�ˮ�����ǰһ��  
%update20201020 ע�����û��ȥ����ˮ�����е�nodata value��Ϊmswepû�л���
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
MSWEP_idir='E:\DATA\MSWEP\MSWEP\MSWEP\';
odir='F:\GlobalTCMask1\NODEduplicate_sum\MSWEPoriginalextend\';
Sodir='F:\GlobalTCMask1\NODEduplicate_sum\MSWEPoriginalextend\SeasonOri\';

G_SUM_PRE =zeros(360,720); 
UntilLastyearG =G_SUM_PRE;
sensonprebase=0;
S1=0;S2=0;S3=0;S4=0;
Y=1979;  

[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

for NO=8899:9018%length(raster_files)% 3593
    
time=ISO_TIME(:,:,NO);
lat=LAT(:,NO);       %���ϵ����ǲ��Ƕ�Ӧ̨���ʱ���ǰ�����Ƿ������ģ�����
lat(isnan(lat) )=[];
lon=LON(:,NO);
lon(isnan(lon) )=[]; 

%һ��̨��ÿ��λ(ÿ��Сʱ)���ճ߶�ʱ�����ֵ��ʽget
year=cell(length(lat));%Ԥ�����ڴ�
month=cell(length(lat));
day=cell(length(lat));

  for i=1:length(lat)
    year(i)=cellstr((time(1:4,i))');%���ַ�����ת��Ϊcell���ַ�������:�ַ�����S�е�ÿ�зָ��Ϊcellϸ��Ԫ��C��һ��Ԫ��,��ɾ��S��ÿ��β���ո�
    month(i)=cellstr((time(6:7,i))');
    day(i)=cellstr((time(9:10,i))');
  end
  
%��ÿһ����ܽ�ˮ��40������彵ˮ���
  if str2double(char(year(1))) > Y       %���Y��������һ��̨������,��ֿ϶���Խ��Խ��,������󳡿���̨��Ľ�ˮ�����ǰһ��
  out = G_SUM_PRE - UntilLastyearG;
  %OUTPUT(odir,['MSWEP_G_SUM_PRE_',num2str(Y)],header,out);
  UntilLastyearG = G_SUM_PRE;
      if str2double(char(year(1)))==2020 %40���ܼƽ�ˮ��� ����ֻ��37�������ϵ�����2020�����Ҫ���һ��37���G_SUM_PRE
         %OUTPUT(odir,'MSWEP_G_SUM_PRE_1980_2019_40YEARS',header,G_SUM_PRE); 
      end
  Y = str2double(char(year(1)));
  S1=0;S2=0;S3=0;S4=0;%ÿ�긴λ���ڿ���
  end
  
%SeasonOri Calculate
if str2double(char(month(1))) == 6 && S1==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
   OUTPUT(Sodir,['MSWEP_GSUMTCP_Origin',char(year(1)),'_S1'],header,SoutO); 
   sensonprebase= G_SUM_PRE;
   S1=1;
end
if str2double(char(month(1))) == 9 && S2==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
   OUTPUT(Sodir,['MSWEP_GSUMTCP_Origin',char(year(1)),'_S2'],header,SoutO);
   sensonprebase= G_SUM_PRE;
   S2=1; 
end
if str2double(char(month(1))) == 12 && S3==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
   OUTPUT(Sodir,['MSWEP_GSUMTCP_Origin',char(year(1)),'_S3'],header,SoutO);
   sensonprebase= G_SUM_PRE;
   S3=1;    
end
if str2double(char(month(1))) == 3 && S4==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
   OUTPUT(Sodir,['MSWEP_GSUMTCP_Origin',num2str(Y-1),'_S4'],header,SoutO);
   sensonprebase= G_SUM_PRE;
   S4=1;    
end
            
%��acrtxtͷ�ļ�����Ĥ %header   %1--ncols  %2--nrows  %3--xllcorner  %4--yllcorner  %5--cellsize  %6--NODATA_value
[~,inputheader,ras]= read_ARCascii ([files_idir,raster_files(NO-5305).name]); %1950����

%process 1 globalize the imcomplete ras(ter) ���������Сarctxt�����ȫ���
raster = globalize_the_imcomplete_raster(inputheader,ras);
TCpre =zeros(180/inputheader(5),360/inputheader(5)); %ÿһ�ζ�Ҫ���¸�ֵ,Ҫ������꽵ˮ�Ļ�����-999��ȫ��Ϳ�����

%process 2 calculate time of each grid from TCRoutedata
[weizhi,passyear,passmonth,passday,locationR,locationC] = raster_grid_TC_pass_time(inputheader,raster,year,month,day,lon,lat);
%ע������������Ǻ���0-360������ֱ�ӽ����ݼ������ȡ��������-180-180

%process 3 calculate Precipitation of each grid 
  for j=1:(length(weizhi))
          
      if locationC(j)>180/inputheader(5)    %0��360to-180��180�ع鵽����ԭʼ�ĵ�������
        locationC(j)=locationC(j)-180/inputheader(5);
      else
        locationC(j)=locationC(j)+180/inputheader(5);
      end

      %�����Ѿ����л����ˣ��������ﲻ�����н���(ת��)
       
      %��ȫ��pre1=ncread([MSWEP_idir,num2str(passyear(j)),num2str(passmonth(j)),'.nc'],'precipitation',[1,1,passday],[360/header(5),180/header(5),1]);ò�Ʋ���
      [dayB,monthB,yearB,dayF,monthF,yearF] = daybackdayforward(str2double(char(passday(j))),str2double(char(passmonth(j))),str2double(char(passyear(j))));
      [dayBB,monthBB,yearBB,~,~,~] = daybackdayforward(str2double(dayB),str2double(monthB),str2double(yearB));
      [~,~,~,dayFF,monthFF,yearFF] = daybackdayforward(str2double(dayF),str2double(monthF),str2double(yearF));
   
      pre1=ncread([MSWEP_idir,yearBB,monthBB,'.nc'],'precipitation',[locationC(j),locationR(j),str2double(dayBB)],[1,1,1]);
      pre2=ncread([MSWEP_idir,yearB,monthB,'.nc'],'precipitation',[locationC(j),locationR(j),str2double(dayB)],[1,1,1]);%ע���ļ����������ά��˳����������λ�ý�������Ϊת��֮���б�Ϊ���б�Ϊ��
      pre3=ncread([MSWEP_idir,char(passyear(j)),char(passmonth(j)),'.nc'],'precipitation',[locationC(j),locationR(j),str2double(char(passday(j)))],[1,1,1]);%ע���ļ����������ά��˳��
      pre4=ncread([MSWEP_idir,yearF,monthF,'.nc'],'precipitation',[locationC(j),locationR(j),str2double(dayF)],[1,1,1]);%ע���ļ����������ά��˳��
      pre5=ncread([MSWEP_idir,yearFF,monthFF,'.nc'],'precipitation',[locationC(j),locationR(j),str2double(dayFF)],[1,1,1]);
           
      TCpre(weizhi(j)) = pre1+pre2+pre3+pre4+pre5;
  end

G_SUM_PRE =G_SUM_PRE + TCpre;   
end

%process 4 Final output
%OUTPUT(odir,'MSWEP_G_SUM_PRE_entir',header,G_SUM_PRE);
