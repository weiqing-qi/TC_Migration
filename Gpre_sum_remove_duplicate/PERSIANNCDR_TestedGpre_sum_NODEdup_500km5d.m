%����ȫ��̨��500kmӰ������ˮ���� ǰ��������ӵİ汾  FOR PERSIANN
%update20210602 ע��ȥ����ˮ�����е�nodata ����̨��Ľ�ˮ�����ǰһ��  
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
Sodir='F:\GlobalTCMask1\NODEduplicate_sum\PERSIANNCDRorigin\Season\';%���ڽ�ˮ
h=waitbar(0,'please wait');%progress bar

G_SUM_PRE =zeros(360,720); 
UntilLastyearG =G_SUM_PRE;
sensonprebase=0;
S1=0;S2=0;S3=0;S4=0;

Y=1983;  

[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

for NO=340:length(raster_files)%1983���һ����ʼ
    
time=ISO_TIME(:,:,NO);
lat=LAT(:,NO);       %���ϵ����ǲ��Ƕ�Ӧ̨���ʱ���ǰ�����Ƿ������ģ�����
lat(isnan(lat) )=[];
lon=LON(:,NO);
lon(isnan(lon) )=[]; 

%һ��̨��ÿ��λ(ÿ��Сʱ)���ճ߶�ʱ�����ֵ��ʽget
year=cell(length(lat),1);%Ԥ�����ڴ�
month=cell(length(lat),1);
day=cell(length(lat),1);
numtime=zeros(length(lat),1);

  for i=1:length(lat)
    year(i)=cellstr(time(1:4,i)');%���ַ�����ת��Ϊcell���ַ�������:�ַ�����S�е�ÿ�зָ��Ϊcellϸ��Ԫ��C��һ��Ԫ��,��ɾ��S��ÿ��β���ո�
    month(i)=cellstr(time(6:7,i)');
    day(i)=cellstr(time(9:10,i)');
    numtime(i,1)=str2double([time(1:4,i)',time(6:7,i)',time(9:10,i)']);%��������ʽ����ʱ�䷽��ȥ��
  end

%��ÿһ����ܽ�ˮ��20������彵ˮ���
  if str2double(char(year(1))) > Y       %���Y��������һ��̨������,��ֿ϶���Խ��Խ��,������󳡿���̨��Ľ�ˮ�����ǰһ��
  out = G_SUM_PRE - UntilLastyearG;
  outO=[out(:,(180/0.5)+1:360/0.5),out(:,1:180/0.5)];%���-180-180
  OUTPUT(odir,['PERSIANNCDR_GSUMTCP_Origin',num2str(Y)],header,outO);
  UntilLastyearG = G_SUM_PRE;
  
      if str2double(char(year(1)))==2021 %20���ܼƽ�ˮ��� ���������ϵ�����2021�����Ҫ���һ��20���G_SUM_PRE
         G_O=[G_SUM_PRE(:,(180/0.5)+1:360/0.5),G_SUM_PRE(:,1:180/0.5)];%���-180-180
         %OUTPUT('F:\GlobalTCMask1\NODEduplicate_sum\SUMandAVE\','PERSIANNCDR_AVE_GTCP_Origin_1983_2020_38YEARS',header,G_O/38); %ƽ��38��!!
         break %����ѭ��
      end
  Y = str2double(char(year(1)));
  
  S1=0;S2=0;S3=0;S4=0;%ÿ�긴λ���ڿ���
  end
  
%ÿ���Ƚ�ˮ��� ��̨���һ�γ��ֵ�ʱ����
if str2double(char(month(1))) == 6 && S1==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
   OUTPUT(Sodir,['PERSIANNCDR_GSUMTCP_Origin',char(year(1)),'_S1'],header,SoutO); 
   sensonprebase= G_SUM_PRE;
   S1=1;
end
if str2double(char(month(1))) == 9 && S2==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
   OUTPUT(Sodir,['PERSIANNCDR_GSUMTCP_Origin',char(year(1)),'_S2'],header,SoutO);
   sensonprebase= G_SUM_PRE;
   S2=1; 
end
if str2double(char(month(1))) == 12 && S3==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
   OUTPUT(Sodir,['PERSIANNCDR_GSUMTCP_Origin',char(year(1)),'_S3'],header,SoutO);
   sensonprebase= G_SUM_PRE;
   S3=1;    
end
if str2double(char(month(1))) == 3 && S4==0
   Sout = G_SUM_PRE - sensonprebase;
   SoutO=[Sout(:,(180/0.5)+1:360/0.5),Sout(:,1:180/0.5)];%���-180-180
   OUTPUT(Sodir,['PERSIANNCDR_GSUMTCP_Origin',num2str(Y-1),'_S4'],header,SoutO);
   sensonprebase= G_SUM_PRE;
   S4=1;    
end
            
%��acrtxtͷ�ļ�����Ĥ %header   %1--ncols  %2--nrows  %3--xllcorner  %4--yllcorner  %5--cellsize  %6--NODATA_value
[~,inputheader,ras]= read_ARCascii ([files_idir,raster_files(NO).name]); 

%process 1 globalize the imcomplete ras(ter) ���������Сarctxt�����ȫ���
raster = globalize_the_imcomplete_raster(inputheader,ras);
TCpre =zeros(180/inputheader(5),360/inputheader(5)); %ÿһ�ζ�Ҫ���¸�ֵ,Ҫ������꽵ˮ�Ļ�����-999��ȫ��Ϳ�����

%process 2 calculate time of each grid from TCRoutedata
[weizhi,passyear,passmonth,passday,locationR,locationC] = raster_grid_TC_pass_time(inputheader,raster,year,month,day,lon,lat);
%ע������������Ǻ���0-360������ֱ�ӽ����ݼ������ȡ����-180-180

%process 3 preread all txt--saving time
[DB,MB,YB,~,~,~] = daybackdayforward(str2double(char(day(1))),str2double(char(month(1))),str2double(char(year(1))));%̨���һ��֮ǰ������
[DBB,MBB,YBB,~,~,~] = daybackdayforward(str2double(DB),str2double(MB),str2double(YB));
[~,~,~,DF,MF,YF] = daybackdayforward(str2double(char(day(length(numtime)))),str2double(char(month(length(numtime)))),str2double(char(year(length(numtime)))));%̨�����һ��֮�������
[~,~,~,DFF,MFF,YFF] = daybackdayforward(str2double(DF),str2double(MF),str2double(YF));
allday=[str2double([YBB,MBB,DBB]);str2double([YB,MB,DB]);numtime;str2double([YF,MF,DF]);str2double([YFF,MFF,DFF])];%��Ҫ����2021��1��1�յ�����
allday=unique(allday);%uniqueĬ��ȥ��֮����������
allpre=zeros(360,720,length(allday));%��ʼ��
for k=1:length(allday)
    [~,~,allpre(:,:,k)] = read_ARCascii([PERSIANN_idir,'PERSIANN_CDR_AW05d_',num2str(floor(allday(k)/10000)),YMD_num(allday(k)),'.txt']);%ÿһ��allday��Ӧһ��allpre
end


%process 4 calculate Precipitation of each grid 
  for j=1:(length(weizhi))
          
      if locationC(j)>180/inputheader(5)    %0��360to-180��180�ع鵽����ԭʼ�ĵ������� ��֮�����weizhi0360��location-180180
        locationC(j)=locationC(j)-180/inputheader(5);
      else
        locationC(j)=locationC(j)+180/inputheader(5);
      end
       
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

G_SUM_PRE =G_SUM_PRE + TCpre;   

    str=['������...',num2str(NO/length(raster_files)*100,'%0.2f'),'%'];%progress bar sprintf('%8.2f',a) or '%03d'
    waitbar(NO/length(raster_files),h,str)%progress bar
end

