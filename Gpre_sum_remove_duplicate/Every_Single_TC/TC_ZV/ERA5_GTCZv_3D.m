%ERA5ÿһ��̨����Χ��ֱ���б�ļ���
%update202205010
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
ERA5_idir_200U='E:\DATA\ERA5\ERA5_UCW_200\';
ERA5_idir_200V='E:\DATA\ERA5\ERA5_VCW_200\';
ERA5_idir_850U='E:\DATA\ERA5\ERA5_UCW_850\';
ERA5_idir_850V='E:\DATA\ERA5\ERA5_VCW_850\';
ERAfn_200U='ERA5_U_Component_Of_Wind_on_200hpa_levels_daily';
ERAfn_200V='ERA5_V_Component_Of_Wind_on_200hpa_levels_daily';
ERAfn_850U='ERA5_U_Component_Of_Wind_on_850hpa_levels_daily';
ERAfn_850V='ERA5_V_Component_Of_Wind_on_850hpa_levels_daily';
odir='\';


[lmheader1,lmheader2,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 �����������ݵ�extend�����޸�
[Garea025,~] = Gridarea(0.25);
[Garea05,~] = Gridarea(0.5);
[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

S_TC_Zv=zeros(8196,1);%ÿ��̨��Ľ�ˮ
parfor NO=8899:13476 %ע�����1950��ʼ
    
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
TCZv =ones(180/inputheader(5),360/inputheader(5))*-9999; %ÿһ�ζ�Ҫ���¸�ֵ,Ҫ������꽵ˮ�Ļ�����-999��ȫ��Ϳ�����

%process 2 calculate time of each grid from TCRoutedata
[weizhi,passyear,passmonth,passday,locationR,locationC] = raster_grid_TC_pass_time(inputheader,raster,year,month,day,lon,lat);
%ע������������Ǻ���0-360������ֱ�ӽ����ݼ������ȡ��������-180-180

%process 3 preread all txt--saving time
[DB,MB,YB,~,~,~] = daybackdayforward(str2double(char(day(1))),str2double(char(month(1))),str2double(char(year(1))));%��һ��֮ǰ������
[~,~,~,DF,MF,YF] = daybackdayforward(str2double(char(day(length(numtime)))),str2double(char(month(length(numtime)))),str2double(char(year(length(numtime)))));%��һ��֮�������
allday=[str2double([YB,MB,DB]);numtime;str2double([YF,MF,DF])];%��Ҫ����2021��1��1�յ�����
allday=unique(allday);%uniqueĬ��ȥ��֮����������
allW_200U=zeros(360,720,length(allday));%��ʼ��
allW_200V=zeros(360,720,length(allday));%��ʼ��
allW_850U=zeros(360,720,length(allday));%��ʼ��
allW_850V=zeros(360,720,length(allday));%��ʼ��
All_Zv       =zeros(360,720,length(allday));%��ʼ��
for k=1:length(allday)
    STRdate = num2str(allday(k));
    
    Amonthdata_200U = ncread([ERA5_idir_200U,ERAfn_200U,STRdate(1:6),'.nc'],'u');%ÿһ��allday��Ӧһ��allpre
    Amonthdata_200V = ncread([ERA5_idir_200V,ERAfn_200V,STRdate(1:6),'.nc'],'v');%ÿһ��allday��Ӧһ��allpre
    Amonthdata_850U = ncread([ERA5_idir_850U,ERAfn_850U,STRdate(1:6),'.nc'],'u');%ÿһ��allday��Ӧһ��allpre
    Amonthdata_850V = ncread([ERA5_idir_850V,ERAfn_850V,STRdate(1:6),'.nc'],'v');%ÿһ��allday��Ӧһ��allpre
    
    allW_200U(:,:,k)=ERA5dailyresample( flipud(Amonthdata_200U(:,:,str2double(STRdate(7:8)))') ,Garea025); %��ȡ������ERA5ԭ���ݾ����ز�������׼���� unit:m/s
    allW_200V(:,:,k)=ERA5dailyresample( flipud(Amonthdata_200V(:,:,str2double(STRdate(7:8)))') ,Garea025); %��ȡ������ERA5ԭ���ݾ����ز�������׼���� unit:m/s
    allW_850U(:,:,k)=ERA5dailyresample( flipud(Amonthdata_850U(:,:,str2double(STRdate(7:8)))') ,Garea025); %��ȡ������ERA5ԭ���ݾ����ز�������׼���� unit:m/s
    allW_850V(:,:,k)=ERA5dailyresample( flipud(Amonthdata_850V(:,:,str2double(STRdate(7:8)))') ,Garea025); %��ȡ������ERA5ԭ���ݾ����ز�������׼���� unit:m/s
    
    Diff_UW=allW_200U(:,:,k)-allW_850U(:,:,k);
    Diff_VW=allW_200V(:,:,k)-allW_850V(:,:,k);
    All_Zv(:,:,k)=sqrt( (Diff_UW.^2) + (Diff_VW.^2) );%extra_process calculating |Zv| >0
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
   
     NOdataexclute=zeros(3,1);
      
      Zv1 =All_Zv(:,:, allday==str2double([yearB,monthB,dayB]) );
      NOdataexclute(1,1)=Zv1(locationR(j),locationC(j));
      
      Zv2 =All_Zv(:,:, allday==str2double([char(passyear(j)),char(passmonth(j)),char(passday(j))]) );
      NOdataexclute(2,1)=Zv2(locationR(j),locationC(j));
      
      Zv3 =All_Zv(:,:, allday==str2double([yearF,monthF,dayF]) );
      NOdataexclute(3,1)=Zv3(locationR(j),locationC(j));
      
        
      NOdataexclute(NOdataexclute<0)=[];%��ˮ����-9999ȥ��
      NOdataexclute(isnan(NOdataexclute))=[];
      TCZv(locationR(j),locationC(j)) = mean(NOdataexclute);%-180180��λ�ò��䣬ע��
  end
  
%OUTPUT(odir,['ERA5_Ori_SingleTC_pre','_SID',SID(:,NO)'],header,TCpre);
Zvcurrent=TCZv;
Zvcurrent(Zvcurrent<0)=0;%nonvalue���0�������
Zvcurrent(isnan(Zvcurrent))=0;%nonvalue���0�������
S_TC_Zv(NO-5305)=sum(sum(Zvcurrent.*Garea05))/sum(sum(Garea05(Zvcurrent>0)));%ÿ��̨��Ľ�ˮ�� ʮ��m3
end
