%ERA5ÿһ��̨�����߽�ˮ����ÿ��ķ�λ̨�罵ˮ����ȡ   
%update20211102
clear; clc; close all;
%Initialize Matlab Parallel Computing Enviornment
%p=parpool(10);
%------------Done--------
header=['ncols           720';
        'nrows           360';
        'xllcorner      -180';
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value  -9999'];
    
TC_idir='E:\DATA\IBTrACS\NEW\IBTrACS.ALL.v04r00.nc';
files_idir ='F:\GlobalTCMask1\TCraster500km\';
raster_files = dir([files_idir,'*.txt']);
ERA5_idir='E:\DATA\ERA5-Totalprecip\';
odir='F:\GlobalTCMask1\single_pre\ERA5\';
ERAfilename='ERA5_Total_precipitation_on_single_levels_daily';

[lmheader1,lmheader2,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 �����������ݵ�extend�����޸�
[Garea025,~] = Gridarea(0.25);
[Garea05,~] = Gridarea(0.5);
[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

S_TC_Pre=zeros(8196,1);%ÿ��̨��Ľ�ˮ
S_Land_Pre=zeros(8196,1);%ÿ��̨���½�ؽ�ˮ
MaximumPre=zeros(8196,1);%ÿ��̨������ˮ
MaximumLPre=zeros(8196,1);%ÿ��̨������½�ؽ�ˮ
MaximumPre10grid=zeros(8196,1);%ÿ��̨������ˮ
MaximumLPre10grid=zeros(8196,1);%ÿ��̨������½�ؽ�ˮ

parfor NO=5306:13476 %ע�����1950��ʼ
    
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
    STRdate = num2str(allday(k));
    Amonthdata = ncread([ERA5_idir,ERAfilename,STRdate(1:6),'.nc'],'tprate');%ÿһ��allday��Ӧһ��allpre
    allpre(:,:,k)=ERA5dailyresample( flipud(Amonthdata(:,:,str2double(STRdate(7:8)))')*1000*86400 ,Garea025); %��ȡ������ERA5ԭ���ݾ����ز�������׼����
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
      TCpre(locationR(j),locationC(j)) = sum(NOdataexclute);%-180180��λ�ò��䣬ע��
  end
  
OUTPUT(odir,['ERA5_Ori_SingleTC_pre','_SID',SID(:,NO)'],header,TCpre);

Pcurrent=TCpre;

LPcurrent=TCpre;
LPcurrent(Clandmask<0)=0;%����ȥ��
LPcurrent(LPcurrent<0)=0;%nonvalue���0�������
Pcurrent(Pcurrent<0)=0;%nonvalue���0�������
OnecolumnP=reshape(Pcurrent,[],1);
OnecolumnLP=reshape(LPcurrent,[],1);
OnecolumnP=sort(OnecolumnP,'descend');%��Ϊ���������һ�����Կ��ܴ���Ӱ�죬������Խ�ˮ��Ϊ����
OnecolumnLP=sort(OnecolumnLP,'descend');

S_TC_Pre(NO-5305)=10^(-6)*sum(sum(Pcurrent.*Garea05));%ÿ��̨��Ľ�ˮ�� ʮ��m3
S_Land_Pre(NO-5305)=10^(-6)*sum(sum(LPcurrent.*Garea05));%ÿ��̨���½�ؽ�ˮ�� ʮ��m3
MaximumPre(NO-5305)=max(max(Pcurrent));%ÿ��̨������ˮǿ�� 0.5������mm/5d
MaximumLPre(NO-5305)=max(max(LPcurrent));%ÿ��̨���½�����ˮǿ�� 0.5������mm/5d
MaximumPre10grid(NO-5305)=mean(OnecolumnP(1:10));%ÿ��̨������10��ˮǿ��0.5������mm/5d
MaximumLPre10grid(NO-5305)=mean(OnecolumnLP(1:10));%ÿ��̨������10½�ؽ�ˮǿ��0.5������mm/5d
end
output=[S_TC_Pre,S_Land_Pre,MaximumPre,MaximumLPre,MaximumPre10grid,MaximumLPre10grid];
%delete(p);
