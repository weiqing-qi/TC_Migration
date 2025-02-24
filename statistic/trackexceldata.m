%excelͳ�Ʒ�������С���
IBTrACSidir='E:\DATA\IBTrACS\NEW\IBTrACS.since1980.v04r00.nc';
ALLIBTrACSidir='E:\DATA\IBTrACS\NEW\IBTrACS.ALL.v04r00.nc';
[SID,LAT,LON,ISO_TIME,BASIN,WMO_WIND,WMO_PRES,NAME,DIST2LAND,LANDFALL,STORM_SPEED,STORM_DIR] = IBTrACS_nc_entire_variable_r(ALLIBTrACSidir);
[Garea05,~] = Gridarea(0.5);
%% ____����ÿ��̨���½ʱ�����Сѹ��
data=zeros(3714,1);
for i=1:3714 %�������������Ӧ����3712
    a=LAT(:,i+5305);
    a(isnan(a))=[];
    data(i,1)=length(a)*3;
end
for i=1:3714 %�������������Ӧ����3712 ��½ʱ��
    a=LANDFALL(:,i+5305);
    a(isnan(a))=[];
    data(i,1)=length(find(a==0));
end
data=data*3;
for i=1:3714 %�������������Ӧ����3712 ��Сѹ��
    a=WMO_PRES(:,i+5305);
    data(i,1)=min(a);
end

%% ____����ÿ��̨��ƽ���ƶ��ٶ� �ܾ���km/��ʱ�䣨-3��h
data=zeros(8171,3);%��һ���� �ڶ���½�� �����к���
for i=5306:13476 %50-20 ÿ��̨��
    a=DIST2LAND(:,i);%landfall����һ�����ӣ�������һ���ӽ�½�صľ����Զ�����յ�½��Ӱ�죬0�
    a(isnan(a))=[];%½����Ϣ
    
    lat=LAT(:,i);
    lat(isnan(lat))=[];
    lon=LON(:,i);
    lon(isnan(lon))=[];
    dist=0;%ÿһ����Ҫ��ʼ��̨���ܾ���
    Landdist=0;%���������е�½��Ķ�����½�ؾ�������
    
    dountland=0;%����½���϶���
    for j= 1:length(lat)-1 
        dist=dist+SphereDist2([lon(j),lat(j)],[lon(j+1),lat(j+1)]);%��λKM
        if a(j)==0 || a(j+1)==0 %С��ɶ��Ҳ������
        Landdist=Landdist+SphereDist2([lon(j),lat(j)],[lon(j+1),lat(j+1)]); %��λKM
        dountland=dountland+1;
        end
    end
    data(i-5305,1)=dist/(length(lat)-1)/3; % total ave
    data(i-5305,2)=Landdist/dountland/3;%km/h land
    data(i-5305,3)=(dist-Landdist)/(length(lat)-1-dountland)/3;%km/h sea
end
data=zeros(8171,2);
for i=1:8196
    a=LANDFALL(:,i+5305);
    a(isnan(a))=[];
    b=STORM_SPEED(:,i+5305);
    b(isnan(b))=[];
    data(i,1)=mean(b(a==0));%½���ٶ�
    data(i,2)=mean(b(a~=0));%�����ٶ�
end

%% ____����ÿ��̨���½ʱ�䣬ȥ������3Сʱ���ڵĵ���С½��
data=zeros(3715,1);
for i=1:3715
    a=LANDFALL(:,i+5305);
    a(isnan(a))=[];
    b=find(a==0);
    c=zeros(length(b),1);
    for j=2:length(b)
        c(j,1)=b(j,1)-b(j-1,1);
    end
    data(i,1)=length(find(c==1));
end
data=data*3;
%% ____��ȡ������ʼ�ͽ�������
data=zeros(13501,2);
for i=5306:9019
    L=LAT(:,i);
    L(isnan(L))=[];
    a=LAT(1,i);
    b=LAT(length(L),i);
    data(i-5305,1)=a;
    data(i-5305,2)=b;
end
%% ____��ȡ̨�翪ʼ������ֺ���أ��������
data=zeros(2272,2);  %1950:5306 1980:9018 3712
for i= 2398:4669
    time=ISO_TIME(:,:,i);
    L=LAT(:,i);
    L(isnan(L))=[];
    a=(time(1:4,1))';
    b=(time(1:4,length(L)))';
    data(i-2397,1)=str2double(a);
    data(i-2397,2)=str2double(b);
end
data=cell(13501,2);
for i= 5306:9019
    time=BASIN(:,:,i);
    L=LAT(:,i);
    L(isnan(L))=[];
    data(i-5305,1)=cellstr((time(1:2,1))');
    data(i-5305,2)=cellstr((time(1:2,length(L)))');
end

data=zeros(13501,1);
for i= 5306:13501
    basin=BASIN(:,:,i);
    L=LAT(:,i);
    L(isnan(L))=[];
    basinnum=zeros(length(L),1);
    for j=1:length(L)%���е���ر�־��ȡ %SI:1 SP:2 SA:3 NI:4 WP:5 EP:6 NA:7
        switch basin(1:2,j)'
                 case 'SI'
                basinnum(j)=1;
                 case 'SP'
                basinnum(j)=2;
                 case 'SA'
                basinnum(j)=3;
                 case 'NI'
                basinnum(j)=4;
                 case 'WP'
                basinnum(j)=5;
                 case 'EP'
                basinnum(j)=6;
                 case 'NA'
                basinnum(j)=7;     
        end
    end
    data(i-5305)=mode(basinnum);
end

%% ____��������ͳ��
%��excel�е���data
data=zeros(8171,1);
year=zeros(8171,1);
out=zeros(71,1);
for o=1:1
  for i=1:71
    a=data(year==i+1949,o);
    a(isnan(a))=[];
    out(i,o)=mean(a);
  end
end
%��excel�е���basin��Ϣ
startLON=zeros(1);
startLAT=zeros(1);
startBASIN=cell(4484,1);
endLON=zeros(1);
endLAT=zeros(1);
endBASIN=cell(4484,1);

out=zeros(71,7);
for i=1:71
    CYData=endLON(year==i+1949); %ÿ�θ����
    CYlocation=startBASIN(year==i+1949);
    out(i,1)=mean(CYData(strcmp(CYlocation,'NI')));
    out(i,2)=mean(CYData(strcmp(CYlocation,'SI')));
    out(i,3)=mean(CYData(strcmp(CYlocation,'SP')));
    out(i,4)=mean(CYData(strcmp(CYlocation,'WP')));
    out(i,5)=mean(CYData(strcmp(CYlocation,'EP')));
    out(i,6)=mean(CYData(strcmp(CYlocation,'NA')));
    out(i,7)=mean(CYData(strcmp(CYlocation,'SA')));
end

for i=1:71
    CYData=startLON(year==i+1949); %ÿ�θ����
    CYlocation=startBASIN(year==i+1949);
    A=CYData(strcmp(CYlocation,'SP'));
    A(A<0)=A(A<0)+360;
    out(i,3)=mean(A);

end

%% ____��ȡ̨�������·���ظ���
datay=zeros(13501,1); 
datam=zeros(13501,1);%1950:5306 1980:9018 3712 
for i= 5306:13501
    time=ISO_TIME(:,:,i);
    L=LAT(:,i);
    L(isnan(L))=[];
    a=(time(1:4,1))';%year
    b=(time(6:7,1))';%month
    datay(i-5305,1)=str2double(a);
    datam(i-5305,1)=str2double(b);
end

statistic=zeros(71,5);%per S1 S2 S3 S4 year
for y=1950:2020
statistic(y-1949,5)=length(find(datay==y));
MON=datam(datay==y);
statistic(y-1949,1)=length(find(MON==3))+length(find(MON==4))+length(find(MON==5));
statistic(y-1949,2)=length(find(MON==6))+length(find(MON==7))+length(find(MON==8));
statistic(y-1949,3)=length(find(MON==9))+length(find(MON==10))+length(find(MON==11));
statistic(y-1949,4)=length(find(MON==12))+length(find(MON==1))+length(find(MON==2));

end

basincount=zeros(71,7);
for i=1:71%�����ⲻ����
    a=startBASIN(datay==Y+1949);
    basincount(i,1)=length(find(a=='NI'));
    basincount(i,2)=length(find(a=='SI'));
    basincount(i,3)=length(find(a=='SP'));
    basincount(i,4)=length(find(a=='WP'));
    basincount(i,5)=length(find(a=='EP'));
    basincount(i,6)=length(find(a=='NA'));
    basincount(i,7)=length(find(a=='SA'));

end

%% ��ȡ̨��½�غͺ���Ӱ�����
idir='F:\GlobalTCMask1\TCraster500km\';
raster_files = dir([files_idir,'*.txt']);
[lmheader1,lmheader2,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���0360
%Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 �����������ݵ�extend�����޸�
p=parpool(10);

areaout1=zeros(length(raster_files),1);
areaout2=zeros(length(raster_files),1);
parfor i=1:length(raster_files)
[~,inputheader,ras]= read_ARCascii ([idir,raster_files(i).name]);
raster = globalize_the_imcomplete_raster(inputheader,ras);%0-360   
[Garea05,~] = Gridarea(0.5);

Garea05(raster<=0)=0;%����̨����״
areaout1(i)=sum(sum(Garea05));%��Ӱ�����
Garea05(Clandmask<0)=0;%����̨��LAND��״
areaout2(i)=sum(sum(Garea05));%��½��Ӱ�����

end
delete (p)
%% ERA5��ȡ̨��½�غͺ���Ӱ���������ǿ
idir='F:\GlobalTCMask1\single_pre\ERA5\';
raster_files = dir([idir,'*.txt']);%-180-180
[lmheader1,lmheader2,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 �����������ݵ�extend�����޸�
all_in_one = Basinmasks(0.5);
all_in_one=[all_in_one(:,361:720),all_in_one(:,1:360)];%-180-180 �����������ݵ�extend�����޸�


areaout1=zeros(length(raster_files),1);%ȫ��
areaout2=zeros(length(raster_files),1);%������ǿ��5��30mmΪ�ֽ�
areaout1L=zeros(length(raster_files),1);%ȫ��
areaout2L=zeros(length(raster_files),1);%������ǿ��5��30mmΪ�ֽ�
areaout1B=zeros(length(raster_files),6);%ȫ��
areaout2B=zeros(length(raster_files),6);%������ǿ��5��30mmΪ�ֽ�
areaout1BL=zeros(length(raster_files),6);%ȫ��
areaout2BL=zeros(length(raster_files),6);%������ǿ��5��30mmΪ�ֽ�
for i=1:8171
[~,inputheader,TCpre]= read_ARCascii ([idir,raster_files(i).name]);
%------ȫ����ˮ
[Garea05,~] = Gridarea(0.5);
Garea05(TCpre<=0.1)=0;%����̨����״
areaout1(i)=sum(sum(Garea05));%�����

for j=1:6
areaout1B(i,j)=sum(sum(Garea05(all_in_one==j)));%������
end

Garea05(Clandmask<0)=0;%����̨��LAND��״
areaout1L(i)=sum(sum(Garea05));%��½��Ӱ����

for j=1:6
areaout1BL(i,j)=sum(sum(Garea05(all_in_one==j)));%���½�����
end

%-------------������ˮ
[Garea05,~] = Gridarea(0.5);
Garea05(TCpre<=30)=0;%����̨����״
areaout2(i)=sum(sum(Garea05));

for j=1:6
areaout2B(i,j)=sum(sum(Garea05(all_in_one==j)));%������
end

Garea05(Clandmask<0)=0;%����̨��LAND��״
areaout2L(i)=sum(sum(Garea05));%��½��Ӱ�����

for j=1:6
areaout2BL(i,j)=sum(sum(Garea05(all_in_one==j)));%���½�����
end

end

%% ����ͼ������ȡ
input=zeros(1);%excel���濽��
year=zeros(1);%��ʼ��
output=ones(200,71)*-999;%û����һ�곬��200̨�磬����ɾ��-999����
for y=1950:2020
  middle=input(year==y); 
  middle(middle==0)=[];  
  output(1:length(middle),y-1949)=middle;     
end

%% ��ȡÿ������һЩ
input=zeros(1);%excel���濽��
year=zeros(1);%��ʼ��
c=zeros(71,1);
for i=1950:2020
    a=input(year==i);
    b=sort(a,'descend');
    d=round(length(b)*0.5);
    c(i-1949,1)=sum(sum( b( 60:80 ,1 )));
end
plot(c);
%% ____����ÿ��̨������������
data=zeros(8200,2);
for i=5306:13500 
    a=WMO_WIND(:,i);
    b=str2double(ISO_TIME(1:4,1,i));
    data(i-5305,1)=max(a);
    data(i-5305,2)=b;
end






