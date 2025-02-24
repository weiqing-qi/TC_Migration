%excel统计分析辅助小组件
IBTrACSidir='E:\DATA\IBTrACS\NEW\IBTrACS.since1980.v04r00.nc';
ALLIBTrACSidir='E:\DATA\IBTrACS\NEW\IBTrACS.ALL.v04r00.nc';
[SID,LAT,LON,ISO_TIME,BASIN,WMO_WIND,WMO_PRES,NAME,DIST2LAND,LANDFALL,STORM_SPEED,STORM_DIR] = IBTrACS_nc_entire_variable_r(ALLIBTrACSidir);
[Garea05,~] = Gridarea(0.5);
%% ____计算每场台风登陆时间和最小压力
data=zeros(3714,1);
for i=1:3714 %多出来两个检验应该是3712
    a=LAT(:,i+5305);
    a(isnan(a))=[];
    data(i,1)=length(a)*3;
end
for i=1:3714 %多出来两个检验应该是3712 登陆时间
    a=LANDFALL(:,i+5305);
    a(isnan(a))=[];
    data(i,1)=length(find(a==0));
end
data=data*3;
for i=1:3714 %多出来两个检验应该是3712 最小压力
    a=WMO_PRES(:,i+5305);
    data(i,1)=min(a);
end

%% ____计算每场台风平均移动速度 总距离km/总时间（-3）h
data=zeros(8171,3);%第一列总 第二列陆地 第三列海洋
for i=5306:13476 %50-20 每场台风
    a=DIST2LAND(:,i);%landfall会少一个格子，这里有一个接近陆地的距离多远算是收到陆地影响，0嘛？
    a(isnan(a))=[];%陆地信息
    
    lat=LAT(:,i);
    lat(isnan(lat))=[];
    lon=LON(:,i);
    lon(isnan(lon))=[];
    dist=0;%每一场都要初始化台风总距离
    Landdist=0;%两个点中有登陆点的都算在陆地距离里面
    
    dountland=0;%计算陆地上段数
    for j= 1:length(lat)-1 
        dist=dist+SphereDist2([lon(j),lat(j)],[lon(j+1),lat(j+1)]);%单位KM
        if a(j)==0 || a(j+1)==0 %小岛啥的也算上了
        Landdist=Landdist+SphereDist2([lon(j),lat(j)],[lon(j+1),lat(j+1)]); %单位KM
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
    data(i,1)=mean(b(a==0));%陆地速度
    data(i,2)=mean(b(a~=0));%海洋速度
end

%% ____计算每场台风登陆时间，去除单个3小时以内的岛屿小陆地
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
%% ____提取场均开始和结束经度
data=zeros(13501,2);
for i=5306:9019
    L=LAT(:,i);
    L(isnan(L))=[];
    a=LAT(1,i);
    b=LAT(length(L),i);
    data(i-5305,1)=a;
    data(i-5305,2)=b;
end
%% ____提取台风开始结束年分和盆地，众数盆地
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
    for j=1:length(L)%所有的盆地标志提取 %SI:1 SP:2 SA:3 NI:4 WP:5 EP:6 NA:7
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

%% ____分年分盆地统计
%从excel中导入data
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
%从excel中导入basin信息
startLON=zeros(1);
startLAT=zeros(1);
startBASIN=cell(4484,1);
endLON=zeros(1);
endLAT=zeros(1);
endBASIN=cell(4484,1);

out=zeros(71,7);
for i=1:71
    CYData=endLON(year==i+1949); %每次改这个
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
    CYData=startLON(year==i+1949); %每次改这个
    CYlocation=startBASIN(year==i+1949);
    A=CYData(strcmp(CYlocation,'SP'));
    A(A<0)=A(A<0)+360;
    out(i,3)=mean(A);

end

%% ____提取台风分年分月分盆地个数
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
for i=1:71%有问题不好用
    a=startBASIN(datay==Y+1949);
    basincount(i,1)=length(find(a=='NI'));
    basincount(i,2)=length(find(a=='SI'));
    basincount(i,3)=length(find(a=='SP'));
    basincount(i,4)=length(find(a=='WP'));
    basincount(i,5)=length(find(a=='EP'));
    basincount(i,6)=length(find(a=='NA'));
    basincount(i,7)=length(find(a=='SA'));

end

%% 提取台风陆地和海洋影响面积
idir='F:\GlobalTCMask1\TCraster500km\';
raster_files = dir([files_idir,'*.txt']);
[lmheader1,lmheader2,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家0360
%Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 根据输入数据的extend进行修改
p=parpool(10);

areaout1=zeros(length(raster_files),1);
areaout2=zeros(length(raster_files),1);
parfor i=1:length(raster_files)
[~,inputheader,ras]= read_ARCascii ([idir,raster_files(i).name]);
raster = globalize_the_imcomplete_raster(inputheader,ras);%0-360   
[Garea05,~] = Gridarea(0.5);

Garea05(raster<=0)=0;%剪裁台风形状
areaout1(i)=sum(sum(Garea05));%总影响面积
Garea05(Clandmask<0)=0;%剪裁台风LAND形状
areaout2(i)=sum(sum(Garea05));%总陆地影响面积

end
delete (p)
%% ERA5提取台风陆地和海洋影响面积分雨强
idir='F:\GlobalTCMask1\single_pre\ERA5\';
raster_files = dir([idir,'*.txt']);%-180-180
[lmheader1,lmheader2,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 根据输入数据的extend进行修改
all_in_one = Basinmasks(0.5);
all_in_one=[all_in_one(:,361:720),all_in_one(:,1:360)];%-180-180 根据输入数据的extend进行修改


areaout1=zeros(length(raster_files),1);%全部
areaout2=zeros(length(raster_files),1);%大中雨强，5天30mm为分界
areaout1L=zeros(length(raster_files),1);%全部
areaout2L=zeros(length(raster_files),1);%大中雨强，5天30mm为分界
areaout1B=zeros(length(raster_files),6);%全部
areaout2B=zeros(length(raster_files),6);%大中雨强，5天30mm为分界
areaout1BL=zeros(length(raster_files),6);%全部
areaout2BL=zeros(length(raster_files),6);%大中雨强，5天30mm为分界
for i=1:8171
[~,inputheader,TCpre]= read_ARCascii ([idir,raster_files(i).name]);
%------全部降水
[Garea05,~] = Gridarea(0.5);
Garea05(TCpre<=0.1)=0;%剪裁台风形状
areaout1(i)=sum(sum(Garea05));%总面积

for j=1:6
areaout1B(i,j)=sum(sum(Garea05(all_in_one==j)));%盆地面积
end

Garea05(Clandmask<0)=0;%剪裁台风LAND形状
areaout1L(i)=sum(sum(Garea05));%总陆地影响面

for j=1:6
areaout1BL(i,j)=sum(sum(Garea05(all_in_one==j)));%盆地陆地面积
end

%-------------显著降水
[Garea05,~] = Gridarea(0.5);
Garea05(TCpre<=30)=0;%剪裁台风形状
areaout2(i)=sum(sum(Garea05));

for j=1:6
areaout2B(i,j)=sum(sum(Garea05(all_in_one==j)));%盆地面积
end

Garea05(Clandmask<0)=0;%剪裁台风LAND形状
areaout2L(i)=sum(sum(Garea05));%总陆地影响面积

for j=1:6
areaout2BL(i,j)=sum(sum(Garea05(all_in_one==j)));%盆地陆地面积
end

end

%% 箱型图按年提取
input=zeros(1);%excel里面拷贝
year=zeros(1);%开始年
output=ones(200,71)*-999;%没有那一年超过200台风，后面删除-999方便
for y=1950:2020
  middle=input(year==y); 
  middle(middle==0)=[];  
  output(1:length(middle),y-1949)=middle;     
end

%% 提取每年最大的一些
input=zeros(1);%excel里面拷贝
year=zeros(1);%开始年
c=zeros(71,1);
for i=1950:2020
    a=input(year==i);
    b=sort(a,'descend');
    d=round(length(b)*0.5);
    c(i-1949,1)=sum(sum( b( 60:80 ,1 )));
end
plot(c);
%% ____计算每场台风最大持续风速
data=zeros(8200,2);
for i=5306:13500 
    a=WMO_WIND(:,i);
    b=str2double(ISO_TIME(1:4,1,i));
    data(i-5305,1)=max(a);
    data(i-5305,2)=b;
end






