%计算各种多年统计结果0360
%台风时间以第一次发现为准，跨年的依然算第一年的台风，可能导致一些误差，假设mswep没有坏点
clear; clc; close all;
header1=['ncols           720';
        'nrows           360';
        'xllcorner      -180'; %注意-180-180
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value  -9999'];
    
TC_idir='E:\DATA\IBTrACS\IBTrACS.since1980.v04r00.nc';
GTCP_idir ='F:\GlobalTCMask1\NODEduplicate_sum\MERRA2Ori\';
GTCPfiles = dir([GTCP_idir,'*.txt']);      %GTCP降水文件路径
STCP_idir ='F:\GlobalTCMask1\single_pre\';
STCPfiles = dir([STCP_idir,'*.txt']);      %STCP降水文件路径
STCras_idir ='F:\GlobalTCMask1\TCraster500km\';
STCras_files = dir([STCras_idir,'*.txt']);      %STCP降水文件路径
odir='';

cs=0.5;
[area,~] = Gridarea(cs);                    %陆地面积 单位：平方千米
[~,~,~,~,~,~,~,Basinmask] = Basinmasks(cs); %大洋块分区SI:1 SP:2 SA:3 NI:4 WP:5 EP&NA:6 0360
Basinmask=[Basinmask(:,(180/cs)+1:360/cs),Basinmask(:,1:180/cs)];%-180-180 根据输入数据的extend进行修改
[lmheader1,lmheader2,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家
Clandmask=[Clandmask(:,(180/cs)+1:360/cs),Clandmask(:,1:180/cs)];%-180-180 根据输入数据的extend进行修改

%% -----------------------basinTC降水量序列---------------------------------
%---------------结果数组初始化------------
Basinpre=zeros(length(GTCPfiles),6);        %SI:1 SP:2 SA:3 NI:4 WP:5 EP&NA:6 初始化每年每basin降水组
Basinlandpre=Basinpre;

for i1 = 1:length(GTCPfiles)
   [ph1,ph2,pre1]= read_ARCascii([GTCP_idir,GTCPfiles(i1).name]);%！！！这里要求pre除了降水以外的点值都为零，不要有nodata value

   Bpre=zeros(180/cs,360/cs,6);             %每年重新赋值 1=SI 2=SP 3=SA 4=NI 5=WP 6=EP=NA 
   
   %basin内总降水量--Basinpre
   for j1 = 1:6                             %每盆地内的所有降水放入Bpre
   cur_p=zeros(180/cs,360/cs);              %初始化一个二维空白
   cur_p(Basinmask==j1)=pre1(Basinmask==j1); %每个盆地放进去
   Bpre(:,:,j1)=cur_p;                      %存储到Bpre方便后面陆地剪裁
   %cur_p(cur_p<0)=0;
   Basinpre(i1,j1)=10^(-6)*sum(sum(cur_p.*area));          %当前盆地每年TC总降水统计值 单位：十亿m3
   end
   
   %basin内陆地上的总降水量--Basinlandpre
   BLpre=Bpre;                              %把这一年每块的降水要过来继续
   for j2 = 1:6                             %每盆地内的所有陆地上降水放入Basinlandpre
   cur_lp=BLpre(:,:,j2);                    %当前盆地内总降水分布
   cur_lp(Clandmask==-9999)=0;              %减去所有非陆地部分 注意这里可能去掉了一些很小很小的岛屿
   BLpre(:,:,j2)=cur_lp;                    %直接覆盖掉原来没剪裁的 BLpre
   %cur_lp(cur_p<0)=0;
   Basinlandpre(i1,j2)=10^(-6)*sum(sum(cur_lp.*area));     %当前盆地每年陆地上TC总降水统计值 单位：十亿m3  ！！！这样计算损失了一些很小很小的岛，而且海陆交界带也不清楚会产生多大误差
   end  
    
end

%% -----------------------每年单台风降水性质---------------------------------
%---------------结果数组初始化------------                    %SI:1 SP:2 SA:3 NI:4 WP:5 EP&NA:6 1980-2016 年的数据
BSTCYave_p=zeros(length(STCPfiles),6);                       %每个盆地每年每场台风的平均降水  
BLSTCYave_p=Basinpre;                                        %每个盆地每年每场登陆台风为陆地带来的的平均降水BASIN LAND SINGLE TC YEARLY AVERAGE
BSTCave_p=0;                                                 %每个盆地多年每场台风的平均降水  
BLSTCave_p=0;                                                %每个盆地多年每场登陆台风为陆地带来的的平均降水
GLSTCave_Lp=0;                                               %全球每场登陆台风平均陆地上降水
GLSTCave_Lprate=0;                                           %全球每场登陆台风平均陆地上降水占总降水之比
 
%% ---全球每场台风平均降水GSTCave_p
GSTC_p=zeros(length(STCPfiles),1);                           %每场台风降水量
GSTCave_p=0;                                                 %全球每场台风平均降水
for i2 = 1:length(STCPfiles)
   [ph3,ph4,pre2]= read_ARCascii([STCP_idir,STCPfiles(i2).name]); %！！！这里nodata value -999 
   GSTC_p(i2)=10^(-6)*sum(pre2(pre2>0).*area(pre2>0));       %注意去除错误值
   GSTCave_p=GSTCave_p+GSTC_p(i2);                           %单位   十亿m3            
end
GSTCave_p=GSTCave_p/length(STCPfiles);
%% 
%% ---每个网格的台风经过次数1980-1999，2000-2019
gridcount8099 = zeros(360,720);
gridcount0019 = zeros(360,720);
for rasNO = 1:2284 %80-99
    [~,inputheader,ras] = read_ARCascii([STCras_idir,STCras_files(rasNO).name]);
    raster = globalize_the_imcomplete_raster(inputheader,ras);
    raster=[raster(:,(180/0.5)+1:360/0.5),raster(:,1:180/0.5)]; %-180-180
    raster(raster <= 0) = 0;
    raster(raster == 500) = 1;
    gridcount8099 = gridcount8099 + raster;
end
for rasNO = 2285:4338 %00-19
    [~,inputheader,ras] = read_ARCascii([STCras_idir,STCras_files(rasNO).name]);
    raster = globalize_the_imcomplete_raster(inputheader,ras);
    raster=[raster(:,(180/0.5)+1:360/0.5),raster(:,1:180/0.5)]; %-180-180
    raster(raster <= 0) = 0;
    raster(raster == 500) = 1;
    gridcount0019 = gridcount0019 + raster;
end
ave_gridcount1980_2019 = (gridcount0019+gridcount8099)/40;
ave_gridcount1980_1999 = gridcount8099/20;
ave_gridcount2000_2019 = gridcount0019/20;
change_gridcount = ave_gridcount2000_2019 - ave_gridcount1980_1999;

gridcountodir = 'C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Statistics\Gridcount\';
OUTPUT(gridcountodir, 'AVE_TCgridcount1980_1999', header1, ave_gridcount1980_1999)
OUTPUT(gridcountodir, 'AVE_TCgridcount2000_2019', header1, ave_gridcount2000_2019)
OUTPUT(gridcountodir, 'AVE_TCgridcount1980_2019', header1, ave_gridcount1980_2019)
OUTPUT(gridcountodir, 'AVE_CH_TCgridcount0019_8099', header1, change_gridcount)
[~,~,CH_PRE] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Statistics\Spatialchanges\spa_changesMSWEP37Y_IMERG3Y.txt');
