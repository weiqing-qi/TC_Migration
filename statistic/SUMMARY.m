%������ֶ���ͳ�ƽ��0360
%̨��ʱ���Ե�һ�η���Ϊ׼���������Ȼ���һ���̨�磬���ܵ���һЩ������mswepû�л���
clear; clc; close all;
header1=['ncols           720';
        'nrows           360';
        'xllcorner      -180'; %ע��-180-180
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value  -9999'];
    
TC_idir='E:\DATA\IBTrACS\IBTrACS.since1980.v04r00.nc';
GTCP_idir ='F:\GlobalTCMask1\NODEduplicate_sum\MERRA2Ori\';
GTCPfiles = dir([GTCP_idir,'*.txt']);      %GTCP��ˮ�ļ�·��
STCP_idir ='F:\GlobalTCMask1\single_pre\';
STCPfiles = dir([STCP_idir,'*.txt']);      %STCP��ˮ�ļ�·��
STCras_idir ='F:\GlobalTCMask1\TCraster500km\';
STCras_files = dir([STCras_idir,'*.txt']);      %STCP��ˮ�ļ�·��
odir='';

cs=0.5;
[area,~] = Gridarea(cs);                    %½����� ��λ��ƽ��ǧ��
[~,~,~,~,~,~,~,Basinmask] = Basinmasks(cs); %��������SI:1 SP:2 SA:3 NI:4 WP:5 EP&NA:6 0360
Basinmask=[Basinmask(:,(180/cs)+1:360/cs),Basinmask(:,1:180/cs)];%-180-180 �����������ݵ�extend�����޸�
[lmheader1,lmheader2,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���
Clandmask=[Clandmask(:,(180/cs)+1:360/cs),Clandmask(:,1:180/cs)];%-180-180 �����������ݵ�extend�����޸�

%% -----------------------basinTC��ˮ������---------------------------------
%---------------��������ʼ��------------
Basinpre=zeros(length(GTCPfiles),6);        %SI:1 SP:2 SA:3 NI:4 WP:5 EP&NA:6 ��ʼ��ÿ��ÿbasin��ˮ��
Basinlandpre=Basinpre;

for i1 = 1:length(GTCPfiles)
   [ph1,ph2,pre1]= read_ARCascii([GTCP_idir,GTCPfiles(i1).name]);%����������Ҫ��pre���˽�ˮ����ĵ�ֵ��Ϊ�㣬��Ҫ��nodata value

   Bpre=zeros(180/cs,360/cs,6);             %ÿ�����¸�ֵ 1=SI 2=SP 3=SA 4=NI 5=WP 6=EP=NA 
   
   %basin���ܽ�ˮ��--Basinpre
   for j1 = 1:6                             %ÿ����ڵ����н�ˮ����Bpre
   cur_p=zeros(180/cs,360/cs);              %��ʼ��һ����ά�հ�
   cur_p(Basinmask==j1)=pre1(Basinmask==j1); %ÿ����طŽ�ȥ
   Bpre(:,:,j1)=cur_p;                      %�洢��Bpre�������½�ؼ���
   %cur_p(cur_p<0)=0;
   Basinpre(i1,j1)=10^(-6)*sum(sum(cur_p.*area));          %��ǰ���ÿ��TC�ܽ�ˮͳ��ֵ ��λ��ʮ��m3
   end
   
   %basin��½���ϵ��ܽ�ˮ��--Basinlandpre
   BLpre=Bpre;                              %����һ��ÿ��Ľ�ˮҪ��������
   for j2 = 1:6                             %ÿ����ڵ�����½���Ͻ�ˮ����Basinlandpre
   cur_lp=BLpre(:,:,j2);                    %��ǰ������ܽ�ˮ�ֲ�
   cur_lp(Clandmask==-9999)=0;              %��ȥ���з�½�ز��� ע���������ȥ����һЩ��С��С�ĵ���
   BLpre(:,:,j2)=cur_lp;                    %ֱ�Ӹ��ǵ�ԭ��û���õ� BLpre
   %cur_lp(cur_p<0)=0;
   Basinlandpre(i1,j2)=10^(-6)*sum(sum(cur_lp.*area));     %��ǰ���ÿ��½����TC�ܽ�ˮͳ��ֵ ��λ��ʮ��m3  ����������������ʧ��һЩ��С��С�ĵ������Һ�½�����Ҳ����������������
   end  
    
end

%% -----------------------ÿ�굥̨�罵ˮ����---------------------------------
%---------------��������ʼ��------------                    %SI:1 SP:2 SA:3 NI:4 WP:5 EP&NA:6 1980-2016 �������
BSTCYave_p=zeros(length(STCPfiles),6);                       %ÿ�����ÿ��ÿ��̨���ƽ����ˮ  
BLSTCYave_p=Basinpre;                                        %ÿ�����ÿ��ÿ����½̨��Ϊ½�ش����ĵ�ƽ����ˮBASIN LAND SINGLE TC YEARLY AVERAGE
BSTCave_p=0;                                                 %ÿ����ض���ÿ��̨���ƽ����ˮ  
BLSTCave_p=0;                                                %ÿ����ض���ÿ����½̨��Ϊ½�ش����ĵ�ƽ����ˮ
GLSTCave_Lp=0;                                               %ȫ��ÿ����½̨��ƽ��½���Ͻ�ˮ
GLSTCave_Lprate=0;                                           %ȫ��ÿ����½̨��ƽ��½���Ͻ�ˮռ�ܽ�ˮ֮��
 
%% ---ȫ��ÿ��̨��ƽ����ˮGSTCave_p
GSTC_p=zeros(length(STCPfiles),1);                           %ÿ��̨�罵ˮ��
GSTCave_p=0;                                                 %ȫ��ÿ��̨��ƽ����ˮ
for i2 = 1:length(STCPfiles)
   [ph3,ph4,pre2]= read_ARCascii([STCP_idir,STCPfiles(i2).name]); %����������nodata value -999 
   GSTC_p(i2)=10^(-6)*sum(pre2(pre2>0).*area(pre2>0));       %ע��ȥ������ֵ
   GSTCave_p=GSTCave_p+GSTC_p(i2);                           %��λ   ʮ��m3            
end
GSTCave_p=GSTCave_p/length(STCPfiles);
%% 
%% ---ÿ�������̨�羭������1980-1999��2000-2019
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
