%��ȡ��һ������ˮ�ռ�ֲ�
clc; clear;close all;
counmaskindir='C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt';%43
[H1,H2,counmask]=read_ARCascii(counmaskindir);%0-360
prefilepath='F:\GlobalTCMask1\NODEduplicate_sum\originalextend\';%������-180-180
prefiles=dir([prefilepath,'*.txt']);
odir='C:\Users\Dell\Desktop\CHINACX\';

header=['ncols           720';
        'nrows           360';
        'xllcorner      -180';
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value   -999'];

cs=0.5; counmaskOri=[counmask(:,(180/cs)+1:360/cs),counmask(:,1:180/cs)];%��Ҫת��Ϊһ�µ�-180-180
[area,~] = Gridarea(cs);
[r,c]=size(counmask); ChinaLandTCpre=zeros(r,c);%��ʼ��
for i=1:length(prefiles)
[PH1,PH2,GTCP]=read_ARCascii([prefilepath,prefiles(i).name]);%��ˮ-180-180
ChinaLandTCpre(counmaskOri==43)=GTCP(counmaskOri==43);%-180-180
ChinaLandTCpre=(ChinaLandTCpre.*area)*0.00001;%��m?
OUTPUT(odir,['ChinaTCLP',prefiles(i).name(17:20)],header,ChinaLandTCpre)
end
