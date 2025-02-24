%读取一年的所有TC经过ERA5，并取平均值，除网格台风经过次数,每个点每年平均每场经过台风的降水强度mm/h
clear; clc; close all;
header=['ncols           720';
        'nrows           360';
        'xllcorner      -180';
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value  -9999'];
    

file_idir='F:\GlobalTCMask1\single_pre\ERA5_Single_Precip3d\';
odir='F:\GlobalTCMask1\AVEPricip\ERA5_Grid_PreTC_AVEPrecip_Ori_Ummd\';
ERAfilename='ERA5_Ori_SingleTC_3dPrecip_Ummh_SID';

countF = dir(['F:\GlobalTCMask1\TCraster500km\','*.txt']);%所有数据
dataindir='C:\Users\Dell\Desktop\Research_new\COMMENT.xlsx';
ALLy=xlsread(dataindir,1,'A2:A8172'); %开始年1950因为下面的文件也是,因为有一个7191号年份不匹配，excel改为1967年raster的数据可能老了和后面新数据的有个别时间不匹配

for y=2018:2020
    files = dir([file_idir,ERAfilename,num2str(y),'*.txt']);%每年数据读一下
    countfiles=countF(ALLy==y);
    
    data=zeros(360,720);
    count=zeros(360,720);%计数
    for i=1:length(files)
        [~,inputheader,ras]= read_ARCascii ([file_idir,files(i).name]);
        ras(ras<0)=0;
        data=data+ras;
        
        [~,inputheader2,Cras]= read_ARCascii (['F:\GlobalTCMask1\TCraster500km\',countfiles(i).name]);
        Cras = globalize_the_imcomplete_raster(inputheader2,Cras);%注意这个输出0-360 
        Cras=[Cras(:,361:720),Cras(:,1:360)];%-180-180
        count=count+Cras;   
    end
    count(count==0)=10000;
    avesst=data./count;
    avesst(avesst<0.00001)=-9999;
    OUTPUT4(odir,['ERA5_Grid_PreTC_AVEPrecip_Ori_Ummd',num2str(y)],header,avesst);%可能太小了必须提升精度，忽略了台风强度每场都不一样
end
