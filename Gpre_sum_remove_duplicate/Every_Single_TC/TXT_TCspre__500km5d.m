%读取逐年全球台风500km影响区降水总量 前后五天相加的版本  FOR TXT
%update20201020 注意这个没有去除降水数据中的nodata value因为mswep没有坏点  跨年台风的降水会算进前一年  
clear; clc; close all;
files_idir ='F:\GlobalTCMask1\single_pre\MSWEP\';
raster_files = dir([files_idir,'*.txt']);
[Garea05,~] = Gridarea(0.5);
S_TC_Pre=zeros(length(raster_files),1);%每场台风的降水
parfor NO=1:length(raster_files)
    [~,inputheader,pre]= read_ARCascii ([files_idir,raster_files(NO).name]);
    P=[pre(:,361:720),pre(:,1:360)];
    P(P<0)=0;
    S_TC_Pre(NO)=10^(-6)*sum(sum(P.*Garea05));
end