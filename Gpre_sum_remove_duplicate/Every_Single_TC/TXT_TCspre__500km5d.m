%��ȡ����ȫ��̨��500kmӰ������ˮ���� ǰ��������ӵİ汾  FOR TXT
%update20201020 ע�����û��ȥ����ˮ�����е�nodata value��Ϊmswepû�л���  ����̨��Ľ�ˮ�����ǰһ��  
clear; clc; close all;
files_idir ='F:\GlobalTCMask1\single_pre\MSWEP\';
raster_files = dir([files_idir,'*.txt']);
[Garea05,~] = Gridarea(0.5);
S_TC_Pre=zeros(length(raster_files),1);%ÿ��̨��Ľ�ˮ
parfor NO=1:length(raster_files)
    [~,inputheader,pre]= read_ARCascii ([files_idir,raster_files(NO).name]);
    P=[pre(:,361:720),pre(:,1:360)];
    P(P<0)=0;
    S_TC_Pre(NO)=10^(-6)*sum(sum(P.*Garea05));
end