clear; clc; close all;
Preindir='F:\GlobalTCMask1\single_pre\ERA5_Single_Precip3d_mmd\'; %��
Prefile=dir([Preindir,'*.txt']);
MaxPre=zeros(length(Prefile),1);

parfor i=1:length(Prefile)
[~,~,prerate] = read_ARCascii([Preindir,Prefile(i).name]);%��λmm/h
MaxPre(i,1)=max(max(prerate));
end
