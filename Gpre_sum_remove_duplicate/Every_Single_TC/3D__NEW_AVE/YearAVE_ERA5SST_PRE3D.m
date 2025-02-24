%��ȡһ�������TC����ERA5sst����ȡƽ��ֵ��������̨�羭������,̨����ÿ����ƽ�������ĺ���
clear; clc; close all;
header=['ncols           720';
        'nrows           360';
        'xllcorner      -180';
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value  -9999'];
    

file_idir='F:\GlobalTCMask1\single_SST\ERA5_Pre3d_sst\';
odir='F:\GlobalTCMask1\AVESST\ERA5_GSUMTC_AVESST_Pre3d_Ori\';
ERAfilename='ERA5_PRE3d_Ori_SingleTC_SST_SID';
for y=1966:2020
    files = dir([file_idir,ERAfilename,num2str(y),'*.txt']);%ÿ�����ݶ�һ��
    data=zeros(360,720);
    count=zeros(360,720);%����
    for i=1:length(files)
        [~,inputheader,ras]= read_ARCascii ([file_idir,files(i).name]);
        ras(ras<0)=0;
        data=data+ras;
        ras(ras>0)=1;%����˵����Ӧ����500km��Ĥ�����Ǻ���û�е�������������಻���ǣ���ˮҪ����
        count=count+ras;   
    end
    count(count==0)=10000;
    avesst=data./count;
    avesst(avesst<0.5)=-9999;
    OUTPUT(odir,['ERA5_GSUMTC_AVESST_Pre3d_Ori',num2str(y)],header,avesst);
end
