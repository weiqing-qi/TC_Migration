%��½�������밶��ÿ�ٹ��ﻷ��������ˮ��Χ
clear; clc; close all;
tic;
file_idir='C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%ע���ǲ��ǰ�˳�򴢴��
%½�ع�������
[~,~,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 �����������ݵ�extend�����޸�
[GA05,~] = Gridarea(0.5);
GA05Land=zeros(360,720);
GA05Land(Clandmask>0)=GA05(Clandmask>0);%ȫ��½���ϵ����������Ϊ��
%% ������뻷�ͺ����߽��� 
CRing_sea=zeros(360,720,20);
CRing_land=zeros(360,720,20);
Cbuffer_sea=zeros(360,720,20);
Cbuffer_land=zeros(360,720,20);
CR_BS=zeros(360,720,20);

    [lmheader1,lmheader2,CR_BS1] = read_ARCascii([file_idir,filename,'100km.txt']);
    CR_BS1(CR_BS1~=1)=0;CR_BS1(CR_BS1==1)=10000;
    CR_BS(:,:,1)=CR_BS1;
for i=200:100:2000
    [~,~,CB_Near] = read_ARCascii([file_idir,filename,num2str(i-100),'km.txt']);
    [~,~,CB_Far] = read_ARCascii([file_idir,filename,num2str(i),'km.txt']);
    CR_BS(:,:,i*0.01)=CB_Far-CB_Near;%��ֵӦ����1-��-9999��=10000
end
    CR_BS(1:60,:,:)=0; CR_BS(301:360,:,:)=0;%ȥ��60ns��� 
    

for j=1:20%100:100:2000
    Csea=zeros(360,720);
    Cland=zeros(360,720);
    
    BS=CR_BS(:,:,j);
    Csea(Clandmask<0)=BS(Clandmask<0);
    Cland(Clandmask>0)=BS(Clandmask>0);
    CRing_sea(:,:,j)=Csea;
    CRing_land(:,:,j)=Cland;
    Cbuffer_sea(:,:,j)=sum(CRing_sea,3);
    Cbuffer_land(:,:,j)=sum(CRing_land,3);
end
%% ����һ�������ˮǿ�ȷֲ�
indir='F:\GlobalTCMask1\single_pre\ERA5_Single_Precip3d_mmd\'; %��
MaxPre=zeros(360,720,41);
parfor y=1980:2020
    files=dir([indir,'ERA5_Ori_SingleTC_3dPrecip_Ummh_SID',num2str(y),'*.txt']);%MM/Dע�ⵥλ %��
    for i=1:length(files)
    [~,~,pre2] = read_ARCascii([indir,files(i).name]);%��λmm/D
    pre1=MaxPre(:,:,y-1979);
    pre1(pre2>pre1)=pre2(pre2>pre1);
    MaxPre(:,:,y-1979)=pre1;
    end
end
%% ����仯
Area=zeros(41,1);%ȫ�������
Brea=zeros(41,7);%��������
R_Area=zeros(41,8);%RINGȫ�������
R_Brea=zeros(41,8,7);%RING��������
% A=sort(MaxPre(MaxPre>0.1));
% PreBar=A(round(length(A)*0.75)); %����75��λ
PreBar=30;
BasMas = Basinmasks_EPNA(0.5,1);
for y2=1980:2020
MP=MaxPre(:,:,y2-1979);
Area(y2-1979,1)=sum(GA05Land(MP>=PreBar));%ȫ�������

    for bas=1:7
        GA05LB=zeros(360,720);%�������½���������Ϊ��
        GA05LB(BasMas==bas)=GA05Land(BasMas==bas);
        Brea(y2-1979,bas)=sum(GA05LB(MP>=PreBar));%��������
    %------ÿ����
        for d=1:8    %800-700...100-0km to coast 
            ZoneLR=CRing_land(:,:,9-d);
            GA05LBR=zeros(360,720);
            GA05LBR(ZoneLR==10000)=GA05LB(ZoneLR==10000);
            R_Brea(y2-1979,d,bas)=sum(GA05LBR(MP>=PreBar)); 
            
            if bas==1
                GA05LR=zeros(360,720);
                GA05LR(ZoneLR==10000)=GA05(ZoneLR==10000);%�������������������Ϊ��
                R_Area(y2-1979,d)=sum(GA05LR(MP>=PreBar));%RINGȫ�������
            end
        end
    end
end

GAR_Regress=zeros(8,4);
BAR_Regress=zeros(8,4,7);
for i=1:8
[b,bint,r,rint,stats]=regress(R_Area(:,i),[ones(2020-1980+1,1),(1980:2020)'],0.05);
GAR_Regress(i,1)=b(2);
GAR_Regress(i,2)=bint(2,1);%ע��
GAR_Regress(i,3)=bint(2,2);
GAR_Regress(i,4)=stats(3); 
    for bas=1:7
    [b,bint,r,rint,stats]=regress(R_Brea(:,i,bas),[ones(2020-1980+1,1),(1980:2020)'],0.05);
    BAR_Regress(i,1,bas)=b(2);
    BAR_Regress(i,2,bas)=bint(2,1);%ע��
    BAR_Regress(i,3,bas)=bint(2,2);
    BAR_Regress(i,4,bas)=stats(3); 
    end
end
out1=[Area,R_Area,Brea];
out2=[];
for bas=1:7
    for i=1:3
    out2=[out2,BAR_Regress(:,i,bas)];
    end
end
out2(:,7:9)=[];
toc;