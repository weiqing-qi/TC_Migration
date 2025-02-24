%��½�������밶��ÿ�ٹ��ﻷ�ھ���̨�罵ˮ���˿�
%ÿ������ÿ�����Ľ�ˮǿ��
%ɸѡ����50mm/d�Ľ�ˮ���µ��˿�
%�û����еõ�ÿ��ÿ���˿���
clear; clc; close all;
file_idir='C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%ע���ǲ��ǰ�˳�򴢴��
%½�ع�������
[~,~,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 �����������ݵ�extend�����޸�
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
%% ����һ�������ˮǿ�ȷֲ����˿ڷֲ�
indir='F:\GlobalTCMask1\single_pre\PERSIANN_Single_Precip3d_mmd\'; %��
MaxPre=zeros(360,720,38);
parfor y=1983:2020
    files=dir([indir,'PERSIANNCDR_Ori_SingleTC_3dPrecip_Ummd_SID',num2str(y),'*.txt']);%MM/Dע�ⵥλ %��
    for i=1:length(files)
    [~,~,pre2] = read_ARCascii([indir,files(i).name]);%��λmm/D
    pre1=MaxPre(:,:,y-1982);
    pre1(pre2>pre1)=pre2(pre2>pre1);
    MaxPre(:,:,y-1982)=pre1;
    end
end
%% ------------------�����˿�
Popindir='E:\DATA\WORLD_Population\Landscan\Landscan05_txt\';
PopDistri_TCP=zeros(360,720,21);%ÿ�꾭��̨�罵ˮ����
PopAll_TCP=zeros(21,1);%ÿ��������
PopRing_TCP=zeros(21,8);   %ÿ��800-700...100-0km to coast
Pop_TCP_Basin=zeros(21,6);   %
PopRing_TCP_Basin=zeros(21,8,6);   %ÿ������800-700...100-0km to coast
A=sort(MaxPre(MaxPre>0.1));
% PreBar=A(round(length(A)*0.75)); %����75��λ  ��00000 0000 00000 0000000000 0000 000
PreBar=0.1;
all_in_one = Basinmasks(0.5); all_in_one=[all_in_one(:,361:720),all_in_one(:,1:360)];
parfor y2=2000:2020
POPTCP=zeros(360,720);

[~,~,Pop] = read_ARCascii([Popindir,'Resample05-landscan-global-',num2str(y2),'.txt']);%ע��������ݼ�������ֵ������
MP=MaxPre(:,:,y2-1982);
POPTCP(MP>PreBar)=Pop(MP>PreBar); %ɸѡ���ڱ�׼��TC��ˮ�����˿ڷֲ�
PopDistri_TCP(:,:,y2-1999)=POPTCP;
PopAll_TCP(y2-1999)=sum(sum(POPTCP));
for basin1=1:6
Pop_TCP_Basin(y2-1999,basin1)=sum(POPTCP(all_in_one==basin1));
end
%------ÿ�����������
    for d=1:8    %800-700...100-0km to coast 
        ZoneLR=CRing_land(:,:,9-d);
        PopRing_TCP(y2-1999,d)=sum(POPTCP(ZoneLR==10000));   
        %-------�ָ����all_in_one
        LandRingPop=zeros(360,720);
        LandRingPop(ZoneLR==10000)=POPTCP(ZoneLR==10000);
        for basin2=1:6
            PopRing_TCP_Basin(y2-1999,d,basin2)=sum(LandRingPop(all_in_one==basin2)); 
        end
    end
end

PR_Regress=zeros(8,4);
for i=1:8
[b,bint,r,rint,stats]=regress(PopRing_TCP(:,i),[ones(2020-2000+1,1),(2000:2020)'],0.05);
    PR_Regress(i,1)=b(2);
    PR_Regress(i,2)=bint(2,1);%ע��
    PR_Regress(i,3)=bint(2,2);
    PR_Regress(i,4)=stats(3); 
end
