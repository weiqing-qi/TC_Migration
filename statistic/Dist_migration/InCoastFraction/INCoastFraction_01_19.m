%������뺣����ͬ���뻷�����ˮռ�������ı����仯
%����ÿ�������ڵĽ�ˮ��ǿ�ȱ仯����Ϊ���ں��ܾ�����
%ע��������60NS
clear; clc; close all;
file_idir='D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%ע���ǲ��ǰ�˳�򴢴��
%Filter
Ylim=[2001 2019];
TC_idir='E:\DATA\5.TCdata\IBTrACS\NEW\IBTrACS.ALL.v04r00.nc';
TC_idir='D:\DATA\Segmentation map tutorial\IBTrACS.since1980.v04r00.240321.nc';
dist=xlsread('D:\Desktop2\Global_cyclone_project\DATA2.xlsx',2,'I2401:I4340');
dist=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted.xlsx',1,'AR330:AR2269');
dist=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted.xlsx',1,'AR3:AR1943');
[~,~,~,ISO_TIME,USA_WIND,~,BASIN,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);
[KEEP_NO,KEEP_Y,KEEP_B,KEEP_i] = FilterNO(Ylim,ISO_TIME,USA_WIND,BASIN,dist);
%��ˮ����
% patchdir='F:\GlobalTCMask1\single_pre\IMERG_Single_Precip3d_mmd\'; %��
% patchfile=dir([patchdir,'*.txt']);
% patchfile=patchfile(1:1940);%��2001��ʼ!!!!!!!!!!!!IMERG

% Preindir='D:\DATA\TC_spatial_data\TMPA_PRE025\'; %��
% Preindir='D:\DATA\TC_spatial_data\IMERG_FV7B_PRE01\'; %��
Preindir='E:\DATA\TCMigrationPJ\TMPA_Single_Precip3d_mmd\'; %��
Prefile=dir([Preindir,'*.txt']);
Prefile=Prefile(1:1941);
Prefile=Prefile(328:end);
Prefile=Prefile(KEEP_i);
% Prefile=Prefile(328:end);%��2001��ʼ!!!!!!!!!!!!TMPA

% Prefile=Prefile(2399:4338);%��2001-2019��ʼ!!!!!!!!!!!!1980-2020
%½�ع�������
[~,~,Clandmask] = read_ARCascii('D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���0360
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
    CR_BS(:,:,i*0.01)=CB_Far-CB_Near;%��ֵӦ����1-��-9999��10000
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
% %% ������뺣����ͬ���뻷�����ˮռ�������ı����仯
% % data=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA2.xlsx',2,'C3:D4461');
% % LMI=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA2.xlsx',1,'D3715:D8173');
% year=xlsread('G:\����\Global_cyclone_project\DATA2.xlsx',2,'C2401:C4461');%��
% MaxPreDist=xlsread('G:\����\Global_cyclone_project\DATA2.xlsx',2,'G2401:G4461');%��
% % year=year(LMI>34);
% % MaxPreDist=MaxPreDist(LMI>34);
% %% ��������͸���
% trend_mpfraction=zeros(28,4);%-600km-2000km,CIUp,CIDown��P
% MKtest=zeros(28,4);
% Z=zeros(28,1);
% for d=-800:100:1900
%     FractPY_R=zeros(20,1);%����
%     CountPY_R=zeros(20,1);%����
%     FractPY_B=zeros(20,1);%����
%     for y=2001:2020
%         PerY_MP=MaxPreDist(year==y);%�õ�ÿ������
%         FractPY_R(y-2000)=length(find(PerY_MP>d & PerY_MP<=d+100))/length(PerY_MP);
%         CountPY_R(y-2000)=length(find(PerY_MP>d & PerY_MP<=d+100));
%         if d<0
%         FractPY_B(y-2000)=length(find(PerY_MP>d & PerY_MP<=0))/length(PerY_MP);
%         else
%         FractPY_B(y-2000)=length(find(PerY_MP>0 & PerY_MP<d+100))/length(PerY_MP);    
%         end
%     end
%  
%     [b,bint,r,rint,stats]=regress(FractPY_B,[ones(2020-2001+1,1),(2001:2020)'],0.05);
%     
%     trend_mpfraction((d+900)*0.01,1)=b(2);
%     trend_mpfraction((d+900)*0.01,2)=bint(2,1);
%     trend_mpfraction((d+900)*0.01,3)=bint(2,2);
%     trend_mpfraction((d+900)*0.01,4)=stats(3); 
%     
%     [Zs, p_value, UFk, UBk2, beta, beta_CI, corr_val]= MKtrend((2001:2020)',FractPY_B);
%     MKtest((d+900)*0.01,1)=p_value;
%     MKtest((d+900)*0.01,2)=beta;
%     MKtest((d+900)*0.01,3)=beta_CI(1);
%     MKtest((d+900)*0.01,4)=beta_CI(2);
%     
%     Z((d+900)*0.01,1)=length(find(FractPY_B==0));
% 
% end

%% ƽ����ˮ�ʽ���
R_PR=zeros(length(Prefile),28);%-800��2000km
B_PR=zeros(length(Prefile),28);
ATCP=zeros(length(Prefile),1);
%ÿ��TCÿһRorB�е�ƽ��ǿ��
parfor i=1:length(Prefile)
    disp(i)
[~,H,prerate] = read_ARCascii([Preindir,Prefile(i).name]);%��λmm/h
% prerate = FillView2Global(H,prerate);
prerate = GRID_resample(prerate, 0.1, 0.5);
prerate(isnan(prerate))=0; prerate(prerate<0)=0;
%-----------------------����50-60NS--------------------------
% [~,~,patch] = read_ARCascii([patchdir,patchfile(i).name]);%��λmm/h
% patch(isnan(patch))=0; patch(patch<0)=0;
% prerate(61:80,:)=patch(61:80,:);%��IMERG������Թص�
% prerate(281:300,:)=patch(281:300,:);%��IMERG����
%-----------------------------------------------------
ATCP(i)=mean(prerate(prerate>0.1));
    for d=1:28%100:2000km
        if d<=8                               %Land
        ZoneLR=CRing_land(:,:,9-d);          %ע�⽵ˮת��
        preLR=prerate(ZoneLR==10000);         %������������н�ˮ
        R_PR(i,d)=mean(preLR(preLR>0.1));     %0.1���ż���ƽ����ǿ
        
        ZoneLB=Cbuffer_land(:,:,9-d);
        preLB=prerate(ZoneLB==10000);
        B_PR(i,d)=mean(preLB(preLB>0.1));
        else                                  %Sea
        ZoneSR=CRing_sea(:,:,d-8);
        preSR=prerate(ZoneSR==10000);
        R_PR(i,d)=mean(preSR(preSR>0.1)); 
        
        ZoneSB=Cbuffer_sea(:,:,d-8);
        preSB=prerate(ZoneSB==10000);
        B_PR(i,d)=mean(preSB(preSB>0.1));
        end
    end
end
R_PR(isnan(R_PR))=0;
B_PR(isnan(B_PR))=0;
ATCP(isnan(ATCP))=0;
%������ƽ�������Իع�
year1=KEEP_Y;   %xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA2.xlsx',2,'C2401:C4340');

R_PR_YAVE=zeros(19,28);%-800��2000km
B_PR_YAVE=zeros(19,28);
parfor yi=2001:2019
    for DR=1:28
       R1C=R_PR(:,DR); 
       B1C=B_PR(:,DR);
       R1C=R1C(year1==yi);
       B1C=B1C(year1==yi);
       R_PR_YAVE(yi-2000,DR)=mean(R1C);
       B_PR_YAVE(yi-2000,DR)=mean(B1C); 
    end
end

PR_Regress=zeros(28,9);
PR_MK=zeros(28,4);
for DInt=1:28
[b,bint,r,rint,stats]=regress(R_PR_YAVE(:,DInt),[ones(2019-2001+1,1),(2001:2019)'],0.05);
    PR_Regress(DInt,1)=b(2);
    PR_Regress(DInt,2)=bint(2,1);%ע��
    PR_Regress(DInt,3)=bint(2,2);
    PR_Regress(DInt,4)=stats(3); 
    
[ci_lower,ci_upper,pvals,b_Test,DW_p] = NeweyWestAdjust(r,(2001:2019)',R_PR_YAVE(:,DInt),1,1);
    PR_Regress(DInt,5)=b_Test(2);
    PR_Regress(DInt,6)=ci_lower(2);%ע��
    PR_Regress(DInt,7)=ci_upper(2);
    PR_Regress(DInt,8)=pvals(2);
    PR_Regress(DInt,9)=DW_p;    
    
 [Zs, p_value, UFk, UBk2, beta, beta_CI, corr_val]= MKtrend((2001:2019)',R_PR_YAVE(:,DInt));
    PR_MK(DInt,1)=p_value;
    PR_MK(DInt,2)=beta;
    PR_MK(DInt,3)=beta_CI(1);
    PR_MK(DInt,4)=beta_CI(2);
end
%DW����ؼ����Լ�Nw����
for i=1:28
    for j=2:4
        if ~isnan(PR_Regress(i,j+4))
            PR_Regress(i,j)=PR_Regress(i,j+4);
        end
    end
end
plot(PR_Regress(:,1))
% plot(mean(R_PR_YAVE))
R_distribution=(mean(R_PR_YAVE))';
