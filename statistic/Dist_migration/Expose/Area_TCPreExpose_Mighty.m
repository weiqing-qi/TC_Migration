%��½�������밶��ÿ�ٹ��ﻷ��������ˮ��Χ
clear; clc; close all;
tic;
file_idir='D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%ע���ǲ��ǰ�˳�򴢴��
%½�ع�������
[~,~,Clandmask] = read_ARCascii('D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 �����������ݵ�extend�����޸�
[GA05,~] = Gridarea(0.5);
GA05Land=zeros(360,720);
GA05Land(Clandmask>0)=GA05(Clandmask>0);%ȫ��½���ϵ����������Ϊ��
num_Land_TC=load('num_Land_TC.mat').num_Land_TC;
% num_Land_TC(num_Land_TC==0)=1;

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
% patchdir='F:\GlobalTCMask1\single_pre\IMERG_Single_Precip3d_mmd\'; %��
indir='D:\DATA\TC_spatial_data\PreRate\CMIP6\BCC-CSM2-MR\'; %��

AVE_TC_SWITCH = false;
StartY=1980;%��
EndY=2014;%��
MaxPre=zeros(360,720,EndY-StartY+1);
for y=StartY:EndY
    files=dir([indir,'*',num2str(y),'*.txt']);%MM/Dע�ⵥλ %��    
    % patchfile=dir([patchdir,'IMERGV6_Ori_SingleTC_3dPrecip_Ummd_SID',num2str(y),'*.txt']);%MM/Dע�ⵥλ %��
  
    for i=1:length(files)
    [~,H,pre2] = read_ARCascii([indir,files(i).name]);%��λmm/D
    %Ԥ����
    pre2 = FillView2Global(H,pre2);
    % pre2 = areaweight_downscale(pre2, 0.25, 0.5);
    pre2 = resample_cmip6_data_slice(pre2, 180, 360, 0.5, 0.5, 'nearest');
    pre2(isnan(pre2))=0; pre2(pre2<0.1)=0;
    %-----------------------����50-60NS-------------------------- ��
    % [~,~,patch] = read_ARCascii([patchdir,patchfile(i).name]);%��λmm/h
    % patch(isnan(patch))=0; patch(patch<0)=0;
    % pre2(61:80,:)=patch(61:80,:);%��IMERG������Թص�
    % pre2(281:300,:)=patch(281:300,:);%��IMERG����
    %-----------------------------------------------------
    pre1=MaxPre(:,:,y-1979);%��
    pre1(pre2>pre1)=pre2(pre2>pre1);
    MaxPre(:,:,y-1979)=pre1;%��
    end
end
%% ����仯
Area=zeros(EndY-StartY+1,1);%ȫ�������
Brea=zeros(EndY-StartY+1,7);%��������
R_Area=zeros(EndY-StartY+1,8);%RINGȫ�������
R_Brea=zeros(EndY-StartY+1,8,7);%RING��������

A=sort(MaxPre(MaxPre>0.1));
PreBar=A(round(length(A)*0.75)); PreBar%����75��λ
PreBar=quantile(A,0.05); PreBar%����75��λ
% PreBar=30;  %��

BasMas = Basinmasks_EPNA(0.5,1);
for y2=StartY:EndY
    MP=MaxPre(:,:,y2-StartY+1);
    Area(y2-StartY+1,1)=sum(GA05Land(MP>=PreBar));%ȫ�������

    for bas=1:7
        GA05LB=zeros(360,720);%�������½���������Ϊ��
        GA05LB(BasMas==bas)=GA05Land(BasMas==bas);
        Brea(y2-StartY+1,bas)=sum(GA05LB(MP>=PreBar));%��������
    %------ÿ����
        for d=1:8    %800-700...100-0km to coast 
            ZoneLR=CRing_land(:,:,9-d);
            GA05LBR=zeros(360,720);
            GA05LBR(ZoneLR==10000)=GA05LB(ZoneLR==10000);
            R_Brea(y2-StartY+1,d,bas)=sum(GA05LBR(MP>=PreBar)); 
            
            if bas==1
                GA05LR=zeros(360,720);
                GA05LR(ZoneLR==10000)=GA05(ZoneLR==10000);%�������������������Ϊ��
                R_Area(y2-StartY+1,d)=sum(GA05LR(MP>=PreBar));%RINGȫ�������
            end
        end
    end
end
% %--------�Ƿ�ƽ��������½��������,���ƽ���Ļ���Ҫ����������������ص�
% if AVE_TC_SWITCH
%     for y2=StartY:EndY
%         MP=MaxPre(:,:,y2-StartY+1);
%         Area(y2-StartY+1,1)=sum(GA05Land(MP>=PreBar));%ȫ�������
% 
%         for bas=1:7
%             GA05LB=zeros(360,720);%�������½���������Ϊ��
%             GA05LB(BasMas==bas)=GA05Land(BasMas==bas);
%             Brea(y2-StartY+1,bas)=sum(GA05LB(MP>=PreBar));%��������
%         %------ÿ����
%             for d=1:8    %800-700...100-0km to coast 
%                 ZoneLR=CRing_land(:,:,9-d);
%                 GA05LBR=zeros(360,720);
%                 GA05LBR(ZoneLR==10000)=GA05LB(ZoneLR==10000);
%                 R_Brea(y2-StartY+1,d,bas)=sum(GA05LBR(MP>=PreBar)); 
% 
%                 if bas==1
%                     GA05LR=zeros(360,720);
%                     GA05LR(ZoneLR==10000)=GA05(ZoneLR==10000);%�������������������Ϊ��
%                     R_Area(y2-StartY+1,d)=sum(GA05LR(MP>=PreBar));%RINGȫ�������
%                 end
%             end
%         end
%     end
%     Area=Area./ num_Land_TC(StartY-1979:EndY-1979,8);
%     R_Area=R_Area./ repmat(num_Land_TC(StartY-1979:EndY-1979,8),1,8);
%     Brea=Brea./num_Land_TC(StartY-1979:EndY-1979,1:7);
%     for j = 1:7
%         R_Brea(:,:,j)=R_Brea(:,:,j)./repmat(num_Land_TC(StartY-1979:EndY-1979,j),1,8);
%     end
% Area(Area<=0)=missing;
% Brea(Brea<=0)=missing;
% R_Area(R_Area<=0)=missing;
% R_Brea(R_Brea<=0)=missing;
% end
%----------------
GAR_Regress=zeros(8,9);
BAR_Regress=zeros(8,9,7);
for i=1:8
[b,bint,r1,rint,stats]=regress(R_Area(:,i),[ones(EndY-StartY+1,1),(StartY:EndY)'],0.05);
GAR_Regress(i,1)=b(2);
GAR_Regress(i,2)=bint(2,1);%ע��
GAR_Regress(i,3)=bint(2,2);
GAR_Regress(i,4)=stats(3); 

[Gci_lower,Gci_upper,Gpvals,Gb_Test,GDW_p] = NeweyWestAdjust(r1,(StartY:EndY)',R_Area(:,i),1,1);
    GAR_Regress(i,5)=Gb_Test(2);
    GAR_Regress(i,6)=Gci_lower(2);%ע��
    GAR_Regress(i,7)=Gci_upper(2);
    GAR_Regress(i,8)=Gpvals(2);
    GAR_Regress(i,9)=GDW_p;
    for bas=1:7
    [b,bint,r2,rint,stats]=regress(R_Brea(:,i,bas),[ones(EndY-StartY+1,1),(StartY:EndY)'],0.05);
    BAR_Regress(i,1,bas)=b(2);
    BAR_Regress(i,2,bas)=bint(2,1);%ע��
    BAR_Regress(i,3,bas)=bint(2,2);
    BAR_Regress(i,4,bas)=stats(3); 
    
    [Bci_lower,Bci_upper,Bpvals,Bb_Test,BDW_p] = NeweyWestAdjust(r2,(StartY:EndY)',R_Brea(:,i,bas),1,1);
    BAR_Regress(i,5,bas)=Bb_Test(2);
    BAR_Regress(i,6,bas)=Bci_lower(2);%ע��
    BAR_Regress(i,7,bas)=Bci_upper(2);
    BAR_Regress(i,8,bas)=Bpvals(2);
    BAR_Regress(i,9,bas)=BDW_p;
    end
end
out1=[Area,R_Area,Brea];

out2=[];out2DT=[];out2_P=[];out2_DWP=[];out2_Padj=[];
for bas=1:7
    for i=1:3
    out2=[out2,BAR_Regress(:,i,bas)];
    out2DT=[out2DT,BAR_Regress(:,i+4,bas)];
    end
    out2_P=[out2_P,BAR_Regress(:,4,bas)];
    out2_DWP=[out2_DWP,BAR_Regress(:,9,bas)];
    out2_Padj=[out2_Padj,BAR_Regress(:,8,bas)];
end
out2_P(:,3)=[];out2_Padj(:,3)=[];out2_DWP(:,3)=[];
out2_P(out2_DWP>0 & out2_DWP<0.05)=out2_Padj(out2_DWP>0 & out2_DWP<0.05);%�滻���޸ĵ�ֵ

out2(:,7:9)=[];out2DT(:,7:9)=[];
out2(out2DT>0 | out2DT<0)=out2DT(out2DT>0 | out2DT<0);%�滻���޸ĵ�ֵ

%DW����ؼ����Լ�Nw����
for i=1:8
    for j=2:4
        if ~isnan(GAR_Regress(i,j+4))
            GAR_Regress(i,j)=GAR_Regress(i,j+4);
        end
    end
end
out3=GAR_Regress;

out4=mean(R_Brea,1);%������ƽ���ķֲ�
out4=permute(out4,[2,3,1]);out4(:,3)=[];

toc;