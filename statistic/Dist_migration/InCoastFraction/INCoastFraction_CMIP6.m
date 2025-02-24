%������뺣����ͬ���뻷�����ˮռ�������ı����仯
%����ÿ�������ڵĽ�ˮ��ǿ�ȱ仯����Ϊ���ں��ܾ�����
%ע��������60NS
clear; clc; close all;
% parpool(15);
file_idir='D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%ע���ǲ��ǰ�˳�򴢴��

%½�ع�������
[~,~,Clandmask] = read_ARCascii('D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 �����������ݵ�extend�����޸�

% ������뻷�ͺ����߽��� 
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

%% Filter
Ylim=[1980 2014];
allmodel_name={'BCC-CSM2-MR','CAS-FGOALS-f3-L','CMCC-ESM2','CNRM-CM6-1-HR','EC-Earth3', ...
            'EC-Earth3-AerChem','EC-Earth3-CC','EC-Earth3-Veg','NOAA-GFDL-CM4','SNU-SAM0-UNICON'};
MOD=10;
model_name = allmodel_name(MOD);
% dist=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted_CMIP6.xlsx',1,'B3:B5000');
% dist=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted.xlsx','PAmount','C3:C4461');%������Ҫû����filted���ܶ�Ӧ��
KEEP_Y=[];
% KEEP_Y=xlsread('D:\Desktop2\Global_cyclone_project\DATA_CMIP6.xlsx','P1G',[num2abc2(MOD*5-4),'3:',num2abc2(MOD*5-4),'5000']);
% TC_idir='D:\DATA\Segmentation map tutorial\NEW\IBTrACS.ALL.v04r00.nc';
% [~,~,~,ISO_TIME,USA_WIND,~,BASIN,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

% [KEEP_NO,KEEP_Y,KEEP_B,KEEP_i] = FilterNO(Ylim,ISO_TIME,USA_WIND,BASIN,dist);
%��ˮ����
%Preindir='E:\DATA\TCMigrationPJ\JRA55_Single_Precip3d_mmd\'; %��
Preindir=['D:\DATA\TC_spatial_data\PreRate\CMIP6\',model_name{:},'\']; %��
Prefile=dir([Preindir,'*.txt']);
% Prefile=Prefile(2060:3999);%��1980��ʼ1955ERA5 2898JRA55
% Prefile=Prefile(KEEP_i);



%% ƽ����ˮ�ʽ���
R_PR=zeros(length(Prefile),28);%-800��2000km
B_PR=zeros(length(Prefile),28);
COUNT=zeros(length(Prefile),28);%���������Ľ�ˮ��
%ÿ��TCÿһRorB�е�ƽ��ǿ��
parfor i=1:length(Prefile)
    i
    [~,H,prerate] = read_ARCascii([Preindir,Prefile(i).name]);%��λmm/D
    prerate = FillView2Global(H,prerate);
    [rows,cols]=size(prerate);
    prerate = resample_cmip6_data_slice(prerate, rows, cols, 0.5, 0.5, 'nearest');
    % prerate = areaweight_downscale(prerate, 0.25, 0.5);
    prerate(isnan(prerate))=0; prerate(prerate<2.4)=0;%��ֵ
    % sstprerate=prerate;sstprerate(Clandmask>0)=0 %�㺣��ȥ��½��
    % ATCP(i)=mean(sstprerate(sstprerate>0.1));
    %ATCP(i)=mean(prerate(prerate>0.1));
        for d=1:28%100:2000km
            if d<=8                               %Land
            ZoneLR=CRing_land(:,:,9-d);          %ע�⽵ˮת��
            preLR=prerate(ZoneLR==10000);         %������������н�ˮ
            R_PR(i,d)=mean(preLR(preLR>0.1));     %0.1���ż���ƽ����ǿ
            COUNT(i,d) = length(preLR(preLR>0.1));
    
            % ZoneLB=Cbuffer_land(:,:,9-d);
            % preLB=prerate(ZoneLB==10000);
            % B_PR(i,d)=mean(preLB(preLB>0.1));
            else                                  %Sea
            ZoneSR=CRing_sea(:,:,d-8);
            preSR=prerate(ZoneSR==10000);
            R_PR(i,d)=mean(preSR(preSR>0.1)); 
            COUNT(i,d) = length(preLR(preLR>0.1));
    
            % ZoneSB=Cbuffer_sea(:,:,d-8);
            % preSB=prerate(ZoneSB==10000);
            % B_PR(i,d)=mean(preSB(preSB>0.1));
            end
        end
end
R_PR(isnan(R_PR))=0;%û�о���0
B_PR(isnan(B_PR))=0;
COUNT(isnan(COUNT))=0;
%������ƽ�������Իع�
year1=KEEP_Y;%xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA3TMPA60NS.xlsx',2,'C3:C4461');

R_PR_YAVE=zeros(Ylim(2)-Ylim(1)+1,28);%-800��2000km
B_PR_YAVE=zeros(Ylim(2)-Ylim(1)+1,28);
for yi=Ylim(1):Ylim(2)
    for DR=1:28
       R1C=R_PR(:,DR); 
       B1C=B_PR(:,DR);
       R1C=R1C(year1==yi);
       B1C=B1C(year1==yi);
       R_PR_YAVE(yi-Ylim(1)+1,DR)=mean(R1C);
       B_PR_YAVE(yi-Ylim(1)+1,DR)=mean(B1C); 
    end
end

PR_Regress=zeros(28,9);
PR_MK=zeros(28,4);
for DInt=1:28
[b,bint,r,rint,stats]=regress(R_PR_YAVE(:,DInt),[ones(Ylim(2)-Ylim(1)+1,1),(Ylim(1):Ylim(2))'],0.05);
    PR_Regress(DInt,1)=b(2);
    PR_Regress(DInt,2)=bint(2,1);%ע��
    PR_Regress(DInt,3)=bint(2,2);
    PR_Regress(DInt,4)=stats(3); 
    
[ci_lower,ci_upper,pvals,b_Test,DW_p] = NeweyWestAdjust(r,(Ylim(1):Ylim(2))',R_PR_YAVE(:,DInt),1,1);
    PR_Regress(DInt,5)=b_Test(2);
    PR_Regress(DInt,6)=ci_lower(2);%ע��
    PR_Regress(DInt,7)=ci_upper(2);
    PR_Regress(DInt,8)=pvals(2);
    PR_Regress(DInt,9)=DW_p;
    
 [Zs, p_value, UFk, UBk2, beta, beta_CI, corr_val]= MKtrend((Ylim(1):Ylim(2))',R_PR_YAVE(:,DInt));
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

% %---------------------200km 
% R_PR_200=zeros(length(Prefile),14);%-800��2000km
% for i=1:14
%     R_PR_200(:,i)=(R_PR(:,2*i-1).*COUNT(:,2*i-1) + R_PR(:,2*i).*COUNT(:,2*i)) ./ sum(COUNT(:,2*i-1:2*i),2);
% end
% R_PR_200(isnan(R_PR_200))=0;
% R_PR_YAVE_200=zeros(Ylim(2)-Ylim(1)+1,14);%-800��2000km
% for yi=Ylim(1):Ylim(2)
%     for DR=1:14
%        R1C=R_PR_200(:,DR); 
%        R1C=R1C(year1==yi);
%        R_PR_YAVE_200(yi-Ylim(1)+1,DR)=mean(R1C);
%     end
% end
% 
% PR_Regress_200=zeros(14,9);
% PR_MK_200=zeros(14,4);
% for DInt=1:14
% [b,bint,r,rint,stats]=regress(R_PR_YAVE_200(:,DInt),[ones(Ylim(2)-Ylim(1)+1,1),(Ylim(1):Ylim(2))'],0.05);
%     PR_Regress_200(DInt,1)=b(2);
%     PR_Regress_200(DInt,2)=bint(2,1);%ע��
%     PR_Regress_200(DInt,3)=bint(2,2);
%     PR_Regress_200(DInt,4)=stats(3); 
% 
% [ci_lower,ci_upper,pvals,b_Test,DW_p] = NeweyWestAdjust(r,(Ylim(1):Ylim(2))',R_PR_YAVE_200(:,DInt),1,1);
%     PR_Regress_200(DInt,5)=b_Test(2);
%     PR_Regress_200(DInt,6)=ci_lower(2);%ע��
%     PR_Regress_200(DInt,7)=ci_upper(2);
%     PR_Regress_200(DInt,8)=pvals(2);
%     PR_Regress_200(DInt,9)=DW_p;
% 
%  [Zs, p_value, UFk, UBk2, beta, beta_CI, corr_val]= MKtrend((Ylim(1):Ylim(2))',R_PR_YAVE_200(:,DInt));
%     PR_MK_200(DInt,1)=p_value;
%     PR_MK_200(DInt,2)=beta;
%     PR_MK_200(DInt,3)=beta_CI(1);
%     PR_MK_200(DInt,4)=beta_CI(2);
% end
% %DW����ؼ����Լ�Nw����
% for i=1:14
%     for j=2:4
%         if ~isnan(PR_Regress_200(i,j+4))
%             PR_Regress_200(i,j)=PR_Regress_200(i,j+4);
%         end
%     end
% end
% figure
% plot(PR_Regress_200(:,1))
% % plot(mean(R_PR_YAVE))
% R_distribution_200=(mean(R_PR_YAVE_200))';
