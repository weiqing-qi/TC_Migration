%������뺣����ͬ���뻷��ȫ���º��ܿɽ�ˮ���仯
%ע��������60NS
clear; clc; close all;
file_idir='D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%ע���ǲ��ǰ�˳�򴢴��
%�������������յ��ַ�������
Ylim=[1980 2020];
% ������������
dates = datetime(Ylim(1), 1, 1) : datetime(Ylim(2), 12, 31);
% ����������ת��Ϊ�ַ������飬��ʽΪ yyyymmdd
daycell = cellstr(datestr(dates, 'yyyymmdd'));
monthcell = unique(cellstr(datestr(dates, 'yyyymm')));


dataname='JRA55TPW';%��
%��ˮ����
if strcmp(dataname,'ERA5SST')
    indir='E:\DATA\1.Reanalysis\EAR5\ERA5-Sea_surface_temperature\'; 
    nametext = 'ERA5_Sea_Surface_Temperature_on_single_levels_daily';
    Yearnum = floor(str2double(monthcell)/100);
elseif strcmp(dataname,'ERA5TPW')
    indir='E:\DATA\1.Reanalysis\EAR5\ERA5-Total-Column-Water\'; 
    nametext = 'ERA5_Total_Column_Water_on_single_levels_daily';
    Yearnum = floor(str2double(monthcell)/100);
elseif strcmp(dataname,'JRA55SST')
    indir='E:\DATA\1.Reanalysis\JRA55\JRA55Bilinear05SST\'; 
    nametext = 'JRA55_Daily_Bilinear_Ori05_SST';
    Yearnum =  floor(str2double(daycell)/10000);
elseif strcmp(dataname,'JRA55TPW')
    indir='E:\DATA\1.Reanalysis\JRA55\JRA55TPW_Bilinear05\'; 
    nametext = 'JRA55_Daily_Bilinear_Ori05_TPW';
    Yearnum =  floor(str2double(daycell)/10000);
end
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

%% ƽ����ˮ�ʽ���


%ÿ��TCÿһRorB�е�ƽ��ǿ��
if strncmp(dataname,'ERA5',4)
    lengthfile=length(monthcell);
    R_PR=zeros(lengthfile,28);%-800��2000km

elseif strncmp(dataname,'JRA5',4)
    lengthfile=length(daycell);
    R_PR=zeros(lengthfile,28);%-800��2000km

end

parfor i=1:lengthfile
    i
    % if strcmp(dataname,'ERA5SST')
    %     timecell=monthcell(i);
    %     Amonthdata = permute(ncread([indir,nametext,timecell{:},'.nc'],'tos'),[2 1 3]);
    %     Amonthdata = flip(Amonthdata(1:360,:,:)-273.15 , 1);
    %     prerate = mean(Amonthdata,3);%���϶ȣ�90--90 -180-180
    %     prerate(Clandmask>0)=missing;
    % 
    % elseif strcmp(dataname,'ERA5TPW')
    %     timecell=monthcell(i);
    %     Amonthdata = permute(ncread([indir,nametext,timecell{:},'.nc'],'tcw'),[2 1 3]);
    %     Amonthdata = flip(Amonthdata(1:360,:,:) , 1);
    %     prerate = mean(Amonthdata,3);%MM��90--90 -180-180
    % 
    % else
    if strncmp(dataname,'JRA55',5)
        timecell=daycell(i);
        prerate =  imread([indir,nametext,timecell{:},'.tif']);%���϶ȣ�90--90 -180-180
        if strcmp(dataname,'JRA55SST')
            prerate(Clandmask>0)=missing;
        end
    end
    prerate(prerate>100)=missing; prerate(prerate<-50)=missing;

    for d=1:28%100:2000km
        if d<=8                               %Land
        ZoneLR=CRing_land(:,:,9-d);          %ע�⽵ˮת��
        preLR=prerate(ZoneLR==10000);         %������������н�ˮ
        R_PR(i,d)=mean(preLR,"omitnan");     %0.1���ż���ƽ����ǿ

        else                                  %Sea
        ZoneSR=CRing_sea(:,:,d-8);
        preSR=prerate(ZoneSR==10000);
        R_PR(i,d)=mean(preSR,"omitnan");
        end
    end
end

R_PR_YAVE=zeros(Ylim(2)-Ylim(1)+1,28);%-800��2000km
for yi=Ylim(1):Ylim(2)
    for DR=1:28
       R1C=R_PR(:,DR); 
       R1C=R1C(Yearnum==yi);
       R_PR_YAVE(yi-Ylim(1)+1,DR)=mean(R1C);

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
