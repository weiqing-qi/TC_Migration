%计算进入海岸不同距离环的最大降水占年总数的比例变化 仅计算basinl的变化
%计算每个环以内的降水率强度变化，分为环内和总距离内
%注意区域是60NS 
clear; clc; close all;

file_idir='D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%注意是不是按顺序储存的
%降水数据
Preindir='D:\DATA\TC_spatial_data\PreRate\PERSIANN\'; %改
Prefile=dir([Preindir,'*.txt']);
Prefile=Prefile(2286:4225);%改1980开始1955ERA5 2898JRA55
Ylim=[2001 2019];

year1=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted.xlsx',1,'AR2288:AR4227');
%陆地国家数据
[~,~,Clandmask] = read_ARCascii('D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 根据输入数据的extend进行修改
BasMas = Basinmasks_EPNA(0.5,1);
%% 计算距离环和海岸边界区 
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
    CR_BS(:,:,i*0.01)=CB_Far-CB_Near;%真值应该是1-（-9999）10000
end
    CR_BS(1:60,:,:)=0; CR_BS(301:360,:,:)=0;%去掉60ns外的 
    
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
%% 平均降水率进入
R_PR=zeros(length(Prefile),28,7);%-800：2000km  %SImask:1 SPmask:2 SAmask:3 NImask:4 WPmask:5 EPmask:6 NAmask:7 
tic;
%每个TC每一RorB中的平均强度
parfor i=1:length(Prefile)
[~,H,prerate] = read_ARCascii([Preindir,Prefile(i).name]);%单位mm/D
prerate = FillView2Global(H,prerate);
% prerate = areaweight_downscale(prerate, 0.25, 0.5)
prerate(isnan(prerate))=0; prerate(prerate<2.4)=0;
    for bas=1:7
        for d=1:28%100:2000km
            if d<=8                               %Land
            ZoneLR=CRing_land(:,:,9-d);           %注意降水转换
            ZoneLR(BasMas~=bas)=0;                %计算的盆地留出来
            preLR=prerate(ZoneLR==10000);         %这个区域内所有降水
            R_PR(i,d,bas)=mean(preLR(preLR>0.1));     %0.1的门槛，平均雨强

            else                                  %Sea
            ZoneSR=CRing_sea(:,:,d-8);
            ZoneSR(BasMas~=bas)=0;                %计算的盆地留出来
            preSR=prerate(ZoneSR==10000);
            R_PR(i,d,bas)=mean(preSR(preSR>0.1)); 
            end
        end
    end
end
toc;
R_PR(isnan(R_PR))=0;
JudgeB=permute(sum(R_PR,2),[1 3 2]);%所有距离区间加在一起，并去掉距离区间的数组维度 %每个盆地的分母是每个盆地

%按年求平均再线性回归

R_PR_YAVE=zeros(Ylim(2)-Ylim(1)+1,28,7);%-800：2000km
for bas=1:7
    for yi=Ylim(1):Ylim(2)
        Jud=JudgeB(:,bas);
        Jud=Jud(year1==yi);%判断每年这个BAS里有没有降水

        for DR=1:28
           R1C=R_PR(:,DR,bas); 
           R1C=R1C(year1==yi); %R1C没有空值
           R1C=R1C(Jud>0);%只要这个盆地里面出现过降水的数据没出现的不要%开关 决定了平均数量
           R_PR_YAVE(yi-Ylim(1)+1,DR,bas)=mean(R1C);
        end
    end
end

PR_Regress=zeros(28,9,7);
PR_MK=zeros(28,4,7);
for bas=1:7
    for DInt=1:28
    [b,bint,r,rint,stats]=regress(R_PR_YAVE(:,DInt,bas),[ones(Ylim(2)-Ylim(1)+1,1),(Ylim(1):Ylim(2))'],0.05);
        PR_Regress(DInt,1,bas)=b(2);
        PR_Regress(DInt,2,bas)=bint(2,1);%注意
        PR_Regress(DInt,3,bas)=bint(2,2);
        PR_Regress(DInt,4,bas)=stats(3); 
        
    [ci_lower,ci_upper,pvals,b_Test,DW_p] = NeweyWestAdjust(r,(Ylim(1):Ylim(2))',R_PR_YAVE(:,DInt,bas),1,1);
        PR_Regress(DInt,5,bas)=b_Test(2);
        PR_Regress(DInt,6,bas)=ci_lower(2);%注意
        PR_Regress(DInt,7,bas)=ci_upper(2);
        PR_Regress(DInt,8,bas)=pvals(2);
        PR_Regress(DInt,9,bas)=DW_p;

     [Zs, p_value, UFk, UBk2, beta, beta_CI, corr_val]= MKtrend((Ylim(1):Ylim(2))',R_PR_YAVE(:,DInt,bas));
        PR_MK(DInt,1,bas)=p_value;
        PR_MK(DInt,2,bas)=beta;
        PR_MK(DInt,3,bas)=beta_CI(1);
        PR_MK(DInt,4,bas)=beta_CI(2);
    end
end
%DW自相关检验以及Nw修正
for bas=1:7
    for i=1:28
        for j=2:4
            if ~isnan(PR_Regress(i,j+4,bas))
                PR_Regress(i,j)=PR_Regress(i,j+4,bas);
            end
        end
    end
end


plot(PR_Regress(:,1,5))
% plot(permute(PR_Regress(:,1,:),[1 3 2]))
% plot(mean(R_PR_YAVE(:,:,5)))
R_distribution=(mean(R_PR_YAVE,1,"omitmissing"));
R_distribution=(permute(R_distribution(1,:,:),[3 2 1]))';
R_distribution(:,3)=[];
OUT=zeros(28,21);
for j=1:7
    OUT(:,3*j-2:3*j)=PR_Regress(:,1:3,j);
end
OUT(:,7:9)=[];
imshow(OUT)

SIG=zeros(28,7);%显著性
for j=1:7
    SIG(:,j)=PR_Regress(:,4,j);
end
SIG(:,3)=[];
toc;
%------------------------DW OUT PUT----------------------------不需要了
% DW=zeros(28,7);%显著性
% for j=1:7
%     DW(:,j)=PR_Regress(:,9,j);
% end
% DW(:,3)=[];
% DWTEST=find(DW < 0.05 & DW > 0)
% DWOUT=zeros(28,28);
% for j=1:7
%     DWOUT(:,4*j-3:4*j)=PR_Regress(:,5:8,j);
% end
% DWOUT(:,9:12)=[];
