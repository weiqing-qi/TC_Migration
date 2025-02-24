%内陆地区距离岸线每百公里环内气旋降水范围
clear; clc; close all;
tic;
file_idir='D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%注意是不是按顺序储存的
%陆地国家数据
[~,~,Clandmask] = read_ARCascii('D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 根据输入数据的extend进行修改
[GA05,~] = Gridarea(0.5);
GA05Land=zeros(360,720);
GA05Land(Clandmask>0)=GA05(Clandmask>0);%全部陆地上的面积，海洋为零
num_Land_TC=load('num_Land_TC.mat').num_Land_TC;
% num_Land_TC(num_Land_TC==0)=1;

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
    CR_BS(:,:,i*0.01)=CB_Far-CB_Near;%真值应该是1-（-9999）=10000
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
%% 计算一年中最大降水强度分布
% patchdir='F:\GlobalTCMask1\single_pre\IMERG_Single_Precip3d_mmd\'; %改
indir='D:\DATA\TC_spatial_data\PreRate\CMIP6\BCC-CSM2-MR\'; %改

AVE_TC_SWITCH = false;
StartY=1980;%改
EndY=2014;%改
MaxPre=zeros(360,720,EndY-StartY+1);
for y=StartY:EndY
    files=dir([indir,'*',num2str(y),'*.txt']);%MM/D注意单位 %改    
    % patchfile=dir([patchdir,'IMERGV6_Ori_SingleTC_3dPrecip_Ummd_SID',num2str(y),'*.txt']);%MM/D注意单位 %改
  
    for i=1:length(files)
    [~,H,pre2] = read_ARCascii([indir,files(i).name]);%单位mm/D
    %预处理
    pre2 = FillView2Global(H,pre2);
    % pre2 = areaweight_downscale(pre2, 0.25, 0.5);
    pre2 = resample_cmip6_data_slice(pre2, 180, 360, 0.5, 0.5, 'nearest');
    pre2(isnan(pre2))=0; pre2(pre2<0.1)=0;
    %-----------------------补充50-60NS-------------------------- 改
    % [~,~,patch] = read_ARCascii([patchdir,patchfile(i).name]);%单位mm/h
    % patch(isnan(patch))=0; patch(patch<0)=0;
    % pre2(61:80,:)=patch(61:80,:);%用IMERG补充可以关掉
    % pre2(281:300,:)=patch(281:300,:);%用IMERG补充
    %-----------------------------------------------------
    pre1=MaxPre(:,:,y-1979);%改
    pre1(pre2>pre1)=pre2(pre2>pre1);
    MaxPre(:,:,y-1979)=pre1;%改
    end
end
%% 面积变化
Area=zeros(EndY-StartY+1,1);%全球总面积
Brea=zeros(EndY-StartY+1,7);%盆地总面积
R_Area=zeros(EndY-StartY+1,8);%RING全球总面积
R_Brea=zeros(EndY-StartY+1,8,7);%RING盆地总面积

A=sort(MaxPre(MaxPre>0.1));
PreBar=A(round(length(A)*0.75)); PreBar%大于75分位
PreBar=quantile(A,0.05); PreBar%大于75分位
% PreBar=30;  %改

BasMas = Basinmasks_EPNA(0.5,1);
for y2=StartY:EndY
    MP=MaxPre(:,:,y2-StartY+1);
    Area(y2-StartY+1,1)=sum(GA05Land(MP>=PreBar));%全球总面积

    for bas=1:7
        GA05LB=zeros(360,720);%单个盆地陆地面积，海为零
        GA05LB(BasMas==bas)=GA05Land(BasMas==bas);
        Brea(y2-StartY+1,bas)=sum(GA05LB(MP>=PreBar));%盆地总面积
    %------每个环
        for d=1:8    %800-700...100-0km to coast 
            ZoneLR=CRing_land(:,:,9-d);
            GA05LBR=zeros(360,720);
            GA05LBR(ZoneLR==10000)=GA05LB(ZoneLR==10000);
            R_Brea(y2-StartY+1,d,bas)=sum(GA05LBR(MP>=PreBar)); 
            
            if bas==1
                GA05LR=zeros(360,720);
                GA05LR(ZoneLR==10000)=GA05(ZoneLR==10000);%单个环的面积网格，其他为零
                R_Area(y2-StartY+1,d)=sum(GA05LR(MP>=PreBar));%RING全球总面积
            end
        end
    end
end
% %--------是否平均气旋登陆场次数量,如果平均的化需要允许气旋存在面积重叠
% if AVE_TC_SWITCH
%     for y2=StartY:EndY
%         MP=MaxPre(:,:,y2-StartY+1);
%         Area(y2-StartY+1,1)=sum(GA05Land(MP>=PreBar));%全球总面积
% 
%         for bas=1:7
%             GA05LB=zeros(360,720);%单个盆地陆地面积，海为零
%             GA05LB(BasMas==bas)=GA05Land(BasMas==bas);
%             Brea(y2-StartY+1,bas)=sum(GA05LB(MP>=PreBar));%盆地总面积
%         %------每个环
%             for d=1:8    %800-700...100-0km to coast 
%                 ZoneLR=CRing_land(:,:,9-d);
%                 GA05LBR=zeros(360,720);
%                 GA05LBR(ZoneLR==10000)=GA05LB(ZoneLR==10000);
%                 R_Brea(y2-StartY+1,d,bas)=sum(GA05LBR(MP>=PreBar)); 
% 
%                 if bas==1
%                     GA05LR=zeros(360,720);
%                     GA05LR(ZoneLR==10000)=GA05(ZoneLR==10000);%单个环的面积网格，其他为零
%                     R_Area(y2-StartY+1,d)=sum(GA05LR(MP>=PreBar));%RING全球总面积
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
GAR_Regress(i,2)=bint(2,1);%注意
GAR_Regress(i,3)=bint(2,2);
GAR_Regress(i,4)=stats(3); 

[Gci_lower,Gci_upper,Gpvals,Gb_Test,GDW_p] = NeweyWestAdjust(r1,(StartY:EndY)',R_Area(:,i),1,1);
    GAR_Regress(i,5)=Gb_Test(2);
    GAR_Regress(i,6)=Gci_lower(2);%注意
    GAR_Regress(i,7)=Gci_upper(2);
    GAR_Regress(i,8)=Gpvals(2);
    GAR_Regress(i,9)=GDW_p;
    for bas=1:7
    [b,bint,r2,rint,stats]=regress(R_Brea(:,i,bas),[ones(EndY-StartY+1,1),(StartY:EndY)'],0.05);
    BAR_Regress(i,1,bas)=b(2);
    BAR_Regress(i,2,bas)=bint(2,1);%注意
    BAR_Regress(i,3,bas)=bint(2,2);
    BAR_Regress(i,4,bas)=stats(3); 
    
    [Bci_lower,Bci_upper,Bpvals,Bb_Test,BDW_p] = NeweyWestAdjust(r2,(StartY:EndY)',R_Brea(:,i,bas),1,1);
    BAR_Regress(i,5,bas)=Bb_Test(2);
    BAR_Regress(i,6,bas)=Bci_lower(2);%注意
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
out2_P(out2_DWP>0 & out2_DWP<0.05)=out2_Padj(out2_DWP>0 & out2_DWP<0.05);%替换了修改的值

out2(:,7:9)=[];out2DT(:,7:9)=[];
out2(out2DT>0 | out2DT<0)=out2DT(out2DT>0 | out2DT<0);%替换了修改的值

%DW自相关检验以及Nw修正
for i=1:8
    for j=2:4
        if ~isnan(GAR_Regress(i,j+4))
            GAR_Regress(i,j)=GAR_Regress(i,j+4);
        end
    end
end
out3=GAR_Regress;

out4=mean(R_Brea,1);%所有年平均的分布
out4=permute(out4,[2,3,1]);out4(:,3)=[];

toc;