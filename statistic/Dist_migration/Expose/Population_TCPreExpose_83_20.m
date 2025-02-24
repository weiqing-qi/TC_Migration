%内陆地区距离岸线每百公里环内经历台风降水的人口
%每个网格每年最大的降水强度
%筛选大于50mm/d的降水量下的人口
%用环剪切得到每年每环人口数
clear; clc; close all;
file_idir='C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%注意是不是按顺序储存的
%陆地国家数据
[~,~,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 根据输入数据的extend进行修改
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
%% 计算一年中最大降水强度分布和人口分布
indir='F:\GlobalTCMask1\single_pre\PERSIANN_Single_Precip3d_mmd\'; %改
MaxPre=zeros(360,720,38);
parfor y=1983:2020
    files=dir([indir,'PERSIANNCDR_Ori_SingleTC_3dPrecip_Ummd_SID',num2str(y),'*.txt']);%MM/D注意单位 %改
    for i=1:length(files)
    [~,~,pre2] = read_ARCascii([indir,files(i).name]);%单位mm/D
    pre1=MaxPre(:,:,y-1982);
    pre1(pre2>pre1)=pre2(pre2>pre1);
    MaxPre(:,:,y-1982)=pre1;
    end
end
%% ------------------计算人口
Popindir='E:\DATA\WORLD_Population\Landscan\Landscan05_txt\';
PopDistri_TCP=zeros(360,720,21);%每年经历台风降水的人
PopAll_TCP=zeros(21,1);%每年总人数
PopRing_TCP=zeros(21,8);   %每环800-700...100-0km to coast
Pop_TCP_Basin=zeros(21,6);   %
PopRing_TCP_Basin=zeros(21,8,6);   %每环大洋800-700...100-0km to coast
A=sort(MaxPre(MaxPre>0.1));
% PreBar=A(round(length(A)*0.75)); %大于75分位  改00000 0000 00000 0000000000 0000 000
PreBar=0.1;
all_in_one = Basinmasks(0.5); all_in_one=[all_in_one(:,361:720),all_in_one(:,1:360)];
parfor y2=2000:2020
POPTCP=zeros(360,720);

[~,~,Pop] = read_ARCascii([Popindir,'Resample05-landscan-global-',num2str(y2),'.txt']);%注意这个数据集除了正值都是零
MP=MaxPre(:,:,y2-1982);
POPTCP(MP>PreBar)=Pop(MP>PreBar); %筛选大于标准的TC降水经历人口分布
PopDistri_TCP(:,:,y2-1999)=POPTCP;
PopAll_TCP(y2-1999)=sum(sum(POPTCP));
for basin1=1:6
Pop_TCP_Basin(y2-1999,basin1)=sum(POPTCP(all_in_one==basin1));
end
%------每个环里的人数
    for d=1:8    %800-700...100-0km to coast 
        ZoneLR=CRing_land(:,:,9-d);
        PopRing_TCP(y2-1999,d)=sum(POPTCP(ZoneLR==10000));   
        %-------分割大洋all_in_one
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
    PR_Regress(i,2)=bint(2,1);%注意
    PR_Regress(i,3)=bint(2,2);
    PR_Regress(i,4)=stats(3); 
end
