%判断气旋是否进入陆地距离范围内
%% 计算距离环和海岸边界区 
file_idir='C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%注意是不是按顺序储存的
%陆地国家数据
[~,~,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 根据输入数据的extend进行修改

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
    CR_BS(:,:,i*0.01)=CB_Far-CB_Near;%真值应该是1-(-9999)=10000
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
imshow(CRing_land(:,:,1))
%% 判断
Preindir='F:\GlobalTCMask1\single_pre\ERA5_Single_Precip3d_mmd\'; %改
Prefile=dir([Preindir,'*.txt']);
switch length(Prefile)
    case 6413 %ERA5
    Prefile=Prefile(1955:end);
    case 7356 %JRA55
    Prefile=Prefile(2898:end);
    case 2061 %IMERG
    Prefile=Prefile(1:1940);
    case 2267 %TMPA
    Prefile=Prefile(328:end);    
end
judgeM=zeros(length(Prefile),1);
tic;
parfor f=1:length(Prefile)
    [~,~,TCpre] = read_ARCascii([Preindir,Prefile(f).name]);
    TCpre(TCpre>=0)=1;TCpre(TCpre<0)=0;
    J=TCpre.*CRing_land(:,:,3);   %改
    if sum(sum(J))>0%是否近陆500km
        judgeM(f)=1;
    else
       judgeM(f)=0; 
    end
end
toc;
%% excel 处理按年求平均再线性回归
%-----1980-2020

% DATA=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA2.xlsx',2,'E3:E4461');%D3:D4461
% year1=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA2.xlsx',2,'C3:C4461');
% 
% year1(DATA>2000)=[];DATA(DATA>2000)=[];%大于2000的不考虑
% year1(DATA<-50)=[];DATA(DATA<-50)=[];%已经确定最大降水登陆的不考虑
% 
% R_PR_YAVE=zeros(41,1);%就是固定的
% for yi=1980:2020  %改这个就可以2001:2019 1980:2020
%     C=DATA(year1==yi);
%     R_PR_YAVE(yi-1979)=mean(C);
% end
%----2001-2019
DATA=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA2.xlsx',2,'I2401:I4340');
year1=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA2.xlsx',2,'C2401:C4340');
year1(DATA>2000)=[];
DATA(DATA>2000)=[];%大于2000的不考虑
year1(DATA<-50)=[];
DATA(DATA<-50)=[];%已经确定最大降水登陆的不考虑

R_PR_YAVE=zeros(19,1);
for yi=2001:2019  
    C=DATA(year1==yi);
    R_PR_YAVE(yi-2000)=mean(C);
end
%-----------显著性回归
[Zs, p_value, UFk, UBk2, beta, beta_CI, corr_val]= MKtrend((2001:2019)',R_PR_YAVE);
[b,bint,r,rint,stats]=regress(R_PR_YAVE,[ones(2019-2001+1,1),(2001:2019)'],0.05);
% [b,bint,r,rint,stats]=regress(R_PR_YAVE,[ones(2020-1980+1,1),(1980:2020)'],0.05);
% [Zs, p_value, UFk, UBk2, beta, beta_CI, corr_val]= MKtrend((1980:2020)',R_PR_YAVE);

[b(2),bint(2,1),bint(2,2),stats(3),b(1),p_value]


