%对比两套降水数据的误差
%% 注意区域是60NS
clear; clc; close all;
file_idir='D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%注意是不是按顺序储存的
%陆地国家数据
[~,~,Clandmask] = read_ARCascii('D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 根据输入数据的extend进行修改
% 计算距离环和海岸边界区 
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
    CR_BS(1:100,:,:)=0; CR_BS(261:360,:,:)=0;%去掉40ns外的 
    
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


%%
%----气旋数据
TC_idir='E:\DATA\5.TCdata\IBTrACS\NEW\IBTrACS.ALL.v04r00.nc';
[SID,LAT,LON,ISO_TIME,USA_WIND,USA_RMW]=IBTrACS_nc_entire_variable_r(TC_idir);

Ylim=[2001 2019];
NOlim=lim_Y2NO(Ylim,ISO_TIME);  

Globaldata1dir='E:\DATA\2.Satellite\PERSIANN\';%原始全球降水
Globaldata2dir='E:\DATA\2.Satellite\PERSIANN-CDR-daily\BIN\';

TCfiles1dir='D:\DATA\TC_spatial_data\PreRate\PERSIANN\';%提取的降水
TCfiles1=dir([TCfiles1dir,'*.txt']);
TCfiles2dir='D:\DATA\TC_spatial_data\PreRate\PERSIANNCDR_PRE025\';%不用变
TCfiles2=[];
for i=Ylim(1):Ylim(2)
    TCfiles2=[TCfiles2;dir([TCfiles2dir,'*SID',num2str(i),'*.txt'])];
end

%% ----------------------------------------------------气旋周围的对比
TC_BIAS=zeros(NOlim(3),28);
TC_RMSE=zeros(NOlim(3),28);
TC_CC=zeros(NOlim(3),28);
TC_ME=zeros(NOlim(3),28);
tic;
parfor i=1:NOlim(3)
    i
    %-------prerate1
    [~,H1,prerate1] = read_ARCascii([TCfiles1dir,TCfiles1(i).name]);%单位mm/D  用作实测
    prerate1 = FillView2Global(H1,prerate1);
    % prerate1 = areaweight_downscale(prerate1, 0.25, 0.5);
    prerate1(isnan(prerate1))=0; prerate1(prerate1<0)=0;              %阈值

    %-------prerate2 PERSIANNCDR
    [~,H2,prerate2] = read_ARCascii([TCfiles2dir,TCfiles2(i).name]);%单位mm/D  PERSIANNCDR 用作估测
    prerate2 = FillView2Global(H2,prerate2);
    % prerate2 = areaweight_downscale(prerate2, 0.25, 0.5);
    prerate2(isnan(prerate2))=0; prerate2(prerate2<0)=0;              %阈值    

    for d=1:28%100:2000km
        if d<=8                               %Land
        ZoneLR=CRing_land(:,:,9-d);          %注意降水转换
        preLR1=prerate1(ZoneLR==10000);         %这个区域内所有降水
        preLR2=prerate2(ZoneLR==10000);         %这个区域内所有降水
        [TC_CC(i,d),TC_RMSE(i,d),TC_ME(i,d),~,TC_BIAS(i,d),~,~,~,~]=CC_RMSE_ME_MAE_BIAS_ABIAS([preLR1,preLR2]);

        else                                  %Sea
        ZoneSR=CRing_sea(:,:,d-8);
        preSR1=prerate1(ZoneSR==10000);
        preSR2=prerate2(ZoneSR==10000);
        [TC_CC(i,d),TC_RMSE(i,d),TC_ME(i,d),~,TC_BIAS(i,d),~,~,~,~]=CC_RMSE_ME_MAE_BIAS_ABIAS([preSR1,preSR2]);

        end
    end
end
toc;
%TC_BIAS(TC_BIAS<=-100000 | TC_BIAS>=100000 | TC_BIAS==0)=missing;
%TC_RMSE(TC_RMSE<=-100000 | TC_RMSE>=100000 | TC_RMSE==0)=missing;
%TC_CC(TC_CC<=-100000 | TC_CC>=100000 | TC_CC==0)=missing;
%TC_ME(TC_ME<=-100000 | TC_ME>=100000 | TC_ME==0)=missing;
%TC_BIAS(TC_BIAS<=-100000 | TC_BIAS>=100000 | TC_BIAS==0)=missing;
%TC_RMSE(TC_RMSE<=-100000 | TC_RMSE>=100000 | TC_RMSE==0)=missing;
%TC_CC(TC_CC<=-100000 | TC_CC>=100000 | TC_CC==0)=missing;
%TC_ME(TC_ME<=-100000 | TC_ME>=100000 | TC_ME==0)=missing;
TC_BIAS(TC_BIAS<=-100000 | TC_BIAS>=100000)=missing;
TC_RMSE(TC_RMSE<=-100000 | TC_RMSE>=100000)=missing;
TC_CC(TC_CC<=-100000 | TC_CC>=100000)=missing;
TC_ME(TC_ME<=-100000 | TC_ME>=100000)=missing;
TC_BIAS(TC_BIAS<=-100000 | TC_BIAS>=100000)=missing;
TC_RMSE(TC_RMSE<=-100000 | TC_RMSE>=100000)=missing;
TC_CC(TC_CC<=-100000 | TC_CC>=100000)=missing;
TC_ME(TC_ME<=-100000 | TC_ME>=100000)=missing;
TC_BIAS_mean=mean(TC_BIAS,"omitmissing");
TC_RMSE_mean=mean(TC_RMSE,"omitmissing");
TC_CC_mean=mean(TC_CC,"omitmissing");
TC_ME_mean=mean(TC_ME,"omitmissing");

OUT_TC=[TC_CC_mean;TC_RMSE_mean;TC_ME_mean;TC_BIAS_mean];
plot(TC_ME_mean)

%% --------------------------------整体数据的误差变化
start_date_vec = datevec('2018.01.01', 'yyyy.mm.dd');
end_date_vec = datevec('2019.12.31', 'yyyy.mm.dd');
num_days = datenum(end_date_vec) - datenum(start_date_vec);% 计算起始和结束日期之间的天数

G_BIAS=zeros(num_days+1,28);
G_RMSE=zeros(num_days+1,28);
G_CC=zeros(num_days+1,28);
G_ME=zeros(num_days+1,28);
cs=0.25;
tic
parfor i = 0:num_days
    i
    % 当前日期
    current_date = addtodate(datenum(start_date_vec), i, 'day');
    c_d_vec = datevec(current_date);
    daynum=c_d_vec(1)*10000 + c_d_vec(2)*100 + c_d_vec(3);
    STRD=num2str(YMD_num(c_d_vec(1),c_d_vec(2),c_d_vec(3)),'%03d');
    %%%prerate1
    if ismember(daynum,[20131230,20131231,20140101,20140102,20140103,20140104])  %这四天有问题
        prerate1 = read_PERCDRbin([Globaldata1dir,'\errordata20131230-20140104\','aB1_d',num2str(mod(c_d_vec(1),100),'%02d'),STRD,'.bin'],cs);
    else
        prerate1 = read_PERSIANNbin([Globaldata1dir,'ms6s4_d',num2str(mod(c_d_vec(1),100),'%02d'),STRD,'.bin'],cs);
    end
    prerate1(isnan(prerate1))=0; prerate1(prerate1<0)=0;
    %%%PRECDR
    prerate2 = read_PERCDRbin([Globaldata2dir,'aB1_d',num2str(mod(c_d_vec(1),100),'%02d'),STRD,'.bin'],cs);
    prerate2(isnan(prerate2))=0; prerate2(prerate2<0)=0;

    for d=1:28%100:2000km
        if d<=8                               %Land
        ZoneLR=CRing_land(:,:,9-d);          %注意降水转换
        preLR1=prerate1(ZoneLR==10000);         %这个区域内所有降水
        preLR2=prerate2(ZoneLR==10000);         %这个区域内所有降水
        [G_CC(i+1,d),G_RMSE(i+1,d),G_ME(i+1,d),~,G_BIAS(i+1,d),~,~,~,~]=CC_RMSE_ME_MAE_BIAS_ABIAS([preLR1,preLR2]);

        else                                  %Sea
        ZoneSR=CRing_sea(:,:,d-8);
        preSR1=prerate1(ZoneSR==10000);
        preSR2=prerate2(ZoneSR==10000);
        [G_CC(i+1,d),G_RMSE(i+1,d),G_ME(i+1,d),~,G_BIAS(i+1,d),~,~,~,~]=CC_RMSE_ME_MAE_BIAS_ABIAS([preSR1,preSR2]);

        end
    end
end
toc
G_BIAS_mean=mean(G_BIAS);
G_RMSE_mean=mean(G_RMSE);
G_CC_mean=mean(G_CC);
G_ME_mean=mean(G_ME);
OUT_G=[G_CC_mean;G_RMSE_mean;G_ME_mean;G_BIAS_mean];
