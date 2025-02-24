%提取每一年前25%的台风，台风的总降水，最高降水量，每年的分位台风降水量提取   
%update20211102
clear; clc; close all;
%Initialize Matlab Parallel Computing Enviornment
%p=parpool(10);
%------------Done--------
header=['ncols           720';
        'nrows           360';
        'xllcorner      -180';
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value  -9999'];
    
TC_idir='E:\DATA\IBTrACS\NEW\IBTrACS.ALL.v04r00.nc';
[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);

files_idir ='F:\GlobalTCMask1\TCraster500km\';
raster_files = dir([files_idir,'*.txt']);

ERA5_Pre_idir='F:\GlobalTCMask1\single_pre\ERA5\';
odir='F:\GlobalTCMask1\NODEduplicate_sum\ERA5origin\MaxLandPreTop25\';
ERAPre_FN='ERA5_Ori_SingleTC_pre_SID';

Startyear=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA.xlsx',5,'F3:F8173');%开始年
Totalpre=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA.xlsx',5,'AC3:AC8173');%总降水 这个可以决定分位时排序的条件

[lmheader1,lmheader2,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 根据输入数据的extend进行修改

Eachyear_QuanPre_Sequence=zeros(50,71);%每一年的前Q分位降水序号
AVE_top25_pre=zeros(360,720);

Q=0.25;%提取的前多少分位

for y=1950:2020 %改初始年份记得改5305
TCpre=zeros(360,720);
YearPre=Totalpre(Startyear==y);  
%-----计算得到每一年开始的那个编号方便后面找位置
if y==1950
    startnumber=0;
else
    startnumber=length(find(Startyear<y))+1;  
end

%-----先找出前25%的编号  
[SortedPre,WZ]=sort(YearPre,'descend');%对每一列排序这里只有一列
Sequence_Q=WZ(1:round(Q * length(WZ)),1);%降序的原位置序列取前当前分位值的数据（0.25）
Sequence_Q=Sequence_Q + startnumber + 5305;%带回生成原来的序号 从50年开始要加上5305

Eachyear_QuanPre_Sequence(1:length(Sequence_Q),y-1949)=Sequence_Q;%存下每一年的前Q分位降水序号

%-----按编号读取每年的前25%加到一起 求平均
    for i=1:length(Sequence_Q)
        [~,inputheader,currentPre]= read_ARCascii ([ ERA5_Pre_idir , ERAPre_FN , SID(:,Sequence_Q(i))' ,'.txt']);
        currentPre(currentPre==-9999)=0;%注意去掉-9999不然结果会错
        TCpre=TCpre+currentPre;
    end
    
OUTPUT(odir,['ERA5_Ori_LandMax_TCpre_Top25_',num2str(y)],header,TCpre);

if y>=1966
AVE_top25_pre=AVE_top25_pre+TCpre;
end

end
AVE=AVE_top25_pre/(2020-1966+1);
OUTPUT('F:\GlobalTCMask1\NODEduplicate_sum\SUMandAVE\','ERA5_AVE_Ori_LandMax_TCpre_Top25_1966_2020',header,AVE); %好像有问题



