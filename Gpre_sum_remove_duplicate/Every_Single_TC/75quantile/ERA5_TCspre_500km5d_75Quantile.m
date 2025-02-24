%��ȡÿһ��ǰ25%��̨�磬̨����ܽ�ˮ����߽�ˮ����ÿ��ķ�λ̨�罵ˮ����ȡ   
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

Startyear=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA.xlsx',5,'F3:F8173');%��ʼ��
Totalpre=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA.xlsx',5,'AC3:AC8173');%�ܽ�ˮ ������Ծ�����λʱ���������

[lmheader1,lmheader2,Clandmask] = read_ARCascii('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 �����������ݵ�extend�����޸�

Eachyear_QuanPre_Sequence=zeros(50,71);%ÿһ���ǰQ��λ��ˮ���
AVE_top25_pre=zeros(360,720);

Q=0.25;%��ȡ��ǰ���ٷ�λ

for y=1950:2020 %�ĳ�ʼ��ݼǵø�5305
TCpre=zeros(360,720);
YearPre=Totalpre(Startyear==y);  
%-----����õ�ÿһ�꿪ʼ���Ǹ���ŷ��������λ��
if y==1950
    startnumber=0;
else
    startnumber=length(find(Startyear<y))+1;  
end

%-----���ҳ�ǰ25%�ı��  
[SortedPre,WZ]=sort(YearPre,'descend');%��ÿһ����������ֻ��һ��
Sequence_Q=WZ(1:round(Q * length(WZ)),1);%�����ԭλ������ȡǰ��ǰ��λֵ�����ݣ�0.25��
Sequence_Q=Sequence_Q + startnumber + 5305;%��������ԭ������� ��50�꿪ʼҪ����5305

Eachyear_QuanPre_Sequence(1:length(Sequence_Q),y-1949)=Sequence_Q;%����ÿһ���ǰQ��λ��ˮ���

%-----����Ŷ�ȡÿ���ǰ25%�ӵ�һ�� ��ƽ��
    for i=1:length(Sequence_Q)
        [~,inputheader,currentPre]= read_ARCascii ([ ERA5_Pre_idir , ERAPre_FN , SID(:,Sequence_Q(i))' ,'.txt']);
        currentPre(currentPre==-9999)=0;%ע��ȥ��-9999��Ȼ������
        TCpre=TCpre+currentPre;
    end
    
OUTPUT(odir,['ERA5_Ori_LandMax_TCpre_Top25_',num2str(y)],header,TCpre);

if y>=1966
AVE_top25_pre=AVE_top25_pre+TCpre;
end

end
AVE=AVE_top25_pre/(2020-1966+1);
OUTPUT('F:\GlobalTCMask1\NODEduplicate_sum\SUMandAVE\','ERA5_AVE_Ori_LandMax_TCpre_Top25_1966_2020',header,AVE); %����������



