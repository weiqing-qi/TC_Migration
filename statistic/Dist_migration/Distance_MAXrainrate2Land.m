%����������㵽դ��½�ذ��ߵľ���
%------------��ȡ������0��1դ��
cs=0.5;
[lmheader1,lmheader2,Clandmask] = read_ARCascii( ...
    'D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\gadm36_0fidmasks01_180.txt');%½�غ͹���
Clandmask = LandM01_resample(Clandmask, 0.1, cs);
Clandmask(Clandmask~=-9999)=0;
Clandmask(Clandmask==-9999)=100;
Landedge=bwperim(Clandmask,8);         %For numeric input, any nonzero pixels are considered to be 1 (true)
Landedge=double(Landedge);
Landedge(:,1)=0;Landedge(:,360/cs)=0;
Landedge(1,:)=0;Landedge(180/cs,:)=0;     %ע��ȥ����Χ�ĺ����ߣ�
%------------ÿһ���������ˮ�ʵ������ľ���
%���������ĳ�ʼֵ�����ˮ��λ���Ѿ�����
% LON=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted.xlsx',1,'AU3:AU4348');
% LAT=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted.xlsx',1,'AV3:AV4348');

%-------------����ֱ����TXT�����ģ��
%���иı࣬�����Ѿ��������TXT����ֱ������
Preindir='D:\DATA\TC_spatial_data\pre_amount\IMERG\'; %��
Prefile=dir([Preindir,'*.txt']);
% Prefile=Prefile(1955:end);
% Prefile=Prefile(1955:end);%��1980��ʼ1955ERA5 2898JRA55

LON=zeros(length(Prefile),1);%ÿ��̨������ˮǿ��---����λ�� 
LAT=zeros(length(Prefile),1);%ÿ��̨�������ս�ˮǿ��
[LonCenter,LatCenter] = GridCenterLocation(0.25);%���Ҫ������H��5)��� %��
tic;
parfor f=1:length(Prefile)
    % if mod(f,100)==0 
    %     disp(f) 
    % end
    [~,H,TCpre] = read_ARCascii([Preindir,Prefile(f).name]);
    [Full_V] = FillView2Global(H,TCpre);
    if max(Full_V,[],"all")<=0 %û�н�ˮ�����
        LON(f)=missing;
        LAT(f)=missing;
        continue
    end
    WZ=find(Full_V==max(Full_V,[],"all"))
    LON(f)=LonCenter(WZ(1));%�п����ҵ����
    LAT(f)=LatCenter(WZ(1));
    % [~,I1] = sort(reshape(TCpre,[],1),'descend');%���һ��֮��˳���ǲ����
    % I1(I1<0.1)=[];
    % LON(f)=mean(LonCenter(I1(1:round(length(I1)*0.25))));%ǰ�ٷ�֮75 ����û��
    % LAT(f)=mean(LatCenter(I1(1:round(length(I1)*0.25))));
end
toc;


coast=find(Landedge==1);%�����������λ������

D=zeros(length(coast),1,'gpuArray');
Dist=zeros(length(LAT),1);
CLON=zeros(length(coast),1,'gpuArray'); %ALL Grid size initialization
CLAT=zeros(length(coast),1,'gpuArray');


[LonCT_gpu,LatCT_gpu] = GridCenterLocation_GPU(cs);

for i=1:length(LAT)
    if isnan(LON(i))
        Dist(i)=missing;
        continue
    end
    CLON(:)=LON(i); 
    CLAT(:)=LAT(i); 
    D=SphereDist2_Matrix([CLON,CLAT],[LonCT_gpu(coast),LatCT_gpu(coast)]);
    Dist(i)=min(D);
    [R,C]=latlon2grid(LON(i),LAT(i),cs);
    if Clandmask(R,C)==0%ע������Clandmask����ֵ
      Dist(i)=(-1)*Dist(i);%��½���Ͼ���Ϊ��ֵ
    end
end
out=[Dist,LON,LAT];
%imshowpair(Clandmask,Landedge,'montage')

