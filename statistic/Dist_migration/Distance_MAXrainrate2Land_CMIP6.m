%����������㵽դ��½�ذ��ߵľ���
allmodel_name={'BCC-CSM2-MR','CAS-FGOALS-f3-L','CMCC-ESM2','CNRM-CM6-1-HR','EC-Earth3', ...
            'EC-Earth3-AerChem','EC-Earth3-CC','EC-Earth3-Veg','NOAA-GFDL-CM4','SNU-SAM0-UNICON'};
ALL_Results = NaN(5000,5*length(allmodel_name));
for MOD = 1:length(allmodel_name)
    model_name = allmodel_name(MOD)

    if ismember(model_name, {''})%������ģ��
        continue
    end

    if ismember(model_name,{'BCC-CSM2-MR','CAS-FGOALS-f3-L','CMCC-ESM2', 'NOAA-GFDL-CM4', 'SNU-SAM0-UNICON'})
        cs = 1;                                                                 
    elseif ismember(model_name,{'EC-Earth3','EC-Earth3-AerChem','EC-Earth3-CC','EC-Earth3-Veg' })
        cs = 180/256;                                                           
    elseif ismember(model_name,{'CNRM-CM6-1-HR'})
        cs = 0.5;                                                               
    end
    %------------��ȡ������0��1դ��
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
    YLONLAT=xlsread(['D:\DATA\CMIP6\',model_name{:},'\',model_name{:},'_TCAVE_PR_Day_500km_Ummd_100500.xlsx'],'Entire','AR3:AT5000');
    YEAR = YLONLAT(:,1);
    LON = YLONLAT(:,2);
    LAT = YLONLAT(:,3);
    
    % %-------------����ֱ����TXT�����ģ��
    % %���иı࣬�����Ѿ��������TXT����ֱ������
    % Preindir='D:\DATA\TC_spatial_data\pre_amount\IMERG\'; %��
    % Prefile=dir([Preindir,'*.txt']);
    % % Prefile=Prefile(1955:end);
    % % Prefile=Prefile(1955:end);%��1980��ʼ1955ERA5 2898JRA55
    % 
    % LON=zeros(length(Prefile),1);%ÿ��̨������ˮǿ��---����λ�� 
    % LAT=zeros(length(Prefile),1);%ÿ��̨�������ս�ˮǿ��
    % [LonCenter,LatCenter] = GridCenterLocation(0.25);%���Ҫ������H��5)��� %��
    % tic;
    % parfor f=1:length(Prefile)
    %     % if mod(f,100)==0 
    %     %     disp(f) 
    %     % end
    %     [~,H,TCpre] = read_ARCascii([Preindir,Prefile(f).name]);
    %     [Full_V] = FillView2Global(H,TCpre);
    %     if max(Full_V,[],"all")<=0 %û�н�ˮ�����
    %         LON(f)=missing;
    %         LAT(f)=missing;
    %         continue
    %     end
    %     WZ=find(Full_V==max(Full_V,[],"all"))
    %     LON(f)=LonCenter(WZ(1));%�п����ҵ����
    %     LAT(f)=LatCenter(WZ(1));
    %     % [~,I1] = sort(reshape(TCpre,[],1),'descend');%���һ��֮��˳���ǲ����
    %     % I1(I1<0.1)=[];
    %     % LON(f)=mean(LonCenter(I1(1:round(length(I1)*0.25))));%ǰ�ٷ�֮75 ����û��
    %     % LAT(f)=mean(LatCenter(I1(1:round(length(I1)*0.25))));
    % end
    % toc;
    
    
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
    ALL_Results(1:length(LON),MOD*5-4:MOD*5)=[YEAR,Dist,LON,LAT,NaN(length(LON),1)];
    %imshowpair(Clandmask,Landedge,'montage')
end
