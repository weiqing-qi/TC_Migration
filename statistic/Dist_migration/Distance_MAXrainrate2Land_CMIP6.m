%计算盆地里，格点到栅格化陆地岸线的距离
allmodel_name={'BCC-CSM2-MR','CAS-FGOALS-f3-L','CMCC-ESM2','CNRM-CM6-1-HR','EC-Earth3', ...
            'EC-Earth3-AerChem','EC-Earth3-CC','EC-Earth3-Veg','NOAA-GFDL-CM4','SNU-SAM0-UNICON'};
ALL_Results = NaN(5000,5*length(allmodel_name));
for MOD = 1:length(allmodel_name)
    model_name = allmodel_name(MOD)

    if ismember(model_name, {''})%跳过的模型
        continue
    end

    if ismember(model_name,{'BCC-CSM2-MR','CAS-FGOALS-f3-L','CMCC-ESM2', 'NOAA-GFDL-CM4', 'SNU-SAM0-UNICON'})
        cs = 1;                                                                 
    elseif ismember(model_name,{'EC-Earth3','EC-Earth3-AerChem','EC-Earth3-CC','EC-Earth3-Veg' })
        cs = 180/256;                                                           
    elseif ismember(model_name,{'CNRM-CM6-1-HR'})
        cs = 0.5;                                                               
    end
    %------------提取海岸线0。1栅格
    [lmheader1,lmheader2,Clandmask] = read_ARCascii( ...
        'D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\gadm36_0fidmasks01_180.txt');%陆地和国家
    Clandmask = LandM01_resample(Clandmask, 0.1, cs);
    Clandmask(Clandmask~=-9999)=0;
    Clandmask(Clandmask==-9999)=100;
    Landedge=bwperim(Clandmask,8);         %For numeric input, any nonzero pixels are considered to be 1 (true)
    Landedge=double(Landedge);
    Landedge(:,1)=0;Landedge(:,360/cs)=0;
    Landedge(1,:)=0;Landedge(180/cs,:)=0;     %注意去掉周围的海岸线；
    %------------每一场计算最大降水率到海岸的距离
    %输入想计算的初始值，最大降水率位置已经有了
    YLONLAT=xlsread(['D:\DATA\CMIP6\',model_name{:},'\',model_name{:},'_TCAVE_PR_Day_500km_Ummd_100500.xlsx'],'Entire','AR3:AT5000');
    YEAR = YLONLAT(:,1);
    LON = YLONLAT(:,2);
    LAT = YLONLAT(:,3);
    
    % %-------------可以直接用TXT计算的模块
    % %进行改编，利用已经算出来的TXT数据直接运行
    % Preindir='D:\DATA\TC_spatial_data\pre_amount\IMERG\'; %改
    % Prefile=dir([Preindir,'*.txt']);
    % % Prefile=Prefile(1955:end);
    % % Prefile=Prefile(1955:end);%改1980开始1955ERA5 2898JRA55
    % 
    % LON=zeros(length(Prefile),1);%每场台风的最大降水强度---都是位置 
    % LAT=zeros(length(Prefile),1);%每场台风的最大日降水强度
    % [LonCenter,LatCenter] = GridCenterLocation(0.25);%这个要与下面H（5)相等 %改
    % tic;
    % parfor f=1:length(Prefile)
    %     % if mod(f,100)==0 
    %     %     disp(f) 
    %     % end
    %     [~,H,TCpre] = read_ARCascii([Preindir,Prefile(f).name]);
    %     [Full_V] = FillView2Global(H,TCpre);
    %     if max(Full_V,[],"all")<=0 %没有降水的情况
    %         LON(f)=missing;
    %         LAT(f)=missing;
    %         continue
    %     end
    %     WZ=find(Full_V==max(Full_V,[],"all"))
    %     LON(f)=LonCenter(WZ(1));%有可能找到多个
    %     LAT(f)=LatCenter(WZ(1));
    %     % [~,I1] = sort(reshape(TCpre,[],1),'descend');%变成一列之后顺序还是不变的
    %     % I1(I1<0.1)=[];
    %     % LON(f)=mean(LonCenter(I1(1:round(length(I1)*0.25))));%前百分之75 但是没用
    %     % LAT(f)=mean(LatCenter(I1(1:round(length(I1)*0.25))));
    % end
    % toc;
    
    
    coast=find(Landedge==1);%岸线网格里的位置索引
    
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
        if Clandmask(R,C)==0%注意上面Clandmask赋的值
          Dist(i)=(-1)*Dist(i);%在陆地上距离为负值
        end
    end
    ALL_Results(1:length(LON),MOD*5-4:MOD*5)=[YEAR,Dist,LON,LAT,NaN(length(LON),1)];
    %imshowpair(Clandmask,Landedge,'montage')
end
