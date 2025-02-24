%气旋经过网格周围气旋变量的值
clear; clc; close all;
%Initialize Matlab Parallel Computing Enviornment
% p=parpool(2);
%----------------------数据预设---------------------------------------------
%----降水数据
ERASETS={'Horizontal_Divergence','Relative_Humidity','Relative_Vorticity', ...
    'Temperature','U_Component_Of_Wind','V_Component_Of_Wind','Vertical_Velocity'};%和文件名里的一样
%----陆地掩膜
[lmheader1,lmheader2,Clandmask] = read_ARCascii( ...
    'D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%陆地和国家
Clandmask_L=Clandmask;
Clandmask_S=Clandmask;
Clandmask_L(Clandmask>=0)=1;
Clandmask_L(Clandmask==-9999)=0;
Clandmask_S(Clandmask>=0)=0;
Clandmask_S(Clandmask==-9999)=1;
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];                    %-180-180还是0360要确认好

%----气旋数据
TC_idir='D:\DATA\Segmentation map tutorial\IBTrACS.since1980.v04r00.240321.nc';
[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);
load('RMRW_AIRadius01_22.mat');%RMRW_AIRadius

%----参数预设
cs = 0.5;                                                               %改
RMRflag=0;                                                              %改
NOlim=lim_Y2NO([2001,2022],ISO_TIME);                                   %改
NODATA_value=-9999;                                                     %改
FN=50;                                                                  %改
[LonCT_gpu,LatCT_gpu] = GridCenLocat_gpu(cs);
[~,Garea2D] = Gridarea3D(cs);
[~,Garea2D_2] = Gridarea3D(cs/2);%ERA5重新采样用的
%----构建NC存储文件
NC_ECfile='D:\Desktop2\Attri-project\Attri_EC_daily_Data_UnetRMRW.nc';
if ~exist(NC_ECfile)
nccreate(NC_ECfile,'record','Dimensions',{'record',NOlim(3)},'FillValue',NaN);
nccreate(NC_ECfile,'hPa','Dimensions',{'hPa',18},'FillValue',NaN);
nccreate(NC_ECfile,'time','Dimensions',{'time',FN},'FillValue',NaN);
nccreate(NC_ECfile,'o','Dimensions',{'o',5},'FillValue',NaN);
nccreate(NC_ECfile,'HD_mean',         'Dimensions',{'record',NOlim(3),'time',FN,'o',5,'hPa',18},'FillValue',NaN);
nccreate(NC_ECfile,'HD_abs_mean',     'Dimensions',{'record',NOlim(3),'time',FN,'o',5,'hPa',18},'FillValue',NaN);
nccreate(NC_ECfile,'RH_mean',         'Dimensions',{'record',NOlim(3),'time',FN,'o',5,'hPa',18},'FillValue',NaN);
nccreate(NC_ECfile,'RV_mean',         'Dimensions',{'record',NOlim(3),'time',FN,'o',5,'hPa',18},'FillValue',NaN);
nccreate(NC_ECfile,'RV_abs_mean',     'Dimensions',{'record',NOlim(3),'time',FN,'o',5,'hPa',18},'FillValue',NaN);
nccreate(NC_ECfile,'T_mean',          'Dimensions',{'record',NOlim(3),'time',FN,'o',5,'hPa',18},'FillValue',NaN);
nccreate(NC_ECfile,'UVW_mean',        'Dimensions',{'record',NOlim(3),'time',FN,'o',5,'hPa',18},'FillValue',NaN);
nccreate(NC_ECfile,'UVW_abs_mean',    'Dimensions',{'record',NOlim(3),'time',FN,'o',5,'hPa',18},'FillValue',NaN);
nccreate(NC_ECfile,'UVW_mean_Azimuth','Dimensions',{'record',NOlim(3),'time',FN,'o',5,'hPa',18},'FillValue',NaN);
nccreate(NC_ECfile,'VV_mean',         'Dimensions',{'record',NOlim(3),'time',FN,'o',5,'hPa',18},'FillValue',NaN);
nccreate(NC_ECfile,'VV_abs_mean',     'Dimensions',{'record',NOlim(3),'time',FN,'o',5,'hPa',18},'FillValue',NaN);

ncwrite(NC_ECfile,'record',(NOlim(1):NOlim(2)));
ncwrite(NC_ECfile,'hPa',(950:-50:100));
ncwrite(NC_ECfile,'time',(1:FN));
ncwrite(NC_ECfile,'o',(1:5));%1'O' 2'I' 3'OR' 4'OL' 5'OS'
end
%-------------------------------------------------------------------------
for ECDS=1:length(ERASETS)
    
    if ECDS == 6                                                        %UV风速要合在一起算所以要跳过一个 但是UV风一定要是第五个和第六个ERASETS
        continue;
    elseif ECDS == 5                                                    %改 注意这里面的5和6是不是UV WIND
        char(ERASETS(ECDS))
        [ERA5_idir_U,Excel_file] = ERA5_switch_dir(char(ERASETS(ECDS)));
        [ERA5_idir_V] = ERA5_switch_dir(char(ERASETS(ECDS+1)));
    else
        char(ERASETS(ECDS))                                             %改
        [ERA5_idir,Excel_file] = ERA5_switch_dir(char(ERASETS(ECDS)));
    end

    for hpa=950:-50:100
        hpa
        %-----初始化结果矩阵
        %360*3/24 最大45天 平均7天左右  时间分辨率需调整"50"
        %每日: 维度: O: Full without RB; I: inner; OR:OutRing; OL: Land result; 5 OS: Sea Result; 6 RB: Rainbelt; FWRB: Full with RB 
        
        TCD_EC=ones(NOlim(3),FN,5)*NODATA_value;        %每场气旋每日环境条件平均值
        if ECDS == 1 || ECDS == 3 || ECDS == 7
        TCD_ABEC=ones(NOlim(3),FN,5)*NODATA_value;      %每场气旋每日环境条件绝对平均值
        elseif ECDS == 5                                                                             %每场气旋平均风向角
        TCD_W=ones(NOlim(3),FN,5)*NODATA_value; 
        end          

        % 列 1 O Full without RB;2 inner;3 OutRing;4 O Land result; 5 O Sea Result; 6 Rainbelt; 7 Full with RB
        Arr3=repmat(NODATA_value,NOlim(3),1);                                  %改 满足TC_MP_Coor
        TC_EC_O=Arr3;   TC_EC_I=Arr3;   TC_EC_OR=Arr3;   TC_EC_OL=Arr3;   TC_EC_OS=Arr3;             %每场气旋每日环境条件平均值
        if ECDS == 1 || ECDS == 3 || ECDS == 7
        TC_ABEC_O=Arr3; TC_ABEC_I=Arr3; TC_ABEC_OR=Arr3; TC_ABEC_OL=Arr3; TC_ABEC_OS=Arr3;           %每场气旋每日环境条件绝对平均值
        elseif ECDS == 5                                                                             %每场气旋平均风向角
        TC_W_O=Arr3;    TC_W_I=Arr3;    TC_W_OR=Arr3;    TC_W_OL=Arr3;    TC_W_OS=Arr3; 
        end   

        % --------------------开始运算-----------------------------------------------
        % 注意parfor中被赋值的参数索引必须是固定的, 同一数组只存一个结果(同一数组存储时只能出现一种索引方法)
        % 除了循环index: NO 之外不能出现任何数组中内容，否则就不能运行parfor
        for NO=1:NOlim(3)
            if rem(NO,100)==0 
                [char(ERASETS(ECDS)),'-',num2str(hpa),'hpa-NO',num2str(NO)] 
            end
            time=ISO_TIME(:,:,NO + NOlim(1)-1);
            lat=LAT(:,NO + NOlim(1)-1); 
            lat(isnan(lat) )=[];
            lon=LON(:,NO + NOlim(1)-1);
            lon(isnan(lon) )=[]; 
        
            %-------每场气旋(每三小时)的日尺度时间点数值形式get-----------------------
            year   =str2num(time(1:4,1:length(lat))');                          %str2num可以对多行字符串使用变成一个数组
            month  =str2num(time(6:7,1:length(lat))');                          %str2double只能对一个值使用
            day    =str2num(time(9:10,1:length(lat))');
            numtime=year*10^4+month*10^2+day;
        
            %计算加减一天
            % [DB,MB,YB,~,~,~] = daybackdayforward(day(1),month(1),year(1));      %第一天之前的两天
            % [~,~,~,DF,MF,YF] = daybackdayforward(day(end),month(end),year(end));%第一天之后的两天
            % allday=[str2double([YB,MB,DB]);numtime;str2double([YF,MF,DF])];     %需要补充跨年前一天和后一天
            % allday=unique(allday);                                              %unique默认去重之后升序排序
            %不计算加减一天
            allday=unique(numtime);    
            allday_S=num2str(allday);
            allmon=unique(year*10^2+month);
            allmon_S=num2str(allmon);
       
            %-------初始化储存矩阵--------------------------------------------------
            [~,Garea3D] = Gridarea3D(cs,length(allday));                        %维度对应的网格面积矩阵
            Clandmask3D_L=repmat(Clandmask_L,[1 1 length(allday)]);             %维度对应的网格陆地掩膜
            Clandmask3D_S=repmat(Clandmask_S,[1 1 length(allday)]);             %维度对应的网格陆地掩膜
        
            allI_Mask=zeros(180/cs,360/cs,length(allday),'gpuArray');                      %每天内核降水掩膜
            allO_Mask=zeros(180/cs,360/cs,length(allday),'gpuArray');                      %每大半径降水掩膜
            allOR_Mask=zeros(180/cs,360/cs,length(allday),'gpuArray');                     %每天外环(雨带是否计算不一样)降水掩膜
            alldata=zeros(180/cs,360/cs,length(allday),'gpuArray');
    
            %-------本场气旋所涉及降水数据读取----------------------------------------
            if ECDS == 5
                [alldata_U] = ERA5_switch_read(char(ERASETS(ECDS)),ERA5_idir_U,hpa,alldata,allmon,allmon_S,allday_S,month,Garea2D_2,cs);
                [alldata_V] = ERA5_switch_read(char(ERASETS(ECDS+1)),ERA5_idir_V,hpa,alldata,allmon,allmon_S,allday_S,month,Garea2D_2,cs);

            else
                [alldata] = ERA5_switch_read(char(ERASETS(ECDS)),ERA5_idir,hpa,alldata,allmon,allmon_S,allday_S,month,Garea2D_2,cs);
            end
            %-------计算每日掩膜----------------------------------------------------
            for d=1:length(allday)                                              %DayBF不要才能和lat数量对上，实际上是unique numtime，时间顺序应该是不会变的所以第一个和最后一个剪掉就行
                %准备每天需要的数据
                DailyIndex=find(numtime==allday(d));
                if d~=length(allday)                                            %补足每一天21点到24点的掩膜,也就是把第二天0点加进来
                    DailyIndex=[DailyIndex;DailyIndex(end)+1];
                end
                Dlon=lon(DailyIndex);                                            
                Dlat=lat(DailyIndex);

                if ECDS == 5                                                    %这个没意义，就是保证L2B能运行 而且固定半径不需要用到dailydata
                    Dailydata_U=alldata_U(:,:,d);
                    Dailydata_V=alldata_V(:,:,d);
                else
                    Dailydata=alldata(:,:,d);
                end
                
                %计算气旋的固定距离掩膜
                [londen,latden] = TCDenseTrackPoint(Dlon,Dlat,cs);              %Dlon Dlat 必须是气旋行进顺序

                [WZ,~] = TCPoint2RLine(londen,latden,cs);                        
                % Oradius=500;                                                    %这里需要一个控制radius的程序!!!!
                % Iradius=100;
                Iradius=RMRW_AIRadius(NO,d,1);
                Oradius=RMRW_AIRadius(NO,d,2);

                if ECDS == 5
                    Dailydata=Dailydata_U;                                      %这个没意义，就是保证L2B能运行 而且固定半径不需要用到dailydata
                end
                DRMW=[];%这个没意义，就是保证L2B能运行 而且固定半径不需要用到RMW
                [Buff_OWZ,Buff_OMask,Buff_IWZ,Buff_IMask] = TCLine2Buffer_GPU(WZ,Oradius,Iradius,cs,LonCT_gpu,LatCT_gpu,Dailydata,RMRflag,DRMW);
        
                ORing_Mask=Buff_OMask;
                ORing_Mask(Buff_IWZ)=0;
    
                %每天计算降水然后整体平均，并把整体的空间图存出来
                allI_Mask(:,:,d)=Buff_IMask;
                allO_Mask(:,:,d)=Buff_OMask;
                allOR_Mask(:,:,d)=ORing_Mask;
    
            end
            %如果计算BF1D的话需要将处理好的掩膜时间相互赋值使得每一天的掩膜包含BF1D的位置即可）可以设计个ifBF1D的口令
        
        
            % -----------结果提取------------------------------------------------
            %---------------每日-----------
            fillarray=repmat(NODATA_value,1,FN-length(allday));
            MPPG_O=allO_Mask.*Garea3D;
            MPPG_I=allI_Mask.*Garea3D;
            MPPG_OR=allOR_Mask.*Garea3D;
            MPPG_OL=allO_Mask.*Garea3D.*Clandmask3D_L;
            MPPG_OS=allO_Mask.*Garea3D.*Clandmask3D_S;

            if ECDS == 1 || ECDS == 3 || ECDS == 7 % DV：s⁻¹ RV：s⁻¹ VV：Pa/s
                %环境变量异常值都是NaN 也不需要阈值
                Data_Mask=alldata;
                Data_Mask(~isnan(Data_Mask))=1;
                %每个环境条件每天的平均
                TCD_EC(NO,:,1)  = [permute(sum(alldata.*MPPG_O,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_O,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,2)  = [permute(sum(alldata.*MPPG_I,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_I,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,3) = [permute(sum(alldata.*MPPG_OR,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OR,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,4) = [permute(sum(alldata.*MPPG_OL,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OL,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,5) = [permute(sum(alldata.*MPPG_OS,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OS,[1 2]) ,[3,1,2])',fillarray]; 
            
                %每个环境条件的绝对值平均
                TCD_ABEC(NO,:,1)  = [permute(sum(abs(alldata).*MPPG_O,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_O,[1 2]) ,[3,1,2])',fillarray];
                TCD_ABEC(NO,:,2)  = [permute(sum(abs(alldata).*MPPG_I,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_I,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_ABEC(NO,:,3) = [permute(sum(abs(alldata).*MPPG_OR,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OR,[1 2]) ,[3,1,2])',fillarray];
                TCD_ABEC(NO,:,4) = [permute(sum(abs(alldata).*MPPG_OL,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OL,[1 2]) ,[3,1,2])',fillarray];
                TCD_ABEC(NO,:,5) = [permute(sum(abs(alldata).*MPPG_OS,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OS,[1 2]) ,[3,1,2])',fillarray];
            
                %每个环境条件整体的网格平均
                TC_EC_O(NO,1)  = sum(alldata.*MPPG_O,"all","omitmissing") / sum(Data_Mask.*MPPG_O,"all"); 
                TC_EC_I(NO,1)  = sum(alldata.*MPPG_I,"all","omitmissing") / sum(Data_Mask.*MPPG_I,"all"); 
                TC_EC_OR(NO,1) = sum(alldata.*MPPG_OR,"all","omitmissing") / sum(Data_Mask.*MPPG_OR,"all"); 
                TC_EC_OL(NO,1) = sum(alldata.*MPPG_OL,"all","omitmissing") / sum(Data_Mask.*MPPG_OL,"all"); 
                TC_EC_OS(NO,1) = sum(alldata.*MPPG_OS,"all","omitmissing") / sum(Data_Mask.*MPPG_OS,"all"); 
            
                %每个环境条件整体绝对值的网格平均
                TC_ABEC_O(NO,1)  = sum(abs(alldata).*MPPG_O,"all","omitmissing") / sum(Data_Mask.*MPPG_O,"all");
                TC_ABEC_I(NO,1)  = sum(abs(alldata).*MPPG_I,"all","omitmissing") / sum(Data_Mask.*MPPG_I,"all"); 
                TC_ABEC_OR(NO,1) = sum(abs(alldata).*MPPG_OR,"all","omitmissing") / sum(Data_Mask.*MPPG_OR,"all");
                TC_ABEC_OL(NO,1) = sum(abs(alldata).*MPPG_OL,"all","omitmissing") / sum(Data_Mask.*MPPG_OL,"all");
                TC_ABEC_OS(NO,1) = sum(abs(alldata).*MPPG_OS,"all","omitmissing") / sum(Data_Mask.*MPPG_OS,"all");

            elseif ECDS == 2 || ECDS == 4 % T:℃ RH：%
                %环境变量异常值都是NaN 也不需要阈值
                Data_Mask=alldata;
                Data_Mask(~isnan(Data_Mask))=1;
                %每个环境条件每天的平均
                TCD_EC(NO,:,1)  = [permute(sum(alldata.*MPPG_O,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_O,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,2)  = [permute(sum(alldata.*MPPG_I,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_I,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,3) = [permute(sum(alldata.*MPPG_OR,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OR,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,4) = [permute(sum(alldata.*MPPG_OL,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OL,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,5) = [permute(sum(alldata.*MPPG_OS,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OS,[1 2]) ,[3,1,2])',fillarray]; 
              
                %每个环境条件整体的网格平均
                TC_EC_O(NO,1)  = sum(alldata.*MPPG_O,"all","omitmissing") / sum(Data_Mask.*MPPG_O,"all"); 
                TC_EC_I(NO,1)  = sum(alldata.*MPPG_I,"all","omitmissing") / sum(Data_Mask.*MPPG_I,"all"); 
                TC_EC_OR(NO,1) = sum(alldata.*MPPG_OR,"all","omitmissing") / sum(Data_Mask.*MPPG_OR,"all"); 
                TC_EC_OL(NO,1) = sum(alldata.*MPPG_OL,"all","omitmissing") / sum(Data_Mask.*MPPG_OL,"all"); 
                TC_EC_OS(NO,1) = sum(alldata.*MPPG_OS,"all","omitmissing") / sum(Data_Mask.*MPPG_OS,"all"); 
              
            elseif ECDS == 5  % U V wind：m/s
                %环境变量异常值都是NaN 也不需要阈值
                Data_Mask_U=alldata_U;
                Data_Mask_U(~isnan(Data_Mask_U))=1;
                Data_Mask_V=alldata_V;
                Data_Mask_V(~isnan(Data_Mask_V))=1;
                %每个环境条件每天的平均 U
                TCD_EC_O_U = permute(sum(alldata_U.*MPPG_O,[1 2],"omitmissing") ./ sum(Data_Mask_U.*MPPG_O,[1 2]) ,[3,1,2])'; 
                TCD_EC_I_U = permute(sum(alldata_U.*MPPG_I,[1 2],"omitmissing") ./ sum(Data_Mask_U.*MPPG_I,[1 2]) ,[3,1,2])'; 
                TCD_EC_OR_U = permute(sum(alldata_U.*MPPG_OR,[1 2],"omitmissing") ./ sum(Data_Mask_U.*MPPG_OR,[1 2]) ,[3,1,2])'; 
                TCD_EC_OL_U = permute(sum(alldata_U.*MPPG_OL,[1 2],"omitmissing") ./ sum(Data_Mask_U.*MPPG_OL,[1 2]) ,[3,1,2])'; 
                TCD_EC_OS_U = permute(sum(alldata_U.*MPPG_OS,[1 2],"omitmissing") ./ sum(Data_Mask_U.*MPPG_OS,[1 2]) ,[3,1,2])'; 
                %    V
                TCD_EC_O_V = permute(sum(alldata_V.*MPPG_O,[1 2],"omitmissing") ./ sum(Data_Mask_V.*MPPG_O,[1 2]) ,[3,1,2])'; 
                TCD_EC_I_V = permute(sum(alldata_V.*MPPG_I,[1 2],"omitmissing") ./ sum(Data_Mask_V.*MPPG_I,[1 2]) ,[3,1,2])'; 
                TCD_EC_OR_V = permute(sum(alldata_V.*MPPG_OR,[1 2],"omitmissing") ./ sum(Data_Mask_V.*MPPG_OR,[1 2]) ,[3,1,2])'; 
                TCD_EC_OL_V = permute(sum(alldata_V.*MPPG_OL,[1 2],"omitmissing") ./ sum(Data_Mask_V.*MPPG_OL,[1 2]) ,[3,1,2])'; 
                TCD_EC_OS_V = permute(sum(alldata_V.*MPPG_OS,[1 2],"omitmissing") ./ sum(Data_Mask_V.*MPPG_OS,[1 2]) ,[3,1,2])'; 
                % 每个环境条件每天绝对值的平均 U
                TCD_ABEC_O_U = permute(sum(abs(alldata_U).*MPPG_O,[1 2],"omitmissing") ./ sum(Data_Mask_U.*MPPG_O,[1 2]) ,[3,1,2])'; 
                TCD_ABEC_I_U = permute(sum(abs(alldata_U).*MPPG_I,[1 2],"omitmissing") ./ sum(Data_Mask_U.*MPPG_I,[1 2]) ,[3,1,2])'; 
                TCD_ABEC_OR_U = permute(sum(abs(alldata_U).*MPPG_OR,[1 2],"omitmissing") ./ sum(Data_Mask_U.*MPPG_OR,[1 2]) ,[3,1,2])'; 
                TCD_ABEC_OL_U = permute(sum(abs(alldata_U).*MPPG_OL,[1 2],"omitmissing") ./ sum(Data_Mask_U.*MPPG_OL,[1 2]) ,[3,1,2])'; 
                TCD_ABEC_OS_U = permute(sum(abs(alldata_U).*MPPG_OS,[1 2],"omitmissing") ./ sum(Data_Mask_U.*MPPG_OS,[1 2]) ,[3,1,2])'; 
                %    V
                TCD_ABEC_O_V = permute(sum(abs(alldata_V).*MPPG_O,[1 2],"omitmissing") ./ sum(Data_Mask_V.*MPPG_O,[1 2]) ,[3,1,2])'; 
                TCD_ABEC_I_V = permute(sum(abs(alldata_V).*MPPG_I,[1 2],"omitmissing") ./ sum(Data_Mask_V.*MPPG_I,[1 2]) ,[3,1,2])'; 
                TCD_ABEC_OR_V = permute(sum(abs(alldata_V).*MPPG_OR,[1 2],"omitmissing") ./ sum(Data_Mask_V.*MPPG_OR,[1 2]) ,[3,1,2])'; 
                TCD_ABEC_OL_V = permute(sum(abs(alldata_V).*MPPG_OL,[1 2],"omitmissing") ./ sum(Data_Mask_V.*MPPG_OL,[1 2]) ,[3,1,2])'; 
                TCD_ABEC_OS_V = permute(sum(abs(alldata_V).*MPPG_OS,[1 2],"omitmissing") ./ sum(Data_Mask_V.*MPPG_OS,[1 2]) ,[3,1,2])'; 
                %-----------------------每日结果---------------------------------
                %每日平均风方位角 输出时候表头注意 atan2d(U,V) N 0°, E 90°, W -90°, S 180°
                TCD_W(NO,:,1)  = [atan2d(TCD_EC_O_U,TCD_EC_O_V),fillarray]; %U: degree
                TCD_W(NO,:,2)  = [atan2d(TCD_EC_I_U,TCD_EC_I_V),fillarray]; 
                TCD_W(NO,:,3) = [atan2d(TCD_EC_OR_U,TCD_EC_OR_V),fillarray]; 
                TCD_W(NO,:,4) = [atan2d(TCD_EC_OL_U,TCD_EC_OL_V),fillarray]; 
                TCD_W(NO,:,5) = [atan2d(TCD_EC_OS_U,TCD_EC_OS_V),fillarray]; 
            
                %每日平均风的模
                TCD_EC(NO,:,1)  = [sqrt(TCD_EC_O_U.^2+TCD_EC_O_V.^2),fillarray];
                TCD_EC(NO,:,2)  = [sqrt(TCD_EC_I_U.^2+TCD_EC_I_V.^2),fillarray]; 
                TCD_EC(NO,:,3) = [sqrt(TCD_EC_OR_U.^2+TCD_EC_OR_V.^2),fillarray];
                TCD_EC(NO,:,4) = [sqrt(TCD_EC_OL_U.^2+TCD_EC_OL_V.^2),fillarray];
                TCD_EC(NO,:,5) = [sqrt(TCD_EC_OS_U.^2+TCD_EC_OS_V.^2),fillarray];

                %每日平均绝对风的模
                TCD_ABEC(NO,:,1)  = [sqrt(TCD_ABEC_O_U.^2+TCD_ABEC_O_V.^2),fillarray];
                TCD_ABEC(NO,:,2)  = [sqrt(TCD_ABEC_I_U.^2+TCD_ABEC_I_V.^2),fillarray]; 
                TCD_ABEC(NO,:,3) = [sqrt(TCD_ABEC_OR_U.^2+TCD_ABEC_OR_V.^2),fillarray];
                TCD_ABEC(NO,:,4) = [sqrt(TCD_ABEC_OL_U.^2+TCD_ABEC_OL_V.^2),fillarray];
                TCD_ABEC(NO,:,5) = [sqrt(TCD_ABEC_OS_U.^2+TCD_ABEC_OS_V.^2),fillarray];
                %-------------------------每场结果------------------------------------
                %每个环境条件每场的平均 U
                TC_EC_O_U = sum(alldata_U.*MPPG_O,"all","omitmissing") / sum(Data_Mask_U.*MPPG_O,"all"); 
                TC_EC_I_U = sum(alldata_U.*MPPG_I,"all","omitmissing") / sum(Data_Mask_U.*MPPG_I,"all"); 
                TC_EC_OR_U = sum(alldata_U.*MPPG_OR,"all","omitmissing") / sum(Data_Mask_U.*MPPG_OR,"all"); 
                TC_EC_OL_U = sum(alldata_U.*MPPG_OL,"all","omitmissing") / sum(Data_Mask_U.*MPPG_OL,"all"); 
                TC_EC_OS_U = sum(alldata_U.*MPPG_OS,"all","omitmissing") / sum(Data_Mask_U.*MPPG_OS,"all"); 
                %    V
                TC_EC_O_V = sum(alldata_V.*MPPG_O,"all","omitmissing") / sum(Data_Mask_V.*MPPG_O,"all"); 
                TC_EC_I_V = sum(alldata_V.*MPPG_I,"all","omitmissing") / sum(Data_Mask_V.*MPPG_I,"all"); 
                TC_EC_OR_V = sum(alldata_V.*MPPG_OR,"all","omitmissing") / sum(Data_Mask_V.*MPPG_OR,"all"); 
                TC_EC_OL_V = sum(alldata_V.*MPPG_OL,"all","omitmissing") / sum(Data_Mask_V.*MPPG_OL,"all"); 
                TC_EC_OS_V = sum(alldata_V.*MPPG_OS,"all","omitmissing") / sum(Data_Mask_V.*MPPG_OS,"all"); 
                % 每个环境条件每场绝对值的平均 U
                TC_ABEC_O_U = sum(abs(alldata_U).*MPPG_O,"all","omitmissing") / sum(Data_Mask_U.*MPPG_O,"all"); 
                TC_ABEC_I_U = sum(abs(alldata_U).*MPPG_I,"all","omitmissing") / sum(Data_Mask_U.*MPPG_I,"all"); 
                TC_ABEC_OR_U = sum(abs(alldata_U).*MPPG_OR,"all","omitmissing") / sum(Data_Mask_U.*MPPG_OR,"all"); 
                TC_ABEC_OL_U = sum(abs(alldata_U).*MPPG_OL,"all","omitmissing") / sum(Data_Mask_U.*MPPG_OL,"all"); 
                TC_ABEC_OS_U = sum(abs(alldata_U).*MPPG_OS,"all","omitmissing") / sum(Data_Mask_U.*MPPG_OS,"all"); 
                %    V
                TC_ABEC_O_V = sum(abs(alldata_V).*MPPG_O,"all","omitmissing") / sum(Data_Mask_V.*MPPG_O,"all"); 
                TC_ABEC_I_V = sum(abs(alldata_V).*MPPG_I,"all","omitmissing") / sum(Data_Mask_V.*MPPG_I,"all"); 
                TC_ABEC_OR_V = sum(abs(alldata_V).*MPPG_OR,"all","omitmissing") / sum(Data_Mask_V.*MPPG_OR,"all"); 
                TC_ABEC_OL_V = sum(abs(alldata_V).*MPPG_OL,"all","omitmissing") / sum(Data_Mask_V.*MPPG_OL,"all"); 
                TC_ABEC_OS_V = sum(abs(alldata_V).*MPPG_OS,"all","omitmissing") / sum(Data_Mask_V.*MPPG_OS,"all"); 
                %-----------------------每场结果---------------------------------
                %每日平均风方位角 输出时候表头注意 atan2d(U,V) N 0°, E 90°, W -90°, S 180°
                TC_W_O(NO,1)  = atan2d(TC_EC_O_U,TC_EC_O_V); %U: degree
                TC_W_I(NO,1)  = atan2d(TC_EC_I_U,TC_EC_I_V); 
                TC_W_OR(NO,1) = atan2d(TC_EC_OR_U,TC_EC_OR_V); 
                TC_W_OL(NO,1) = atan2d(TC_EC_OL_U,TC_EC_OL_V); 
                TC_W_OS(NO,1) = atan2d(TC_EC_OS_U,TC_EC_OS_V); 
            
                %每日平均风的模
                TC_EC_O(NO,1)  = sqrt(TC_EC_O_U^2+TC_EC_O_V^2);
                TC_EC_I(NO,1)  = sqrt(TC_EC_I_U^2+TC_EC_I_V^2); 
                TC_EC_OR(NO,1) = sqrt(TC_EC_OR_U^2+TC_EC_OR_V^2);
                TC_EC_OL(NO,1) = sqrt(TC_EC_OL_U^2+TC_EC_OL_V^2);
                TC_EC_OS(NO,1) = sqrt(TC_EC_OS_U^2+TC_EC_OS_V^2);

                %每日平均绝对风的模
                TC_ABEC_O(NO,1)  = sqrt(TC_ABEC_O_U^2+TC_ABEC_O_V^2);
                TC_ABEC_I(NO,1)  = sqrt(TC_ABEC_I_U^2+TC_ABEC_I_V^2); 
                TC_ABEC_OR(NO,1) = sqrt(TC_ABEC_OR_U^2+TC_ABEC_OR_V^2);
                TC_ABEC_OL(NO,1) = sqrt(TC_ABEC_OL_U^2+TC_ABEC_OL_V^2);
                TC_ABEC_OS(NO,1) = sqrt(TC_ABEC_OS_U^2+TC_ABEC_OS_V^2);

            end
            % %-----------空间输出---------暂时不要
            % %整个气旋期间的降水强度空间分布  TCPass期间零降水也考虑进去，与上面的整体数值和日数值计算不同
            % Spa_Data=alldata;
            % Spa_Data(allO_Mask==0)=missing;
            % Spa_PR=mean(Spa_Data,3,"omitmissing");
            % %Spa_PAm=sum(Spa_Data.*Garea3D,3,"omitmissing")*10^-6;
            % 
            % %剪裁输出空间降水图
            % ZipView_ARCtxt('cut&save',Spa_PR,odir,[Outname,'_SID',SID(:,NO + NOlim(1)-1)'],cs,LatCT_gpu,LonCT_gpu,NODATA_value,2);
        
        end
              
        %% ------------------结果输出-------------------------------------------
        if ECDS == 7 %vu风速合在一起输出点的话就会提前一个sheet
            sheet= ECDS-1;
        else
            sheet= ECDS;
        end
        %----------------结果输出到excel表格
            OUT_EC=      [TC_EC_O   ;NaN(10,1);   TC_EC_I ;NaN(10,1);   TC_EC_OR ;NaN(10,1);   TC_EC_OL ;NaN(10,1);   TC_EC_OS]; 
            OUT_EC(OUT_EC   ==NODATA_value)=missing;
            xlswrite(Excel_file,OUT_EC    ,1,    [num2abc2(5+(950-hpa)/50 + (sheet-1)*18*3),num2str(3),':',num2abc2(5+(950-hpa)/50 + (sheet-1)*18*3),num2str(3+ NOlim(3)*5+40 -1)]);%改   
        %----------------
        if ECDS == 1 || ECDS == 3 || ECDS == 5 || ECDS == 7
            OUT_ABEC=    [TC_ABEC_O ;NaN(10,1); TC_ABEC_I ;NaN(10,1); TC_ABEC_OR ;NaN(10,1); TC_ABEC_OL ;NaN(10,1); TC_ABEC_OS];
            OUT_ABEC(OUT_ABEC   ==NODATA_value)=missing;
            OUT_ABEC=gather(OUT_ABEC);%有一次gpuarray死活输不进去excel
            %                     一个sheet一个变量    第5列开始 每个气压层 同一个平均方法 之前一个平均方法的总数     
            xlswrite(Excel_file,OUT_ABEC    ,1,   [num2abc2(5+(950-hpa)/50 +18 + (sheet-1)*18*3) , num2str(3),':',num2abc2(5+(950-hpa)/50 +18 + (sheet-1)*18*3),num2str(3+ NOlim(3)*5+40 -1)]);%改
        %----------------
            if ECDS == 5   
                OUT_W=       [TC_W_O    ;NaN(10,1);    TC_W_I ;NaN(10,1);    TC_W_OR ;NaN(10,1);    TC_W_OL ;NaN(10,1);    TC_W_OS];  
                OUT_W(OUT_W   ==NODATA_value)=missing;
                %                     一个sheet一个变量    第5列开始 每个气压层 同一个平均方法 之前一个平均方法的总数     
                xlswrite(Excel_file,OUT_W ,1,         [num2abc2(5+(950-hpa)/50 +18*2 + (sheet-1)*18*3) , num2str(3),':',num2abc2(5+(950-hpa)/50 +18*2 + (sheet-1)*18*3),num2str(3+ NOlim(3)*5+40 -1)]); %改
            end
        end

        if strcmp(char(ERASETS(ECDS)),'Horizontal_Divergence')
            ncwrite(NC_ECfile,'HD_mean'          ,TCD_EC,[1 1 1 (950-hpa)/50+1]);
            ncwrite(NC_ECfile,'HD_abs_mean'      ,TCD_ABEC,[1 1 1 (950-hpa)/50+1]);
        elseif strcmp(char(ERASETS(ECDS)),'Relative_Humidity')
            ncwrite(NC_ECfile,'RH_mean'          ,TCD_EC,[1 1 1 (950-hpa)/50+1]);
        elseif strcmp(char(ERASETS(ECDS)),'Relative_Vorticity')
            ncwrite(NC_ECfile,'RV_mean'          ,TCD_EC,[1 1 1 (950-hpa)/50+1]);
            ncwrite(NC_ECfile,'RV_abs_mean'      ,TCD_ABEC,[1 1 1 (950-hpa)/50+1]);
        elseif strcmp(char(ERASETS(ECDS)),'Temperature')
            ncwrite(NC_ECfile,'T_mean'           ,TCD_EC,[1 1 1 (950-hpa)/50+1]);
        elseif strcmp(char(ERASETS(ECDS)),'U_Component_Of_Wind')
            ncwrite(NC_ECfile,'UVW_mean'         ,TCD_EC,[1 1 1 (950-hpa)/50+1]);
            TCD_ABEC=gather(TCD_ABEC);
            ncwrite(NC_ECfile,'UVW_abs_mean'     ,TCD_ABEC,[1 1 1 (950-hpa)/50+1]);
            ncwrite(NC_ECfile,'UVW_mean_Azimuth' ,TCD_W,[1 1 1 (950-hpa)/50+1]);
        elseif strcmp(char(ERASETS(ECDS)),'Vertical_Velocity')
            ncwrite(NC_ECfile,'VV_mean'          ,TCD_EC,[1 1 1 (950-hpa)/50+1]);
            ncwrite(NC_ECfile,'VV_abs_mean'      ,TCD_ABEC,[1 1 1 (950-hpa)/50+1]);
        end
        
    end
end
%----------------输出表头----------------------------------------
O=char({'_O','_I','_OR','_OL','_OS'});
for j1 =1:5
    % %----------------Horizontal_Divergence
    % THead_HD=[{'年份','编号',['MeanDaily_HD',O(j1,:)]},repmat({NaN},1,18*FN-1),{['ABS_MeanDaily_HD',O(j1,:)]},repmat({NaN},1,18*FN-1) ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18*2)];
    % for j2=3:FN:3+FN*18-1
    %     THead_HD(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    %     THead_HD(2,j2+FN*18)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_HD,1,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18*2-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    % %----------------Relative_Humidity
    % THead_RH=[{'年份','编号',['MeanDaily_RH',O(j1,:)]},repmat({NaN},1,18*FN-1) ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18)];
    % for j2=3:FN:3+FN*18-1
    %     THead_RH(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_RH,2,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    % %----------------Relative_Vorticity
    % THead_RV=[{'年份','编号',['MeanDaily_RV',O(j1,:)]},repmat({NaN},1,18*FN-1),{['ABS_MeanDaily_RV',O(j1,:)]},repmat({NaN},1,18*FN-1) ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18*2)];
    % for j2=3:FN:3+FN*18-1
    %     THead_RV(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    %     THead_RV(2,j2+FN*18)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_RV,3,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18*2-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    % %----------------Temperature
    % THead_T=[{'年份','编号',['MeanDaily_T',O(j1,:)]},repmat({NaN},1,18*FN-1) ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18)];
    % for j2=3:FN:3+FN*18-1
    %     THead_T(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_T,4,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    % %----------------UV_Component_Of_Wind
    % THead_UVW=[{'年份','编号',['MeanDaily_UVW',O(j1,:)]},repmat({NaN},1,18*FN-1),{['ABS_MeanDaily_UVW',O(j1,:)]},repmat({NaN},1,18*FN-1),{['MeanDaily_Azimuth',O(j1,:)]},repmat({NaN},1,18*FN-1) ...
    %     ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18*3)];
    % for j2=3:FN:3+FN*18-1
    %     THead_UVW(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    %     THead_UVW(2,j2+FN*18)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    %     THead_UVW(2,j2+FN*18*2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_UVW,5,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18*3-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    % %----------------Vertical_Velocity
    % THead_VV=[{'年份','编号',['MeanDaily_VV',O(j1,:)]},repmat({NaN},1,18*FN-1),{['ABS_MeanDaily_VV',O(j1,:)]},repmat({NaN},1,18*FN-1) ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18*2)];
    % for j2=3:FN:3+FN*18-1
    %     THead_VV(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    %     THead_VV(2,j2+FN*18)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_VV,6,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18*2-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    %----------------Entire
    THead_Entire=[{'年份','编号',['Mean_HD',O(j1,:),' s⁻¹']},repmat({NaN},1,18-1),{['ABS_Mean_HD',O(j1,:),' s⁻¹']},repmat({NaN},1,18*2-1),{['Mean_RH',O(j1,:),' %']},repmat({NaN},1,18*3-1)...
       ,{['Mean_RV',O(j1,:),' s⁻¹']},repmat({NaN},1,18-1),{['ABS_Mean_RV',O(j1,:),' s⁻¹']},repmat({NaN},1,18*2-1),{['Mean_T',O(j1,:),' ℃']},repmat({NaN},1,18*3-1) ...
       ,{['Mean_UVW',O(j1,:),' m/s']},repmat({NaN},1,18-1),{['ABS_Mean_UVW',O(j1,:),' m/s']},repmat({NaN},1,18-1),{['Mean_WindAzimuth',O(j1,:),' degree N 0°, E 90°, W -90°, S 180°']},repmat({NaN},1,18-1)...
       ,{['Mean_VV',O(j1,:),' Pa/s']},repmat({NaN},1,18-1),{['abs_Mean_VV',O(j1,:),' Pa/s']},repmat({NaN},1,18*2-1)...
       ; repmat({NaN},1,2),repmat(num2cell(950:-50:100),1,3*6)];%为了方便存，每个变量都是按照三种平均3*18来预留的位置，但是只有有数据的地方会命名
    xlswrite(Excel_file,THead_Entire,1,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+18*3*6-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
end

%% -------------------------Test & visualization--------------------------
% imshow(Spa_PR,[],Colormap=CustomColormap)
% P=polybuffer([lon,lat],'lines',5);
% plot(P)
% A=allO_Mask;
% A(allO_Mask==0)=missing;
% A=mean(A,3,"omitmissing");

% MaskO=mean(allO_Mask,3);
% hold on
% [R,C]=find(MaskO>0);
% plot(C,R,"o",'MarkerEdgeColor','r',"MarkerSize",1)

%判断Lon lat是不是都是相同的符号，因为区间并不是-180-180
% FH=zeros(4757,1);
% for j =1:4757
%     lon=LON(:,j);
%     if length(find(lon>180))>0
%         FH(j)= length(find(lon<0));
%     elseif length(find(lon<-180))>0
%         FH(j)= length(find(lon>0));
%     end
% end
% sum(FH)

