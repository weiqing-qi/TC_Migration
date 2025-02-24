%��������������Χ����������ֵ
clear; clc; close all;
%Initialize Matlab Parallel Computing Enviornment
% p=parpool(2);
%----------------------����Ԥ��---------------------------------------------
%----��ˮ����
ERASETS={'Horizontal_Divergence','Relative_Humidity','Relative_Vorticity', ...
    'Temperature','U_Component_Of_Wind','V_Component_Of_Wind','Vertical_Velocity'};%���ļ������һ��
%----½����Ĥ
[lmheader1,lmheader2,Clandmask] = read_ARCascii( ...
    'D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���
Clandmask_L=Clandmask;
Clandmask_S=Clandmask;
Clandmask_L(Clandmask>=0)=1;
Clandmask_L(Clandmask==-9999)=0;
Clandmask_S(Clandmask>=0)=0;
Clandmask_S(Clandmask==-9999)=1;
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];                    %-180-180����0360Ҫȷ�Ϻ�

%----��������
TC_idir='D:\DATA\Segmentation map tutorial\IBTrACS.since1980.v04r00.240321.nc';
[SID,LAT,LON,ISO_TIME,~,~,~,~,~,~,~,~,~]=IBTrACS_nc_entire_variable_r(TC_idir);
load('RMRW_AIRadius01_22.mat');%RMRW_AIRadius

%----����Ԥ��
cs = 0.5;                                                               %��
RMRflag=0;                                                              %��
NOlim=lim_Y2NO([2001,2022],ISO_TIME);                                   %��
NODATA_value=-9999;                                                     %��
FN=50;                                                                  %��
[LonCT_gpu,LatCT_gpu] = GridCenLocat_gpu(cs);
[~,Garea2D] = Gridarea3D(cs);
[~,Garea2D_2] = Gridarea3D(cs/2);%ERA5���²����õ�
%----����NC�洢�ļ�
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
    
    if ECDS == 6                                                        %UV����Ҫ����һ��������Ҫ����һ�� ����UV��һ��Ҫ�ǵ�����͵�����ERASETS
        continue;
    elseif ECDS == 5                                                    %�� ע���������5��6�ǲ���UV WIND
        char(ERASETS(ECDS))
        [ERA5_idir_U,Excel_file] = ERA5_switch_dir(char(ERASETS(ECDS)));
        [ERA5_idir_V] = ERA5_switch_dir(char(ERASETS(ECDS+1)));
    else
        char(ERASETS(ECDS))                                             %��
        [ERA5_idir,Excel_file] = ERA5_switch_dir(char(ERASETS(ECDS)));
    end

    for hpa=950:-50:100
        hpa
        %-----��ʼ���������
        %360*3/24 ���45�� ƽ��7������  ʱ��ֱ��������"50"
        %ÿ��: ά��: O: Full without RB; I: inner; OR:OutRing; OL: Land result; 5 OS: Sea Result; 6 RB: Rainbelt; FWRB: Full with RB 
        
        TCD_EC=ones(NOlim(3),FN,5)*NODATA_value;        %ÿ������ÿ�ջ�������ƽ��ֵ
        if ECDS == 1 || ECDS == 3 || ECDS == 7
        TCD_ABEC=ones(NOlim(3),FN,5)*NODATA_value;      %ÿ������ÿ�ջ�����������ƽ��ֵ
        elseif ECDS == 5                                                                             %ÿ������ƽ�������
        TCD_W=ones(NOlim(3),FN,5)*NODATA_value; 
        end          

        % �� 1 O Full without RB;2 inner;3 OutRing;4 O Land result; 5 O Sea Result; 6 Rainbelt; 7 Full with RB
        Arr3=repmat(NODATA_value,NOlim(3),1);                                  %�� ����TC_MP_Coor
        TC_EC_O=Arr3;   TC_EC_I=Arr3;   TC_EC_OR=Arr3;   TC_EC_OL=Arr3;   TC_EC_OS=Arr3;             %ÿ������ÿ�ջ�������ƽ��ֵ
        if ECDS == 1 || ECDS == 3 || ECDS == 7
        TC_ABEC_O=Arr3; TC_ABEC_I=Arr3; TC_ABEC_OR=Arr3; TC_ABEC_OL=Arr3; TC_ABEC_OS=Arr3;           %ÿ������ÿ�ջ�����������ƽ��ֵ
        elseif ECDS == 5                                                                             %ÿ������ƽ�������
        TC_W_O=Arr3;    TC_W_I=Arr3;    TC_W_OR=Arr3;    TC_W_OL=Arr3;    TC_W_OS=Arr3; 
        end   

        % --------------------��ʼ����-----------------------------------------------
        % ע��parfor�б���ֵ�Ĳ������������ǹ̶���, ͬһ����ֻ��һ�����(ͬһ����洢ʱֻ�ܳ���һ����������)
        % ����ѭ��index: NO ֮�ⲻ�ܳ����κ����������ݣ�����Ͳ�������parfor
        for NO=1:NOlim(3)
            if rem(NO,100)==0 
                [char(ERASETS(ECDS)),'-',num2str(hpa),'hpa-NO',num2str(NO)] 
            end
            time=ISO_TIME(:,:,NO + NOlim(1)-1);
            lat=LAT(:,NO + NOlim(1)-1); 
            lat(isnan(lat) )=[];
            lon=LON(:,NO + NOlim(1)-1);
            lon(isnan(lon) )=[]; 
        
            %-------ÿ������(ÿ��Сʱ)���ճ߶�ʱ�����ֵ��ʽget-----------------------
            year   =str2num(time(1:4,1:length(lat))');                          %str2num���ԶԶ����ַ���ʹ�ñ��һ������
            month  =str2num(time(6:7,1:length(lat))');                          %str2doubleֻ�ܶ�һ��ֵʹ��
            day    =str2num(time(9:10,1:length(lat))');
            numtime=year*10^4+month*10^2+day;
        
            %����Ӽ�һ��
            % [DB,MB,YB,~,~,~] = daybackdayforward(day(1),month(1),year(1));      %��һ��֮ǰ������
            % [~,~,~,DF,MF,YF] = daybackdayforward(day(end),month(end),year(end));%��һ��֮�������
            % allday=[str2double([YB,MB,DB]);numtime;str2double([YF,MF,DF])];     %��Ҫ�������ǰһ��ͺ�һ��
            % allday=unique(allday);                                              %uniqueĬ��ȥ��֮����������
            %������Ӽ�һ��
            allday=unique(numtime);    
            allday_S=num2str(allday);
            allmon=unique(year*10^2+month);
            allmon_S=num2str(allmon);
       
            %-------��ʼ���������--------------------------------------------------
            [~,Garea3D] = Gridarea3D(cs,length(allday));                        %ά�ȶ�Ӧ�������������
            Clandmask3D_L=repmat(Clandmask_L,[1 1 length(allday)]);             %ά�ȶ�Ӧ������½����Ĥ
            Clandmask3D_S=repmat(Clandmask_S,[1 1 length(allday)]);             %ά�ȶ�Ӧ������½����Ĥ
        
            allI_Mask=zeros(180/cs,360/cs,length(allday),'gpuArray');                      %ÿ���ں˽�ˮ��Ĥ
            allO_Mask=zeros(180/cs,360/cs,length(allday),'gpuArray');                      %ÿ��뾶��ˮ��Ĥ
            allOR_Mask=zeros(180/cs,360/cs,length(allday),'gpuArray');                     %ÿ���⻷(����Ƿ���㲻һ��)��ˮ��Ĥ
            alldata=zeros(180/cs,360/cs,length(allday),'gpuArray');
    
            %-------�����������漰��ˮ���ݶ�ȡ----------------------------------------
            if ECDS == 5
                [alldata_U] = ERA5_switch_read(char(ERASETS(ECDS)),ERA5_idir_U,hpa,alldata,allmon,allmon_S,allday_S,month,Garea2D_2,cs);
                [alldata_V] = ERA5_switch_read(char(ERASETS(ECDS+1)),ERA5_idir_V,hpa,alldata,allmon,allmon_S,allday_S,month,Garea2D_2,cs);

            else
                [alldata] = ERA5_switch_read(char(ERASETS(ECDS)),ERA5_idir,hpa,alldata,allmon,allmon_S,allday_S,month,Garea2D_2,cs);
            end
            %-------����ÿ����Ĥ----------------------------------------------------
            for d=1:length(allday)                                              %DayBF��Ҫ���ܺ�lat�������ϣ�ʵ������unique numtime��ʱ��˳��Ӧ���ǲ��������Ե�һ�������һ����������
                %׼��ÿ����Ҫ������
                DailyIndex=find(numtime==allday(d));
                if d~=length(allday)                                            %����ÿһ��21�㵽24�����Ĥ,Ҳ���ǰѵڶ���0��ӽ���
                    DailyIndex=[DailyIndex;DailyIndex(end)+1];
                end
                Dlon=lon(DailyIndex);                                            
                Dlat=lat(DailyIndex);

                if ECDS == 5                                                    %���û���壬���Ǳ�֤L2B������ ���ҹ̶��뾶����Ҫ�õ�dailydata
                    Dailydata_U=alldata_U(:,:,d);
                    Dailydata_V=alldata_V(:,:,d);
                else
                    Dailydata=alldata(:,:,d);
                end
                
                %���������Ĺ̶�������Ĥ
                [londen,latden] = TCDenseTrackPoint(Dlon,Dlat,cs);              %Dlon Dlat �����������н�˳��

                [WZ,~] = TCPoint2RLine(londen,latden,cs);                        
                % Oradius=500;                                                    %������Ҫһ������radius�ĳ���!!!!
                % Iradius=100;
                Iradius=RMRW_AIRadius(NO,d,1);
                Oradius=RMRW_AIRadius(NO,d,2);

                if ECDS == 5
                    Dailydata=Dailydata_U;                                      %���û���壬���Ǳ�֤L2B������ ���ҹ̶��뾶����Ҫ�õ�dailydata
                end
                DRMW=[];%���û���壬���Ǳ�֤L2B������ ���ҹ̶��뾶����Ҫ�õ�RMW
                [Buff_OWZ,Buff_OMask,Buff_IWZ,Buff_IMask] = TCLine2Buffer_GPU(WZ,Oradius,Iradius,cs,LonCT_gpu,LatCT_gpu,Dailydata,RMRflag,DRMW);
        
                ORing_Mask=Buff_OMask;
                ORing_Mask(Buff_IWZ)=0;
    
                %ÿ����㽵ˮȻ������ƽ������������Ŀռ�ͼ�����
                allI_Mask(:,:,d)=Buff_IMask;
                allO_Mask(:,:,d)=Buff_OMask;
                allOR_Mask(:,:,d)=ORing_Mask;
    
            end
            %�������BF1D�Ļ���Ҫ������õ���Ĥʱ���໥��ֵʹ��ÿһ�����Ĥ����BF1D��λ�ü��ɣ�������Ƹ�ifBF1D�Ŀ���
        
        
            % -----------�����ȡ------------------------------------------------
            %---------------ÿ��-----------
            fillarray=repmat(NODATA_value,1,FN-length(allday));
            MPPG_O=allO_Mask.*Garea3D;
            MPPG_I=allI_Mask.*Garea3D;
            MPPG_OR=allOR_Mask.*Garea3D;
            MPPG_OL=allO_Mask.*Garea3D.*Clandmask3D_L;
            MPPG_OS=allO_Mask.*Garea3D.*Clandmask3D_S;

            if ECDS == 1 || ECDS == 3 || ECDS == 7 % DV��s�6�3�0�1 RV��s�6�3�0�1 VV��Pa/s
                %���������쳣ֵ����NaN Ҳ����Ҫ��ֵ
                Data_Mask=alldata;
                Data_Mask(~isnan(Data_Mask))=1;
                %ÿ����������ÿ���ƽ��
                TCD_EC(NO,:,1)  = [permute(sum(alldata.*MPPG_O,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_O,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,2)  = [permute(sum(alldata.*MPPG_I,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_I,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,3) = [permute(sum(alldata.*MPPG_OR,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OR,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,4) = [permute(sum(alldata.*MPPG_OL,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OL,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,5) = [permute(sum(alldata.*MPPG_OS,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OS,[1 2]) ,[3,1,2])',fillarray]; 
            
                %ÿ�����������ľ���ֵƽ��
                TCD_ABEC(NO,:,1)  = [permute(sum(abs(alldata).*MPPG_O,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_O,[1 2]) ,[3,1,2])',fillarray];
                TCD_ABEC(NO,:,2)  = [permute(sum(abs(alldata).*MPPG_I,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_I,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_ABEC(NO,:,3) = [permute(sum(abs(alldata).*MPPG_OR,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OR,[1 2]) ,[3,1,2])',fillarray];
                TCD_ABEC(NO,:,4) = [permute(sum(abs(alldata).*MPPG_OL,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OL,[1 2]) ,[3,1,2])',fillarray];
                TCD_ABEC(NO,:,5) = [permute(sum(abs(alldata).*MPPG_OS,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OS,[1 2]) ,[3,1,2])',fillarray];
            
                %ÿ�������������������ƽ��
                TC_EC_O(NO,1)  = sum(alldata.*MPPG_O,"all","omitmissing") / sum(Data_Mask.*MPPG_O,"all"); 
                TC_EC_I(NO,1)  = sum(alldata.*MPPG_I,"all","omitmissing") / sum(Data_Mask.*MPPG_I,"all"); 
                TC_EC_OR(NO,1) = sum(alldata.*MPPG_OR,"all","omitmissing") / sum(Data_Mask.*MPPG_OR,"all"); 
                TC_EC_OL(NO,1) = sum(alldata.*MPPG_OL,"all","omitmissing") / sum(Data_Mask.*MPPG_OL,"all"); 
                TC_EC_OS(NO,1) = sum(alldata.*MPPG_OS,"all","omitmissing") / sum(Data_Mask.*MPPG_OS,"all"); 
            
                %ÿ�����������������ֵ������ƽ��
                TC_ABEC_O(NO,1)  = sum(abs(alldata).*MPPG_O,"all","omitmissing") / sum(Data_Mask.*MPPG_O,"all");
                TC_ABEC_I(NO,1)  = sum(abs(alldata).*MPPG_I,"all","omitmissing") / sum(Data_Mask.*MPPG_I,"all"); 
                TC_ABEC_OR(NO,1) = sum(abs(alldata).*MPPG_OR,"all","omitmissing") / sum(Data_Mask.*MPPG_OR,"all");
                TC_ABEC_OL(NO,1) = sum(abs(alldata).*MPPG_OL,"all","omitmissing") / sum(Data_Mask.*MPPG_OL,"all");
                TC_ABEC_OS(NO,1) = sum(abs(alldata).*MPPG_OS,"all","omitmissing") / sum(Data_Mask.*MPPG_OS,"all");

            elseif ECDS == 2 || ECDS == 4 % T:�� RH��%
                %���������쳣ֵ����NaN Ҳ����Ҫ��ֵ
                Data_Mask=alldata;
                Data_Mask(~isnan(Data_Mask))=1;
                %ÿ����������ÿ���ƽ��
                TCD_EC(NO,:,1)  = [permute(sum(alldata.*MPPG_O,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_O,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,2)  = [permute(sum(alldata.*MPPG_I,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_I,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,3) = [permute(sum(alldata.*MPPG_OR,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OR,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,4) = [permute(sum(alldata.*MPPG_OL,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OL,[1 2]) ,[3,1,2])',fillarray]; 
                TCD_EC(NO,:,5) = [permute(sum(alldata.*MPPG_OS,[1 2],"omitmissing") ./ sum(Data_Mask.*MPPG_OS,[1 2]) ,[3,1,2])',fillarray]; 
              
                %ÿ�������������������ƽ��
                TC_EC_O(NO,1)  = sum(alldata.*MPPG_O,"all","omitmissing") / sum(Data_Mask.*MPPG_O,"all"); 
                TC_EC_I(NO,1)  = sum(alldata.*MPPG_I,"all","omitmissing") / sum(Data_Mask.*MPPG_I,"all"); 
                TC_EC_OR(NO,1) = sum(alldata.*MPPG_OR,"all","omitmissing") / sum(Data_Mask.*MPPG_OR,"all"); 
                TC_EC_OL(NO,1) = sum(alldata.*MPPG_OL,"all","omitmissing") / sum(Data_Mask.*MPPG_OL,"all"); 
                TC_EC_OS(NO,1) = sum(alldata.*MPPG_OS,"all","omitmissing") / sum(Data_Mask.*MPPG_OS,"all"); 
              
            elseif ECDS == 5  % U V wind��m/s
                %���������쳣ֵ����NaN Ҳ����Ҫ��ֵ
                Data_Mask_U=alldata_U;
                Data_Mask_U(~isnan(Data_Mask_U))=1;
                Data_Mask_V=alldata_V;
                Data_Mask_V(~isnan(Data_Mask_V))=1;
                %ÿ����������ÿ���ƽ�� U
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
                % ÿ����������ÿ�����ֵ��ƽ�� U
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
                %-----------------------ÿ�ս��---------------------------------
                %ÿ��ƽ���緽λ�� ���ʱ���ͷע�� atan2d(U,V) N 0��, E 90��, W -90��, S 180��
                TCD_W(NO,:,1)  = [atan2d(TCD_EC_O_U,TCD_EC_O_V),fillarray]; %U: degree
                TCD_W(NO,:,2)  = [atan2d(TCD_EC_I_U,TCD_EC_I_V),fillarray]; 
                TCD_W(NO,:,3) = [atan2d(TCD_EC_OR_U,TCD_EC_OR_V),fillarray]; 
                TCD_W(NO,:,4) = [atan2d(TCD_EC_OL_U,TCD_EC_OL_V),fillarray]; 
                TCD_W(NO,:,5) = [atan2d(TCD_EC_OS_U,TCD_EC_OS_V),fillarray]; 
            
                %ÿ��ƽ�����ģ
                TCD_EC(NO,:,1)  = [sqrt(TCD_EC_O_U.^2+TCD_EC_O_V.^2),fillarray];
                TCD_EC(NO,:,2)  = [sqrt(TCD_EC_I_U.^2+TCD_EC_I_V.^2),fillarray]; 
                TCD_EC(NO,:,3) = [sqrt(TCD_EC_OR_U.^2+TCD_EC_OR_V.^2),fillarray];
                TCD_EC(NO,:,4) = [sqrt(TCD_EC_OL_U.^2+TCD_EC_OL_V.^2),fillarray];
                TCD_EC(NO,:,5) = [sqrt(TCD_EC_OS_U.^2+TCD_EC_OS_V.^2),fillarray];

                %ÿ��ƽ�����Է��ģ
                TCD_ABEC(NO,:,1)  = [sqrt(TCD_ABEC_O_U.^2+TCD_ABEC_O_V.^2),fillarray];
                TCD_ABEC(NO,:,2)  = [sqrt(TCD_ABEC_I_U.^2+TCD_ABEC_I_V.^2),fillarray]; 
                TCD_ABEC(NO,:,3) = [sqrt(TCD_ABEC_OR_U.^2+TCD_ABEC_OR_V.^2),fillarray];
                TCD_ABEC(NO,:,4) = [sqrt(TCD_ABEC_OL_U.^2+TCD_ABEC_OL_V.^2),fillarray];
                TCD_ABEC(NO,:,5) = [sqrt(TCD_ABEC_OS_U.^2+TCD_ABEC_OS_V.^2),fillarray];
                %-------------------------ÿ�����------------------------------------
                %ÿ����������ÿ����ƽ�� U
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
                % ÿ����������ÿ������ֵ��ƽ�� U
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
                %-----------------------ÿ�����---------------------------------
                %ÿ��ƽ���緽λ�� ���ʱ���ͷע�� atan2d(U,V) N 0��, E 90��, W -90��, S 180��
                TC_W_O(NO,1)  = atan2d(TC_EC_O_U,TC_EC_O_V); %U: degree
                TC_W_I(NO,1)  = atan2d(TC_EC_I_U,TC_EC_I_V); 
                TC_W_OR(NO,1) = atan2d(TC_EC_OR_U,TC_EC_OR_V); 
                TC_W_OL(NO,1) = atan2d(TC_EC_OL_U,TC_EC_OL_V); 
                TC_W_OS(NO,1) = atan2d(TC_EC_OS_U,TC_EC_OS_V); 
            
                %ÿ��ƽ�����ģ
                TC_EC_O(NO,1)  = sqrt(TC_EC_O_U^2+TC_EC_O_V^2);
                TC_EC_I(NO,1)  = sqrt(TC_EC_I_U^2+TC_EC_I_V^2); 
                TC_EC_OR(NO,1) = sqrt(TC_EC_OR_U^2+TC_EC_OR_V^2);
                TC_EC_OL(NO,1) = sqrt(TC_EC_OL_U^2+TC_EC_OL_V^2);
                TC_EC_OS(NO,1) = sqrt(TC_EC_OS_U^2+TC_EC_OS_V^2);

                %ÿ��ƽ�����Է��ģ
                TC_ABEC_O(NO,1)  = sqrt(TC_ABEC_O_U^2+TC_ABEC_O_V^2);
                TC_ABEC_I(NO,1)  = sqrt(TC_ABEC_I_U^2+TC_ABEC_I_V^2); 
                TC_ABEC_OR(NO,1) = sqrt(TC_ABEC_OR_U^2+TC_ABEC_OR_V^2);
                TC_ABEC_OL(NO,1) = sqrt(TC_ABEC_OL_U^2+TC_ABEC_OL_V^2);
                TC_ABEC_OS(NO,1) = sqrt(TC_ABEC_OS_U^2+TC_ABEC_OS_V^2);

            end
            % %-----------�ռ����---------��ʱ��Ҫ
            % %���������ڼ�Ľ�ˮǿ�ȿռ�ֲ�  TCPass�ڼ��㽵ˮҲ���ǽ�ȥ���������������ֵ������ֵ���㲻ͬ
            % Spa_Data=alldata;
            % Spa_Data(allO_Mask==0)=missing;
            % Spa_PR=mean(Spa_Data,3,"omitmissing");
            % %Spa_PAm=sum(Spa_Data.*Garea3D,3,"omitmissing")*10^-6;
            % 
            % %��������ռ併ˮͼ
            % ZipView_ARCtxt('cut&save',Spa_PR,odir,[Outname,'_SID',SID(:,NO + NOlim(1)-1)'],cs,LatCT_gpu,LonCT_gpu,NODATA_value,2);
        
        end
              
        %% ------------------������-------------------------------------------
        if ECDS == 7 %vu���ٺ���һ�������Ļ��ͻ���ǰһ��sheet
            sheet= ECDS-1;
        else
            sheet= ECDS;
        end
        %----------------��������excel���
            OUT_EC=      [TC_EC_O   ;NaN(10,1);   TC_EC_I ;NaN(10,1);   TC_EC_OR ;NaN(10,1);   TC_EC_OL ;NaN(10,1);   TC_EC_OS]; 
            OUT_EC(OUT_EC   ==NODATA_value)=missing;
            xlswrite(Excel_file,OUT_EC    ,1,    [num2abc2(5+(950-hpa)/50 + (sheet-1)*18*3),num2str(3),':',num2abc2(5+(950-hpa)/50 + (sheet-1)*18*3),num2str(3+ NOlim(3)*5+40 -1)]);%��   
        %----------------
        if ECDS == 1 || ECDS == 3 || ECDS == 5 || ECDS == 7
            OUT_ABEC=    [TC_ABEC_O ;NaN(10,1); TC_ABEC_I ;NaN(10,1); TC_ABEC_OR ;NaN(10,1); TC_ABEC_OL ;NaN(10,1); TC_ABEC_OS];
            OUT_ABEC(OUT_ABEC   ==NODATA_value)=missing;
            OUT_ABEC=gather(OUT_ABEC);%��һ��gpuarray�����䲻��ȥexcel
            %                     һ��sheetһ������    ��5�п�ʼ ÿ����ѹ�� ͬһ��ƽ������ ֮ǰһ��ƽ������������     
            xlswrite(Excel_file,OUT_ABEC    ,1,   [num2abc2(5+(950-hpa)/50 +18 + (sheet-1)*18*3) , num2str(3),':',num2abc2(5+(950-hpa)/50 +18 + (sheet-1)*18*3),num2str(3+ NOlim(3)*5+40 -1)]);%��
        %----------------
            if ECDS == 5   
                OUT_W=       [TC_W_O    ;NaN(10,1);    TC_W_I ;NaN(10,1);    TC_W_OR ;NaN(10,1);    TC_W_OL ;NaN(10,1);    TC_W_OS];  
                OUT_W(OUT_W   ==NODATA_value)=missing;
                %                     һ��sheetһ������    ��5�п�ʼ ÿ����ѹ�� ͬһ��ƽ������ ֮ǰһ��ƽ������������     
                xlswrite(Excel_file,OUT_W ,1,         [num2abc2(5+(950-hpa)/50 +18*2 + (sheet-1)*18*3) , num2str(3),':',num2abc2(5+(950-hpa)/50 +18*2 + (sheet-1)*18*3),num2str(3+ NOlim(3)*5+40 -1)]); %��
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
%----------------�����ͷ----------------------------------------
O=char({'_O','_I','_OR','_OL','_OS'});
for j1 =1:5
    % %----------------Horizontal_Divergence
    % THead_HD=[{'���','���',['MeanDaily_HD',O(j1,:)]},repmat({NaN},1,18*FN-1),{['ABS_MeanDaily_HD',O(j1,:)]},repmat({NaN},1,18*FN-1) ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18*2)];
    % for j2=3:FN:3+FN*18-1
    %     THead_HD(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    %     THead_HD(2,j2+FN*18)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_HD,1,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18*2-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    % %----------------Relative_Humidity
    % THead_RH=[{'���','���',['MeanDaily_RH',O(j1,:)]},repmat({NaN},1,18*FN-1) ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18)];
    % for j2=3:FN:3+FN*18-1
    %     THead_RH(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_RH,2,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    % %----------------Relative_Vorticity
    % THead_RV=[{'���','���',['MeanDaily_RV',O(j1,:)]},repmat({NaN},1,18*FN-1),{['ABS_MeanDaily_RV',O(j1,:)]},repmat({NaN},1,18*FN-1) ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18*2)];
    % for j2=3:FN:3+FN*18-1
    %     THead_RV(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    %     THead_RV(2,j2+FN*18)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_RV,3,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18*2-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    % %----------------Temperature
    % THead_T=[{'���','���',['MeanDaily_T',O(j1,:)]},repmat({NaN},1,18*FN-1) ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18)];
    % for j2=3:FN:3+FN*18-1
    %     THead_T(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_T,4,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    % %----------------UV_Component_Of_Wind
    % THead_UVW=[{'���','���',['MeanDaily_UVW',O(j1,:)]},repmat({NaN},1,18*FN-1),{['ABS_MeanDaily_UVW',O(j1,:)]},repmat({NaN},1,18*FN-1),{['MeanDaily_Azimuth',O(j1,:)]},repmat({NaN},1,18*FN-1) ...
    %     ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18*3)];
    % for j2=3:FN:3+FN*18-1
    %     THead_UVW(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    %     THead_UVW(2,j2+FN*18)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    %     THead_UVW(2,j2+FN*18*2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_UVW,5,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18*3-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    % %----------------Vertical_Velocity
    % THead_VV=[{'���','���',['MeanDaily_VV',O(j1,:)]},repmat({NaN},1,18*FN-1),{['ABS_MeanDaily_VV',O(j1,:)]},repmat({NaN},1,18*FN-1) ; repmat({NaN},1,2),repmat(num2cell(1:FN),1,18*2)];
    % for j2=3:FN:3+FN*18-1
    %     THead_VV(2,j2)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    %     THead_VV(2,j2+FN*18)={['hPa:',num2str(950-(j2-3)/FN*50),' Day:1']};
    % end
    % xlswrite(Excel_file,THead_VV,6,[num2abc2(3),num2str(1+(j1-1)*(NOlim(3)+10)),':',num2abc2(5+FN*18*2-1),num2str(2+(j1-1)*(NOlim(3)+10))]);
    %----------------Entire
    THead_Entire=[{'���','���',['Mean_HD',O(j1,:),' s�6�3�0�1']},repmat({NaN},1,18-1),{['ABS_Mean_HD',O(j1,:),' s�6�3�0�1']},repmat({NaN},1,18*2-1),{['Mean_RH',O(j1,:),' %']},repmat({NaN},1,18*3-1)...
       ,{['Mean_RV',O(j1,:),' s�6�3�0�1']},repmat({NaN},1,18-1),{['ABS_Mean_RV',O(j1,:),' s�6�3�0�1']},repmat({NaN},1,18*2-1),{['Mean_T',O(j1,:),' ��']},repmat({NaN},1,18*3-1) ...
       ,{['Mean_UVW',O(j1,:),' m/s']},repmat({NaN},1,18-1),{['ABS_Mean_UVW',O(j1,:),' m/s']},repmat({NaN},1,18-1),{['Mean_WindAzimuth',O(j1,:),' degree N 0��, E 90��, W -90��, S 180��']},repmat({NaN},1,18-1)...
       ,{['Mean_VV',O(j1,:),' Pa/s']},repmat({NaN},1,18-1),{['abs_Mean_VV',O(j1,:),' Pa/s']},repmat({NaN},1,18*2-1)...
       ; repmat({NaN},1,2),repmat(num2cell(950:-50:100),1,3*6)];%Ϊ�˷���棬ÿ���������ǰ�������ƽ��3*18��Ԥ����λ�ã�����ֻ�������ݵĵط�������
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

%�ж�Lon lat�ǲ��Ƕ�����ͬ�ķ��ţ���Ϊ���䲢����-180-180
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

