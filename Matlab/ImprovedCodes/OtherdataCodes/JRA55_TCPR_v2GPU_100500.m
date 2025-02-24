%�������������(-1&+1day)ƽ����ˮ���ܽ�ˮ������Josh Unet AI ��Ĥ����
clear; clc; close all;
%Initialize Matlab Parallel Computing Enviornment
% p=parpool(2);
%----------------------����Ԥ��---------------------------------------------
%----��ˮ����

odir='D:\DATA\TC_spatial_data\pre_amount\JRA55\';
Outname='JRA55_TC_PreAmount_Day05_500km_Um3';
%Outname='JRA55_TC_AVE_PreRate_Day05_500km_Ummd';
% Excel_file='D:\Desktop2\Attri-project\Results\Otherdatasets\JRA55_100500.xlsx';
ErrorDate=load('E:\DATA\doc\JRA-55\JRA55_TCerror_DATE_NA.txt');
ERA5_idir='E:\DATA\1.Reanalysis\EAR5\ERA5-Totalprecip\';
JRA55_idir='E:\DATA\1.Reanalysis\JRA55\JRA55Daily05\';
JRAfilename='JRA55_Daily_Bilinear_Ori05_TotalPre';
ERAfilename='ERA5_Total_precipitation_on_single_levels_daily';
%AI_Mask_dir= 'D:\DATA\AI_IRcloud_TC_mask_daily.nc4';
%outDtxtName='IMERGFV7B_TC_DailyPR_Unet_mmd_240321IBtrack_44NS_USAW34_1Index';%��������ÿ�����Ĥ
%odir2='D:\DATA\TC_spatial_data\IMERG_FV7B_PRE01_Unet_DailyView\';

%----��������
TC_idir='E:\DATA\5.TCdata\IBTrACS\NEW\IBTrACS.ALL.v04r00.nc';
[SID,LAT,LON,ISO_TIME,USA_WIND,USA_RMW,BASIN]=IBTrACS_nc_entire_variable_r(TC_idir);

%----����Ԥ��
cs = 0.5;                                                              %��
thres=2.4;                                                              %�� 0.1 mm/h 2.4 mm/d ��̬��ֵ���а�
RMRflag=0;                                                              %��
NOlim=lim_Y2NO([1980,2020],ISO_TIME);                                   %��
NODATA_value=-9999;                                                     %��
[LonCT,LatCT] = GridCenterLocation_GPU(cs);
[~,Garea2D] = Gridarea3D(cs);
[~,Garea025] = Gridarea3D(0.25);
%----½����Ĥ
NAmask = Basinmasks(0.5);  %EPmask&NAmask:6 
[lmheader1,lmheader2,Clandmask] = read_ARCascii( ...
    'D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\gadm36_0fidmasks01_180.txt');%½�غ͹���
Clandmask = LandM01_resample(Clandmask, 0.1, cs);
Clandmask_L=Clandmask;
Clandmask_S=Clandmask;
Clandmask_L(Clandmask>=0)=1;
Clandmask_L(Clandmask==-9999)=0;
Clandmask_S(Clandmask>=0)=0;
Clandmask_S(Clandmask==-9999)=1;
% Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];                  %-180-180����0360Ҫȷ�Ϻ�

%-----��ʼ���������
%360*3/24 ���45�� ƽ��7������  ʱ��ֱ��������"50"
%ÿ��: ά��: O: Full without RB; I: inner; OR:OutRing; OL: Land result; 5 OS: Sea Result; 6 RB: Rainbelt; FWRB: Full with RB 
ColuScal=50;                                                            %��

Arr1=repmat(NODATA_value,NOlim(3),ColuScal);
TCD_PR_O=Arr1;   TCD_PR_I=Arr1;   TCD_PR_OR=Arr1;   TCD_PR_OL=Arr1;   TCD_PR_OS=Arr1;        %ÿ������ÿ�ս�ˮǿ�� 
TCD_PAm_O=Arr1;  TCD_PAm_I=Arr1;  TCD_PAm_OR=Arr1;  TCD_PAm_OL=Arr1;  TCD_PAm_OS=Arr1;       %ÿ������ÿ�ս�ˮ��
TCD_MaxP_O=Arr1; TCD_MaxP_I=Arr1; TCD_MaxP_OR=Arr1; TCD_MaxP_OL=Arr1; TCD_MaxP_OS=Arr1;      %ÿ������ÿ�����ˮǿ��
TCD_PA_O=Arr1;   TCD_PA_I=Arr1;   TCD_PA_OR=Arr1;   TCD_PA_OL=Arr1;   TCD_PA_OS=Arr1;        %ÿ������ÿ�ս�ˮ��� 

%ÿ�����ˮǿ�Ⱦ�γ�� ���������lon x coluScal; lat x coluScal
Arr2=repmat(NODATA_value,NOlim(3),ColuScal*2);
TCD_MPCo_O=Arr2; TCD_MPCo_I=Arr2; TCD_MPCo_OR=Arr2; TCD_MPCo_OL=Arr2; TCD_MPCo_OS=Arr2;  %ÿ������ÿ�����ˮǿ�Ⱦ�γ��

% �� 1 O Full without RB;2 inner;3 OutRing;4 O Land result; 5 O Sea Result; 6 Rainbelt; 7 Full with RB
Arr3=repmat(NODATA_value,NOlim(3),10);                                  %�� ����TC_MP_Coor
TC_PR=Arr3; TC_PAm=Arr3; TC_Max_P=Arr3; TC_PA=Arr3; TC_MP_Coor=Arr3;    %һ��lonһ��lat����

Arr4=repmat(NODATA_value,NOlim(3),1); 
TCR_PA_O=Arr4; TCR_PA_L=Arr4; TCR_PA_S=Arr4;%��ʵӰ�����

% --------------------��ʼ����-----------------------------------------------
% ע��parfor�б���ֵ�Ĳ������������ǹ̶���, ͬһ����ֻ��һ�����(ͬһ����洢ʱֻ�ܳ���һ����������)
% ����ѭ��index: NO ֮�ⲻ�ܳ����κ����������ݣ�����Ͳ�������parfor
ST=1;
EN=NOlim(3);
% tic
for NO=ST:EN
    NO
    time=ISO_TIME(:,:,NO + NOlim(1)-1);
    lat=LAT(:,NO + NOlim(1)-1); 
    lat(isnan(lat) )=[];
    lon=LON(:,NO + NOlim(1)-1);
    lon(isnan(lon) )=[]; 
    USA_WIND_1NO=USA_WIND(:,NO + NOlim(1)-1);
    USA_WIND_1NO(isnan(lat))=[]; 
    USA_RMW_1NO=USA_RMW(:,NO + NOlim(1)-1);
    USA_RMW_1NO(isnan(lat))=[];%���������latע��
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
    allyear_S=allday_S(:,1:4);
    allyear=str2num(allyear_S);
    allmon=unique(year*10^2+month);
    allmon_S=num2str(allmon);
    %-------��ʼ���������--------------------------------------------------
    [~,Garea3D] = Gridarea3D(cs,length(allday));                        %ά�ȶ�Ӧ�������������
    Clandmask3D_L=repmat(Clandmask_L,[1 1 length(allday)]);             %ά�ȶ�Ӧ������½����Ĥ
    Clandmask3D_S=repmat(Clandmask_S,[1 1 length(allday)]);             %ά�ȶ�Ӧ������½����Ĥ

    allI_Mask=zeros(180/cs,360/cs,length(allday),'gpuArray');                      %ÿ���ں˽�ˮ��Ĥ
    allO_Mask=zeros(180/cs,360/cs,length(allday),'gpuArray');                      %ÿ��뾶��ˮ��Ĥ
    allOR_Mask=zeros(180/cs,360/cs,length(allday),'gpuArray');                     %ÿ���⻷(����Ƿ���㲻һ��)��ˮ��Ĥ
    alldata=zeros(180/cs,360/cs,length(allday),'gpuArray');                        %bwconncomp��֧��GPU����

    %-------�����������漰��ˮ���ݶ�ȡ--------------------------------
    for k=1:length(allday)
        iserror = ismember(allday(k),ErrorDate);
        STRdate = num2str(allday(k));
        if iserror == 1
            Amonthdata = ncread([ERA5_idir,ERAfilename,STRdate(1:6),'.nc'],'tprate');%ÿһ��allday��Ӧһ��allpre
            patch=ERA5dailyresample( flipud(Amonthdata(:,:,str2double(STRdate(7:8)))')*1000*86400 ,Garea025); %��ȡ������ERA5ԭ���ݾ����ز�������׼����
            errorpre=imread([JRA55_idir,JRAfilename,STRdate,'.tif'])/8;%��λmm/8d
            errorpre(NAmask==6)=patch(NAmask==6); %����ERA5����JRA55
            alldata(:,:,k)=errorpre;
        else 
            alldata(:,:,k)=imread([JRA55_idir,JRAfilename,STRdate,'.tif'])/8;%��λmm/8d
        end
    end
    alldata(alldata<0)=0;
    alldata(isnan(alldata))=0;

    %-------����ÿ����Ĥ----------------------------------------------------
    for d=1:length(allday)                                              %DayBF��Ҫ���ܺ�lat�������ϣ�ʵ������unique numtime��ʱ��˳��Ӧ���ǲ��������Ե�һ�������һ����������
        %׼��ÿ����Ҫ������

        DailyIndex=find(numtime==allday(d));                            %������day��Ϊ���ڳ���30���������ʧЧ
        if d~=length(allday)                                            %����ÿһ��21�㵽24�����Ĥ,Ҳ���ǰѵڶ���0��ӽ���
            DailyIndex=[DailyIndex;DailyIndex(end)+1];
        end
        Dlon=lon(DailyIndex);                                            
        Dlat=lat(DailyIndex);
        DRMW=USA_RMW_1NO(DailyIndex);
        DRMW(isnan(DRMW))=[];
        DWIND=USA_WIND_1NO(DailyIndex);
        DWIND(isnan(DWIND))=[];
        Dailydata=alldata(:,:,d);
        %���������Ĺ̶�������Ĥ
        [londen,latden] = TCDenseTrackPoint(Dlon,Dlat,cs);              %Dlon Dlat �����������н�˳��

        [WZ,~] = TCPoint2RLine(londen,latden,cs);                        
        Oradius=500;                                                   
        Iradius=100;
        [Buff_OWZ,Buff_OMask,Buff_IWZ,Buff_IMask,RMR] = TCLine2Buffer_GPU(WZ,Oradius,Iradius,cs,LonCT,LatCT,Dailydata,RMRflag,DRMW);
        
        %[ORing_Mask,Buff_OMask,~] = OuterRainBelt(thres,Garea2D,Dailydata,Buff_OMask,Buff_OWZ,Buff_IWZ,cs);
                
        ORing_Mask=Buff_OMask;
        ORing_Mask(Buff_IWZ)=0;

        allI_Mask(:,:,d)=Buff_IMask;
        allO_Mask(:,:,d)=Buff_OMask;
        allOR_Mask(:,:,d)=ORing_Mask;

        %--------��usawind>34,44NS,ʱ���������������������--------
        % if ~isempty(DWIND)
        %     if max(DWIND)>34 && max(Dlat)<=44 && min(Dlat)>=-44
        %         DOUT=AI_OMask.*Dailydata;
        %         if sum(DOUT,"all","omitmissing")~=0  %���Ϊ�վͲ�Ҫ�����
        %             DOUT(DOUT==0)=missing;
        %             ZipView_ARCtxt('cut&save',DOUT,odir2, ...
        %                 [outDtxtName,num2str(NO+NOlim(1)-1,'%04d'),'D',num2str(d,'%02d'),'_SID',SID(:,NO + NOlim(1)-1)'], ...
        %                 cs,LatCT_gpu,LonCT_gpu,0,2);
        %         end
        %     end
        % end
    end
    %�������BF1D�Ļ���Ҫ������õ���Ĥʱ���໥��ֵʹ��ÿһ�����Ĥ����BF1D��λ�ü��ɣ�������Ƹ�ifBF1D�Ŀ���


    % -----------�����ȡ-------------------------------------------------
    alldata(alldata<thres)=0;
    Data_Mask=alldata;
    Data_Mask(Data_Mask>0)=1;
    %---------------ÿ��-----------
    % fillarray=repmat(NODATA_value,1,ColuScal-length(allday));
    % 
    % %ÿ����������ƽ��ÿ�ս�ˮǿ�� MM/D
    % TCD_PR_O(NO,:)  = [permute(sum(alldata.*allO_Mask.*Garea3D,[1 2]) ./ sum(Data_Mask.*allO_Mask.*Garea3D,[1 2]) ,[3,1,2])',fillarray]; 
    % TCD_PR_I(NO,:)  = [permute(sum(alldata.*allI_Mask.*Garea3D,[1 2]) ./ sum(Data_Mask.*allI_Mask.*Garea3D,[1 2]) ,[3,1,2])',fillarray]; 
    % TCD_PR_OR(NO,:) = [permute(sum(alldata.*allOR_Mask.*Garea3D,[1 2]) ./ sum(Data_Mask.*allOR_Mask.*Garea3D,[1 2]) ,[3,1,2])',fillarray]; 
    % TCD_PR_OL(NO,:) = [permute(sum(alldata.*allO_Mask.*Garea3D.*Clandmask3D_L,[1 2]) ./ sum(Data_Mask.*allO_Mask.*Garea3D.*Clandmask3D_L,[1 2]) ,[3,1,2])',fillarray]; 
    % TCD_PR_OS(NO,:) = [permute(sum(alldata.*allO_Mask.*Garea3D.*Clandmask3D_S,[1 2]) ./ sum(Data_Mask.*allO_Mask.*Garea3D.*Clandmask3D_S,[1 2]) ,[3,1,2])',fillarray]; 
    % 
    % %ÿ������ÿ�ս�ˮ�� 10^(-3+6-9) 10��m�0�6=����ǧ��
    % TCD_PAm_O(NO,:)  = [permute(sum(alldata.*allO_Mask.*Garea3D,[1 2]) *10^(-6) ,[3,1,2])',fillarray];
    % TCD_PAm_I(NO,:)  = [permute(sum(alldata.*allI_Mask.*Garea3D,[1 2]) *10^(-6) ,[3,1,2])',fillarray]; 
    % TCD_PAm_OR(NO,:) = [permute(sum(alldata.*allOR_Mask.*Garea3D,[1 2]) *10^(-6) ,[3,1,2])',fillarray];
    % TCD_PAm_OL(NO,:) = [permute(sum(alldata.*allO_Mask.*Garea3D.*Clandmask3D_L,[1 2]) *10^(-6) ,[3,1,2])',fillarray];
    % TCD_PAm_OS(NO,:) = [permute(sum(alldata.*allO_Mask.*Garea3D.*Clandmask3D_S,[1 2]) *10^(-6) ,[3,1,2])',fillarray];
    % 
    % %ÿ������ÿ�����ˮǿ�� MM/D �Լ�alldata��������������
    % [M_MPO,I_MPO]  = max(alldata.*allO_Mask,[],[1 2],"linear"); 
    % [M_MPI,I_MPI]  = max(alldata.*allI_Mask,[],[1 2],"linear");
    % [M_MPOR,I_MPOR]= max(alldata.*allOR_Mask,[],[1 2],"linear");
    % [M_MPOL,I_MPOL]= max(alldata.*allO_Mask.*Clandmask3D_L,[],[1 2],"linear");
    % [M_MPOS,I_MPOS]= max(alldata.*allO_Mask.*Clandmask3D_S,[],[1 2],"linear");
    % 
    % TCD_MaxP_O(NO,:)  = [permute(M_MPO , [3,1,2])',fillarray]; 
    % TCD_MaxP_I(NO,:)  = [permute(M_MPI , [3,1,2])',fillarray];
    % TCD_MaxP_OR(NO,:) = [permute(M_MPOR, [3,1,2])',fillarray];
    % TCD_MaxP_OL(NO,:) = [permute(M_MPOL, [3,1,2])',fillarray];
    % TCD_MaxP_OS(NO,:) = [permute(M_MPOS, [3,1,2])',fillarray];
    % 
    % TCD_MPCo_O(NO,:)  = [LonCT(rem(permute(I_MPO , [3,1,2])',180/cs*360/cs)),fillarray, ...
    %     LatCT(rem(permute(I_MPO , [3,1,2])',180/cs*360/cs)),fillarray]; %�ȼۣ�LonCT_gpu3D(permute(I_MPO , [3,1,2]));
    % TCD_MPCo_I(NO,:)  = [LonCT(rem(permute(I_MPI , [3,1,2])',180/cs*360/cs)),fillarray, ...
    %     LatCT(rem(permute(I_MPI , [3,1,2])',180/cs*360/cs)),fillarray];
    % TCD_MPCo_OR(NO,:) = [LonCT(rem(permute(I_MPOR, [3,1,2])',180/cs*360/cs)),fillarray, ...
    %     LatCT(rem(permute(I_MPOR, [3,1,2])',180/cs*360/cs)),fillarray];
    % TCD_MPCo_OL(NO,:) = [LonCT(rem(permute(I_MPOL, [3,1,2])',180/cs*360/cs)),fillarray, ...
    %     LatCT(rem(permute(I_MPOL, [3,1,2])',180/cs*360/cs)),fillarray];
    % TCD_MPCo_OS(NO,:) = [LonCT(rem(permute(I_MPOS, [3,1,2])',180/cs*360/cs)),fillarray, ...
    %     LatCT(rem(permute(I_MPOS, [3,1,2])',180/cs*360/cs)),fillarray];
    % 
    % %ÿ������ÿ�ս�ˮ���
    % TCD_PA_O(NO,:)  = [permute(sum(Data_Mask.*allO_Mask.*Garea3D,[1 2]) ,[3,1,2])',fillarray];
    % TCD_PA_I(NO,:)  = [permute(sum(Data_Mask.*allI_Mask.*Garea3D,[1 2]) ,[3,1,2])',fillarray];
    % TCD_PA_OR(NO,:) = [permute(sum(Data_Mask.*allOR_Mask.*Garea3D,[1 2]) ,[3,1,2])',fillarray];
    % TCD_PA_OL(NO,:) = [permute(sum(Data_Mask.*allO_Mask.*Garea3D.*Clandmask3D_L,[1 2]) ,[3,1,2])',fillarray];
    % TCD_PA_OS(NO,:) = [permute(sum(Data_Mask.*allO_Mask.*Garea3D.*Clandmask3D_S,[1 2]) ,[3,1,2])',fillarray];

    %ÿ��������ʵӰ��������ֿ�O I OR��ûɶ���壬��ͬ����϶��кܶ��غϵ����
    PA=Data_Mask.*allO_Mask;
    PA(PA==0)=missing;
    PA=mean(PA,3,"omitmissing");
    TCR_PA_O(NO,1)=sum(PA.*Garea2D,"all","omitmissing");
    TCR_PA_L(NO,1)=sum(PA.*Clandmask_L.*Garea2D,"all","omitmissing");
    TCR_PA_S(NO,1)=sum(PA.*Clandmask_S.*Garea2D,"all","omitmissing");

    %ÿ������ÿ�ս�ˮǿ�����˼·��for���쵫��parfor�����
    % TCD_PR_O(NO,:)  = [TCD_PAm_O(NO,1:length(allday)) *10^(6) ./ TCD_PA_O(NO,1:length(allday)),fillarray]; 

    %-----------�ռ����---------
    %���������ڼ�Ľ�ˮǿ�ȿռ�ֲ�  TCPass�ڼ��㽵ˮҲ���ǽ�ȥ���������������ֵ������ֵ���㲻ͬ
    Spa_Data=alldata;
    Spa_Data(allO_Mask==0)=missing;
    allO_Mask_FLAT=sum(allO_Mask,3);
    %Spa_PR=mean(Spa_Data,3,"omitmissing");
    Spa_PAm=sum(Spa_Data.*Garea3D,3,"omitmissing");%��λ��Ҫ��10e-6���ܱ��billion������
    Spa_PAm(allO_Mask_FLAT==0)=missing;%�и�����Ϊȱʧ���ݻᵼ�¼������

    %��������ռ併ˮͼ
    %ZipView_ARCtxt('cut&save',Spa_PR,odir,[Outname,'_SID',SID(:,NO + NOlim(1)-1)'],cs,LatCT,LonCT,NODATA_value,2);
    ZipView_ARCtxt('cut&save',Spa_PAm,odir,[Outname,'_SID',SID(:,NO + NOlim(1)-1)'],cs,LatCT,LonCT,NODATA_value,2);

end
% toc
TCD_MPCo_OL(TCD_MPCo_OL== 90-0.5*cs | TCD_MPCo_OL== -180+0.5*cs)=missing;
TCD_MPCo_OS(TCD_MPCo_OS== 90-0.5*cs | TCD_MPCo_OS== -180+0.5*cs)=missing;%���û�г�����½������ô��Ĥȫ����ͻ��ҵ���һ����

%ȫ���ڣ���Ϊ��parfor����̫Ҫ�趨̫��������ֲ��������еĽ�����˷�ʱ��, �ʻ���parfor����ŵ������������ͳһ����
%------------ȫ���ڽ������--------------    
for i = ST:EN
    time=ISO_TIME(:,:,i + NOlim(1)-1);
    lat=LAT(:,i + NOlim(1)-1); 
    lat(isnan(lat) )=[];
    %-------ÿ������(ÿ��Сʱ)���ճ߶�ʱ�����ֵ��ʽget-----------------------
    year   =str2num(time(1:4,1:length(lat))');                          %str2num���ԶԶ����ַ���ʹ�ñ��һ������
    month  =str2num(time(6:7,1:length(lat))');                          %str2doubleֻ�ܶ�һ��ֵʹ��
    day    =str2num(time(9:10,1:length(lat))');
    numtime=year*10^4+month*10^2+day;
    allday=unique(numtime);
    %ÿ��������ˮ��
    TC_PAm(i,1) = sum(TCD_PAm_O(i,1:length(allday)));
    TC_PAm(i,2) = sum(TCD_PAm_I(i,1:length(allday)));
    TC_PAm(i,3) = sum(TCD_PAm_OR(i,1:length(allday)));
    TC_PAm(i,4) = sum(TCD_PAm_OL(i,1:length(allday)));
    TC_PAm(i,5) = sum(TCD_PAm_OS(i,1:length(allday)));
    %ÿ���������ˮǿ��
    [TC_Max_P(i,1),EI_MPO]  = max(TCD_MaxP_O(i,1:length(allday)));
    [TC_Max_P(i,2),EI_MPI]  = max(TCD_MaxP_I(i,1:length(allday)));
    [TC_Max_P(i,3),EI_MPOR] = max(TCD_MaxP_OR(i,1:length(allday)));
    [TC_Max_P(i,4),EI_MPOL] = max(TCD_MaxP_OL(i,1:length(allday)));
    [TC_Max_P(i,5),EI_MPOS] = max(TCD_MaxP_OS(i,1:length(allday)));
    
    TC_MP_Coor(i,1) = TCD_MPCo_O(i,EI_MPO);%һ��lonһ��lat����
    TC_MP_Coor(i,3) = TCD_MPCo_I(i,EI_MPI);
    TC_MP_Coor(i,5) = TCD_MPCo_OR(i,EI_MPOR);
    TC_MP_Coor(i,7) = TCD_MPCo_OL(i,EI_MPOL);
    TC_MP_Coor(i,9) = TCD_MPCo_OS(i,EI_MPOS);
    
    TC_MP_Coor(i,2) = TCD_MPCo_O(i,EI_MPO   +ColuScal);%�պ���coluScal*2
    TC_MP_Coor(i,4) = TCD_MPCo_I(i,EI_MPI   +ColuScal);
    TC_MP_Coor(i,6) = TCD_MPCo_OR(i,EI_MPOR +ColuScal);
    TC_MP_Coor(i,8) = TCD_MPCo_OL(i,EI_MPOL +ColuScal);
    TC_MP_Coor(i,10)= TCD_MPCo_OS(i,EI_MPOS +ColuScal);

    % ÿ��������ˮ��� ע�����ﲻ��������Ӱ�����
    % ֻΪ�˼��㽵ˮǿ��ÿ������ÿ�죬��ˮ�������϶�û����
    % �����е������������е�����������ֱ�ӳ���Ӱ������������ǿ��Ӧ���Ǵ�
    % ����Ϊ�������൱��ÿ�����񶼿�������ͬʱ���Ľ�ˮ
    % �����㲻�� �����������TC_PR��Ҫ��
    TC_PA(i,1) = sum(TCD_PA_O(i,1:length(allday)));
    TC_PA(i,2) = sum(TCD_PA_I(i,1:length(allday)));
    TC_PA(i,3) = sum(TCD_PA_OR(i,1:length(allday)));
    TC_PA(i,4) = sum(TCD_PA_OL(i,1:length(allday)));
    TC_PA(i,5) = sum(TCD_PA_OS(i,1:length(allday)));

    %ÿ��������ˮǿ��
    TC_PR(i,1) = TC_PAm(i,1) *10^(6) / TC_PA(i,1);
    TC_PR(i,2) = TC_PAm(i,2) *10^(6) / TC_PA(i,2);
    TC_PR(i,3) = TC_PAm(i,3) *10^(6) / TC_PA(i,3);
    TC_PR(i,4) = TC_PAm(i,4) *10^(6) / TC_PA(i,4);
    TC_PR(i,5) = TC_PAm(i,5) *10^(6) / TC_PA(i,5);
end

%% ------------------������-------------------------------------------
% TCD_PR_O=Arr1;   TCD_PR_I=Arr1;   TCD_PR_OR=Arr1;   TCD_PR_OL=Arr1;   TCD_PR_OS=Arr1;        %ÿ������ÿ�ս�ˮǿ�� 
% TCD_PAm_O=Arr1;  TCD_PAm_I=Arr1;  TCD_PAm_OR=Arr1;  TCD_PAm_OL=Arr1;  TCD_PAm_OS=Arr1;       %ÿ������ÿ�ս�ˮ��
% TCD_MaxP_O=Arr1; TCD_MaxP_I=Arr1; TCD_MaxP_OR=Arr1; TCD_MaxP_OL=Arr1; TCD_MaxP_OS=Arr1;      %ÿ������ÿ�����ˮǿ��
% TCD_PA_O=Arr1;   TCD_PA_I=Arr1;   TCD_PA_OR=Arr1;   TCD_PA_OL=Arr1;   TCD_PA_OS=Arr1;        %ÿ������ÿ�ս�ˮ��� 
% TCD_MPCo_O=Arr2; TCD_MPCo_I=Arr2; TCD_MPCo_OR=Arr2; TCD_MPCo_OL=Arr2; TCD_MPCo_OS=Arr2;      %ÿ������ÿ�����ˮǿ�Ⱦ�γ��                                              %�� ����TC_MP_Coor
% TC_PR=Arr3; TC_PAm=Arr3; TC_Max_P=Arr3; TC_PA=Arr3; TC_MP_Coor=Arr3;
% TCR_PA_O=Arr4; TCR_PA_L=Arr4; TCR_PA_S=Arr4;%��ʵӰ�����
%----------------��������excel���
TCR_PA=[TCR_PA_O,Arr4,Arr4,TCR_PA_L,TCR_PA_S];

OUT_PR=   [TCD_PR_O ,   TCD_PR_I ,   TCD_PR_OR ,   TCD_PR_OL ,   TCD_PR_OS];
OUT_Am=   [TCD_PAm_O ,  TCD_PAm_I ,  TCD_PAm_OR ,  TCD_PAm_OL ,  TCD_PAm_OS];
OUT_MaxP= [TCD_MaxP_O , TCD_MaxP_I , TCD_MaxP_OR , TCD_MaxP_OL , TCD_MaxP_OS];
OUT_PA=   [TCD_PA_O ,   TCD_PA_I ,   TCD_PA_OR ,   TCD_PA_OL ,   TCD_PA_OS];
OUT_MPCo= [TCD_MPCo_O , TCD_MPCo_I , TCD_MPCo_OR , TCD_MPCo_OL , TCD_MPCo_OS];
OUT_Entir=[TC_PR,       TC_PAm,      TC_Max_P,     TCR_PA,       TC_MP_Coor];

OUT_PR(OUT_PR   ==NODATA_value)=missing;
OUT_Am(OUT_Am   ==NODATA_value)=missing;
OUT_MaxP(OUT_MaxP ==NODATA_value)=missing;
OUT_PA(OUT_PA   ==NODATA_value)=missing;
OUT_MPCo(OUT_MPCo ==NODATA_value)=missing;
OUT_Entir(OUT_Entir==NODATA_value)=missing;

xlswrite(Excel_file,OUT_PR(ST:EN,:)    ,'PR',[num2abc2(5),num2str(ST+2),':',num2abc2(5+250-1),num2str(EN+2)]);%��
xlswrite(Excel_file,OUT_Am(ST:EN,:)    ,'PAm',[num2abc2(5),num2str(ST+2),':',num2abc2(5+250-1),num2str(EN+2)]);%��
xlswrite(Excel_file,OUT_MaxP(ST:EN,:)  ,'MaxP',[num2abc2(5),num2str(ST+2),':',num2abc2(5+250-1),num2str(EN+2)]);%��  
xlswrite(Excel_file,OUT_PA(ST:EN,:)    ,'PA',[num2abc2(5),num2str(ST+2),':',num2abc2(5+250-1),num2str(EN+2)]);%��
xlswrite(Excel_file,OUT_MPCo(ST:EN,:)  ,'MPCo',[num2abc2(5),num2str(ST+2),':',num2abc2(5+500-1),num2str(EN+2)]);%��
xlswrite(Excel_file,OUT_Entir(ST:EN,:) ,'Entire',[num2abc2(5),num2str(ST+2),':',num2abc2(5+50 -1),num2str(EN+2)]);%��


%------------���ɱ�ͷ-----------------
THead_PR=[{'���','���','TCD_PR_O'},repmat({NaN},1,49),{'TCD_PR_I'},repmat({NaN},1,49),{'TCD_PR_OR'},repmat({NaN},1,49),{'TCD_PR_OL'},repmat({NaN},1,49),{'TCD_PR_OS'},repmat({NaN},1,49); ...
    repmat({NaN},1,2),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50)];
THead_PAm=[{'���','���','TCD_PAm_O'},repmat({NaN},1,49),{'TCD_PAm_I'},repmat({NaN},1,49),{'TCD_PAm_OR'},repmat({NaN},1,49),{'TCD_PAm_OL'},repmat({NaN},1,49),{'TCD_PAm_OS'},repmat({NaN},1,49); ...
    repmat({NaN},1,2),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50)];
THead_MaxP=[{'���','���','TCD_MaxP_O'},repmat({NaN},1,49),{'TCD_MaxP_I'},repmat({NaN},1,49),{'TCD_MaxP_OR'},repmat({NaN},1,49),{'TCD_MaxP_OL'},repmat({NaN},1,49),{'TCD_MaxP_OS'},repmat({NaN},1,49); ...
    repmat({NaN},1,2),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50)];
THead_PA=[{'���','���','TCD_PA_O'},repmat({NaN},1,49),{'TCD_PA_I'},repmat({NaN},1,49),{'TCD_PA_OR'},repmat({NaN},1,49),{'TCD_PA_OL'},repmat({NaN},1,49),{'TCD_PA_OS'},repmat({NaN},1,49); ...
    repmat({NaN},1,2),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50)];
THead_MPCo=[{'���','���','TCD_MPCo_O_LON'},repmat({NaN},1,49),{'TCD_MPCo_O_LAT'},repmat({NaN},1,49),{'TCD_MPCo_I_LON'},repmat({NaN},1,49),{'TCD_MPCo_I_LAT'},repmat({NaN},1,49),{'TCD_MPCo_OR_LON'},repmat({NaN},1,49),{'TCD_MPCo_OR_LAT'},repmat({NaN},1,49),{'TCD_MPCo_OL_LON'},repmat({NaN},1,49),{'TCD_MPCo_OL_LAT'},repmat({NaN},1,49),{'TCD_MPCo_OS_LON'},repmat({NaN},1,49),{'TCD_MPCo_OS_LAT'},repmat({NaN},1,49); ...
    repmat({NaN},1,2),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50),num2cell(1:50)];
THead_Entir=[{'���','���','������ˮǿ�� mm/d/ƽ��0.1�������� '},repmat({NaN},1,9),{'������ˮ�� 10�������� or ����ǧ��'},repmat({NaN},1,9),{'�������ˮǿ�� mm/d 0.1������'},repmat({NaN},1,9),{'������ˮ������� ƽ��ǧ��'},repmat({NaN},1,9),{'�������ˮ�ľ�γ��'},repmat({NaN},1,14); ...
    repmat({NaN},1,2),{'O','I','OR','OL','OS'},repmat({NaN},1,5),{'O','I','OR','OL','OS'},repmat({NaN},1,5),{'O','I','OR','OL','OS'},repmat({NaN},1,5),{'O','I','OR','OL','OS'},repmat({NaN},1,5),{'O LON','O LAT','I LON','I LAT','OR LON','OR LAT','OL LON','OL LAT','OS LON','OS LAT'},repmat({NaN},1,5)];

xlswrite(Excel_file,THead_PR     ,1,[num2abc2(3),num2str(1),':',num2abc2(5+250-1),num2str(2)]);%��
xlswrite(Excel_file,THead_PAm    ,2,[num2abc2(3),num2str(1),':',num2abc2(5+250-1),num2str(2)]);%��
xlswrite(Excel_file,THead_MaxP   ,3,[num2abc2(3),num2str(1),':',num2abc2(5+250-1),num2str(2)]);%��  
xlswrite(Excel_file,THead_PA     ,4,[num2abc2(3),num2str(1),':',num2abc2(5+250-1),num2str(2)]);%��
xlswrite(Excel_file,THead_MPCo   ,5,[num2abc2(3),num2str(1),':',num2abc2(5+500-1),num2str(2)]);%��
xlswrite(Excel_file,THead_Entir  ,6,[num2abc2(3),num2str(1),':',num2abc2(5+55 -1),num2str(2)]); %��

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

