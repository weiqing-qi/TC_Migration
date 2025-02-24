%计算每天的RMRW和动态oradius加快速度
AI_Mask_dir= 'D:\DATA\AI_IRcloud_TC_mask_daily.nc4';
IMERG_idir='E:\DATA\2.Satellite\IMERG_Final_V7B_Daily\';

TC_idir='D:\DATA\Segmentation map tutorial\IBTrACS.since1980.v04r00.240321.nc';
[SID,LAT,LON,ISO_TIME,USA_WIND,USA_RMW]=IBTrACS_nc_entire_variable_r(TC_idir);
cs=0.1;
NOlim=lim_Y2NO([2001,2022],ISO_TIME); 
[LonCT_gpu,LatCT_gpu] = GridCenLocat_gpu(cs);

ORADIUUS=zeros(NOlim(3),50);
RMRW=zeros(NOlim(3),50);

for NO=1:NOlim(3)
    NO
    time=ISO_TIME(:,:,NO + NOlim(1)-1);
    lat=LAT(:,NO + NOlim(1)-1); 
    lat(isnan(lat) )=[];
    lon=LON(:,NO + NOlim(1)-1);
    lon(isnan(lon) )=[]; 
    USA_RMW_1NO=USA_RMW(:,NO + NOlim(1)-1);
    USA_RMW_1NO(isnan(lat))=[];%这里必须是lat注意
    %-------每场气旋(每三小时)的日尺度时间点数值形式get-----------------------
    year   =str2num(time(1:4,1:length(lat))');                          %str2num可以对多行字符串使用变成一个数组
    month  =str2num(time(6:7,1:length(lat))');                          %str2double只能对一个值使用
    day    =str2num(time(9:10,1:length(lat))');
    numtime=year*10^4+month*10^2+day;

    allday=unique(numtime);    
    allday_S=num2str(allday);
    allmon=unique(year*10^2+month);
    allmon_S=num2str(allmon);

    alldata=zeros(180/cs,360/cs,length(allday),'gpuArray');
    %-------本场气旋所涉及降水数据读取--------------------------------
    for k=1:length(allday)
        data = ncread( ...
            [IMERG_idir,'3B-DAY.MS.MRG.3IMERG.',allday_S(k,:), ...
            '-S000000-E235959.V07B.nc4'],'precipitation');              % 改
        alldata(:,:,k)=flipud(data);                                    %读取出来的IMERG数据需要翻转，具体见文件内
    end
    alldata(alldata<0)=0;
    alldata(isnan(alldata))=0;

    EventMasks=permute(ncread(AI_Mask_dir,'daily_mask',[NO,1,1],[1,inf,45]),[2 3 1]);
    for d=1:length(allday)
        DailyIndex=find(numtime==allday(d));
        if d~=length(allday)                                            %补足每一天21点到24点的掩膜,也就是把第二天0点加进来
            DailyIndex=[DailyIndex;DailyIndex(end)+1];
        end
        Dlon=lon(DailyIndex);                                            
        Dlat=lat(DailyIndex);
        DRMW=USA_RMW_1NO(DailyIndex);
        DRMW(isnan(DRMW))=[];

        %计算气旋的固定距离掩膜
        [londen,latden] = TCDenseTrackPoint(Dlon,Dlat,cs);              %Dlon Dlat 必须是气旋行进顺序
        londen(londen>180)=londen(londen>180)-360;                      % 注意最佳路径数据中的经纬度并不是-180-180，当跨越180经度的时候符号会保持不变
        latden(londen<-180)=londen(londen<-180)+360;                    %确保输入TCPoint2RLine为-180-180                
        [WZ,~] = TCPoint2RLine(londen,latden,cs); 
        
        Oradius=500;                                                    %这里需要一个控制radius的程序!!!!
        Iradius=100;
        %与固定距离不同这里只能导入IMERG来计算RMR
        Dailydata=alldata(:,:,d);
        [Buff_OWZ,Buff_OMask,Buff_IWZ,Buff_IMask,RMR] = TCLine2Buffer_GPU(WZ,Oradius,Iradius,cs,LonCT_gpu,LatCT_gpu,Dailydata,1,DRMW);
        
        %提取AI_Mask   
        AI_OMask=zeros(180/cs,360/cs);
        DayMask=EventMasks(:,d);
        DayMask(DayMask==0)=[]; %MISSING VALUE=0
        DayMask(isnan(DayMask))=[];
        AI_OMask(DayMask)=1;
        %--------------remove flash--------------
        AIcc = bwconncomp(AI_OMask,4);                                  %计算连通分量，用的是4连通，默认8连通
        LableM = labelmatrix(AIcc);                                     %用标签矩阵标注连通分量 背景值：0                                                                     
        LableLink2Buff=unique(LableM(Buff_OWZ));
        LableLink2Buff(LableLink2Buff==0)=[];                           %500km中心掩膜必须相连通
        AI_OMask=ismember(LableM,LableLink2Buff);                      %逻辑矩阵在数值上下文中被直接当作数据
        AI_OMask=double(AI_OMask);%unit8不能和double.*
        %-------最终的AImask被用来计算一个新的Oradius
        % 考虑到RB的不稳定以及可能并无变化，稳一点还是不用RB计算半径了，只考虑AImask
        AIWZ=find(AI_OMask==1);
        TGCLON=zeros(length(londen),1,'gpuArray'); %ALL Grid size initialization
        TGCLAT=zeros(length(londen),1,'gpuArray');
        if ~isempty(AIWZ)
            D=zeros(length(AIWZ),1);
            for i=1:10:length(AIWZ)%太慢了，缩小10倍
                TGCLON(:)=LonCT_gpu(AIWZ(i)); 
                TGCLAT(:)=LatCT_gpu(AIWZ(i));
                CRadi = SphereDist2_Matrix([TGCLON,TGCLAT],[londen,latden]); %每一列对应一个denlon位置
                D(i)= min(CRadi);%先求每个掩膜点到路径的最短距离然后把最大的D作为半径   
            end
            OradiusNEW=max(D);
        else 
            OradiusNEW=500;
        end
        ORADIUUS(NO,d)= OradiusNEW;
        RMRW(NO,d)= RMR;
    end
end