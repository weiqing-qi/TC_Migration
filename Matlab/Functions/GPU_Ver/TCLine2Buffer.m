function [Buff_OWZ,Buff_OMask,Buff_IWZ,Buff_IMask,RMR] = TCLine2Buffer(WZ,Oradius,Iradius,cs,LonCT_gpu,LatCT_gpu,Dailydata,RMRflag,DRMW)
%气旋轨迹点格计算buffer input from TCPoint2RLine -180-180 
%cs=0.5时13个点耗时0.106″ Unit radius: km 输入两个值，分别是气旋外范围和内核半径

%[LonCT_gpu,LatCT_gpu] = GridCenLocat_gpu(cs);%格点中心经纬度索引表
TGCLON=zeros(180/cs*360/cs,1); %ALL Grid size initialization
TGCLAT=zeros(180/cs*360/cs,1);
Buff_OMask=zeros(180/cs,360/cs);
Buff_IMask=zeros(180/cs,360/cs);

if RMRflag==1
    D=zeros(180/cs*360/cs,length(WZ));

    for i=1:length(WZ)
        %一次计算当前点与所有网格点的距离
        TGCLON(:)=LonCT_gpu(WZ(i)); 
        TGCLAT(:)=LatCT_gpu(WZ(i));
        D(:,i) = SphereDist2_Matrix([TGCLON,TGCLAT],[reshape(LonCT_gpu,[],1),reshape(LatCT_gpu,[],1)]); %每一列对应一个denlon位置
    end

    minD=min(D,[],2);           %每个网格到路径的最小距离(近似)
    %如何计算RMR
    if ~isempty(DRMW) %有RMW用RMW
       RMR=max(DRMW);
       if RMR<cs*100/2  
           RMR=cs*100/2;                 %内核不能太小，小于网格点
       end
    else    
        M_amv=0;
        for dist=0:cs*100:500
            C_amv=mean(Dailydata(minD>=dist & minD<dist+cs*100));
            if C_amv>M_amv
                M_amv=C_amv;
                RMR=dist+cs*100;
            end
        end
        if RMR>100  
            RMR=100;                %内核不能超过200km 感觉和RMW差的还有点多存在异常值
        end
        if RMR<cs*100/2  
            RMR=cs*100/2;                 %内核不能太小，小于网格点
        end
    end
    Buff_IMask(minD<=RMR*2)=1;                                           %这里尽可能采用较大的RMR        
    Buff_OMask(minD<=Oradius)=1;                                         %这里注意BM要和D大小一样

    Buff_OWZ=find(Buff_OMask==1);
    Buff_IWZ=find(Buff_IMask==1);
    
elseif RMRflag==0
    RMR=[];
    for i=1:length(WZ)
        %一次计算当前点与所有网格点的距离
        TGCLON(:)=LonCT_gpu(WZ(i)); 
        TGCLAT(:)=LatCT_gpu(WZ(i));
        D = SphereDist2_Matrix([TGCLON,TGCLAT],[reshape(LonCT_gpu,[],1),reshape(LatCT_gpu,[],1)]);
        Buff_OMask(D<=Oradius)=1;%这里注意BM要和D大小一样
        Buff_IMask(D<=Iradius)=1;
    end
    Buff_OWZ=find(Buff_OMask==1);
    Buff_IWZ=find(Buff_IMask==1);
end
end