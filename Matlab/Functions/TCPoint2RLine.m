function [WZ,Grid_track] = TCPoint2RLine(londen, latden,cs)
%气旋加密轨迹点作为轨迹线映射到网格尺度 -180-180 NX1 input from TCDenseTrackPoint
%240120TESTED
Longrid=ceil((180+londen)/cs);%第几列
Longrid(Longrid==0)=1;
Latgrid=ceil((90-latden)/cs);%第几行
Latgrid(Latgrid==0)=1;
%去重得到位置 nx1 注意这里之后WZ无法与读取的最佳路径中的时间数组对应
WZ=unique((Longrid-1)*(180/cs)+Latgrid);
%制作轨迹线掩膜
Grid_track=zeros(180/cs,360/cs);
Grid_track(WZ)=1;
end

