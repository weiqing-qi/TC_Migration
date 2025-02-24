function [WZ,Grid_track] = TCDenPoint2Raster(londen, latden,cs)
%气旋加密轨迹点作为轨迹线映射到网格尺度 -180-180
Longrid=ceil((180+londen)/cs);
Longrid(Longrid==0)=1;
Latgrid=ceil((90-latden)/cs);
Latgrid(Latgrid==0)=1;
%去重得到位置
WZ=unique((Longrid-1)*(360/cs)+Latgrid);
%制作轨迹线掩膜
Grid_track=zeros(180/cs,360/cs);
Grid_track(WZ)=1;
end

