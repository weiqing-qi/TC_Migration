function [londen, latden] = TCDenseTrackPoint(lon,lat,cs)
%输入nx1经纬度进行加密度处理，使得轨迹上每个网格内都有经纬度点
%lon lat 顺序是气旋时间顺序 输入无所谓，输出为-180-180
%240119 TESTED
londen=lon;
latden=lat;

D = Dist_i_ip1(lon, lat);%每两点之间的距离
MAXD=GL_i_ip1(lat);      %每两点之间的最大距离
InterpWZ=find(D >= MAXD);

while InterpWZ
    for i=1:length(InterpWZ)
        %---对每一个超过距离需要插入中点的位置计算中点经纬度
        [lonMid, latMid] = Midpoint(lon(InterpWZ(i)), lat(InterpWZ(i)), ...
                                    lon(InterpWZ(i)+1), lat(InterpWZ(i)+1));
        %插入值到数组中间
        londen=[londen(1:InterpWZ(i)+i-1);lonMid;londen(InterpWZ(i)+i:end)];%用i考虑之前已经插入过的几个变量
        latden=[latden(1:InterpWZ(i)+i-1);latMid;latden(InterpWZ(i)+i:end)];%InterpWZ(i)+(i-1)+1
    end
    %更新InterpWZ，lon，lat
    lon=londen;
    lat=latden;
    InterpWZ=find( Dist_i_ip1(londen, latden) >= GL_i_ip1(latden) );
end

londen(londen>180)=londen(londen>180)-360;                      % 注意最佳路径数据中的经纬度并不是-180-180，当跨越180经度的时候符号会保持不变
londen(londen<-180)=londen(londen<-180)+360;                    %确保输入TCPoint2RLine为-180-180

%-------------------------------function-----------------------------------
function[distance]=Dist_i_ip1(loncolume, latcolume)
    %计算列，nx1经，nx1纬，每两个顺序相邻点之间的距离
    distance=zeros(length(loncolume)-1,1);
    for j=1:length(loncolume)-1
        distance(j)=SphereDist2_Matrix([loncolume(j),latcolume(j)],[loncolume(j+1),latcolume(j+1)]);
    end
end

function[GridLength]=GL_i_ip1(latcolume)
    %计算列，nx1经，nx1纬，每两个顺序相邻点纬度更高的点对应的格网边长
    GridLength=zeros(length(latcolume)-1,1);
    for j=1:length(latcolume)-1
        highgridlat=ceil(max(abs( [latcolume(j),latcolume(j+1)] ))/cs)*cs;%得到两点间最大的纬度所在网格最短的一条边的纬度
        GridLength(j)=SphereDist2_Matrix([0,highgridlat],[cs,highgridlat]);       %计算边长在上述最短网格边长以内就不会出现间隔格网                                     
    end
end

end