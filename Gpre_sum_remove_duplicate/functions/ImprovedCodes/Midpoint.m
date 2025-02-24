function [lonMid, latMid] = Midpoint(lon1, lat1, lon2, lat2)
% midpoint of two lat long cord on a sphere, all units are deg
%-180-180 0-360都可以但是结果不一样 
Bx = cosd(lat2) * cosd(lon2-lon1);
By = cosd(lat2) * sind(lon2-lon1);
latMid = atan2d(sind(lat1) + sind(lat2), ...
               sqrt( (cosd(lat1)+Bx)*(cosd(lat1)+Bx) + By*By ) );
lonMid = lon1 + atan2d(By, cosd(lat1) + Bx);
end