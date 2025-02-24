function [LonCenTab,LatCenTab] = GridCenterLocation(cs)
%GRIDCENTERLOCATION INPUT grid size return center lon lat table -180-180
%ио90об-90вС-180ср180

i=(1:180/cs)';
LatCenTab=repmat((90-i*cs+0.5*cs),[1 360/cs]);

j=(1:360/cs);
LonCenTab=repmat((-180+j*cs-0.5*cs),[180/cs 1]);


end
