function [LonCenTab,LatCenTab] = GridCenterLocation_GPU(cs)
%GRIDCENTERLOCATION INPUT grid size return center lon lat table -180-180
%��90��-90��-180��180

i=gpuArray((1:180/cs)');
LatCenTab=repmat((90-i*cs+0.5*cs),[1 360/cs]);

j=gpuArray((1:360/cs));
LonCenTab=repmat((-180+j*cs-0.5*cs),[180/cs 1]);

end