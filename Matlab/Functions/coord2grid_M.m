function [WZ] = coord2grid_M(LON,LAT,cs)
%�������⾭γ�ȷ����������� INPUT -180-180 Nx1֧�� 
%��90��-90��-180��180
GC=ceil((LON+180)/cs);
GC(GC==0)=1;

GR=ceil((90-LAT)/cs);
GR(GR==0)=1;
%����λ��
WZ=(GC-1)*(180/cs)+GR;
end
