function [GR,GC] = latlon2grid(LON,LAT,size)
%�������⾭γ�ȷ�����������LON-180-180 
%��90��-90��-180��180
GC=ceil((LON+180)/size);
if GC==0 
    GC=1; 
end

GR=ceil((90-LAT)/size);
if GR==0 
    GR=1; 
end
%�����������������
end
