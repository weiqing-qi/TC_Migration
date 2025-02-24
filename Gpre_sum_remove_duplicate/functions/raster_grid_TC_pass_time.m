function [weizhi,passyear,passmonth,passday,locationR,locationC] = raster_grid_TC_pass_time(header,raster,year,month,day,lon,lat)
%calculate time of each grid from TCRoutedata
%update20201020

[locationR,locationC]= find(raster==1);
weizhi=find(raster==1);
passyear=cell(length(weizhi),1);
passmonth=cell(length(weizhi),1);
passday=cell(length(weizhi),1);

for i= 1:length(weizhi) %ÿһ����Ĥ����
    GCLON=locationC(i)*header(5)-0.5*header(5);%locationC��0-360����������
    GCLAT=90-locationR(i)*header(5)+0.5*header(5);
    minDist=40076;%����ܳ�
    
    for j=1:length(lon)   %ÿһ��̨��λ�� ����������Ӧ��ʱ��
        Dist=SphereDist([GCLON,GCLAT],[lon(j),lat(j)]); %��λKM lon��-180-180��������������Һ������һ��
        if Dist<minDist
            minDist=Dist;
            Syear=year(j);    %����������Ӧ��ʱ��
            Smonth=month(j);
            Sday=day(j);
        end
    end
    
   passyear(i)= Syear;%Ӧ�û�ȱһ�������յ�ÿ��ڼ����ת��  
   passmonth(i)= Smonth;
   passday(i)= Sday;
end
