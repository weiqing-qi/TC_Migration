function [weizhi,passyear,passmonth,passday,locationR,locationC] = raster_grid_TC_pass_period(header,raster,year,month,day,lon,lat,ran)
%This is an update version of calculate time of each grid from TCRoutedata 
%The passing time is calculated rigidly, including every period it appeared
%in 500 km. And we add the former & later 1 day.
%update2023.9.6
%Raster is 0-360 

[locationR,locationC]= find(raster==1);
weizhi=find(raster==1);

Ltrac = length(lon);
LG=length(weizhi);

%ÿһ��һ������ÿ�м�¼������ʱ�䣬���ΪLtrac��ȡ3LtracΪ�˷������BF1day�ν�
%��չ��2Ltrac��֤����BF1dayʱ���Բ��������鳤��
passyear=cell(LG,3*Ltrac);                     
passmonth=cell(LG,3*Ltrac);
passday=cell(LG,3*Ltrac);
                         
for i= 1:LG                                                                %ÿһ����Ĥ����
    %����500km�������ľ�γ��
    GCLON=locationC(i)*header(5)-0.5*header(5);
    GCLAT=90-locationR(i)*header(5)+0.5*header(5);
    MGLON=GCLON*ones(Ltrac,1);                                             %ͳһSphereDist_Matrix�����ά��
    MGLAT=GCLAT*ones(Ltrac,1);

    %����ÿ���켣�㵽��ǰ�������ľ���(Spheredist�������0-360��-180-180�������ý������)
    Dist=SphereDist_Matrix([MGLON,MGLAT],[lon,lat]);                       %��λKM ������
    
    
    Lpass=length(find(Dist<=ran));
    passyear (i,1:Lpass) = year (Dist<=ran)';                              %�п�����Ȼ�Ҳ����κ�ֵ������ʱ�����Ϊ��
    passmonth(i,1:Lpass) = month(Dist<=ran)';                              %Lpass=0ʱ������и�ֵ
    passday  (i,1:Lpass) = day  (Dist<=ran)'; 
 
    %ÿ��ÿ��ʱ��ϵ�
    Index=find(Dist<=ran);                                                 %�õ�TC���������λ�����ڵ�ȫ������
    IndexB=unique(Index-(1:Lpass)');                                       %�����ֵlike:000333555���������������Ȳ�1,ȥ��Ȼ�����������Ѱ
    Ls=length(IndexB);
    
    BF1D = zeros(1,2*Ls);                                                  %ÿһ�������������˶�Ҫ��һ�죬����������
    BFM  = zeros(1,2*Ls);
    BFY  = zeros(1,2*Ls);

    for j= 1:Ls
        %-------����ϵ��Backdayǰһ��
        Ord = find(index==IndexB(j),1,'first');                            
        [BF1D(2*j-1),BFM(2*j-1),BFY(2*j-1),~,~,~]=daybackdayforward( ...   
            str2double(char(day(Ord))), ...
            str2double(char(month(Ord))), ...
            str2double(char(year(Ord))) );
        %-------����ϵ��Forwardday��һ��
        inOrd= find(index==IndexB(j),1,'last' ); 
        [~,~,~,BF1D(2*j),BFM(2*j),BFY(2*j)]=daybackdayforward( ...         
            str2double(char(day(inOrd))), ...
            str2double(char(month(inOrd))), ...
            str2double(char(year(inOrd))) );
    end

    %��BF1day���������ӵ�passdmy������ⲿ��������ȥ������Ͳ�ȥ����
    passyear(i,Ltrac+1:Ltrac+2*Ls)=BFY;  %Lpass<Ltrac                   
    passmonth(i,Ltrac+1:Ltrac+2*Ls)=BFM;
    passday(i,Ltrac+1:Ltrac+2*Ls)=BF1D;
end
end
