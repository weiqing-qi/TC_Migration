%�������
clear; clc;
S='3:';
E='10000';

Year=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted_CMIP6.xlsx',1,['A',S,'A',E]);
Distance=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted_CMIP6.xlsx',1,['B',S,'B',E]);%��
LON=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted_CMIP6.xlsx',1,['C',S,'C',E]);%��
LAT=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted_CMIP6.xlsx',1,['D',S,'D',E]);%��

[num2str(length(LAT)),'=',num2str(length(LON)),'=',num2str(length(Year)),'=',num2str(length(Distance))]
Bas=zeros(length(LON),1);
[LonCenter,LatCenter] = GridCenterLocation(0.5);
Basinmasks = Basinmasks_EPNA(0.5,1);

for i=1:length(LON)
R=ceil((90-LAT(i))*2);
C=ceil((180+LON(i))*2);
Bas(i)=Basinmasks(R,C);%Bas����
end

AveYB=zeros(Year(length(Year))-Year(1)+1,7);
for Y=Year(1):Year(length(Year))        
    Ydis=Distance(Year==Y);%ÿ��ÿ��TC�ľ���
    YB=Bas(Year==Y);%ÿ��ÿ��TC�����
    YB(Ydis<-50 | Ydis>2000)=[];%ȥ��-50��2000
    Ydis(Ydis<-50 | Ydis>2000)=[];
    
    for b=1:7
        if length(Ydis)==length(YB)
            AveYB(Y-Year(1)+1,b)=mean(Ydis(YB==b));
        else
            disp('FALSE')
            break
        end
    end  
end

