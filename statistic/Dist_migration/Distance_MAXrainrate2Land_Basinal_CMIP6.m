%计算盆地
clear; clc; close all;
allmodel_name={'BCC-CSM2-MR','CAS-FGOALS-f3-L','CMCC-ESM2','CNRM-CM6-1-HR','EC-Earth3', ...
            'EC-Earth3-AerChem','EC-Earth3-CC','EC-Earth3-Veg','NOAA-GFDL-CM4','SNU-SAM0-UNICON'};
ALL_Results = NaN(2014-1980+1,9*length(allmodel_name));
for MOD = 1:length(allmodel_name)

    model_name = allmodel_name(MOD)

    if ismember(model_name, {''})%跳过的模型
        continue
    end

    % if ismember(model_name,{'BCC-CSM2-MR','CAS-FGOALS-f3-L','CMCC-ESM2', 'NOAA-GFDL-CM4', 'SNU-SAM0-UNICON'})
    %     cs = 1;                                                                 
    % elseif ismember(model_name,{'EC-Earth3','EC-Earth3-AerChem','EC-Earth3-CC','EC-Earth3-Veg' })
    %     cs = 180/256;                                                           
    % elseif ismember(model_name,{'CNRM-CM6-1-HR'})
    %     cs = 0.5;                                                               
    % end
    cs=0.5;

    YDLL=xlsread('D:\Desktop2\Global_cyclone_project\DATA_CMIP6.xlsx',1, ...
        [num2abc2(MOD*5-4),'3:',num2abc2(MOD*5),'5000']);
    
    Year=YDLL(:,1);
    Distance=YDLL(:,2);
    LON=YDLL(:,3);
    LAT=YDLL(:,4);
    
    [num2str(length(LAT)),'=',num2str(length(LON)),'=',num2str(length(Year)),'=',num2str(length(Distance))]
    Bas=zeros(length(LON),1);
    [LonCenter,LatCenter] = GridCenterLocation(cs);
    Basinmasks = Basinmasks_EPNA(cs,1);
    
    for i=1:length(LON)
    R=ceil((90-LAT(i))/cs);
    C=ceil((180+LON(i))/cs);
    Bas(i)=Basinmasks(R,C);%Bas数字
    end
    
    AveYB=zeros(Year(length(Year))-Year(1)+1,7);
    for Y=Year(1):Year(length(Year))        
        Ydis=Distance(Year==Y);%每年每场TC的距离
        YB=Bas(Year==Y);%每年每场TC的盆地
        YB(Ydis<-50 | Ydis>2000)=[];%去掉-50和2000
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
    ALL_Results(:,MOD*9-8:MOD*9) = [(1980:2014)',AveYB,NaN(35,1)];
end
