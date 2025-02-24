%计算观测和在分析之间的相关关系表，DATA3TMPA60NS
%% -------------DIST---------------
%输入4列，每个4列中前两列为再分析，后两列为观测，计算观测和在分析之间的关系
clear
data=zeros(1940,4); %ERA5 JRA55 IMERG TRMM
data(:,1)=xlsread('C:\Users\29585\Desktop\Desktop2\Global_cyclone_project\DATA3TMPA60NS.xlsx',1,'H6113:H8052');
data(:,2)=xlsread('C:\Users\29585\Desktop\Desktop2\Global_cyclone_project\DATA3TMPA60NS.xlsx',1,'X6113:X8052');
data(:,3)=xlsread('C:\Users\29585\Desktop\Desktop2\Global_cyclone_project\DATA3TMPA60NS.xlsx',1,'AY6113:AY8052');
data(:,4)=xlsread('C:\Users\29585\Desktop\Desktop2\Global_cyclone_project\DATA3TMPA60NS.xlsx',1,'BG6113:BG8052');

P_val=zeros(2,2);
corr_val=zeros(2,2);
test=zeros(1,4);

for re=1:2
    for ob=3:4
        [corr_val(re,ob-2),P_val(re,ob-2)]=corr(data(:,re),data(:,ob),'Type','Pearson');%Spearman_rho
    end
end
for i =1:4
    test(i)=lillietest(data(:,i));
end

%% 从migration文件夹里的Correlation表中计算每个距离区间的降水强度关系
%由于都是降水，感觉没有必要计算相关系数，因为肯定相关嘛，所以也可以计算一下CC等
clear;
P_val=zeros(2,2);
corr_val=zeros(2,2);
test=zeros(2,2);
RMSE=zeros(2,2,29);ME=zeros(2,2,29);MAE=zeros(2,2,29);BIAS=zeros(2,2,29);ABIAS=zeros(2,2,29);
P_val_d=zeros(2,2,28);corr_val_d=zeros(2,2,28);%带_d的意思就是距离区间
n=zeros(2,2,29);%第一个都是整体的，后面28个是分盆地的

G=3;            %算单个记录1还是分区年均0还是单个记录加分区3
Count0=0;       %是否考虑都是零的值？1考虑0否
tic;
for RE=1:2      %再分析
    if G==1
        D1=xlsread('C:\Users\29585\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',RE,'B3:B4461');%降水
    elseif G==0
        D1=xlsread('C:\Users\29585\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',RE,'AF3:BG43');%降水C3:AD4461 AF3:BG43
    elseif G==3
        D1=xlsread('C:\Users\29585\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',RE,'C3:AD4461');%降水
    end
    for OB=3:4  %观测值
        
        if G==1
            D2=xlsread('C:\Users\29585\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',OB,'B3:B4461');%降水
        elseif G==0
            D2=xlsread('C:\Users\29585\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',OB,'AF3:BG43');%影响因素  
        elseif G==3
            D2=xlsread('C:\Users\29585\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',OB,'C3:AD4461');%降水
        end

        if G==0 || G==3
        %-----------------计算分距离区间的误差
        for disint=1:28
            co1_d=D1(:,disint);
            co2_d=D2(:,disint);

            if Count0==0 co1_d(co2_d==0 & co1_d==0)=-9999; end%都等于零的话就不要了

            co2_d(co1_d==-9999)=[];co1_d(co1_d==-9999)=[];
            co1_d(co2_d==-9999)=[];co2_d(co2_d==-9999)=[];
            

            [~,RMSE(RE,OB-2,disint+1),ME(RE,OB-2,disint+1),MAE(RE,OB-2,disint+1),BIAS(RE,OB-2,disint+1),ABIAS(RE,OB-2,disint+1),~,~,~]=CC_RMSE_ME_MAE_BIAS_ABIAS([co2_d,co1_d],1,1,1); 
            [corr_val_d(RE,OB-2,disint),P_val_d(RE,OB-2,disint)]=corr(co1_d,co2_d,'Type','Spearman') ;%Spearman_rho
            n(RE,OB-2,disint+1)=length(co1_d);
        end
        %-----------------基本处理
            co1=reshape(D1,[28*length(D1),1]);
            co2=reshape(D2,[28*length(D1),1]);
        elseif G==1
            co1=D1;
            co2=D2;
        end

        if Count0==0 co1(co2==0 & co1==0)=-9999; end%都等于零的话就不要了

        co2(co1==-9999)=[];co1(co1==-9999)=[];      %适用于年份不足的时候
        co1(co2==-9999)=[];co2(co2==-9999)=[];      %适用于海温的情况，excel里面原始是没有nan只有空值会有nan

        
        %-----------------计算不区分distance intervals
        [~,RMSE(RE,OB-2,1),ME(RE,OB-2,1),MAE(RE,OB-2,1),BIAS(RE,OB-2,1),ABIAS(RE,OB-2,1),~,~,~]=CC_RMSE_ME_MAE_BIAS_ABIAS([co2,co1],1,1,1);
        [corr_val(RE,OB-2),P_val(RE,OB-2)]=corr(co1,co2,'Type','Spearman') ;%Spearman_rho
        lillietest(co1)
        lillietest(co2)
        n(RE,OB-2,1)=(length(co1));
    end
end











