%% 计算相关性TC-周围影响因素，这种用28个数据的样本太少，应该不合适，还是主要用下面 的方法
% corr_valT=zeros(4,3*4);
% corr_valD=zeros(4,3*4);
% P_valD=zeros(4,3*4);
% P_valT=zeros(4,3*4);
% 
% D=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA3TMPA60NS.xlsx',5,'M97:X124');%改影响因素L97:W124 L128:W155
% d1=D(:,7:10);%改 影响因素
% d2=D(:,11:12);%改 PRE数据
% d1=d1(1:28,:);%改
% d2=d2(1:28,:);%改 
% 
% for j=3:4%影响因素 改
%     for i=3:4 %PRE数据
%     [corr_valT(j,i*3-2),P_valT(j,i*3-2)]=corr(d1(:,j),d2(:,i-2),'Type','Pearson');%Pearson_COE
%     [corr_valT(j,i*3-1),P_valT(j,i*3-1)]=corr(d1(:,j),d2(:,i-2),'Type','Kendall');%kendall_tau
%     [corr_valT(j,i*3  ),P_valT(j,i*3  )]=corr(d1(:,j),d2(:,i-2),'Type','Spearman');%Spearman_rho
%     end
% end
% 
% out=[corr_valT;corr_valD];
%------------------年平均
% Data=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',5,'C3:AD4461');%降水
% year1=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',5,'A3:A4461');%降水
% R_PR_YAVE=zeros(41,28);%-800：2000km
% for yi=1980:2020
%     for DR=1:28
%        R1C=Data(:,DR); 
%        R1C=R1C(year1==yi);
% %        R1C(R1C<0.1 & R1C>=0)=[];
%        R_PR_YAVE(yi-1979,DR)=mean(R1C);
%     end
% end
%----------------------------------------------------------------------
%% 计算每场台风每一个距离区间的降水
P_val=zeros(4,3*4);
corr_val=zeros(4,3*4);
test=zeros(4,3*4);
G=0;%算全球1还是分区0改
U=1;%改
D=0;%改
tic;
for pre=1:4 %全球
    if G==1
        D1=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',pre,'B3:B4461');%降水
    elseif G==0
        D1=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',pre,'AF3:BG43');%降水C3:AD4461 AF3:BG43
    end
    for influ=5:8 %盆地内
        
        if G==1
            D2=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',influ,'B3:B4461');%降水
        elseif G==0
            D2=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',influ,'AF3:BG43');%影响因素    
        end
        %-----------------基本处理
        if G==0
            co1=reshape(D1,[28*length(D1),1]);
            co2=reshape(D2,[28*length(D1),1]);
        elseif G==1
            co1=D1;
            co2=D2;
        end
        co2(co1==-9999)=[];co1(co1==-9999)=[];%适用于年份不足的时候
        co1(co2==-9999)=[];co2(co2==-9999)=[];%适用于海温的情况，excel里面原始是没有nan只有空值会有nan
        co2(co1<=0.1)=[];co1(co1<=0.1)=[];%降水0不要，只有降水会小于0.1
        %-----------------限制条件,选择大小门槛中的数据计算关系
        thre=co1;
        I1=sort(co1);%只用降水做门槛 升序
        threU=I1(round(length(I1)*U));%75%百分位数
        if D==0
            threD=0.1;
        else
            threD=I1(round(length(I1)*D));%5%百分位数
        end
        co1=co1(thre > threD & thre <= threU);
        co2=co2(thre > threD & thre <= threU);
        %-----------------计算
        [corr_val(influ-4,pre*3-2),P_val(influ-4,pre*3-2)]=corr(co1,co2,'Type','Pearson');%Pearson_COE
        [corr_val(influ-4,pre*3-1),P_val(influ-4,pre*3-1)]= ;%Spearman_rho
        [corr_val(influ-4,pre*3  ),P_val(influ-4,pre*3  )]=corr(co1,co2,'Type','Kendall');%kendall_tau
        test(influ-4,pre*3-2  )=lillietest(co1);test(influ-4,pre*3-1  )=lillietest(co2);
        n=(length(co1))
    end
end
toc;
P_val;
[G D U]
corr_val











