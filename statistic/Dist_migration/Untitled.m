%从Migration表格中读取并进行线性回归
clear; clc;
Str1='CFILORUX';
Str2='AAADAGAJAMAP';
jud=1; %改

if jud==1%land
    
out1=zeros(8,3); 
for i=1:8
Data=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\MigrationTMPA60andDWtest.xlsx',7,[Str1(i),'5:',Str1(i),'12']);%改
[b,bint,r,rint,stats]=regress(Data(:,1),[ones(length(Data(:,1)),1),(50-100*length(Data(:,1)):100:-50)'],0.05);
out1(i,1)=b(2);
out1(i,2)=b(2)-bint(2,1);%注意
out1(i,3)=stats(3); 
end
out2=zeros(6,3); 
for i=1:6
Data=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\MigrationTMPA60andDWtest.xlsx',7,[Str2(2*i-1:2*i),'5:',Str2(2*i-1:2*i),'12']);%改
[b,bint,r,rint,stats]=regress(Data(:,1),[ones(length(Data(:,1)),1),(50-100*length(Data(:,1)):100:-50)'],0.05);
out2(i,1)=b(2);
out2(i,2)=b(2)-bint(2,1);%注意
out2(i,3)=stats(3); 
end

elseif jud==2%sea
    
out1=zeros(8,3); 
for i=1:8

if i==7 || i==8%NI有两行没有
    Data=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\MigrationTMPA60andDWtest.xlsx',8,[Str1(i),'13:',Str1(i),'30']);%改
    [b,bint,r,rint,stats]=regress(Data,[ones(18,1),(50:100:1750)'],0.05);    
else
    Data=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\MigrationTMPA60andDWtest.xlsx',8,[Str1(i),'13:',Str1(i),'32']);%改
    [b,bint,r,rint,stats]=regress(Data,[ones(20,1),(50:100:1950)'],0.05);
end
out1(i,1)=b(2);
out1(i,2)=b(2)-bint(2,1);%注意
out1(i,3)=stats(3); 
end
out2=zeros(6,3); 
for i=1:6
Data=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\MigrationTMPA60andDWtest.xlsx',8,[Str2(2*i-1:2*i),'13:',Str2(2*i-1:2*i),'32']);%改
[b,bint,r,rint,stats]=regress(Data,[ones(20,1),(50:100:1950)'],0.05);
out2(i,1)=b(2);
out2(i,2)=b(2)-bint(2,1);%注意
out2(i,3)=stats(3); 
end    
out2(:,1)=out2(:,1)*-1;%注意趋势反了
out1(:,1)=out1(:,1)*-1;
end
% %% data2 P1表格盆地
%  %Data=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA3TMPA60NS.xlsx',3,'AI48:AL88');%AI48:AL88
% out=zeros(4,6);
% out2=zeros(4,6);
% for i=1:4
% d=Data(:,i);d(isnan(d))=[];
% 
% if length(d)==41
% [Zs, p_value, UFk, UBk2, beta, beta_CI, corr_val]= MKtrend((1980:2020)',d);
% [b,bint,r,rint,stats]=regress(d,[ones(2020-1980+1,1),(1980:2020)'],0.05);
% %====DW检验NW校正=====
% [ci_lower,ci_upper,pvals,b_Test,DW_p] = NeweyWestAdjust(r,(1980:2020)',d,1,1);
% 
% elseif length(d)==19
% [Zs, p_value, UFk, UBk2, beta, beta_CI, corr_val]= MKtrend((2001:2019)',d);
% [b,bint,r,rint,stats]=regress(d,[ones(2019-2001+1,1),(2001:2019)'],0.05);
% %====DW检验NW校正=====
% [ci_lower,ci_upper,pvals,b_Test,DW_p] = NeweyWestAdjust(r,(2001:2019)',d,1,1);
% end
% 
% b_Test-b%看检验的和计算的是不是一样的序列一样的结果
% out(i,:)=[b_Test(2),ci_lower(2),ci_upper(2),pvals(2),b_Test(1),DW_p];
% out2(i,:)=[b(2),bint(2,1),bint(2,2),stats(3),b(1),p_value];
% end


