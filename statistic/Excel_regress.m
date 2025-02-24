%回归，返回画图用的两个回归端点还有回归斜率和斜率上下限
DATA=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted.xlsx',3,'AT3:AX43');%改
drawdata=zeros(2,10);%五套数据 两点一线
PR_Regress=zeros(5,9);%五套数据 四个值
for i=1:5%五套数据 ERA5 JRA55 IMERG TRMM PERSIANN
X=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted.xlsx',3,'A3:A43');%改
Y=DATA(:,i);
X(isnan(Y))=[];%长度一致
Y(isnan(Y))=[];
if i == 4
    X=X(3:end);%长度一致
    Y=Y(3:end);
elseif i==3
    X=X(1:end-1);%长度一致
    Y=Y(1:end-1);   
end
if length(Y)==length(X)
    [b,bint,r,rint,stats]=regress(Y,[ones(length(X),1),X],0.05);
    PR_Regress(i,1)=b(2);
    PR_Regress(i,2)=bint(2,1);%注意
    PR_Regress(i,3)=bint(2,2);
    PR_Regress(i,4)=stats(3);

    [ci_lower,ci_upper,pvals,b_Test,DW_p] = NeweyWestAdjust(r,X,Y,1,1);
    PR_Regress(i,5)=b_Test(2);
    PR_Regress(i,6)=ci_lower(2);%注意
    PR_Regress(i,7)=ci_upper(2);
    PR_Regress(i,8)=pvals(2);
    PR_Regress(i,9)=DW_p;  
    switch i
        case 1
        drawdata(1,2*i-1)=1980;
        drawdata(2,2*i-1)=2020;
        drawdata(1,2*i)=b(2)*drawdata(1,2*i-1)+b(1);
        drawdata(2,2*i)=b(2)*drawdata(2,2*i-1)+b(1);
        case 2
        drawdata(1,2*i-1)=1980;
        drawdata(2,2*i-1)=2020;
        drawdata(1,2*i)=b(2)*drawdata(1,2*i-1)+b(1);
        drawdata(2,2*i)=b(2)*drawdata(2,2*i-1)+b(1);
        case 3
        drawdata(1,2*i-1)=2001;%改
        drawdata(2,2*i-1)=2019;%改
        drawdata(1,2*i)=b(2)*drawdata(1,2*i-1)+b(1);
        drawdata(2,2*i)=b(2)*drawdata(2,2*i-1)+b(1);         
        case 4
        drawdata(1,2*i-1)=2001;%改
        drawdata(2,2*i-1)=2019;%改
        drawdata(1,2*i)=b(2)*drawdata(1,2*i-1)+b(1);
        drawdata(2,2*i)=b(2)*drawdata(2,2*i-1)+b(1);          
        case 5    
        drawdata(1,2*i-1)=1983;%改
        drawdata(2,2*i-1)=2020;%改
        drawdata(1,2*i)=b(2)*drawdata(1,2*i-1)+b(1);
        drawdata(2,2*i)=b(2)*drawdata(2,2*i-1)+b(1);
    end
else
    disp('FALSE')
    break
end
end
for i=1:5
    for j=2:4
        if ~isnan(PR_Regress(i,j+4))
            PR_Regress(i,j)=PR_Regress(i,j+4);
        end
    end
end
%% FOR Migration.xlsx
DATA=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Migration.xlsx',1,'D152:AM159');%改
D=zeros(8,12);
for i=1:12
    D(:,i)=DATA(:,3*i-2);
end %读取所有28R所有数据按照basin分

X=(-750:100:1950);
D=D(1:8,:);%land or sea
X=X(1:8)';
Out=zeros(3,12);
for i=1:12
    [b,bint,r,rint,stats]=regress(D(:,i),[ones(length(D(:,i)),1),X],0.05);
    Out(1,i)=b(2);
    Out(2,i)=b(2)-bint(2,1);%注意
    Out(3,i)=bint(2,2)-b(2);%注意
end
