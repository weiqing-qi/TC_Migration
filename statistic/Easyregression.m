%A为任意一个单列或多列的年序列(NAN填充)
A_in=[];
A_in(A_in==0)=missing;
[~,C]=size(A_in);
PR_Regress=zeros(C,10);
PR_MK=zeros(C,4);
for i=1:C
    A=A_in(:,i);
    A(isnan(A))=[];
    R=length(A);
[b,bint,r,rint,stats]=regress(A,[ones(R,1),(1:R)'],0.05);
    PR_Regress(i,1)=b(2);
    PR_Regress(i,2)=bint(2,1);%注意
    PR_Regress(i,3)=bint(2,2);
    PR_Regress(i,4)=stats(3); 
    
[ci_lower,ci_upper,pvals,b_Test,DW_p] = NeweyWestAdjust(r,(1:R)',A,1,1);
    PR_Regress(i,5)=b_Test(2);
    PR_Regress(i,6)=ci_lower(2);%注意
    PR_Regress(i,7)=ci_upper(2);
    PR_Regress(i,8)=pvals(2);
    PR_Regress(i,9)=DW_p;
    
 [Zs, p_value, UFk, UBk2, beta, beta_CI, corr_val]= MKtrend((1:R)',A);
    PR_MK(i,1)=p_value;
    PR_MK(i,2)=beta;
    PR_MK(i,3)=beta_CI(1);
    PR_MK(i,4)=beta_CI(2);
end

%DW自相关检验以及Nw修正
for i=1:C
    for j=2:4
        if ~isnan(PR_Regress(i,j+4))
            PR_Regress(i,j)=PR_Regress(i,j+4);
        end
    end
end
stats(3)
%筛选之后进行每年平均
Y_D=[];
A=[];
for MOD=1:10
    Y=Y_D(:,MOD*2-1);
    D=Y_D(:,MOD*2);

    D(Y==0)=[];
    Y(Y==0)=[];

    Ylim=[1980 2014];
    yAVE=[];
    for y=Ylim(1):Ylim(2)
        YM=mean(D(Y==y),"omitmissing");
        yAVE=[yAVE;YM];
    end
    A=[A,yAVE(~isnan(yAVE))];
end
