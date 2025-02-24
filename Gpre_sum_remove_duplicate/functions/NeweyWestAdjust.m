function [ci_lower,ci_upper,pvals,b,DW_p] = NeweyWestAdjust(e,X,y,L,constant)
% 利用NeweyWest调整线性回归的多阶自相关性，得到可靠的线性回归参数
% 基于NeweyWest函数得到的经过调整的标准误进一步计算
% 参考GPT和https://stackoverflow.com/questions/36876141/newey-west-p-values-in-matlab
% https://www.stata.com/manuals13/tsnewey.pdf 中的值进行验证
% https://zhuanlan.zhihu.com/p/54913149
%---------------------------------------------------
% where: b = k MLR coefficients 在不知道异方差和自相关矩阵是采用OLR
%        e = T x 1 vector of model residuals
%        X = T x k matrix of independant variables
%        L = lag length to use (Default: Newey-West(1994) plug-in
%        procedure)

%        constant = 0: no constant to be added;
%                 = 1: constant term to be added (Default = 1)
%
%        nwse = Newey-West standard errors
%------------------------
%-----------Durbin-Watson test是不是具有AR（1）过程--------------
mdl = fitlm(X,y);
DW_p = dwtest(mdl);
b=mdl.Coefficients.Estimate;
if DW_p < 0.05 %具有lag-1自相关
    % -------------------计算 NW-adjusted 统计量-----------------
    nwse = NeweyWest(e,X,L,constant); % 经过调整的标准误，注意标准误和标准差的区别

    indexxx = sum(isnan(X),2)==0;%这个是判断每一行是不是有空值，X只选择所有列都存在的行进行运算
    X = X(indexxx,:);
    [n,k] = size(X);

    t = b ./ nwse; %b是参数估计，t是t统计量
    pvals = NaN(size(b));  %这个是NW方法得到的P值
    df_r = n - k;
    pvals(t >= 0) = 2 * (1 - tcdf(t(t>=0), df_r));%df_r是残差自由度，df_r= n（样本数）- k（计算的统计参数量）
    pvals(t < 0)  = 2 * tcdf(t(t<0), df_r);
    alpha = 0.05; % 置信水平
    t_critical = tinv(1 - alpha/2, df_r); % 自由度为 n - 2 的 t 分布的临界值，对应df_r
    ci_lower = b - t_critical * nwse;
    ci_upper = b + t_critical * nwse;
else
    ci_lower = NaN(size(b));
    ci_upper = NaN(size(b));
    pvals = NaN(size(b));
end

