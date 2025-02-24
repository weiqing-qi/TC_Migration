function nwse = NeweyWest(e,X,L,constant)
% PURPOSE: computes Newey-West adjusted heteroscedastic-serial
%          consistent standard errors
%
%        Adapted by Guillaume Nolin from the original code by Ian Gow
%---------------------------------------------------
% where: e = T x 1 vector of model residuals
%        X = T x k matrix of independant variables
%        L = lag length to use (Default: Newey-West(1994) plug-in
%        procedure)

%        constant = 0: no constant to be added;
%                 = 1: constant term to be added (Default = 1)
%
%        nwse = Newey-West standard errors
%---------------------------------------------------

%% Variables

if nargin < 4 || constant ~= 0
    constant = 1;
end

if ~exist('X','var') && constant ~= 0
    X=ones(size(e,1),1);
end

indexxx = sum(isnan(X),2)==0;%这个是判断每一行是不是有空值，X只选择所有列都存在的行进行运算
X = X(indexxx,:);
e = e(indexxx,:);

[N,k] = size(X);

if nargin < 3 || L < 0
    % Newey-West (1994) plug-in procedure
    L = floor(4*((N/100)^(2/9)));
end
    
if any(all(X==1,1),2)
    constant=0;
end

if constant == 1
    k = k+1;
    X = [ones(N,1),X];
end

%% Computation

Q = 0;
for l = 0:L
    w_l = 1-l/(L+1);
    for t = l+1:N
        if (l==0)   % This calculates the S_0 portion
            Q = Q  + e(t) ^2 * X(t, :)' * X(t,:);
        else        % This calculates the off-diagonal terms
            Q = Q + w_l * e(t) * e(t-l)* ...
                (X(t, :)' * X(t-l,:) + X(t-l, :)' * X(t,:));
        end
    end
end
Q = (1/(N-k)) .*Q;

nwse = sqrt(diag(N.*((X'*X)\Q/(X'*X))));%这个是NW方法修正后的标准误，注意标准误和标准差不一样
% 标准误是协方差矩阵的对角元素开方，每个数字对应一个输入的x参数，注意还有一个常数项

end
%参考GPT和https://stackoverflow.com/questions/36876141/newey-west-p-values-in-matlab
% https://www.stata.com/manuals13/tsnewey.pdf 中的值进行验证
% https://zhuanlan.zhihu.com/p/54913149
%-------------------计算 NW-adjusted 统计量-----------------
% t = b ./ nwse; %b是参数估计，t是t统计量
% pvals = NaN(size(b));  %这个是NW方法得到的P值
% df_r = n - 2;
% pvals(t >= 0) = 2 * (1 - tcdf(t(t>=0), df_r));%df_r是残差自由度，df_r= n（样本数）- k（计算的统计参数量）
% pvals(t < 0)  = 2 * tcdf(t(t<0), df_r);
% alpha = 0.05; % 置信水平
% t_critical = tinv(1 - alpha/2, df_r); % 自由度为 n - 2 的 t 分布的临界值，对应df_r
% ci_lower = b - t_critical * nwse;
% ci_upper = b + t_critical * nwse;