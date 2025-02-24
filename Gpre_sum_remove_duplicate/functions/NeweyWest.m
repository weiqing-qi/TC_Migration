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

indexxx = sum(isnan(X),2)==0;%������ж�ÿһ���ǲ����п�ֵ��Xֻѡ�������ж����ڵ��н�������
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

nwse = sqrt(diag(N.*((X'*X)\Q/(X'*X))));%�����NW����������ı�׼��ע���׼��ͱ�׼�һ��
% ��׼����Э�������ĶԽ�Ԫ�ؿ�����ÿ�����ֶ�Ӧһ�������x������ע�⻹��һ��������

end
%�ο�GPT��https://stackoverflow.com/questions/36876141/newey-west-p-values-in-matlab
% https://www.stata.com/manuals13/tsnewey.pdf �е�ֵ������֤
% https://zhuanlan.zhihu.com/p/54913149
%-------------------���� NW-adjusted ͳ����-----------------
% t = b ./ nwse; %b�ǲ������ƣ�t��tͳ����
% pvals = NaN(size(b));  %�����NW�����õ���Pֵ
% df_r = n - 2;
% pvals(t >= 0) = 2 * (1 - tcdf(t(t>=0), df_r));%df_r�ǲв����ɶȣ�df_r= n����������- k�������ͳ�Ʋ�������
% pvals(t < 0)  = 2 * tcdf(t(t<0), df_r);
% alpha = 0.05; % ����ˮƽ
% t_critical = tinv(1 - alpha/2, df_r); % ���ɶ�Ϊ n - 2 �� t �ֲ����ٽ�ֵ����Ӧdf_r
% ci_lower = b - t_critical * nwse;
% ci_upper = b + t_critical * nwse;