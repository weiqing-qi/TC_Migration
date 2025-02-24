function [ci_lower,ci_upper,pvals,b,DW_p] = NeweyWestAdjust(e,X,y,L,constant)
% ����NeweyWest�������Իع�Ķ��������ԣ��õ��ɿ������Իع����
% ����NeweyWest�����õ��ľ��������ı�׼���һ������
% �ο�GPT��https://stackoverflow.com/questions/36876141/newey-west-p-values-in-matlab
% https://www.stata.com/manuals13/tsnewey.pdf �е�ֵ������֤
% https://zhuanlan.zhihu.com/p/54913149
%---------------------------------------------------
% where: b = k MLR coefficients �ڲ�֪���췽�������ؾ����ǲ���OLR
%        e = T x 1 vector of model residuals
%        X = T x k matrix of independant variables
%        L = lag length to use (Default: Newey-West(1994) plug-in
%        procedure)

%        constant = 0: no constant to be added;
%                 = 1: constant term to be added (Default = 1)
%
%        nwse = Newey-West standard errors
%------------------------
%-----------Durbin-Watson test�ǲ��Ǿ���AR��1������--------------
mdl = fitlm(X,y);
DW_p = dwtest(mdl);
b=mdl.Coefficients.Estimate;
if DW_p < 0.05 %����lag-1�����
    % -------------------���� NW-adjusted ͳ����-----------------
    nwse = NeweyWest(e,X,L,constant); % ���������ı�׼��ע���׼��ͱ�׼�������

    indexxx = sum(isnan(X),2)==0;%������ж�ÿһ���ǲ����п�ֵ��Xֻѡ�������ж����ڵ��н�������
    X = X(indexxx,:);
    [n,k] = size(X);

    t = b ./ nwse; %b�ǲ������ƣ�t��tͳ����
    pvals = NaN(size(b));  %�����NW�����õ���Pֵ
    df_r = n - k;
    pvals(t >= 0) = 2 * (1 - tcdf(t(t>=0), df_r));%df_r�ǲв����ɶȣ�df_r= n����������- k�������ͳ�Ʋ�������
    pvals(t < 0)  = 2 * tcdf(t(t<0), df_r);
    alpha = 0.05; % ����ˮƽ
    t_critical = tinv(1 - alpha/2, df_r); % ���ɶ�Ϊ n - 2 �� t �ֲ����ٽ�ֵ����Ӧdf_r
    ci_lower = b - t_critical * nwse;
    ci_upper = b + t_critical * nwse;
else
    ci_lower = NaN(size(b));
    ci_upper = NaN(size(b));
    pvals = NaN(size(b));
end

