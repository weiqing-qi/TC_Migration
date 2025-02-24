function [CC,RMSE,ME,MAE,BIAS,ABIAS,POD,FAR,CSI]=CC_RMSE_ME_MAE_BIAS_ABIAS(INPUT,F,M,H)
%输入一个两列的矩阵，x在第一列是地面数据y在第二列是卫星数据
if nargin < 2
    F=0.001; M=0.001; H=0.001;
end

GI=INPUT(:,1);
SI=INPUT(:,2);
GBAR=mean(GI,1);
SBAR=mean(SI,1);

S_G=zeros(length(GI),1);
for i=1:length(GI)
S_G(i)=SI(i)-GI(i);
end

SG_ABS=zeros(length(GI),1);
for i=1:length(GI)
SG_ABS(i)=abs(SI(i)-GI(i));
end

GI_GB=zeros(length(GI),1);
for i=1:length(GI)
GI_GB(i)=GI(i)-GBAR;
end

SI_SB=zeros(length(GI),1);
for i=1:length(GI)
SI_SB(i)=SI(i)-SBAR;
end
%--------------------------------------------------------------------------
%CC
G_GBXS_SB=zeros(length(GI),1);
for i=1:length(GI)
G_GBXS_SB(i)=GI_GB(i)*SI_SB(i);
end

ccu=sum(G_GBXS_SB);

GI_GBsqrt=zeros(length(GI),1);
for i=1:length(GI)
GI_GBsqrt(i)=(GI_GB(i))^2;
end

SI_SBsqrt=zeros(length(GI),1);
for i=1:length(GI)
SI_SBsqrt(i)=(SI_SB(i))^2;
end

ccd=(sum(GI_GBsqrt))^0.5*(sum(SI_SBsqrt))^0.5;
CC=ccu/ccd;
%--------------------------------------------------------------------------
%RMSE
S_Gsqrt=zeros(length(GI),1);
for i=1:length(GI)
S_Gsqrt(i)=(S_G(i))^2;
end
RMSE=(mean(S_Gsqrt))^0.5;
%--------------------------------------------------------------------------
%ME&MAE
ME=mean(S_G);
MAE=mean(SG_ABS);
%--------------------------------------------------------------------------
%BIAS&ABIAS
BIAS=(sum(S_G))/(sum(GI));
ABIAS=(sum(SG_ABS))/(sum(GI));
%--------------------------------------------------------------------------
POD=H/(H+M);
FAR=F/(H+F);
CSI=H/(H+M+F);
