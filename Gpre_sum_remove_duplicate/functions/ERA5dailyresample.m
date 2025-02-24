%ERA5数据重新采样到标准网格上
%需要原始一帧ERA5，加0.25全球网格面积数据
function [STD_pre] = ERA5dailyresample(Ori_pre,Garea)
Ori_pre=[Ori_pre,Ori_pre(:,1)]; %ERA5每个网格中心在我们要的标准网格四周，所以是周围四个网格插值计算，考虑到
STD_pre=zeros(360,720);
for i=1:360
    for j=1:720
        STD_pre(i,j)=(Ori_pre(i,j)*Garea(2*i-1,2*j-1) + Ori_pre(i+1,j)*Garea(2*i,2*j-1) + Ori_pre(i,j+1)*Garea(2*i-1,2*j) + Ori_pre(i+1,j+1)*Garea(2*i,2*j))/sum(sum(Garea(2*i-1:2*i,2*j-1:2*j)));
    end
end
end
