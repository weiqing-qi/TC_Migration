%ERA5数据重新采样到标准网格上
%需要原始一帧ERA5，加0.5*CS全球网格面积数据
function [STD_pre] = ERA5FixResample(Ori_Grid,Garea2D_2,cs)
Ori_Grid=[Ori_Grid,Ori_Grid(:,1)]; %ERA5每个网格中心在我们要的标准网格四周，所以是周围四个网格插值计算，在右侧补上180°

STD_pre=zeros(180/cs,360/cs);
for i=1:180/cs
    for j=1:360/cs
        STD_pre(i,j)=(Ori_Grid(i,j)*Garea2D_2(2*i-1,2*j-1) + Ori_Grid(i+1,j)*Garea2D_2(2*i,2*j-1) + Ori_Grid(i,j+1)*Garea2D_2(2*i-1,2*j) + Ori_Grid(i+1,j+1)*Garea2D_2(2*i,2*j))/sum(Garea2D_2(2*i-1:2*i,2*j-1:2*j),"all");
    end
end

