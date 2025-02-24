%ERA5�������²�������׼������
%��Ҫԭʼһ֡ERA5����0.5*CSȫ�������������
function [STD_pre] = ERA5FixResample(Ori_Grid,Garea2D_2,cs)
Ori_Grid=[Ori_Grid,Ori_Grid(:,1)]; %ERA5ÿ����������������Ҫ�ı�׼�������ܣ���������Χ�ĸ������ֵ���㣬���Ҳಹ��180��

STD_pre=zeros(180/cs,360/cs);
for i=1:180/cs
    for j=1:360/cs
        STD_pre(i,j)=(Ori_Grid(i,j)*Garea2D_2(2*i-1,2*j-1) + Ori_Grid(i+1,j)*Garea2D_2(2*i,2*j-1) + Ori_Grid(i,j+1)*Garea2D_2(2*i-1,2*j) + Ori_Grid(i+1,j+1)*Garea2D_2(2*i,2*j))/sum(Garea2D_2(2*i-1:2*i,2*j-1:2*j),"all");
    end
end

