%ERA5�������²�������׼������
%��Ҫԭʼһ֡ERA5����0.25ȫ�������������
function [STD_pre] = ERA5dailyresample(Ori_pre,Garea)
Ori_pre=[Ori_pre,Ori_pre(:,1)]; %ERA5ÿ����������������Ҫ�ı�׼�������ܣ���������Χ�ĸ������ֵ���㣬���ǵ�
STD_pre=zeros(360,720);
for i=1:360
    for j=1:720
        STD_pre(i,j)=(Ori_pre(i,j)*Garea(2*i-1,2*j-1) + Ori_pre(i+1,j)*Garea(2*i,2*j-1) + Ori_pre(i,j+1)*Garea(2*i-1,2*j) + Ori_pre(i+1,j+1)*Garea(2*i,2*j))/sum(sum(Garea(2*i-1:2*i,2*j-1:2*j)));
    end
end
end
