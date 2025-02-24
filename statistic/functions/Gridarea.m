function [area,C] = Gridarea(cs)
%地球网格面积计算 0360 km2
%输入网格大小，输出全球格网面积
area=zeros(180/cs,360/cs);
C=zeros(180/cs,1);

   for i=1:180/cs
       C(i)= 6371^2 * cs * (pi()/180) * (sin( (90-(i-1)*cs) * pi()/180) - sin( (90-(i)*cs) * pi()/180));
   end

   for j=1:360/cs
       area(:,j)=C;
   end

end

