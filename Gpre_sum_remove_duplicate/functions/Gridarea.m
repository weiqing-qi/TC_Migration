function [area,C] = Gridarea(cs)
%��������������� 0360 km2
%���������С�����ȫ��������
area=zeros(180/cs,360/cs);
C=zeros(180/cs,1);

   for i=1:180/cs
       C(i)= 6371^2 * cs * (pi()/180) * (sin( (90-(i-1)*cs) * pi()/180) - sin( (90-(i)*cs) * pi()/180));
   end

   for j=1:360/cs
       area(:,j)=C;
   end

end

