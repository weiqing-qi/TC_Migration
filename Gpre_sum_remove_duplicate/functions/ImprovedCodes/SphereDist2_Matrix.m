function D = SphereDist2_Matrix(x,y,R)
%SphereDist1����������ܶ�ʱ,������༸���׵�����,�ᵼ�½ϴ���������,Haversine �����ܹ���������⡣��λKM
%��������ľ�γ�ȼ����Բ����(�����������ҹ�ʽ)
%xΪA��[����, γ��], yΪB��[����, γ��],x,y����nx2������

%��лhttps://zhuanlan.zhihu.com/p/42948839?utm_source=qq
%��лhttp://blog.sina.com.cn/s/blog_658a93570101hynw.html
%��лhttps://blog.csdn.net/u011001084/article/details/52980834
if nargin < 3
    R = 6371;
end

x = D2R(x);
y = D2R(y);
h = HaverSin(abs(x(:,2)-y(:,2)))+cos(x(:,2)).*cos(y(:,2)).*HaverSin(abs(x(:,1)-y(:,1)));
D = 2 * R * asin(sqrt(h));

function h = HaverSin(theta)
    h=sin(theta./2).^2;
end

end