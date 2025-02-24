function D = SphereDist2_Matrix(x,y,R)
%SphereDist1当两点间距离很短时,比如相距几百米的两点,会导致较大的舍入误差,Haversine 方法能够避免该问题。单位KM
%根据两点的经纬度计算大圆距离(基于球面余弦公式)
%x为A点[经度, 纬度], y为B点[经度, 纬度],x,y输入nx2的数组

%鸣谢https://zhuanlan.zhihu.com/p/42948839?utm_source=qq
%鸣谢http://blog.sina.com.cn/s/blog_658a93570101hynw.html
%鸣谢https://blog.csdn.net/u011001084/article/details/52980834
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