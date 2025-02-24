function D = SphereDist_Matrix(x,y,R)
%根据两点的经纬度计算大圆距离(基于球面余弦公式), 默认使用地球半径 单位KM
%x为A点[经度, 纬度], y为B点[经度, 纬度];x,y输入nx2的数组
%输入0-360或-180-180都可以

%鸣谢https://zhuanlan.zhihu.com/p/42948839?utm_source=qq
%鸣谢http://blog.sina.com.cn/s/blog_658a93570101hynw.html
%鸣谢https://blog.csdn.net/u011001084/article/details/52980834

if nargin < 3 
    R = 6371;%输入变量少于3个时默认使用地球半径
end

x = D2R(x);
y = D2R(y);
DeltaS = acos( cos(x(:,2)).*cos(y(:,2)).*cos(x(:,1)-y(:,1))+sin(x(:,2)).*sin(y(:,2)) );
D= R*DeltaS;
end

