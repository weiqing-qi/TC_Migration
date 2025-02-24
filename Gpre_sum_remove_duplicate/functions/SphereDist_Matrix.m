function D = SphereDist_Matrix(x,y,R)
%��������ľ�γ�ȼ����Բ����(�����������ҹ�ʽ), Ĭ��ʹ�õ���뾶 ��λKM
%xΪA��[����, γ��], yΪB��[����, γ��];x,y����nx2������
%����0-360��-180-180������

%��лhttps://zhuanlan.zhihu.com/p/42948839?utm_source=qq
%��лhttp://blog.sina.com.cn/s/blog_658a93570101hynw.html
%��лhttps://blog.csdn.net/u011001084/article/details/52980834

if nargin < 3 
    R = 6371;%�����������3��ʱĬ��ʹ�õ���뾶
end

x = D2R(x);
y = D2R(y);
DeltaS = acos( cos(x(:,2)).*cos(y(:,2)).*cos(x(:,1)-y(:,1))+sin(x(:,2)).*sin(y(:,2)) );
D= R*DeltaS;
end

