function [area1D,area3D] = Gridarea3D(cs,LEN)
%地球网格面积计算一列，可以指定维度 km2
%输入网格大小，输出全球格网面积
if nargin == 1
    LEN = 1;
end

CA=(1:180/cs)';
area1D=6371^2 * cs * (pi()/180) * (sin( (90-(CA-1)*cs) * pi()/180) - sin( (90-(CA)*cs) * pi()/180));
area3D=repmat(area1D,[1 360/cs LEN]);

end
