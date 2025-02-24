function [area1D,area3D] = Gridarea3D(cs,LEN)
%���������������һ�У�����ָ��ά�� km2
%���������С�����ȫ��������
if nargin == 1
    LEN = 1;
end

CA=(1:180/cs)';
area1D=6371^2 * cs * (pi()/180) * (sin( (90-(CA-1)*cs) * pi()/180) - sin( (90-(CA)*cs) * pi()/180));
area3D=repmat(area1D,[1 360/cs LEN]);

end
