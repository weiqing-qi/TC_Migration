function [NUM] = YMD_num(Y,M,D)
%YMD_NUM ������ת��Ϊ����ȵڼ���

% �˴���ֱ��ʹ������������Ϊ�����յ����
% Y=floor(numdata/10000);
% M=floor((numdata-Y*10000)/100);
% D=floor(numdata-Y*10000-M*100);
if M>12 || M<1 || D>eomday(Y,M) || D<0 || Y<0
    disp('ERROR DATE')
    return
end

NUM=0;
if M~=1
    for i=1:M-1
        NUM=NUM+eomday(Y,i);
    end
end
NUM = NUM+D;
% STR=num2str(NUM,'%03d');
end
