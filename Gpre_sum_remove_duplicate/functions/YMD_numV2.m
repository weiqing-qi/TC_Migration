function [STR,NUM] = YMD_num(numdata)
%YMD_NUM 年月日转换为第几天
%   此处显示详细说明
Y=floor(numdata/10000);
M=floor((numdata-Y*10000)/100);
D=floor(numdata-Y*10000-M*100);

NUM=0;
if M~=1
    for i=1:M-1
        NUM=NUM+eomday(Y,i);
    end
end
NUM = NUM+D;
STR=num2str(NUM,'%03d');
end
