%��һ�����ڼ���3��ǰ
%��������������������
%�������������������
function [d3_ago] = threedayago(strdata)
Y=str2double(strdata(1:4));
M=str2double(strdata(5:6));
D=str2double(strdata(7:8));

if strcmp(strdata(5:8),'0101')
    d3_ago=[num2str(Y-1),'1229'];
    return
end
if strcmp(strdata(5:8),'0102')
    d3_ago=[num2str(Y-1),'1230'];
    return
end
if strcmp(strdata(5:8),'0103')
    d3_ago=[num2str(Y-1),'1231'];
    return
end
%----not switch year
NUM=0;
if M~=1
    for i=1:M-1
        NUM=NUM+eomday(Y,i);
    end
end
NUM = NUM+D-3;
d3_ago= num2YMD(NUM,Y);




