function [DATE] = num2YMD(numS,SY)
%第几天换为年月日转 numS:第几天:INPUT:DOUBLE
%   numS(年月日的STR)：19990101

normyear=[0,31,59,90,120,151,181,212,243,273,304,334,365];
leapyear=[0,31,60,91,121,152,182,213,244,274,305,335,366];
if eomday(SY,2)==28   
    M = find(normyear>=numS,1,'first')-1;
    D =numS-normyear(M);
else
    M = find(leapyear>=numS,1,'first')-1;
    D =numS-leapyear(M);   
end
DATE=[num2str(SY),num2str(M,'%02d'),num2str(D,'%02d')];






