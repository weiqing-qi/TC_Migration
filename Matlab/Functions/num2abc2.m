function str =num2abc2(num) 
%将 10进制列数转换为26进制字母excel表列名
string={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

str =[];
while num>0
    m=mod(num,26);%返回整数
    if m==0 
        m=26;
    end
    str = [str string{m}];
    num=(num-m)/26;
end
str = fliplr(str);
end


