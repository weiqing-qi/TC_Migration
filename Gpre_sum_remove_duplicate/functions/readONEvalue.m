function [header1,header2,pre] = read_ARCascii(filepath,R,C,headersize)
%--------------读取用于arcgis的txt文档中的一个值------------update20210512
fid=fopen(filepath,'r','l');      %MATLAB跳过头文件读txt数据 不能用fread因为txt是ascii文档而非二进制
fseek(fid,[headersize+((R-1)*总列数+C-1)*]*位数+分隔符大小,'bof');%bof cof eof 文档开头 现在位置 结尾%不可行解

pre = fscanf(fid,'%f',[header2(1) Inf]);%跳过头文件后每次读一行，读到文档末尾
pre=pre';                           %上一步横向读取后竖向写入，需要转置
fclose(fid);

end
