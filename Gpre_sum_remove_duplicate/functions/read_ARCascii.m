function [header1,header2,pre] = read_ARCascii(filepath)
%--------------读取用于arcgis的txt文档-------------------update20201020
[header1,header2]=textread(filepath,'%s %f',6); %delimiter加上反而会出错
       
fid=fopen(filepath,'r','l');      %MATLAB跳过头文件读txt数据 不能用fread因为txt是ascii文档而非二进制

for i=1:6                         %skip first 6 lines
   line=fgetl(fid);
end

pre = fscanf(fid,'%f',[header2(1) Inf]);%跳过头文件后每次读一行，读到文档末尾
pre=pre';                           %上一步横向读取后竖向写入，需要转置
fclose(fid);

end

