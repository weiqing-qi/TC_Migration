function [] = OUTPUT(odir,name,header,P)
%OUTPUT  update20201020

fid=fopen([odir,name,'.txt'],'w');%加个【】括号不然识别不出
%----------------------------------------------------
%写头文件
[c,d]=size(header);            % 得到矩阵的行数和列数 
for a=1:c
  for b=1:d
  fprintf(fid,'%s',header(a,b));
  end
 fprintf(fid,'\r\n'); 
end
%---------------------------------------------------- 
%写文件
[e,f]=size(P);            % 得到矩阵的行数和列数
for m=1:e
  for n=1:f
  fprintf(fid,'%.4f\t',P(m,n));
  end
  fprintf(fid,'\r\n');%\n是换行，英文是New line，表示使光标到行首 \r是回车，英文是Carriage return，表示使光标下移一格。
end
fclose(fid);
end
