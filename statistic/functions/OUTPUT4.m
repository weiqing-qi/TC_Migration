function [] = OUTPUT(odir,name,header,P)
%OUTPUT  update20201020

fid=fopen([odir,name,'.txt'],'w');%�Ӹ��������Ų�Ȼʶ�𲻳�
%----------------------------------------------------
%дͷ�ļ�
[c,d]=size(header);            % �õ���������������� 
for a=1:c
  for b=1:d
  fprintf(fid,'%s',header(a,b));
  end
 fprintf(fid,'\r\n'); 
end
%---------------------------------------------------- 
%д�ļ�
[e,f]=size(P);            % �õ����������������
for m=1:e
  for n=1:f
  fprintf(fid,'%.4f\t',P(m,n));
  end
  fprintf(fid,'\r\n');%\n�ǻ��У�Ӣ����New line����ʾʹ��굽���� \r�ǻس���Ӣ����Carriage return����ʾʹ�������һ��
end
fclose(fid);
end
