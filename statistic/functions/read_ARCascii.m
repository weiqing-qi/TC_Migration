function [header1,header2,pre] = read_ARCascii(filepath)
%--------------��ȡ����arcgis��txt�ĵ�-------------------update20201020
[header1,header2]=textread(filepath,'%s %f',6); %delimiter���Ϸ��������
       
fid=fopen(filepath,'r','l');      %MATLAB����ͷ�ļ���txt���� ������fread��Ϊtxt��ascii�ĵ����Ƕ�����

for i=1:6                         %skip first 6 lines
   line=fgetl(fid);
end

pre = fscanf(fid,'%f',[header2(1) Inf]);%����ͷ�ļ���ÿ�ζ�һ�У������ĵ�ĩβ
pre=pre';                           %��һ�������ȡ������д�룬��Ҫת��
fclose(fid);

end

