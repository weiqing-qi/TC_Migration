function [header1,header2,pre] = read_ARCascii(filepath,R,C,headersize)
%--------------��ȡ����arcgis��txt�ĵ��е�һ��ֵ------------update20210512
fid=fopen(filepath,'r','l');      %MATLAB����ͷ�ļ���txt���� ������fread��Ϊtxt��ascii�ĵ����Ƕ�����
fseek(fid,[headersize+((R-1)*������+C-1)*]*λ��+�ָ�����С,'bof');%bof cof eof �ĵ���ͷ ����λ�� ��β%�����н�

pre = fscanf(fid,'%f',[header2(1) Inf]);%����ͷ�ļ���ÿ�ζ�һ�У������ĵ�ĩβ
pre=pre';                           %��һ�������ȡ������д�룬��Ҫת��
fclose(fid);

end
