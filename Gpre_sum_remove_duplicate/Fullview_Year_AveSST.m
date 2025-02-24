%����ERA5ÿһ��ƽ������
header=['ncols           720';
        'nrows           360';
        'xllcorner      -180';
        'yllcorner       -90';
        'cellsize        0.5';
        'NODATA_value  -9999'];
    
ERA5_idir='E:\DATA\ERA5\ERA5-Sea_surface_temperature\';
ERAfilename='ERA5_Sea_Surface_Temperature_on_single_levels_daily';
[Garea025,~] = Gridarea(0.25);
[Garea05,~] = Gridarea(0.5);
ACCSST=zeros(360,720);%�ۻ�ƽ�����¿ռ�ֲ�
AVENUM=0;%ÿ��ƽ�����¿ռ�ֲ�
BY=1966;
YEAR_AVE_SST=zeros(360,720,2020-BY+1);
for Year=BY:2020
    for mon=1:12
       STRdate=[num2str(Year'),num2str(mon,'%02d')];
       try
       Amonthdata = ncread([ERA5_idir,ERAfilename,STRdate(1:6),'.nc'],'sst');%ÿһ��allday��Ӧһ��allpre
       catch
       Amonthdata = ncread([ERA5_idir,ERAfilename,STRdate(1:6),'.nc'],'tos');%ÿһ��allday��Ӧһ��allpre
       end
       Amonthdata(Amonthdata<0)=missing;%ֻҪ��һ����ֵ�Ͳ�Ҫ��
       
           for day=1:eomday(Year,mon)
               %Tempo_sst=ERA5dailyresample( flipud(Amonthdata(:,:,str2double(STRdate(7:8)))')-273.15 ,Garea025); %��ȡ������ERA5ԭ���ݾ����ز�������׼����
               Tempo_sst=flipud(Amonthdata(:,:,day)')-273.15;%ת��Ϊ���϶�
               ACCSST=ACCSST+Tempo_sst(1:360,:);%ԭʼ���ݶ���һ��
               
               %������ƽ��
               if strcmp([STRdate(5:6),num2str(day,'%02d')],'1231')%�Ƚ��ַ���ÿ�����
                  [~,DofY] = YMD_num(Year*10000+1231);
                  YEAR_AVE_SST(:,:,Year-BY+1)=ACCSST/DofY;
                  ACCSST=zeros(360,720);
               end
               
           end
    end
    
end
a=zeros(360,720);
b=zeros(360,720);
for i=1:15
    a=a+YEAR_AVE_SST(:,:,i);
end
a=a/15;
for j=16:35
    b=b+YEAR_AVE_SST(:,:,j);
end
b=b/20;
OUTPUT('C:\Users\Dell\Desktop\',['ERA5_Ori_SST_0120-6600fullChange'],header,b-a);
