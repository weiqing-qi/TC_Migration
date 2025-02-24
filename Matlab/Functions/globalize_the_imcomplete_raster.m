function [raster] = globalize_the_imcomplete_raster(header,ras)
%globalize the imcomplete ras(ter) ���������Сarctxt�����ȫ�� in-180180out0360 90NS
%update20201020
%header   %1--ncols  %2--nrows  %3--xllcorner  %4--yllcorner  %5--cellsize  %6--NODATA_value
ras(ras ~= header(6))=1;                %����һ�¶���������Ĥ
ras(ras == header(6))=0;                %һ��Ҫ���û�����ֵ���û�nodatavalue����Ȼ���ܺ�����ֵ�ظ�
raster=zeros(180/header(5),360/header(5)); %�ȶ��㵽90����Ҫ�ټ���


if header(3)<0                               %��Χ��-180-180,ת����0-360
  if (header(3)+ header(5)*header(1))>0      %��Խ���
     headeradjust=header(3)+0.001; %��ֹras xll������gridsizeʱ�������һ��
     xSTARTgrid=ceil(headeradjust/header(5));   %������ ������ buffer�������ȱʧ   
     xSTARTgrid=xSTARTgrid + 360/header(5);
     xENDgrid= xSTARTgrid + header(1) - 1- 360/header(5); 
     
     yENDgridgrid= 0-(ceil(header(4)/header(5))-1) + 90/header(5); %���ϵ��� �ӱ�����
     ySTARTgrid= yENDgridgrid- header(2)+ 1;
 
     raster(ySTARTgrid:yENDgridgrid, xSTARTgrid:360/header(5))=ras(:,1:((360/header(5))-xSTARTgrid +1));
     raster(ySTARTgrid:yENDgridgrid, 1:xENDgrid)=ras(:,(((360/header(5))-xSTARTgrid +1)+1):size(ras,2));
     
     %-------Ϊ�������һ�������λ�����,����Χ,  Ҫ�ģ�����λ��������
     weizhi=find(raster==1);
     for i=1:length(find(raster==1))
         if((weizhi(i)+size(raster,1)) <= size(raster,1)*size(raster,2))  %������Χ�ľͲ�������
             raster(weizhi(i)+size(raster,1))=1;
             
         end
         if ((weizhi(i)-size(raster,2)) > 0)
             raster(weizhi(i)-1)=1;
         end
     end
     
  else                                       %����Խ���
      headeradjust=header(3)-0.001; %��ֹras xll������gridsizeʱ�������һ��
     xSTARTgrid=floor(headeradjust/header(5));   %������ ������
     xSTARTgrid=xSTARTgrid + 360/header(5);
     xENDgrid= xSTARTgrid + header(1)-1;
     yENDgridgrid= 0-(ceil(header(4)/header(5))-1) + 90/header(5); %���ϵ��� �ӱ�����
     ySTARTgrid= yENDgridgrid- header(2)+ 1;

     raster(ySTARTgrid:yENDgridgrid,xSTARTgrid:xENDgrid)=ras;
     
     %-------Ϊ�������һ�������λ�����,����Χ,  Ҫ�ģ�����λ��������
     weizhi=find(raster==1);
     for i=1:length(find(raster==1))
         if((weizhi(i)+size(raster,1)) <= size(raster,1)*size(raster,2))  %������Χ�ľͲ�������
             raster(weizhi(i)+size(raster,1))=1;
             
         end
         if ((weizhi(i)-size(raster,2)) > 0)
             raster(weizhi(i)-1)=1;
         end
     end
     
  end
    
else                                    %header(3)>=0 ������0-360����-180-180��һ�� 
headeradjust=header(3)-0.001; %��ֹras xll������gridsizeʱ�������һ��
xSTARTgrid=ceil(headeradjust/header(5));   %������ ������
xENDgrid= xSTARTgrid + header(1) - 1; 

yENDgridgrid= 0-(ceil(header(4)/header(5))-1) + 90/header(5); %���ϵ��� �ӱ�����
ySTARTgrid= yENDgridgrid- header(2)+ 1;

raster(ySTARTgrid:yENDgridgrid,xSTARTgrid:xENDgrid)=ras;

     %-------Ϊ�������һ�������λ�����,����Χ,  Ҫ�ģ�����λ��������
     weizhi=find(raster==1);
     for i=1:length(find(raster==1))
         if((weizhi(i)+size(raster,1)) <= size(raster,1)*size(raster,2))  %������Χ�ľͲ�������
             raster(weizhi(i)+size(raster,1))=1;
             
         end
         if ((weizhi(i)-size(raster,2)) > 0)
             raster(weizhi(i)-1)=1;
         end
     end
end
end

