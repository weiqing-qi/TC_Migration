function [raster] = globalize_the_imcomplete_raster(header,ras)
%globalize the imcomplete ras(ter) 任意格网大小arctxt处理成全球 in-180180out0360 90NS
%update20201020
%header   %1--ncols  %2--nrows  %3--xllcorner  %4--yllcorner  %5--cellsize  %6--NODATA_value
ras(ras ~= header(6))=1;                %整理一下读出来的掩膜
ras(ras == header(6))=0;                %一定要先置换数据值后置换nodatavalue，不然可能和数据值重复
raster=zeros(180/header(5),360/header(5)); %先都搞到90有需要再剪切


if header(3)<0                               %范围从-180-180,转换成0-360
  if (header(3)+ header(5)*header(1))>0      %跨越零度
     headeradjust=header(3)+0.001; %防止ras xll能整除gridsize时候左面多一格
     xSTARTgrid=ceil(headeradjust/header(5));   %从左到右 从西向东 buffer东侧可能缺失   
     xSTARTgrid=xSTARTgrid + 360/header(5);
     xENDgrid= xSTARTgrid + header(1) - 1- 360/header(5); 
     
     yENDgridgrid= 0-(ceil(header(4)/header(5))-1) + 90/header(5); %从上到下 从北向南
     ySTARTgrid= yENDgridgrid- header(2)+ 1;
 
     raster(ySTARTgrid:yENDgridgrid, xSTARTgrid:360/header(5))=ras(:,1:((360/header(5))-xSTARTgrid +1));
     raster(ySTARTgrid:yENDgridgrid, 1:xENDgrid)=ras(:,(((360/header(5))-xSTARTgrid +1)+1):size(ras,2));
     
     %-------为避免最大一个网格的位置误差,扩大范围,  要改！！！位置是竖着
     weizhi=find(raster==1);
     for i=1:length(find(raster==1))
         if((weizhi(i)+size(raster,1)) <= size(raster,1)*size(raster,2))  %超过范围的就不扩充了
             raster(weizhi(i)+size(raster,1))=1;
             
         end
         if ((weizhi(i)-size(raster,2)) > 0)
             raster(weizhi(i)-1)=1;
         end
     end
     
  else                                       %不跨越零度
      headeradjust=header(3)-0.001; %防止ras xll能整除gridsize时候左面多一格
     xSTARTgrid=floor(headeradjust/header(5));   %从左到右 从西向东
     xSTARTgrid=xSTARTgrid + 360/header(5);
     xENDgrid= xSTARTgrid + header(1)-1;
     yENDgridgrid= 0-(ceil(header(4)/header(5))-1) + 90/header(5); %从上到下 从北向南
     ySTARTgrid= yENDgridgrid- header(2)+ 1;

     raster(ySTARTgrid:yENDgridgrid,xSTARTgrid:xENDgrid)=ras;
     
     %-------为避免最大一个网格的位置误差,扩大范围,  要改！！！位置是竖着
     weizhi=find(raster==1);
     for i=1:length(find(raster==1))
         if((weizhi(i)+size(raster,1)) <= size(raster,1)*size(raster,2))  %超过范围的就不扩充了
             raster(weizhi(i)+size(raster,1))=1;
             
         end
         if ((weizhi(i)-size(raster,2)) > 0)
             raster(weizhi(i)-1)=1;
         end
     end
     
  end
    
else                                    %header(3)>=0 不管是0-360还是-180-180都一样 
headeradjust=header(3)-0.001; %防止ras xll能整除gridsize时候左面多一格
xSTARTgrid=ceil(headeradjust/header(5));   %从左到右 从西向东
xENDgrid= xSTARTgrid + header(1) - 1; 

yENDgridgrid= 0-(ceil(header(4)/header(5))-1) + 90/header(5); %从上到下 从北向南
ySTARTgrid= yENDgridgrid- header(2)+ 1;

raster(ySTARTgrid:yENDgridgrid,xSTARTgrid:xENDgrid)=ras;

     %-------为避免最大一个网格的位置误差,扩大范围,  要改！！！位置是竖着
     weizhi=find(raster==1);
     for i=1:length(find(raster==1))
         if((weizhi(i)+size(raster,1)) <= size(raster,1)*size(raster,2))  %超过范围的就不扩充了
             raster(weizhi(i)+size(raster,1))=1;
             
         end
         if ((weizhi(i)-size(raster,2)) > 0)
             raster(weizhi(i)-1)=1;
         end
     end
end
end

