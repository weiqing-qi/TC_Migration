function [KEEP_NO,KEEP_Y,KEEP_B,KEEP_i] = FilterNO(Ylim,ISO_TIME,USA_WIND,BASIN,dist)
%对整场气旋做Filter给出对应年份之间合适的编号
% TC_idir='D:\DATA\Segmentation map tutorial\NEW\IBTrACS.ALL.v04r00.nc';
% [SID,LAT,LON,ISO_TIME,USA_WIND,USA_RMW,BASIN]=IBTrACS_nc_entire_variable_r(TC_idir);
% dist=[];
% Ylim=[1981 2020]; dist=dist0;
NOlim = lim_Y2NO(Ylim,ISO_TIME);
NO = (NOlim(1):NOlim(2))';  %在文件里的位置编号
KEEP_i=[];
Y=zeros(NOlim(3),1);
B=cell(NOlim(3),1);


for i=1:NOlim(3)
    time=ISO_TIME(:,:,i+NOlim(1)-1);
    Y(i)=str2double(time(1:4,1));   
  
    B(i)=cellstr(BASIN(1:2,1,i+NOlim(1)-1)');
    W = USA_WIND(:,i+NOlim(1)-1);
    W(isnan(W))=[];
    MW=max(W);
    %if ~isempty(W) & MW>34 & dist(i)<2000 & dist(i)>0
    if ~isnan(dist(i))
        KEEP_i=[KEEP_i;i];
    end
end
KEEP_NO=NO(KEEP_i);
KEEP_Y=Y(KEEP_i);
KEEP_B=B(KEEP_i);
KEEP_D=dist(KEEP_i);

yAVE=[];
for y=Ylim(1):Ylim(2)
    YM=mean(KEEP_D(KEEP_Y==y));
    yAVE=[yAVE;YM];
end
X=(Ylim(1):Ylim(2))';
[b,bint,r,rint,stats]=regress(yAVE,[ones(length(X),1),X],0.05);
%A=[];随便输入一个序列即可做回归
%[b,bint,r,rint,stats]=regress(A,[ones(length(A),1),(1:length(A))'],0.05);
disp(b(2))
disp(stats(3))

end