function [ORing_Mask,F_Mask,RB_Mask] = OuterRainBelt(thres,Garea2D,Dailydata,Buff_OMask,Buff_OWZ,Buff_IWZ,cs)
%计算外部雨带的掩膜，如果不计算的话S_RB_Buff=0只返回外部环状掩膜ORing_Mask
%按照阈值做一个filter(不然可能连接太多，阈值需要确认)
BW_Ddata=Dailydata;                                             %这里dailydata没有filter，处理的是BW和以后的连通
BW_Ddata(Dailydata< thres)=0;                                 
BW_Ddata(Dailydata>=thres)=1;
%计算与buffer连通的降水带
cc = bwconncomp(BW_Ddata,4);                                    %计算连通分量，用的是4连通，默认8连通
LableM = labelmatrix(cc);                                       %用标签矩阵标注连通分量 背景值：0                                                                     
LableLink2Buff=unique(LableM(Buff_OWZ));
LableLink2Buff(LableLink2Buff==0)=[];                           %筛选与buffer有交集的标签
%如果整个Buff_OMask范围里面没有降水
if isempty(LableLink2Buff)
    F_Mask=Buff_OMask;
    ORing_Mask=F_Mask;
    ORing_Mask(Buff_IWZ)=0;
    RB_Mask=[];
    return
end
%计算质心并删除质心在Buffer以外的降水带
Fullmask=ismember(LableM,LableLink2Buff);
Mmask=labelmatrix(bwconncomp(Fullmask,4));                      %Fullmusk逻辑数组也可以作为二值图像输入bwconncomp
WCen = regionprops( ...                                         %注意这里一定要用bwconncomp4连通，不然regionprops默认8连通会导致不对齐
    bwconncomp(Mmask,4),Garea2D.*Dailydata,'WeightedCentroid');   %改 这里也可以只考虑面积或data        
WCen = cat(1,WCen.WeightedCentroid);                            %第一列：列坐标x，第二列行坐标y
CenWZ=coord2grid_M(WCen(:,1)*cs-180,90-WCen(:,2)*cs,cs);
LWCinB=(1:length(LableLink2Buff))';
LWCinB=LWCinB(ismember(CenWZ,Buff_OWZ));                        %CenWZ 和LWCinB 需要一样
Mmask(~ismember(Mmask,LWCinB))=0;
%生成最终结果
F_Mask=Mmask;
F_Mask(F_Mask>0)=1;
RB_Mask=F_Mask;
F_Mask(Buff_OWZ)=1;
ORing_Mask=F_Mask;
ORing_Mask(Buff_IWZ)=0; 
%% 图片检验
% imshow(label2rgb(F_Mask,@copper,"c","shuffle"))
% hold on
% [R,C]=find(Buff_OMask==1);
% plot(C,R,"o",'MarkerEdgeColor','r',"MarkerSize",1)
%% --------检测-------------
%经过检测,当两者使用同一个C时，regionprops按照L编码顺序返回中心。       
% A=[0,0,1;0,1,0;0,0,0;1,1,1]
% B=[1,1,1;1,1,1;1,1,1;1,1,1]
% C=bwconncomp(A,4);
% L=labelmatrix(C) 
% WC=regionprops( C,B,'WeightedCentroid');
% cat(1,WC.WeightedCentroid)
end

