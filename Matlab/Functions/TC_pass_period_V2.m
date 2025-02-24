%生成每个网格经过日期索引(BF+1)搜索日降水数据buffer外连通的降水雨带
%质心在一定范围内 input -180-180
%function: 1 bwconncomp 2 labelmatrix 3 regionprops bwmorph bwperim
%由于轨迹点经过加密，经纬度和时间数组无法对应
%每日寻找passtime最后再整合为一个位置和时间的索引元胞数组
function [RB_WZ4,passyear,passmonth,passday] = TC_pass_period_V2(alldata,Buff_WZ,Buff_Mask,cs,)
LG=length(Buff_WZ);
先找二维降水的联通













end