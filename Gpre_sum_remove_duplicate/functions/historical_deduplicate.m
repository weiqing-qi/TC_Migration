function [PRE1,PRE2,PRE3,PRE4,PRE5,Hist_val] = historical_deduplicate(PRE1,PRE2,PRE3,PRE4,PRE5,yearBB,monthBB,dayBB,yearB,monthB,dayB,year,month,day,yearF,monthF,dayF,yearFF,monthFF,dayFF,Hist_val,WZ)
%HISTORICAL_DEDUPLICATE 此处显示有关此函数的摘要
%   此处显示详细说明

      L=Hist_val(:,WZ);
      L(L==0)=[];
      if ~isempty(find(L==str2double([yearBB,monthBB,dayBB]), 1)) %在这个位置找到这一天  注意 '~'
          PRE1=0;
      else
          Hist_val(length(L)+1,WZ)=str2double([yearBB,monthBB,dayBB]);
      end
%-------------------------------------------------------- 
      L=Hist_val(:,WZ);
      L(L==0)=[];
      if ~isempty(find(L==str2double([yearB,monthB,dayB]), 1)) %在这个位置找到这一天  注意 '~'
          PRE2=0;
      else
          Hist_val(length(L)+1,WZ)=str2double([yearB,monthB,dayB]);
      end
%--------------------------------------------------------  
      L=Hist_val(:,WZ);
      L(L==0)=[];
      if ~isempty(find(L==str2double([year,month,day]), 1)) %在这个位置找到这一天  注意 '~'
          PRE3=0;
      else
          Hist_val(length(L)+1,WZ)=str2double([year,month,day]);
      end  
%--------------------------------------------------------  
      L=Hist_val(:,WZ);
      L(L==0)=[];
      if ~isempty(find(L==str2double([yearF,monthF,dayF]), 1)) %在这个位置找到这一天  注意 '~'
          PRE4=0;
      else
          Hist_val(length(L)+1,WZ)=str2double([yearF,monthF,dayF]);
      end 
%--------------------------------------------------------  
      L=Hist_val(:,WZ);
      L(L==0)=[];
      if ~isempty(find(L==str2double([yearFF,monthFF,dayFF]), 1)) %在这个位置找到这一天  注意 '~'
          PRE5=0;
      else
          Hist_val(length(L)+1,WZ)=str2double([yearFF,monthFF,dayFF]);
      end
end

