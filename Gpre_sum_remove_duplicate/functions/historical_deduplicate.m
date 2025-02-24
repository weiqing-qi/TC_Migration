function [PRE1,PRE2,PRE3,PRE4,PRE5,Hist_val] = historical_deduplicate(PRE1,PRE2,PRE3,PRE4,PRE5,yearBB,monthBB,dayBB,yearB,monthB,dayB,year,month,day,yearF,monthF,dayF,yearFF,monthFF,dayFF,Hist_val,WZ)
%HISTORICAL_DEDUPLICATE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

      L=Hist_val(:,WZ);
      L(L==0)=[];
      if ~isempty(find(L==str2double([yearBB,monthBB,dayBB]), 1)) %�����λ���ҵ���һ��  ע�� '~'
          PRE1=0;
      else
          Hist_val(length(L)+1,WZ)=str2double([yearBB,monthBB,dayBB]);
      end
%-------------------------------------------------------- 
      L=Hist_val(:,WZ);
      L(L==0)=[];
      if ~isempty(find(L==str2double([yearB,monthB,dayB]), 1)) %�����λ���ҵ���һ��  ע�� '~'
          PRE2=0;
      else
          Hist_val(length(L)+1,WZ)=str2double([yearB,monthB,dayB]);
      end
%--------------------------------------------------------  
      L=Hist_val(:,WZ);
      L(L==0)=[];
      if ~isempty(find(L==str2double([year,month,day]), 1)) %�����λ���ҵ���һ��  ע�� '~'
          PRE3=0;
      else
          Hist_val(length(L)+1,WZ)=str2double([year,month,day]);
      end  
%--------------------------------------------------------  
      L=Hist_val(:,WZ);
      L(L==0)=[];
      if ~isempty(find(L==str2double([yearF,monthF,dayF]), 1)) %�����λ���ҵ���һ��  ע�� '~'
          PRE4=0;
      else
          Hist_val(length(L)+1,WZ)=str2double([yearF,monthF,dayF]);
      end 
%--------------------------------------------------------  
      L=Hist_val(:,WZ);
      L(L==0)=[];
      if ~isempty(find(L==str2double([yearFF,monthFF,dayFF]), 1)) %�����λ���ҵ���һ��  ע�� '~'
          PRE5=0;
      else
          Hist_val(length(L)+1,WZ)=str2double([yearFF,monthFF,dayFF]);
      end
end

