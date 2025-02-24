function [dayB,monthB,yearB,dayF,monthF,yearF] = daybackdayforward(day,month,year) 
%YDAY2MDAY 计算前一天后一天的日期  update20201020
%输入double返回str
%返回一个月的天数eomday

%----------------------------------------------跨年
if day ==1 && month==1
        dayB='31';monthB='12';yearB=num2str(year-1);dayF='02';monthF='01';yearF=num2str(year);
        return
end

if day==31 && month==12
    dayB='30';monthB='12';yearB=num2str(year);dayF='01';monthF='01';yearF=num2str(year+1);
        return
end

%----------------------------------------------跨月     
if day ==1 
        dayB=num2str(eomday(year,month-1));
        dayF='02';
        monthB=num2str(month-1); 
        if length(monthB)==1  
            monthB=['0',monthB];
        end
        monthF=num2str(month);
        if length(monthF)==1  
            monthF=['0',monthF];
        end
        yearF=num2str(year);
        yearB=yearF;
        return
end

if day == eomday(year,month)
        dayB=num2str(day-1);
        dayF='01';
        monthB=num2str(month); 
        if length(monthB)==1  
            monthB=['0',monthB];
        end
        monthF=num2str(month+1);
        if length(monthF)==1  
            monthF=['0',monthF];
        end
        yearF=num2str(year);
        yearB=yearF;
        return
end
%----------------------------------------------月中
        dayB=num2str(day-1);
        if length(dayB)==1  
            dayB=['0',dayB];
        end
        dayF=num2str(day+1);
        if length(dayF)==1  
            dayF=['0',dayF];
        end
        monthB=num2str(month); 
        if length(monthB)==1  
            monthB=['0',monthB];
        end
        monthF=monthB;
        yearF=num2str(year);
        yearB=yearF;
        return  
end

