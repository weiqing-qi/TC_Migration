function [NOlim] = lim_Y2NO(Ylim,ISO_TIME)
%input the year scale return the start end NO in IBTrACS
%Ylim:[STARTY ENDY] NOlim:[STARTNO ENDNO COUNT]
FirstY = str2num(permute( ISO_TIME(1:4,1,:),[3 1 2] ));
StartNO=find(FirstY==Ylim(1),1,'first');
EndNO=find(FirstY==Ylim(2),1,'last');
NOlim=[StartNO,EndNO,EndNO-StartNO+1];
end