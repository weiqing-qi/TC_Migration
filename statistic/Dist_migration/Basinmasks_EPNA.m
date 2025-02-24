function [output] = Basinmasks_EPNA(cs,method)
%返回每个Basin的掩模 EP和NA以-105°W分开 Method1：-180-180 2：0-360
%SImask:1 SPmask:2 SAmask:3 NImask:4 WPmask:5 EPmask:6 NAmask:7

SImask=zeros(180/cs,360/cs);
SImask(90/cs+1:180/cs,10/cs+1:135/cs)=1;

SPmask=zeros(180/cs,360/cs);
SPmask(90/cs+1:180/cs,135/cs+1:290/cs)=1;

SAmask=zeros(180/cs,360/cs);
SAmask(90/cs+1:180/cs,1:10/cs)=1;
SAmask(90/cs+1:180/cs,290/cs+1:360/cs)=1;

NImask=zeros(180/cs,360/cs);
NImask(1:90/cs,30/cs+1:100/cs)=1;

WPmask=zeros(180/cs,360/cs);
WPmask(1:90/cs,100/cs+1:180/cs)=1;

EPmask=zeros(180/cs,360/cs);
EPmask(1:90/cs,180/cs+1:255/cs)=1;

NAmask=zeros(180/cs,360/cs);
NAmask(1:90/cs,1:30/cs)=1;
NAmask(1:90/cs,255/cs+1:360/cs)=1;

all_in_one=SImask+SPmask*2+SAmask*3+NImask*4+WPmask*5+EPmask*6+NAmask*7;

if method==2
    output=all_in_one;
end    

if method==1
    output=[all_in_one(:,180/cs+1:360/cs),all_in_one(:,1:180/cs)];  
end

end

