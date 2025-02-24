%TMPA补充60-50纬度带之间的降水值
IMERG_idir='F:\GlobalTCMask1\single_pre\IMERG_Single_Precip3d_mmd\';
IMERGfilename='IMERGV6_Ori_SingleTC_3dPrecip_Ummd_SID';
Pre19=zeros(360,720,19);%总和降水
parfor y=2001:2019
    %每一年的降水全加起来放进Pre41;每年总降水
    files=dir([IMERG_idir,IMERGfilename,num2str(y),'*.txt']);
    Precontent=zeros(360,720);
    countcontent=zeros(360,720);
    for f=1:length(files)
    [~,inputheader,Preci]= read_ARCascii ([IMERG_idir,files(f).name]);
    countcontent(Preci>=0.1)=countcontent(Preci>=0.1)+1;%降水大于0.1mm/d
    Preci(Preci<0)=0; Preci(isnan(Preci))=0;
    Precontent=Precontent+Preci;
    end
    Pre19(:,:,y-2000)=Precontent;
end
ALL=Pre19;
ALL(ALL<0)=0;ALL(isnan(ALL))=0;
N60=sum(sum(sum(ALL(60:300,:,:))));
N50=sum(sum(sum(ALL(80:280,:,:))));
100-N50/N60*100