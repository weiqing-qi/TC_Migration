%������뺣����ͬ���뻷�����ˮռ�������ı����仯 ������basinl�ı仯
%����ÿ�������ڵĽ�ˮ��ǿ�ȱ仯����Ϊ���ں��ܾ�����
%ע��������60NS 
clear; clc; close all;

file_idir='D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\CoastLine110m\Coast_Buffer_txt\';
filename = 'ns60_coast_bothside_buffer';%ע���ǲ��ǰ�˳�򴢴��
%��ˮ����
Preindir='D:\DATA\TC_spatial_data\PreRate\PERSIANN\'; %��
Prefile=dir([Preindir,'*.txt']);
Prefile=Prefile(2286:4225);%��1980��ʼ1955ERA5 2898JRA55
Ylim=[2001 2019];

year1=xlsread('D:\Desktop2\Global_cyclone_project\DATA_Filted.xlsx',1,'AR2288:AR4227');
%½�ع�������
[~,~,Clandmask] = read_ARCascii('D:\Desktop2\Global_cyclone_project\Picture_materials\ARCGIS\Rastermasks\countriesmasks.txt');%½�غ͹���0360
Clandmask=[Clandmask(:,361:720),Clandmask(:,1:360)];%-180-180 �����������ݵ�extend�����޸�
BasMas = Basinmasks_EPNA(0.5,1);
%% ������뻷�ͺ����߽��� 
CRing_sea=zeros(360,720,20);
CRing_land=zeros(360,720,20);
Cbuffer_sea=zeros(360,720,20);
Cbuffer_land=zeros(360,720,20);
CR_BS=zeros(360,720,20);

    [lmheader1,lmheader2,CR_BS1] = read_ARCascii([file_idir,filename,'100km.txt']);
    CR_BS1(CR_BS1~=1)=0;CR_BS1(CR_BS1==1)=10000;
    CR_BS(:,:,1)=CR_BS1;
for i=200:100:2000
    [~,~,CB_Near] = read_ARCascii([file_idir,filename,num2str(i-100),'km.txt']);
    [~,~,CB_Far] = read_ARCascii([file_idir,filename,num2str(i),'km.txt']);
    CR_BS(:,:,i*0.01)=CB_Far-CB_Near;%��ֵӦ����1-��-9999��10000
end
    CR_BS(1:60,:,:)=0; CR_BS(301:360,:,:)=0;%ȥ��60ns��� 
    
for j=1:20%100:100:2000
    Csea=zeros(360,720);
    Cland=zeros(360,720);
    
    BS=CR_BS(:,:,j);
    Csea(Clandmask<0)=BS(Clandmask<0);
    Cland(Clandmask>0)=BS(Clandmask>0);
    CRing_sea(:,:,j)=Csea;
    CRing_land(:,:,j)=Cland;
    Cbuffer_sea(:,:,j)=sum(CRing_sea,3);
    Cbuffer_land(:,:,j)=sum(CRing_land,3);
end
%% ƽ����ˮ�ʽ���
R_PR=zeros(length(Prefile),28,7);%-800��2000km  %SImask:1 SPmask:2 SAmask:3 NImask:4 WPmask:5 EPmask:6 NAmask:7 
tic;
%ÿ��TCÿһRorB�е�ƽ��ǿ��
parfor i=1:length(Prefile)
[~,H,prerate] = read_ARCascii([Preindir,Prefile(i).name]);%��λmm/D
prerate = FillView2Global(H,prerate);
% prerate = areaweight_downscale(prerate, 0.25, 0.5)
prerate(isnan(prerate))=0; prerate(prerate<2.4)=0;
    for bas=1:7
        for d=1:28%100:2000km
            if d<=8                               %Land
            ZoneLR=CRing_land(:,:,9-d);           %ע�⽵ˮת��
            ZoneLR(BasMas~=bas)=0;                %��������������
            preLR=prerate(ZoneLR==10000);         %������������н�ˮ
            R_PR(i,d,bas)=mean(preLR(preLR>0.1));     %0.1���ż���ƽ����ǿ

            else                                  %Sea
            ZoneSR=CRing_sea(:,:,d-8);
            ZoneSR(BasMas~=bas)=0;                %��������������
            preSR=prerate(ZoneSR==10000);
            R_PR(i,d,bas)=mean(preSR(preSR>0.1)); 
            end
        end
    end
end
toc;
R_PR(isnan(R_PR))=0;
JudgeB=permute(sum(R_PR,2),[1 3 2]);%���о����������һ�𣬲�ȥ���������������ά�� %ÿ����صķ�ĸ��ÿ�����

%������ƽ�������Իع�

R_PR_YAVE=zeros(Ylim(2)-Ylim(1)+1,28,7);%-800��2000km
for bas=1:7
    for yi=Ylim(1):Ylim(2)
        Jud=JudgeB(:,bas);
        Jud=Jud(year1==yi);%�ж�ÿ�����BAS����û�н�ˮ

        for DR=1:28
           R1C=R_PR(:,DR,bas); 
           R1C=R1C(year1==yi); %R1Cû�п�ֵ
           R1C=R1C(Jud>0);%ֻҪ������������ֹ���ˮ������û���ֵĲ�Ҫ%���� ������ƽ������
           R_PR_YAVE(yi-Ylim(1)+1,DR,bas)=mean(R1C);
        end
    end
end

PR_Regress=zeros(28,9,7);
PR_MK=zeros(28,4,7);
for bas=1:7
    for DInt=1:28
    [b,bint,r,rint,stats]=regress(R_PR_YAVE(:,DInt,bas),[ones(Ylim(2)-Ylim(1)+1,1),(Ylim(1):Ylim(2))'],0.05);
        PR_Regress(DInt,1,bas)=b(2);
        PR_Regress(DInt,2,bas)=bint(2,1);%ע��
        PR_Regress(DInt,3,bas)=bint(2,2);
        PR_Regress(DInt,4,bas)=stats(3); 
        
    [ci_lower,ci_upper,pvals,b_Test,DW_p] = NeweyWestAdjust(r,(Ylim(1):Ylim(2))',R_PR_YAVE(:,DInt,bas),1,1);
        PR_Regress(DInt,5,bas)=b_Test(2);
        PR_Regress(DInt,6,bas)=ci_lower(2);%ע��
        PR_Regress(DInt,7,bas)=ci_upper(2);
        PR_Regress(DInt,8,bas)=pvals(2);
        PR_Regress(DInt,9,bas)=DW_p;

     [Zs, p_value, UFk, UBk2, beta, beta_CI, corr_val]= MKtrend((Ylim(1):Ylim(2))',R_PR_YAVE(:,DInt,bas));
        PR_MK(DInt,1,bas)=p_value;
        PR_MK(DInt,2,bas)=beta;
        PR_MK(DInt,3,bas)=beta_CI(1);
        PR_MK(DInt,4,bas)=beta_CI(2);
    end
end
%DW����ؼ����Լ�Nw����
for bas=1:7
    for i=1:28
        for j=2:4
            if ~isnan(PR_Regress(i,j+4,bas))
                PR_Regress(i,j)=PR_Regress(i,j+4,bas);
            end
        end
    end
end


plot(PR_Regress(:,1,5))
% plot(permute(PR_Regress(:,1,:),[1 3 2]))
% plot(mean(R_PR_YAVE(:,:,5)))
R_distribution=(mean(R_PR_YAVE,1,"omitmissing"));
R_distribution=(permute(R_distribution(1,:,:),[3 2 1]))';
R_distribution(:,3)=[];
OUT=zeros(28,21);
for j=1:7
    OUT(:,3*j-2:3*j)=PR_Regress(:,1:3,j);
end
OUT(:,7:9)=[];
imshow(OUT)

SIG=zeros(28,7);%������
for j=1:7
    SIG(:,j)=PR_Regress(:,4,j);
end
SIG(:,3)=[];
toc;
%------------------------DW OUT PUT----------------------------����Ҫ��
% DW=zeros(28,7);%������
% for j=1:7
%     DW(:,j)=PR_Regress(:,9,j);
% end
% DW(:,3)=[];
% DWTEST=find(DW < 0.05 & DW > 0)
% DWOUT=zeros(28,28);
% for j=1:7
%     DWOUT(:,4*j-3:4*j)=PR_Regress(:,5:8,j);
% end
% DWOUT(:,9:12)=[];
