%% ���������TC-��ΧӰ�����أ�������28�����ݵ�����̫�٣�Ӧ�ò����ʣ�������Ҫ������ �ķ���
% corr_valT=zeros(4,3*4);
% corr_valD=zeros(4,3*4);
% P_valD=zeros(4,3*4);
% P_valT=zeros(4,3*4);
% 
% D=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\DATA3TMPA60NS.xlsx',5,'M97:X124');%��Ӱ������L97:W124 L128:W155
% d1=D(:,7:10);%�� Ӱ������
% d2=D(:,11:12);%�� PRE����
% d1=d1(1:28,:);%��
% d2=d2(1:28,:);%�� 
% 
% for j=3:4%Ӱ������ ��
%     for i=3:4 %PRE����
%     [corr_valT(j,i*3-2),P_valT(j,i*3-2)]=corr(d1(:,j),d2(:,i-2),'Type','Pearson');%Pearson_COE
%     [corr_valT(j,i*3-1),P_valT(j,i*3-1)]=corr(d1(:,j),d2(:,i-2),'Type','Kendall');%kendall_tau
%     [corr_valT(j,i*3  ),P_valT(j,i*3  )]=corr(d1(:,j),d2(:,i-2),'Type','Spearman');%Spearman_rho
%     end
% end
% 
% out=[corr_valT;corr_valD];
%------------------��ƽ��
% Data=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',5,'C3:AD4461');%��ˮ
% year1=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',5,'A3:A4461');%��ˮ
% R_PR_YAVE=zeros(41,28);%-800��2000km
% for yi=1980:2020
%     for DR=1:28
%        R1C=Data(:,DR); 
%        R1C=R1C(year1==yi);
% %        R1C(R1C<0.1 & R1C>=0)=[];
%        R_PR_YAVE(yi-1979,DR)=mean(R1C);
%     end
% end
%----------------------------------------------------------------------
%% ����ÿ��̨��ÿһ����������Ľ�ˮ
P_val=zeros(4,3*4);
corr_val=zeros(4,3*4);
test=zeros(4,3*4);
G=0;%��ȫ��1���Ƿ���0��
U=1;%��
D=0;%��
tic;
for pre=1:4 %ȫ��
    if G==1
        D1=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',pre,'B3:B4461');%��ˮ
    elseif G==0
        D1=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',pre,'AF3:BG43');%��ˮC3:AD4461 AF3:BG43
    end
    for influ=5:8 %�����
        
        if G==1
            D2=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',influ,'B3:B4461');%��ˮ
        elseif G==0
            D2=xlsread('C:\Users\Dell\Desktop\Desktop2\Global_cyclone_project\PICbin\migration\Correlation.xlsx',influ,'AF3:BG43');%Ӱ������    
        end
        %-----------------��������
        if G==0
            co1=reshape(D1,[28*length(D1),1]);
            co2=reshape(D2,[28*length(D1),1]);
        elseif G==1
            co1=D1;
            co2=D2;
        end
        co2(co1==-9999)=[];co1(co1==-9999)=[];%��������ݲ����ʱ��
        co1(co2==-9999)=[];co2(co2==-9999)=[];%�����ں��µ������excel����ԭʼ��û��nanֻ�п�ֵ����nan
        co2(co1<=0.1)=[];co1(co1<=0.1)=[];%��ˮ0��Ҫ��ֻ�н�ˮ��С��0.1
        %-----------------��������,ѡ���С�ż��е����ݼ����ϵ
        thre=co1;
        I1=sort(co1);%ֻ�ý�ˮ���ż� ����
        threU=I1(round(length(I1)*U));%75%�ٷ�λ��
        if D==0
            threD=0.1;
        else
            threD=I1(round(length(I1)*D));%5%�ٷ�λ��
        end
        co1=co1(thre > threD & thre <= threU);
        co2=co2(thre > threD & thre <= threU);
        %-----------------����
        [corr_val(influ-4,pre*3-2),P_val(influ-4,pre*3-2)]=corr(co1,co2,'Type','Pearson');%Pearson_COE
        [corr_val(influ-4,pre*3-1),P_val(influ-4,pre*3-1)]= ;%Spearman_rho
        [corr_val(influ-4,pre*3  ),P_val(influ-4,pre*3  )]=corr(co1,co2,'Type','Kendall');%kendall_tau
        test(influ-4,pre*3-2  )=lillietest(co1);test(influ-4,pre*3-1  )=lillietest(co2);
        n=(length(co1))
    end
end
toc;
P_val;
[G D U]
corr_val











