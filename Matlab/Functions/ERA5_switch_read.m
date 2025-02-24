function [alldata] = ERA5_switch_read(ERARead,ERA5_idir,hpa,alldata,allmon,allmon_S,allday_S,month,Garea2D_2,cs)
%ERA5_SWITCH_SET 根据输入选择合适的读取方法
switch ERARead
    case 'Horizontal_Divergence'
        add=0;
        for m=1:length(allmon)
            Amonthdata = ncread( ...
                [ERA5_idir,num2str(hpa),'\ERA5_',ERARead,'_on_',num2str(hpa),'hpa_levels_daily',allmon_S(m,:),'.nc'],'d');% 改
            mon_d=find(str2num(allday_S(:,1:6))==allmon(m));
            for k=1+add:length(mon_d)+add
                % alldata(:,:,k)=ERA5FixResample( flipud(Amonthdata(:,:,str2double(allday_S(k,7:8)))') , Garea2D_2, cs );    %Unit: s⁻¹ 
                alldata(:,:,k)=flipud(Amonthdata(:,2:180/cs+1,str2double(allday_S(k,7:8)))');    %Unit: s⁻¹ 
            end                          
            add=add+length(mon_d);
        end
    case 'Relative_Humidity'
        add=0;
        for m=1:length(allmon)
            Amonthdata = ncread( ...
                [ERA5_idir,num2str(hpa),'\ERA5_',ERARead,'_on_',num2str(hpa),'hpa_levels_daily',allmon_S(m,:),'.nc'],'r');% 改
            mon_d=find(str2num(allday_S(:,1:6))==allmon(m));
            for k=1+add:length(mon_d)+add
                % alldata(:,:,k)=ERA5FixResample( flipud(Amonthdata(:,:,str2double(allday_S(k,7:8)))') , Garea2D_2, cs );    %Unit: % 大于零度是水小于零度是冰
                alldata(:,:,k)=flipud(Amonthdata(:,2:180/cs+1,str2double(allday_S(k,7:8)))');    %Unit: % 大于零度是水小于零度是冰
            end                             
            add=add+length(mon_d);
        end
    case 'Relative_Vorticity'
        add=0;
        for m=1:length(allmon)
            Amonthdata = ncread( ...
                [ERA5_idir,num2str(hpa),'\ERA5_',ERARead,'_on_',num2str(hpa),'hpa_levels_daily',allmon_S(m,:),'.nc'],'vo');% 改
            mon_d=find(str2num(allday_S(:,1:6))==allmon(m));
            for k=1+add:length(mon_d)+add
                % alldata(:,:,k)=ERA5FixResample( flipud(Amonthdata(:,:,str2double(allday_S(k,7:8)))') , Garea2D_2, cs );    %Unit: s⁻¹ 
                alldata(:,:,k)=flipud(Amonthdata(:,2:180/cs+1,str2double(allday_S(k,7:8)))');    %Unit: s⁻¹
            end                        
            add=add+length(mon_d);
        end
    case 'Temperature'
        add=0;
        for m=1:length(allmon)
            Amonthdata = ncread( ...
                [ERA5_idir,num2str(hpa),'\ERA5_',ERARead,'_on_',num2str(hpa),'hpa_levels_daily',allmon_S(m,:),'.nc'],'t');% 改
            mon_d=find(str2num(allday_S(:,1:6))==allmon(m));
            for k=1+add:length(mon_d)+add
                % alldata(:,:,k)=ERA5FixResample( flipud(Amonthdata(:,:,str2double(allday_S(k,7:8)))')-273.15 , Garea2D_2, cs );    %Unit: K 
                alldata(:,:,k)=flipud(Amonthdata(:,2:180/cs+1,str2double(allday_S(k,7:8)))')-273.15;    %Unit: ℃ 
            end                      
            add=add+length(mon_d);
        end
    case 'U_Component_Of_Wind'
       add=0;
        for m=1:length(allmon)
            Amonthdata = ncread( ...
                [ERA5_idir,num2str(hpa),'\ERA5_',ERARead,'_on_',num2str(hpa),'hpa_levels_daily',allmon_S(m,:),'.nc'],'u');% 改
            mon_d=find(str2num(allday_S(:,1:6))==allmon(m));
            for k=1+add:length(mon_d)+add
                % alldata(:,:,k)=ERA5FixResample( flipud(Amonthdata(:,:,str2double(allday_S(k,7:8)))') , Garea2D_2, cs );    %Unit: m/s →
                alldata(:,:,k)=flipud(Amonthdata(:,2:180/cs+1,str2double(allday_S(k,7:8)))');    %Unit: m/s →
            end                         
            add=add+length(mon_d);
        end
    case 'V_Component_Of_Wind'
       add=0;
        for m=1:length(allmon)
            Amonthdata = ncread( ...
                [ERA5_idir,num2str(hpa),'\ERA5_',ERARead,'_on_',num2str(hpa),'hpa_levels_daily',allmon_S(m,:),'.nc'],'v');% 改
           mon_d=find(str2num(allday_S(:,1:6))==allmon(m));
            for k=1+add:length(mon_d)+add
                % alldata(:,:,k)=ERA5FixResample( flipud(Amonthdata(:,:,str2double(allday_S(k,7:8)))') , Garea2D_2, cs );    %Unit: m/s ↑
                alldata(:,:,k)=flipud(Amonthdata(:,2:180/cs+1,str2double(allday_S(k,7:8)))');    %Unit: m/s ↑
            end                      
            add=add+length(mon_d);
        end
    case 'Vertical_Velocity'
       add=0;
        for m=1:length(allmon)
            Amonthdata = ncread( ...
                [ERA5_idir,num2str(hpa),'\ERA5_',ERARead,'_on_',num2str(hpa),'hpa_levels_daily',allmon_S(m,:),'.nc'],'w');% 改
            mon_d=find(str2num(allday_S(:,1:6))==allmon(m));
            for k=1+add:length(mon_d)+add
                % alldata(:,:,k)=ERA5FixResample( flipud(Amonthdata(:,:,str2double(allday_S(k,7:8)))') , Garea2D_2, cs );    %Unit: Pa/s 正值下降负值上升
                alldata(:,:,k)=flipud(Amonthdata(:,2:180/cs+1,str2double(allday_S(k,7:8)))');    %Unit: Pa/s 正值下降负值上升
            end                           
            add=add+length(mon_d);
        end
end

end

