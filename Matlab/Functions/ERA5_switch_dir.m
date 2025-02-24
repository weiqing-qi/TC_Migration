function [ERA5_idir,Excel_file,odir,Outname] = ERA5_switch_dir(ERASET)
%ERA5_SWITCH_SET 根据输入选择合适的文件路径
switch ERASET
    case 'Horizontal_Divergence'
        ERA5_idir='E:\DATA\1.Reanalysis\EAR5\ERA5-HDIV\';
        % odir='D:\DATA\TC_spatial_data\EAR5\ERA5-HDIV\';
        % Outname='ERA5_TC_AVE_Horizontal_Divergence_Day05_500KM_Ummd';
        Excel_file='D:\Desktop2\Attri-project\Results\Unet\Attri_EC_Data_UnetRMRW.xlsx';
    case 'Relative_Humidity'
        ERA5_idir='E:\DATA\1.Reanalysis\EAR5\ERA5-RH\';
        % odir='D:\DATA\TC_spatial_data\EAR5\ERA5-RH\';
        % Outname='ERA5_TC_AVE_Relative_Humidity_Day05_500KM_Ummd';
        Excel_file='D:\Desktop2\Attri-project\Results\Unet\Attri_EC_Data_UnetRMRW.xlsx';
    case 'Relative_Vorticity'
        ERA5_idir='E:\DATA\1.Reanalysis\EAR5\ERA5-RV\';
        % odir='D:\DATA\TC_spatial_data\EAR5\ERA5-RV\';
        % Outname='ERA5_TC_AVE_Relative_Vorticity_Day05_500KM_Ummd';
        Excel_file='D:\Desktop2\Attri-project\Results\Unet\Attri_EC_Data_UnetRMRW.xlsx';
    case 'Temperature'
        ERA5_idir='E:\DATA\1.Reanalysis\EAR5\ERA5-Temperature\';
        % odir='D:\DATA\TC_spatial_data\EAR5\ERA5-Temperature\';
        % Outname='ERA5_TC_AVE_Temperature_Day05_500KM_Ummd';
        Excel_file='D:\Desktop2\Attri-project\Results\Unet\Attri_EC_Data_UnetRMRW.xlsx';
    case 'U_Component_Of_Wind'
        ERA5_idir='E:\DATA\1.Reanalysis\EAR5\ERA5-UCW\';
        % odir='D:\DATA\TC_spatial_data\EAR5\ERA5-UCW';
        % Outname='ERA5_TC_AVE_UCW_Day05_500KM_Ummd';
        Excel_file='D:\Desktop2\Attri-project\Results\Unet\Attri_EC_Data_UnetRMRW.xlsx';
    case 'V_Component_Of_Wind'
        ERA5_idir='E:\DATA\1.Reanalysis\EAR5\ERA5-VCW\';
        % odir='D:\DATA\TC_spatial_data\EAR5\ERA5-VCW\';
        % Outname='ERA5_TC_AVE_VCW_Day05_500KM_Ummd';
        Excel_file='D:\Desktop2\Attri-project\Results\Unet\Attri_EC_Data_UnetRMRW.xlsx';
    case 'Vertical_Velocity'
        ERA5_idir='E:\DATA\1.Reanalysis\EAR5\ERA5-Vertical-Velocity\';
        % odir='D:\DATA\TC_spatial_data\EAR5\ERA5-Vertical_Velocity\';
        % Outname='ERA5_TC_AVE_Vertical_Velocity_Day05_500KM_Ummd';
        Excel_file='D:\Desktop2\Attri-project\Results\Unet\Attri_EC_Data_UnetRMRW.xlsx';
end
end

