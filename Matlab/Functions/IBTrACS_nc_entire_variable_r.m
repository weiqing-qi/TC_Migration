%一次读取最佳路径数据集文件_nc   update20201020
function [SID6,LAT6,LON6,ISO_TIME6,USA_WIND6,USA_RMW6,BASIN6,WMO_WIND6,WMO_PRES6,NAME6,DIST2LAND6,LANDFALL6,STORM_SPEED6,STORM_DIR6]= IBTrACS_nc_entire_variable_r(idir)
  SID6=ncread(idir,'sid');
  LAT6=ncread(idir,'lat');
  LON6=ncread(idir,'lon');
  ISO_TIME6=ncread(idir,'iso_time');
  USA_WIND6=ncread(idir,'usa_wind');
  USA_RMW6=ncread(idir,'usa_rmw');
  BASIN6=ncread(idir,'basin');
  WMO_WIND6=ncread(idir,'wmo_wind');  
  WMO_PRES6=ncread(idir,'wmo_pres');
  NAME6=ncread(idir,'name');
  DIST2LAND6=ncread(idir,'dist2land');
  LANDFALL6=ncread(idir,'landfall');
  STORM_SPEED6=ncread(idir,'storm_speed');
  STORM_DIR6=ncread(idir,'storm_dir');
end
