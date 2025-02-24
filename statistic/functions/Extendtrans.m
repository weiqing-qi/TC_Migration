function [] = Extendtrans(inpath,outpath,cs)
%EXTENDCHANGE 0.5¡ã ×ª»»0-360Óë-180-180

header2=['ncols           720';
         'nrows           360';
         'xllcorner      -180';
         'yllcorner       -90';
         'cellsize        0.5';
         'NODATA_value   -999'];
infiles = dir([inpath,'*.txt']);

for NO=1:length(infiles)
    [~,~,pre]= read_ARCascii ([inpath,infiles(NO).name]); 
    P=[pre(:,(180/cs)+1:360/cs),pre(:,1:180/cs)];
    OUTPUT(outpath,[infiles(NO).name(1:20),'original'],header2,P)
end

end
