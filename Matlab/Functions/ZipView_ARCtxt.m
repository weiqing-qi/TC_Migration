function [Cut_IG] = ZipView_ARCtxt(CorF,InputGrid,PR_odir,name,cs,LatCT,LonCT,NODATA_value,precise)
%剪裁到最小方框或填充到指定分辨率全球格网 输入必须包含NaN 并且输入变量全是-180-180 90-90分辨率 
% 输出可能为0-360 或180-180 达到最小但是读取的时候要注意
if strcmp(CorF,'cut&save')

    [ScalR,ScalC]=find(~isnan(InputGrid));

    Min_SR=min(ScalR);
    Max_SR=max(ScalR);
    Min_SC=min(ScalC);
    Max_SC=max(ScalC);

    if Min_SC==1 && Max_SC==360/cs %这个肯定跨180度了
        InputGrid=[InputGrid(:,(180/cs)+1:360/cs),InputGrid(:,1:180/cs)];
        
        [ScalR,ScalC]=find(~isnan(InputGrid));

        Min_SR=min(ScalR);
        Max_SR=max(ScalR);
        Min_SC=min(ScalC);
        Max_SC=max(ScalC);
    
        Cut_IG = InputGrid(Min_SR : Max_SR , Min_SC : Max_SC);
        Cut_IG(Cut_IG<0)=NODATA_value;
        Cut_IG(isnan(Cut_IG))=NODATA_value;
    
        ncols=Max_SC-Min_SC+1;
        nrows=Max_SR-Min_SR+1;

        Min_SC180=Min_SC+180/cs; %此时Min_SC一定是1-180，从0360转换到-180-180索引
        xllcorner=LonCT(Max_SR,Min_SC180)-0.5*cs;
        yllcorner=LatCT(Max_SR,Min_SC180)-0.5*cs;
    
        Cut_header = GistxtHeader(cs,ncols,nrows,xllcorner,yllcorner,NODATA_value); 
    
        if precise==2
            OUTPUT(PR_odir,name,Cut_header,Cut_IG);
        elseif precise==4
            OUTPUT4(PR_odir,name,Cut_header,Cut_IG);   
        end 
    else
        Cut_IG = InputGrid(Min_SR : Max_SR , Min_SC : Max_SC);
        Cut_IG(Cut_IG<0)=NODATA_value;
        Cut_IG(isnan(Cut_IG))=NODATA_value;
    
        ncols=Max_SC-Min_SC+1;
        nrows=Max_SR-Min_SR+1;
    
        xllcorner=LonCT(Max_SR,Min_SC)-0.5*cs;
        yllcorner=LatCT(Max_SR,Min_SC)-0.5*cs;
    
        Cut_header = GistxtHeader(cs,ncols,nrows,xllcorner,yllcorner,NODATA_value); 
    
        if precise==2
            OUTPUT(PR_odir,name,Cut_header,Cut_IG);
        elseif precise==4
            OUTPUT4(PR_odir,name,Cut_header,Cut_IG);   
        end
    end

elseif strcmp(CorF,'cut')

    [ScalR,ScalC]=find(InputGrid);

    Min_SR=min(ScalR);
    Max_SR=max(ScalR);
    Min_SC=min(ScalC);
    Max_SC=max(ScalC);

    Cut_IG = InputGrid(Min_SR : Max_SR , Min_SC : Max_SC);
end
end

