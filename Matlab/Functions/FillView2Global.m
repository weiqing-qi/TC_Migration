function [Full_V] = FillView2Global(H,Input)
%任意剪裁ARCtxt格网填充到-180-180 90--90指定分辨率全球格网 input from ZipView_ARCtxt

L_leftlon=round(H(3)/H(5)) * H(5);
R_rightlon=L_leftlon + H(5)*H(1);

if L_leftlon + H(5)*H(1)>=180                                              %一定是">=" yes-> 0-360 no-> -180-180

    Lfill=repmat( H(6), H(2), round(L_leftlon/H(5)) );                            % fill矩阵可以是空的
    Rfill=repmat( H(6), H(2), round((360-R_rightlon)/H(5)) );
    LRFG=[Lfill,Input,Rfill];
    LRFG=[LRFG(:,(180/H(5))+1:360/H(5)),LRFG(:,1:180/H(5))];
else

    Lfill=repmat( H(6), H(2), round((180+L_leftlon)/H(5)) );
    Rfill=repmat( H(6), H(2), round((180-R_rightlon)/H(5)) );
    LRFG=[Lfill,Input,Rfill];
end

D_bottomlat=round(H(4)/H(5))*H(5);
U_toplat=D_bottomlat+H(5)*H(2);

Ufill=repmat(H(6), round((90-U_toplat)/H(5)), 360/H(5));
Dfill=repmat(H(6), round((90+D_bottomlat)/H(5)), 360/H(5));

LRUDFG=[Ufill;LRFG;Dfill];

[R,C]=size(LRUDFG);

if R==180/H(5) && C==360/H(5)
    Full_V=LRUDFG;
else
    Full_V='error extend';
end

end