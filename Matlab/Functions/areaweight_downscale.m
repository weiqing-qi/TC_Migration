function [out] = areaweight_downscale(in, in_cs, out_cs, inarea, outarea)
% in & OUT & inarea & outarea are globl scale 2D geo-grids
% out_cs > in_cs
% in 中的missing 值应该是-9999
    if nargin <= 3 
        % 计算输入和输出网格的每个格网的地理面积
        [inarea, ~] = Gridarea(in_cs);
        [outarea, ~] = Gridarea(out_cs);
    end
    % 计算重采样比例因子
    factor = out_cs / in_cs;

    % 初始化输出网格
    out_rows = size(in, 1) / factor;
    out_cols = size(in, 2) / factor;

    % 重塑输入数据和面积矩阵以适应输出网格
    sub_in = reshape(in, factor, out_rows, factor, out_cols);
    sub_inarea = reshape(inarea, factor, out_rows, factor, out_cols);

    % 计算每个输出像素中的加权平均
    weighted_sum = sum(sum(sub_in .* sub_inarea, 1), 3);  % 沿着第1和第3维求和
    total_area = reshape(outarea, out_rows, out_cols);

    % 最终输出是加权总和除以每个新网格点的面积
    out = squeeze(weighted_sum) ./ total_area;
    
    % out = zeros(out_rows, out_cols);
    % % 遍历输出网格的每个单元
    % for r = 1:out_rows
    %     for c = 1:out_cols
    %         % 计算输入网格中对应的范围
    %         in_row_start = (r-1) * factor + 1;
    %         in_row_end = r * factor;
    %         in_col_start = (c-1) * factor + 1;
    %         in_col_end = c * factor;
    % 
    %         % 计算子矩阵中的所有输入像素和对应的面积
    %         sub_in = in(in_row_start:in_row_end, in_col_start:in_col_end);
    %         sub_inarea = inarea(in_row_start:in_row_end, in_col_start:in_col_end);
    % 
    %         % 进行面积加权平均
    %         weighted_sum = sum(sub_in .* sub_inarea, "all");
    %         out(r, c) = weighted_sum / outarea(r, c);
    %     end
    % end

    out(out<0)=-9999;
end

