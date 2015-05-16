function [ line_count ] = get_line_count( x, dim )
    % assume 2D lines
    if (~exist('dim','var') )
        dim = 2;
    end
    stride = dim + 1;
    line_count = numel(x) / stride;
end

