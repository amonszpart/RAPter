function [ line ] = get_line( id0, x, dim )
    stride = dim + 1;
    line = x( (id0 * stride)+1 : (id0+1) * stride);
end

