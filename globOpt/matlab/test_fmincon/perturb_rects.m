function [ rot_lines ] = perturb_rects( lines, dim, angs, groups )
    %groups has rows with line ids that should be rotated by the same value
    
    if ( ~exist('dim','var') )
        dim = 2;
    end
    
    stride    = dim + 1;
    nLines    = size(lines,1)/stride;
    rot_lines = lines;
    
    for g = 1 : size( groups, 1 )
        R = rotz( angs(g) );
        for i = 1 : size( groups, 2 )
            if groups(g,i) == 0
                break;
            end
            id = groups(g,i);
            
            n = zeros(3,1);
            for d = 1 : dim
                n(d) = lines( (id-1) * stride + d );
            end
            n
            rot_n = (R * n)';
            rot_n
            rot_n_dim = rot_n(1:dim);
            rot_n_dim = rot_n_dim / norm( rot_n_dim );
            rot_lines( (id-1) * stride + 1 : (id-1) * stride + dim ) = rot_n_dim;
            
        end
    end
    [lines,rot_lines]
end