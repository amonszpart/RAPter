function [ output_args ] = show_stats( path, dim )
    if ( ~exist('dim','var' ) )
        dim = 2;
    end
    x = load([path filesep 'x_solved.mat']);
    x = x.x;
   
    dots = [];
    x'
    line_count = get_line_count( x, dim )
    for lid = 1 : line_count-2
        for lid1 = lid+1 : line_count-1
            line0 = get_line(lid,x,dim);
            line1 = get_line(lid1,x,dim);
            dots(end+1) = dot( line0(1:dim), line1(1:dim) );
        end
    end
    
    subplot(1,2,1);
    plot(dots,'.');
    title('n * n-1 dot products');
    xlabel(sprintf('line-pair id (0-1,0-2...%d-%d)',line_count-2,line_count-1));
    ylabel('dot product of two normals')
    
    subplot(1,2,2);
    abs_dots = sort(real(acos((abs(dots))))) / pi * 180.0
    plot(abs_dots,'.');
    title('SORTED n * n-1 acos(abs(angles))');
    xlabel(sprintf('line-pair id (order not kept)',line_count-2,line_count-1));
    ylabel('absolute angle of two normals in Degrees');
    
    
end

