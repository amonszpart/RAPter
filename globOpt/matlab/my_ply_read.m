function points = my_ply_read( cloud_p, Dim )
% cloud_p cloud path
% Dim how many columns to read

fid = fopen( cloud_p );

tline = fgetl( fid );
start = false;
points = zeros( 0, Dim );
N = 0;
pid = 0;
while ischar( tline ) && ( (N == 0) || (size(points,1) < N) )
    tline = fgetl(fid);
    
    tokens = strsplit( tline );
    
    if ( start )
        for i = 1 : Dim
            pnt(i) = sscanf( tokens{i}, '%f' );
        end
        %fprintf(' pnt: {%f,%f,%f}\n', pnt );
        points(end+1,1:Dim) = pnt(1:Dim);
    end
    
    if ( strcmp(tokens{1}, 'end_header') )
        start = true;
    elseif ( strcmp(tokens{1}, 'element') && strcmp(tokens{2}, 'vertex') )
        N = sscanf( tokens{3}, '%d');
        %fprintf('------------- N: %d----------------\n', N );
    end
    
end

fclose( fid );

end