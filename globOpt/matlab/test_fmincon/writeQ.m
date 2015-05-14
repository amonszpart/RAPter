function writeQ( Q, path )
    path
    fid = fopen( path, 'w', 'n', 'UTF-8' );
    % size
    fprintf(fid, '%d,%d\n', size(Q,1), size(Q,2) );
    fprintf(     '%d,%d\n', size(Q,1), size(Q,2) );
    % elements
    for r = 1 : size(Q,1)
        for c = 1 : size(Q,2)
            if Q(r,c) > 0
                %fprintf(    '%d,%d,%f\n', r, c, Q(r,c) );
                fprintf(fid,'%d,%d,%f\n', r, c, Q(r,c) );
            end
        end
    end
    fclose( fid );
end