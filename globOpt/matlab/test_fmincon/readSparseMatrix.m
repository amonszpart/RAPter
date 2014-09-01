function mx = readSparseMatrix( path, offset )
    fid = fopen( path, 'r' );
    line = fgetl( fid );
    hw = sscanf( line, '%d,%d');
    height = hw(1);
    width  = hw(2);
    
    mx = sparse(height,width);
    while ( ~feof(fid) )
        line = fgetl( fid );
        entry = sscanf(line,'%f,%f,%f');
        mx( entry(1) + offset, entry(2) + offset ) = entry(3);
    end
    fclose(fid);
end

