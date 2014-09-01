function write_files( fit_lines, path )
%WRITE_FILES Summary of this function goes here
%   Detailed explanation goes here
path
fid = fopen( path, 'w', 'n', 'UTF-8' );
for l = 1 : size(fit_lines,1)
    for d = 1 : size(fit_lines,2)
        fprintf( fid, '%f,', fit_lines(l,d) );
    end
    fprintf( fid, '\n');
end
fclose( fid );

end

