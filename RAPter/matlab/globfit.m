clc;
p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/input_20140620_1547';

if ( 0 )
    %% read points
    cloud_p = [p filesep 'cloud.ply'];
    %points = my_ply_read( cloud_p, 3 );
    %fprintf('Read %d points from %s\n', size(points,1), cloud_p );
    
    %% read primitives
    Dim = 6;
    primitives_p = [p filesep 'primitives_00060.txt'];
    text = textread( primitives_p, '%s', 'delimiter','\n');
    primitives = zeros(0,Dim);
    for i = 1 : numel(text)
        disp(text{i});
        tokens = strsplit( text{i}, ',' );
        
        for i = 1 : Dim
            pnt(i) = sscanf( tokens{i}, '%f' );
        end
        primitives(end+1,1:Dim) = pnt(1:Dim);
    end
end

scale = 0.05;
A = zeros( 4 * size(points,1), 4 * size(points,1) );
plane_normal = [0,0,1];
for lid = 1 : size(primitives,1)
    line    = primitives(lid,:);
    normal  = line_normal( line, plane_normal );
    
    
end
