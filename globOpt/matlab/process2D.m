% takes a 2D image, an prepares it for sampling
if 0
    close all;
    figure();
    p = '/home/bontius/workspace/globOpt/data/scenes/floorplan_niloy/floorplan.png'
    H = 2; W = 2; k = 1;

    I = imread(p);
    subplot(H,W,k); k=k+1;
    imshow(I);

    subplot(H,W,k); k=k+1;
    G = rgb2gray( I );
    imshow(G);

    subplot(H,W,k); k=k+1;
    bw = edge(G,'log');
    imshow( bw );

    figure();
    imshow( bw );
end

N = 100000.
[ H, W ] = size( bw );
rat = N / (H * W);
P = zeros(0,2);
for y = 1 : H
    for x = 1 : W
        if ( rand < rat && bw(y,x) > 0)
            P(end+1,:) = [ x/W, y/H ];
        end
    end
end

fid = fopen('out.ply','w','n','UTF-8' );
fprintf( fid, ['ply\n' ...
               'format ascii 1.0\n' ...
                'comment Aron generated\n' ...
                'element vertex %d\n' ...
                'property float x\n' ...
                'property float y\n' ...
                'property float z\n' ...
                'end_header\n'], size(P,1) );
for i = 1 : size(P,1)
    fprintf( fid, '%f %f 0.0\n', P(i,1), P(i,2) );
end
fclose( fid );