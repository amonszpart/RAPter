% takes a 2D image, an prepares it for sampling
close all;
figure();
p = '/home/bontius/workspace/globOpt/data/scenes/floorplan3/Ken_Pl_alterations_plan.jpg'
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