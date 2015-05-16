function [ obj ] = objective( x )
%OBJECTIVE Summary of this function goes here
%   Detailed explanation goes here

%0 nx + 0 ny + 0 d + [ 10 nx ^ 2 + 12 nx * d + 10 nx * ny + 10 ny ^ 2 + 12 ny * d + 4 d ^ 2 ] / 2 
Q = [ 5, 5, 6; 0, 5, 6; 0, 0, 2 ];

obj = x * Q * x';

end

