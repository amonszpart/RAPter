function [ normal ] = line_normal( line, plane_normal )
%LINE_NORMAL Summary of this function goes here
%   Detailed explanation goes here

    dir     = line(4:6);
    perp    = plane_normal * dot( dir,plane_normal );
    par     = dir - perp;
    par     = par / norm(par);
    n       = cross( par, plane_normal );
    normal  = n / norm(n);
end

