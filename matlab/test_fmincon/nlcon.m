function [ c, ceq ] = nlcon( x )
    c = [];
    ceq = x(1)^2 + x(2)^2 - 1;
end

