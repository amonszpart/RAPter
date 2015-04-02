function [outputParameters, initialFittingError, exitFittingError, exitflag] = OptimizePoint(inputParameters, maxIterNum, numVertices, primitiveType, coordX, coordY, coordZ, normalX, normalY, normalZ, confVertices, constraints)

LoadNameMap

numPrimitives = size(numVertices, 1);
numConstraints = size(constraints, 1);
if 0 < numConstraints
    constraints(1 : numConstraints, 2 : 5) = constraints(1 : numConstraints, 2 : 5) + 1;
end

% normalization
for i = 1 : numPrimitives
    if shape.plane == primitiveType(i)
        inputParameters(i, 7) = inputParameters(i, 7) / norm(inputParameters(i, 1 : 3));
    end
    if shape.plane == primitiveType(i) || shape.cylinder == primitiveType(i) || shape.cone == primitiveType(i)
        inputParameters(i, 1 : 3) = inputParameters(i, 1 : 3) / norm(inputParameters(i, 1 : 3));
    end
end

fixedParameters = inputParameters;

% expression expansion
coefficients = zeros(numPrimitives, 35);
for i = 1 : numPrimitives
    if shape.sphere == primitiveType(i)
        % sphere cost: ((a - x)^2 - r^2 + (b - y)^2 + (c - z)^2)^2
        % x^4 + y^4 + z^4
        % 2*(x^2*y^2 + x^2*z^2 + y^2*z^2) +
        % (-4*a)*(x^3 + x*y^2 + x*z^2) +
        % (-4*b)*(x^2*y + y^3 + y*z^2) +
        % (-4*c)*(x^2*z + y^2*z + z^3) +
        % (6*a^2 + 2*b^2 + 2*c^2 - 2*r^2)*x^2 +
        % (2*a^2 + 6*b^2 + 2*c^2 - 2*r^2)*y^2 +
        % (2*a^2 + 2*b^2 + 6*c^2 - 2*r^2)*z^2 +
        % (8*a*b)*x*y +
        % (8*a*c)*x*z +
        % (8*b*c)*y*z +
        % (-4*a*(a^2 + b^2 + c^2 - r^2))*x +
        % (-4*b*(a^2 + b^2 + c^2 - r^2))*y +
        % (-4*c*(a^2 + b^2 + c^2 - r^2))*z +
        % (a^2 + b^2 + c^2 - r^2)^2
        coefficients(i, 1) = sum(confVertices(i, 1 : numVertices(i)));
        coefficients(i, 2) = 2*coefficients(i, 1);
        coefficients(i, 3) = -4*dot(coordX(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 4) = -4*dot(coordY(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 5) = -4*dot(coordZ(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        r = fixedParameters(i, 7);
        coefficients(i, 6) = dot(confVertices(i, 1 : numVertices(i)), 6*coordX(i, 1 : numVertices(i)).^2 + 2*coordY(i, 1 : numVertices(i)).^2 + 2*coordZ(i, 1 : numVertices(i)).^2 - 2*r^2);
        coefficients(i, 7) = dot(confVertices(i, 1 : numVertices(i)), 2*coordX(i, 1 : numVertices(i)).^2 + 6*coordY(i, 1 : numVertices(i)).^2 + 2*coordZ(i, 1 : numVertices(i)).^2 - 2*r^2);
        coefficients(i, 8) = dot(confVertices(i, 1 : numVertices(i)), 2*coordX(i, 1 : numVertices(i)).^2 + 2*coordY(i, 1 : numVertices(i)).^2 + 6*coordZ(i, 1 : numVertices(i)).^2 - 2*r^2);
        coefficients(i, 9) = 8 * dot(coordX(i, 1 : numVertices(i)) .* coordY(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 10)= 8 * dot(coordX(i, 1 : numVertices(i)) .* coordZ(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 11)= 8 * dot(coordY(i, 1 : numVertices(i)) .* coordZ(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        diff = coordX(i, 1 : numVertices(i)).^2 + coordY(i, 1 : numVertices(i)).^2 + coordZ(i, 1 : numVertices(i)).^2 - r^2;
        coefficients(i, 12)= -4 * dot(coordX(i, 1 : numVertices(i)) .* diff, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 13)= -4 * dot(coordY(i, 1 : numVertices(i)) .* diff, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 14)= -4 * dot(coordZ(i, 1 : numVertices(i)) .* diff, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 15) = dot(confVertices(i, 1 : numVertices(i)), diff.^2);
    elseif shape.cylinder == primitiveType(i)
        % cylinder cost: ((a-x)^2+(b-y)^2+(c-z)^2-((a-x)*nx+(b-y)*ny+(c-z)*nz)^2-r^2)^2
        nx = fixedParameters(i, 1);
        ny = fixedParameters(i, 2);
        nz = fixedParameters(i, 3);
        n = numVertices(i);
        confSum = sum(confVertices(i, 1:n));
        a = coordX(i, 1:n);
        b = coordY(i, 1:n);
        c = coordZ(i, 1:n);
        r = fixedParameters(i, 7);
        dotp = a*nx + b*ny + c*nz;
        diff = a.^2 - dotp.^2 + b.^2 + c.^2 - r^2;
        % (nx^2 - 1)^2*x^4 +
        % (ny^2 - 1)^2*y^4 +
        % (nz^2 - 1)^2*z^4 +
        coefficients(i, 1) = (nx^2 - 1)^2*confSum;
        coefficients(i, 2) = (ny^2 - 1)^2*confSum;
        coefficients(i, 3) = (nz^2 - 1)^2*confSum;
        % 
        % (4*nx*ny*(nx^2 - 1))*x^3*y +
        % (4*nx*nz*(nx^2 - 1))*x^3*z +
        % (4*nx*ny*(ny^2 - 1))*x*y^3 +
        % (4*nx*nz*(nz^2 - 1))*x*z^3 +
        % (4*ny*nz*(ny^2 - 1))*y^3*z +
        % (4*ny*nz*(nz^2 - 1))*y*z^3 +
        coefficients(i, 4) = (4*nx*ny*(nx^2 - 1))*confSum;
        coefficients(i, 5) = (4*nx*nz*(nx^2 - 1))*confSum;
        coefficients(i, 6) = (4*nx*ny*(ny^2 - 1))*confSum;
        coefficients(i, 7) = (4*nx*nz*(nz^2 - 1))*confSum;
        coefficients(i, 8) = (4*ny*nz*(ny^2 - 1))*confSum;
        coefficients(i, 9) = (4*ny*nz*(nz^2 - 1))*confSum;
        % 
        % (2*(nx^2 - 1)*(2*a - 2*nx*(a*nx + b*ny + c*nz)))*x^3 +
        % (2*(ny^2 - 1)*(2*b - 2*ny*(a*nx + b*ny + c*nz)))*y^3 +
        % (2*(nz^2 - 1)*(2*c - 2*nz*(a*nx + b*ny + c*nz)))*z^3 +
        coefficients(i, 10) = dot(confVertices(i, 1 : n), (2*(nx^2 - 1)*(2*a - 2*nx*dotp)));
        coefficients(i, 11) = dot(confVertices(i, 1 : n), (2*(ny^2 - 1)*(2*b - 2*ny*dotp)));
        coefficients(i, 12) = dot(confVertices(i, 1 : n), (2*(nz^2 - 1)*(2*c - 2*nz*dotp)));
        % 
        % (2*(nx^2 - 1)*(ny^2 - 1) + 4*nx^2*ny^2)*x^2*y^2 +
        % (2*(nx^2 - 1)*(nz^2 - 1) + 4*nx^2*nz^2)*x^2*z^2 +
        % (2*(ny^2 - 1)*(nz^2 - 1) + 4*ny^2*nz^2)*y^2*z^2 +
        coefficients(i, 13) = (2*(nx^2 - 1)*(ny^2 - 1) + 4*nx^2*ny^2)*confSum;
        coefficients(i, 14) = (2*(nx^2 - 1)*(nz^2 - 1) + 4*nx^2*nz^2)*confSum;
        coefficients(i, 15) = (2*(ny^2 - 1)*(nz^2 - 1) + 4*ny^2*nz^2)*confSum;
        % 
        % (8*nx^2*ny*nz + 4*ny*nz*(nx^2 - 1))*x^2*y*z +
        % (8*nx*ny^2*nz + 4*nx*nz*(ny^2 - 1))*x*y^2*z +
        % (8*nx*ny*nz^2 + 4*nx*ny*(nz^2 - 1))*x*y*z^2 +
        coefficients(i, 16) = (8*nx^2*ny*nz + 4*ny*nz*(nx^2 - 1))*confSum;
        coefficients(i, 17) = (8*nx*ny^2*nz + 4*nx*nz*(ny^2 - 1))*confSum;
        coefficients(i, 18) = (8*nx*ny*nz^2 + 4*nx*ny*(nz^2 - 1))*confSum;
        % 
        % (2*(nx^2 - 1)*(2*b - 2*ny*(a*nx + b*ny + c*nz)) + 4*nx*ny*(2*a - 2*nx*(a*nx + b*ny + c*nz)))*x^2*y +
        % (2*(nx^2 - 1)*(2*c - 2*nz*(a*nx + b*ny + c*nz)) + 4*nx*nz*(2*a - 2*nx*(a*nx + b*ny + c*nz)))*x^2*z +
        % (2*(ny^2 - 1)*(2*a - 2*nx*(a*nx + b*ny + c*nz)) + 4*nx*ny*(2*b - 2*ny*(a*nx + b*ny + c*nz)))*x*y^2 +
        % (2*(nz^2 - 1)*(2*a - 2*nx*(a*nx + b*ny + c*nz)) + 4*nx*nz*(2*c - 2*nz*(a*nx + b*ny + c*nz)))*x*z^2 +
        % (2*(ny^2 - 1)*(2*c - 2*nz*(a*nx + b*ny + c*nz)) + 4*ny*nz*(2*b - 2*ny*(a*nx + b*ny + c*nz)))*y^2*z +
        % (2*(nz^2 - 1)*(2*b - 2*ny*(a*nx + b*ny + c*nz)) + 4*ny*nz*(2*c - 2*nz*(a*nx + b*ny + c*nz)))*y*z^2 +
        coefficients(i, 19) = dot(confVertices(i, 1 : n), (2*(nx^2 - 1)*(2*b - 2*ny*dotp) + 4*nx*ny*(2*a - 2*nx*dotp)));
        coefficients(i, 20) = dot(confVertices(i, 1 : n), (2*(nx^2 - 1)*(2*c - 2*nz*dotp) + 4*nx*nz*(2*a - 2*nx*dotp)));
        coefficients(i, 21) = dot(confVertices(i, 1 : n), (2*(ny^2 - 1)*(2*a - 2*nx*dotp) + 4*nx*ny*(2*b - 2*ny*dotp)));
        coefficients(i, 22) = dot(confVertices(i, 1 : n), (2*(nz^2 - 1)*(2*a - 2*nx*dotp) + 4*nx*nz*(2*c - 2*nz*dotp)));
        coefficients(i, 23) = dot(confVertices(i, 1 : n), (2*(ny^2 - 1)*(2*c - 2*nz*dotp) + 4*ny*nz*(2*b - 2*ny*dotp)));
        coefficients(i, 24) = dot(confVertices(i, 1 : n), (2*(nz^2 - 1)*(2*b - 2*ny*dotp) + 4*ny*nz*(2*c - 2*nz*dotp)));
        % 
        % (- 2*(nx^2 - 1)*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2 - r^2) + (2*a - 2*nx*(a*nx + b*ny + c*nz))^2)*x^2 +
        % (- 2*(ny^2 - 1)*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2 - r^2) + (2*b - 2*ny*(a*nx + b*ny + c*nz))^2)*y^2 +
        % (- 2*(nz^2 - 1)*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2 - r^2) + (2*c - 2*nz*(a*nx + b*ny + c*nz))^2)*z^2 +
        coefficients(i, 25) = dot(confVertices(i, 1 : n), (- 2*(nx^2 - 1)*diff + (2*a - 2*nx*dotp).^2));
        coefficients(i, 26) = dot(confVertices(i, 1 : n), (- 2*(ny^2 - 1)*diff + (2*b - 2*ny*dotp).^2));
        coefficients(i, 27) = dot(confVertices(i, 1 : n), (- 2*(nz^2 - 1)*diff + (2*c - 2*nz*dotp).^2));
        % 
        % (4*ny*nz*(2*a - 2*nx*(a*nx + b*ny + c*nz)) + 4*nx*nz*(2*b - 2*ny*(a*nx + b*ny + c*nz)) + 4*nx*ny*(2*c - 2*nz*(a*nx + b*ny + c*nz)))*x*y*z +
        coefficients(i, 28) = dot(confVertices(i, 1 : n), (4*ny*nz*(2*a - 2*nx*dotp) + 4*nx*nz*(2*b - 2*ny*dotp) + 4*nx*ny*(2*c - 2*nz*dotp)));
        %
        % (2*(2*a - 2*nx*(a*nx + b*ny + c*nz))*(2*b - 2*ny*(a*nx + b*ny + c*nz)) - 4*nx*ny*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2 - r^2))*x*y +
        % (2*(2*a - 2*nx*(a*nx + b*ny + c*nz))*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 4*nx*nz*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2 - r^2))*x*z +
        % (2*(2*b - 2*ny*(a*nx + b*ny + c*nz))*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 4*ny*nz*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2 - r^2))*y*z +
        coefficients(i, 29) = dot(confVertices(i, 1 : n), (2*(2*a - 2*nx*dotp).*(2*b - 2*ny*dotp) - 4*nx*ny*diff));
        coefficients(i, 30) = dot(confVertices(i, 1 : n), (2*(2*a - 2*nx*dotp).*(2*c - 2*nz*dotp) - 4*nx*nz*diff));
        coefficients(i, 31) = dot(confVertices(i, 1 : n), (2*(2*b - 2*ny*dotp).*(2*c - 2*nz*dotp) - 4*ny*nz*diff));
        % 
        % (-2*(2*a - 2*nx*(a*nx + b*ny + c*nz))*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2 - r^2))*x +
        % (-2*(2*b - 2*ny*(a*nx + b*ny + c*nz))*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2 - r^2))*y +
        % (-2*(2*c - 2*nz*(a*nx + b*ny + c*nz))*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2 - r^2))*z +
        coefficients(i, 32) = dot(confVertices(i, 1 : n), (-2*(2*a - 2*nx*dotp).*diff));
        coefficients(i, 33) = dot(confVertices(i, 1 : n), (-2*(2*b - 2*ny*dotp).*diff));
        coefficients(i, 34) = dot(confVertices(i, 1 : n), (-2*(2*c - 2*nz*dotp).*diff));
        % 
        % (a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2 - r^2)^2
        coefficients(i, 35) = dot(confVertices(i, 1 : n), diff.^2);
    elseif shape.cone == primitiveType(i)
        % cone cost: (s*((a-x)*nx+(b-y)*ny+(c-z)*nz)^2 - o*((a - x)^2 - ((a-x)*nx+(b-y)*ny+(c-z)*nz)^2 + (b - y)^2 + (c - z)^2))^2
        nx = fixedParameters(i, 1);
        ny = fixedParameters(i, 2);
        nz = fixedParameters(i, 3);
        n = numVertices(i);
        confSum = sum(confVertices(i, 1:n));
        a = coordX(i, 1:n);
        b = coordY(i, 1:n);
        c = coordZ(i, 1:n);
        o = cos(fixedParameters(i, 7))^2;
        s = sin(fixedParameters(i, 7))^2;
        dotp = a*nx + b*ny + c*nz;
        diff = a.^2 - dotp.^2 + b.^2 + c.^2;
        diff2 = (o*diff - s*dotp.^2);
        % (nx^2*s + o*(nx^2 - 1))^2*x^4 +
        % (ny^2*s + o*(ny^2 - 1))^2*y^4 +
        % (nz^2*s + o*(nz^2 - 1))^2*z^4 +
        coefficients(i, 1) = (nx^2*s + o*(nx^2 - 1))^2*confSum;
        coefficients(i, 2) = (ny^2*s + o*(ny^2 - 1))^2*confSum;
        coefficients(i, 3) = (nz^2*s + o*(nz^2 - 1))^2*confSum;
        % 
        % (2*(2*nx*ny*o + 2*nx*ny*s)*(nx^2*s + o*(nx^2 - 1)))*x^3*y +
        % (2*(2*nx*nz*o + 2*nx*nz*s)*(nx^2*s + o*(nx^2 - 1)))*x^3*z +
        % (2*(2*nx*ny*o + 2*nx*ny*s)*(ny^2*s + o*(ny^2 - 1)))*x*y^3 +
        % (2*(2*nx*nz*o + 2*nx*nz*s)*(nz^2*s + o*(nz^2 - 1)))*x*z^3 +
        % (2*(2*ny*nz*o + 2*ny*nz*s)*(ny^2*s + o*(ny^2 - 1)))*y^3*z +
        % (2*(2*ny*nz*o + 2*ny*nz*s)*(nz^2*s + o*(nz^2 - 1)))*y*z^3 +
        coefficients(i, 4) = (2*(2*nx*ny*o + 2*nx*ny*s)*(nx^2*s + o*(nx^2 - 1)))*confSum;
        coefficients(i, 5) = (2*(2*nx*nz*o + 2*nx*nz*s)*(nx^2*s + o*(nx^2 - 1)))*confSum;
        coefficients(i, 6) = (2*(2*nx*ny*o + 2*nx*ny*s)*(ny^2*s + o*(ny^2 - 1)))*confSum;
        coefficients(i, 7) = (2*(2*nx*nz*o + 2*nx*nz*s)*(nz^2*s + o*(nz^2 - 1)))*confSum;
        coefficients(i, 8) = (2*(2*ny*nz*o + 2*ny*nz*s)*(ny^2*s + o*(ny^2 - 1)))*confSum;
        coefficients(i, 9) = (2*(2*ny*nz*o + 2*ny*nz*s)*(nz^2*s + o*(nz^2 - 1)))*confSum;
        % 
        % (2*(o*(2*a - 2*nx*(a*nx + b*ny + c*nz)) - 2*nx*s*(a*nx + b*ny + c*nz))*(nx^2*s + o*(nx^2 - 1)))*x^3 +
        % (2*(o*(2*b - 2*ny*(a*nx + b*ny + c*nz)) - 2*ny*s*(a*nx + b*ny + c*nz))*(ny^2*s + o*(ny^2 - 1)))*y^3 +
        % (2*(o*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 2*nz*s*(a*nx + b*ny + c*nz))*(nz^2*s + o*(nz^2 - 1)))*z^3 +
        coefficients(i, 10) = dot(confVertices(i, 1 : n), (2*(o*(2*a - 2*nx*dotp) - 2*nx*s*dotp)*(nx^2*s + o*(nx^2 - 1))));
        coefficients(i, 11) = dot(confVertices(i, 1 : n), (2*(o*(2*b - 2*ny*dotp) - 2*ny*s*dotp)*(ny^2*s + o*(ny^2 - 1))));
        coefficients(i, 12) = dot(confVertices(i, 1 : n), (2*(o*(2*c - 2*nz*dotp) - 2*nz*s*dotp)*(nz^2*s + o*(nz^2 - 1))));
        % 
        % ((2*nx*ny*o + 2*nx*ny*s)^2 + 2*(nx^2*s + o*(nx^2 - 1))*(ny^2*s + o*(ny^2 - 1)))*x^2*y^2 +
        % ((2*nx*nz*o + 2*nx*nz*s)^2 + 2*(nx^2*s + o*(nx^2 - 1))*(nz^2*s + o*(nz^2 - 1)))*x^2*z^2 +
        % ((2*ny*nz*o + 2*ny*nz*s)^2 + 2*(ny^2*s + o*(ny^2 - 1))*(nz^2*s + o*(nz^2 - 1)))*y^2*z^2 +
        coefficients(i, 13) = ((2*nx*ny*o + 2*nx*ny*s)^2 + 2*(nx^2*s + o*(nx^2 - 1))*(ny^2*s + o*(ny^2 - 1)))*confSum;
        coefficients(i, 14) = ((2*nx*nz*o + 2*nx*nz*s)^2 + 2*(nx^2*s + o*(nx^2 - 1))*(nz^2*s + o*(nz^2 - 1)))*confSum;
        coefficients(i, 15) = ((2*ny*nz*o + 2*ny*nz*s)^2 + 2*(ny^2*s + o*(ny^2 - 1))*(nz^2*s + o*(nz^2 - 1)))*confSum;
        % 
        % (2*(2*ny*nz*o + 2*ny*nz*s)*(nx^2*s + o*(nx^2 - 1)) + 2*(2*nx*ny*o + 2*nx*ny*s)*(2*nx*nz*o + 2*nx*nz*s))*x^2*y*z +
        % (2*(2*nx*nz*o + 2*nx*nz*s)*(ny^2*s + o*(ny^2 - 1)) + 2*(2*nx*ny*o + 2*nx*ny*s)*(2*ny*nz*o + 2*ny*nz*s))*x*y^2*z +
        % (2*(2*nx*ny*o + 2*nx*ny*s)*(nz^2*s + o*(nz^2 - 1)) + 2*(2*nx*nz*o + 2*nx*nz*s)*(2*ny*nz*o + 2*ny*nz*s))*x*y*z^2 +
        coefficients(i, 16) = (2*(2*ny*nz*o + 2*ny*nz*s)*(nx^2*s + o*(nx^2 - 1)) + 2*(2*nx*ny*o + 2*nx*ny*s)*(2*nx*nz*o + 2*nx*nz*s))*confSum;
        coefficients(i, 17) = (2*(2*nx*nz*o + 2*nx*nz*s)*(ny^2*s + o*(ny^2 - 1)) + 2*(2*nx*ny*o + 2*nx*ny*s)*(2*ny*nz*o + 2*ny*nz*s))*confSum;
        coefficients(i, 18) = (2*(2*nx*ny*o + 2*nx*ny*s)*(nz^2*s + o*(nz^2 - 1)) + 2*(2*nx*nz*o + 2*nx*nz*s)*(2*ny*nz*o + 2*ny*nz*s))*confSum;
        % 
        % (2*(o*(2*a - 2*nx*(a*nx + b*ny + c*nz)) - 2*nx*s*(a*nx + b*ny + c*nz))*(2*nx*ny*o + 2*nx*ny*s) + 2*(o*(2*b - 2*ny*(a*nx + b*ny + c*nz)) - 2*ny*s*(a*nx + b*ny + c*nz))*(nx^2*s + o*(nx^2 - 1)))*x^2*y +
        % (2*(o*(2*a - 2*nx*(a*nx + b*ny + c*nz)) - 2*nx*s*(a*nx + b*ny + c*nz))*(2*nx*nz*o + 2*nx*nz*s) + 2*(o*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 2*nz*s*(a*nx + b*ny + c*nz))*(nx^2*s + o*(nx^2 - 1)))*x^2*z +
        % (2*(o*(2*b - 2*ny*(a*nx + b*ny + c*nz)) - 2*ny*s*(a*nx + b*ny + c*nz))*(2*nx*ny*o + 2*nx*ny*s) + 2*(o*(2*a - 2*nx*(a*nx + b*ny + c*nz)) - 2*nx*s*(a*nx + b*ny + c*nz))*(ny^2*s + o*(ny^2 - 1)))*x*y^2 +
        % (2*(o*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 2*nz*s*(a*nx + b*ny + c*nz))*(2*nx*nz*o + 2*nx*nz*s) + 2*(o*(2*a - 2*nx*(a*nx + b*ny + c*nz)) - 2*nx*s*(a*nx + b*ny + c*nz))*(nz^2*s + o*(nz^2 - 1)))*x*z^2 +
        % (2*(o*(2*b - 2*ny*(a*nx + b*ny + c*nz)) - 2*ny*s*(a*nx + b*ny + c*nz))*(2*ny*nz*o + 2*ny*nz*s) + 2*(o*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 2*nz*s*(a*nx + b*ny + c*nz))*(ny^2*s + o*(ny^2 - 1)))*y^2*z +
        % (2*(o*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 2*nz*s*(a*nx + b*ny + c*nz))*(2*ny*nz*o + 2*ny*nz*s) + 2*(o*(2*b - 2*ny*(a*nx + b*ny + c*nz)) - 2*ny*s*(a*nx + b*ny + c*nz))*(nz^2*s + o*(nz^2 - 1)))*y*z^2 +
        coefficients(i, 19) = dot(confVertices(i, 1 : n), (2*(o*(2*a - 2*nx*dotp) - 2*nx*s*dotp)*(2*nx*ny*o + 2*nx*ny*s) + 2*(o*(2*b - 2*ny*dotp) - 2*ny*s*dotp)*(nx^2*s + o*(nx^2 - 1))));
        coefficients(i, 20) = dot(confVertices(i, 1 : n), (2*(o*(2*a - 2*nx*dotp) - 2*nx*s*dotp)*(2*nx*nz*o + 2*nx*nz*s) + 2*(o*(2*c - 2*nz*dotp) - 2*nz*s*dotp)*(nx^2*s + o*(nx^2 - 1))));
        coefficients(i, 21) = dot(confVertices(i, 1 : n), (2*(o*(2*b - 2*ny*dotp) - 2*ny*s*dotp)*(2*nx*ny*o + 2*nx*ny*s) + 2*(o*(2*a - 2*nx*dotp) - 2*nx*s*dotp)*(ny^2*s + o*(ny^2 - 1))));
        coefficients(i, 22) = dot(confVertices(i, 1 : n), (2*(o*(2*c - 2*nz*dotp) - 2*nz*s*dotp)*(2*nx*nz*o + 2*nx*nz*s) + 2*(o*(2*a - 2*nx*dotp) - 2*nx*s*dotp)*(nz^2*s + o*(nz^2 - 1))));
        coefficients(i, 23) = dot(confVertices(i, 1 : n), (2*(o*(2*b - 2*ny*dotp) - 2*ny*s*dotp)*(2*ny*nz*o + 2*ny*nz*s) + 2*(o*(2*c - 2*nz*dotp) - 2*nz*s*dotp)*(ny^2*s + o*(ny^2 - 1))));
        coefficients(i, 24) = dot(confVertices(i, 1 : n), (2*(o*(2*c - 2*nz*dotp) - 2*nz*s*dotp)*(2*ny*nz*o + 2*ny*nz*s) + 2*(o*(2*b - 2*ny*dotp) - 2*ny*s*dotp)*(nz^2*s + o*(nz^2 - 1))));
        % 
        % ((o*(2*a - 2*nx*(a*nx + b*ny + c*nz)) - 2*nx*s*(a*nx + b*ny + c*nz))^2 - 2*(o*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2) - s*(a*nx + b*ny + c*nz)^2)*(nx^2*s + o*(nx^2 - 1)))*x^2 +
        % ((o*(2*b - 2*ny*(a*nx + b*ny + c*nz)) - 2*ny*s*(a*nx + b*ny + c*nz))^2 - 2*(o*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2) - s*(a*nx + b*ny + c*nz)^2)*(ny^2*s + o*(ny^2 - 1)))*y^2 +
        % ((o*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 2*nz*s*(a*nx + b*ny + c*nz))^2 - 2*(o*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2) - s*(a*nx + b*ny + c*nz)^2)*(nz^2*s + o*(nz^2 - 1)))*z^2 +
        coefficients(i, 25) = dot(confVertices(i, 1 : n), ((o*(2*a - 2*nx*dotp) - 2*nx*s*dotp).^2 - 2*diff2*(nx^2*s + o*(nx^2 - 1))));
        coefficients(i, 26) = dot(confVertices(i, 1 : n), ((o*(2*b - 2*ny*dotp) - 2*ny*s*dotp).^2 - 2*diff2*(ny^2*s + o*(ny^2 - 1))));
        coefficients(i, 27) = dot(confVertices(i, 1 : n), ((o*(2*c - 2*nz*dotp) - 2*nz*s*dotp).^2 - 2*diff2*(nz^2*s + o*(nz^2 - 1))));
        % 
        % (2*(o*(2*a - 2*nx*(a*nx + b*ny + c*nz)) - 2*nx*s*(a*nx + b*ny + c*nz))*(2*ny*nz*o + 2*ny*nz*s) + 2*(o*(2*b - 2*ny*(a*nx + b*ny + c*nz)) - 2*ny*s*(a*nx + b*ny + c*nz))*(2*nx*nz*o + 2*nx*nz*s) + 2*(o*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 2*nz*s*(a*nx + b*ny + c*nz))*(2*nx*ny*o + 2*nx*ny*s))*x*y*z +
        coefficients(i, 28) = dot(confVertices(i, 1 : n), (2*(o*(2*a - 2*nx*dotp) - 2*nx*s*dotp)*(2*ny*nz*o + 2*ny*nz*s) + 2*(o*(2*b - 2*ny*dotp) - 2*ny*s*dotp)*(2*nx*nz*o + 2*nx*nz*s) + 2*(o*(2*c - 2*nz*dotp) - 2*nz*s*dotp)*(2*nx*ny*o + 2*nx*ny*s)));
        % 
        % (- 2*(o*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2) - s*(a*nx + b*ny + c*nz)^2)*(2*nx*ny*o + 2*nx*ny*s) + 2*(o*(2*a - 2*nx*(a*nx + b*ny + c*nz)) - 2*nx*s*(a*nx + b*ny + c*nz))*(o*(2*b - 2*ny*(a*nx + b*ny + c*nz)) - 2*ny*s*(a*nx + b*ny + c*nz)))*x*y +
        % (- 2*(o*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2) - s*(a*nx + b*ny + c*nz)^2)*(2*nx*nz*o + 2*nx*nz*s) + 2*(o*(2*a - 2*nx*(a*nx + b*ny + c*nz)) - 2*nx*s*(a*nx + b*ny + c*nz))*(o*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 2*nz*s*(a*nx + b*ny + c*nz)))*x*z +
        % (- 2*(o*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2) - s*(a*nx + b*ny + c*nz)^2)*(2*ny*nz*o + 2*ny*nz*s) + 2*(o*(2*b - 2*ny*(a*nx + b*ny + c*nz)) - 2*ny*s*(a*nx + b*ny + c*nz))*(o*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 2*nz*s*(a*nx + b*ny + c*nz)))*y*z +
        coefficients(i, 29) = dot(confVertices(i, 1 : n), (- 2*diff2*(2*nx*ny*o + 2*nx*ny*s) + 2*(o*(2*a - 2*nx*dotp) - 2*nx*s*dotp).*(o*(2*b - 2*ny*dotp) - 2*ny*s*dotp)));
        coefficients(i, 30) = dot(confVertices(i, 1 : n), (- 2*diff2*(2*nx*nz*o + 2*nx*nz*s) + 2*(o*(2*a - 2*nx*dotp) - 2*nx*s*dotp).*(o*(2*c - 2*nz*dotp) - 2*nz*s*dotp)));
        coefficients(i, 31) = dot(confVertices(i, 1 : n), (- 2*diff2*(2*ny*nz*o + 2*ny*nz*s) + 2*(o*(2*b - 2*ny*dotp) - 2*ny*s*dotp).*(o*(2*c - 2*nz*dotp) - 2*nz*s*dotp)));
        % 
        % (-2*(o*(2*a - 2*nx*(a*nx + b*ny + c*nz)) - 2*nx*s*(a*nx + b*ny + c*nz))*(o*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2) - s*(a*nx + b*ny + c*nz)^2))*x +
        % (-2*(o*(2*b - 2*ny*(a*nx + b*ny + c*nz)) - 2*ny*s*(a*nx + b*ny + c*nz))*(o*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2) - s*(a*nx + b*ny + c*nz)^2))*y +
        % (-2*(o*(2*c - 2*nz*(a*nx + b*ny + c*nz)) - 2*nz*s*(a*nx + b*ny + c*nz))*(o*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2) - s*(a*nx + b*ny + c*nz)^2))*z +
        coefficients(i, 32) = dot(confVertices(i, 1 : n), (-2*(o*(2*a - 2*nx*dotp) - 2*nx*s*dotp).*diff2));
        coefficients(i, 33) = dot(confVertices(i, 1 : n), (-2*(o*(2*b - 2*ny*dotp) - 2*ny*s*dotp).*diff2));
        coefficients(i, 34) = dot(confVertices(i, 1 : n), (-2*(o*(2*c - 2*nz*dotp) - 2*nz*s*dotp).*diff2));
        % 
        % (o*(a^2 - (a*nx + b*ny + c*nz)^2 + b^2 + c^2) - s*(a*nx + b*ny + c*nz)^2)^2
        coefficients(i, 35) = dot(confVertices(i, 1 : n), diff2.^2);
    end
end


initialFittingError = FittingError(inputParameters, numPrimitives, primitiveType, coefficients);

% optimization
options = optimset('Display', 'notify-detailed', 'Algorithm', 'interior-point', 'MaxFunEvals', Inf, 'MaxIter', maxIterNum, 'DerivativeCheck', 'off');
fobj = @(inputParameters)FittingError(inputParameters, numPrimitives, primitiveType, coefficients);
fcon = @(inputParameters)ConstrainPoint(inputParameters, numPrimitives, constraints, numConstraints, fixedParameters);
tic;
[outputParameters, exitFittingError, exitflag] = fmincon(fobj, inputParameters, [], [], [], [], [], [], fcon, options);
toc;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [energy] = FittingError(inputParameters, numPrimitives, primitiveType, coefficients)

LoadNameMap

% energy
energy = 0.0;

% fitting error
for i = 1 : numPrimitives
    x = inputParameters(i, 4);
    y = inputParameters(i, 5);
    z = inputParameters(i, 6);
    if shape.sphere == primitiveType(i)
        energy = energy + coefficients(i, 1)*(x^4 + y^4 + z^4) ...
            + coefficients(i, 2)*(x^2*y^2 + x^2*z^2 + y^2*z^2) ...
            + coefficients(i, 3)*(x^3 + x*y^2 + x*z^2) ...
            + coefficients(i, 4)*(x^2*y + y^3 + y*z^2) ...
            + coefficients(i, 5)*(x^2*z + y^2*z + z^3) ...
            + coefficients(i, 6)* x^2 ...
            + coefficients(i, 7)* y^2 ...
            + coefficients(i, 8)* z^2 ...
            + coefficients(i, 9)* x*y ...
            + coefficients(i, 10)* x*z ...
            + coefficients(i, 11)* y*z ...
            + coefficients(i, 12)* x...
            + coefficients(i, 13)* y...
            + coefficients(i, 14)* z...
            + coefficients(i, 15);
    elseif  shape.cylinder == primitiveType(i) || shape.cone == primitiveType(i)
        energy = energy + ...
            coefficients(i, 1)*x^4 + ...
            coefficients(i, 2)*y^4 + ...
            coefficients(i, 3)*z^4 + ...
 ...
            coefficients(i, 4)*x^3*y + ...
            coefficients(i, 5)*x^3*z + ...
            coefficients(i, 6)*x*y^3 + ...
            coefficients(i, 7)*x*z^3 + ...
            coefficients(i, 8)*y^3*z + ...
            coefficients(i, 9)*y*z^3 + ...
 ...
            coefficients(i, 10)*x^3 + ...
            coefficients(i, 11)*y^3 + ...
            coefficients(i, 12)*z^3 + ...
 ...
            coefficients(i, 13)*x^2*y^2 + ...
            coefficients(i, 14)*x^2*z^2 + ...
            coefficients(i, 15)*y^2*z^2 + ...
 ...
            coefficients(i, 16)*x^2*y*z + ...
            coefficients(i, 17)*x*y^2*z + ...
            coefficients(i, 18)*x*y*z^2 + ...
 ...
            coefficients(i, 19)*x^2*y + ...
            coefficients(i, 20)*x^2*z + ...
            coefficients(i, 21)*x*y^2 + ...
            coefficients(i, 22)*x*z^2 + ...
            coefficients(i, 23)*y^2*z + ...
            coefficients(i, 24)*y*z^2 + ...
 ...
            coefficients(i, 25)*x^2 + ...
            coefficients(i, 26)*y^2 + ...
            coefficients(i, 27)*z^2 + ...
 ...
            coefficients(i, 28)*x*y*z + ...
 ...
            coefficients(i, 29)*x*y + ...
            coefficients(i, 30)*x*z + ...
            coefficients(i, 31)*y*z + ...
 ...
            coefficients(i, 32)*x + ...
            coefficients(i, 33)*y + ...
            coefficients(i, 34)*z + ...
 ...
            coefficients(i, 35);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c, ceq, gradc, gradceq] = ConstrainPoint(inputParameters, numPrimitives, constraints, numConstraints, fixedParameters)

LoadNameMap

% Nonlinear inequality constraints
c = [];

% Nonlinear equality constraints
ceq = zeros(numConstraints, 1);

for i = 1 : numConstraints
    normal = fixedParameters(constraints(i, 2), 1 : 3);
    if sum(normal) == 0
        normal = fixedParameters(constraints(i, 3), 1 : 3);
    end
    vector = inputParameters(constraints(i, 2), 4 : 6)-inputParameters(constraints(i, 3), 4 : 6);
    ceq(i) = vector*vector'-(dot(vector, normal))^2;
end


gradc = [];

gradceq = zeros(size(inputParameters, 2) * numPrimitives, numConstraints);

end