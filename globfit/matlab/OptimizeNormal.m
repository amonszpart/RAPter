function [outputParameters, initialFittingError, exitFittingError, exitFlag] = OptimizeNormal(inputParameters, maxIterNum, numVertices, primitiveType, coordX, coordY, coordZ, normalX, normalY, normalZ, confVertices, constraints)

LoadNameMap

numPrimitives = size(numVertices, 1);
numConstraints = size(constraints, 1);

fprintf('Num Primitives: %d\n', numPrimitives); 
fprintf('Num Constraints: %d\n', numConstraints); 


if 0 < numConstraints
    constraints(1 : numConstraints, 2 : 5) = constraints(1 : numConstraints, 2 : 5) + 1;
end

fprintf('Computing parallel collapse map\n'); 
% compute parallel collapse map
collapseMap = 1 : 1 : numPrimitives;
for i = 1 : numConstraints
    if relation.parallel == constraints(i, 1)
        collapseMap(constraints(i, 3)) = constraints(i, 2);
    end
end

fprintf('Normalize\n'); 
% normalization
for i = 1 : numPrimitives
    if shape.plane == primitiveType(i)
        inputParameters(i, 7) = inputParameters(i, 7) / norm(inputParameters(i, 1 : 3));
    end
    if (shape.plane == primitiveType(i)) || (shape.cylinder == primitiveType(i)) || (shape.cone == primitiveType(i))
        inputParameters(i, 1 : 3) = inputParameters(i, 1 : 3) / norm(inputParameters(i, 1 : 3));
    end
    
    %  re-orientation for planes
    if shape.plane == primitiveType(i) && dot(inputParameters(i, 1 : 3), inputParameters(collapseMap(i), 1 : 3)) < 0
        inputParameters(i, 1:3) = - inputParameters(i, 1:3);
        inputParameters(i, 7) = - inputParameters(i, 7);
    end
    
    %  re-orientation for cones
    if shape.cone == primitiveType(i) && dot(inputParameters(i, 1 : 3), inputParameters(collapseMap(i), 1 : 3)) < 0
        inputParameters(i, 1:3) = - inputParameters(i, 1:3);
        normalX = - normalX;
        normalY = - normalY;
        normalZ = - normalZ;
    end
end

% expression expansion
coefficients = zeros(numPrimitives, 10);
for i = 1 : numPrimitives
    if shape.plane == primitiveType(i)
        % plane cost: Euclidean distance
        % (ax + by + cz + d)^2: (ax)^2 + (by)^2 + (cz)^2 + 2abxy + 2bcyz + 2caxz + 2adx + 2bdy + 2cdz + d^2
        % d = inputParameters(i, 7);
        coefficients(i, 1) = dot(coordX(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 2) = dot(coordY(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 3) = dot(coordZ(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 4) = 2 * dot(coordX(i, 1 : numVertices(i)) .* coordY(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 5) = 2 * dot(coordY(i, 1 : numVertices(i)) .* coordZ(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 6) = 2 * dot(coordZ(i, 1 : numVertices(i)) .* coordX(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 7) = 2 * dot(coordX(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));% * d;
        coefficients(i, 8) = 2 * dot(coordY(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));% * d;
        coefficients(i, 9) = 2 * dot(coordZ(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));% * d;
        coefficients(i, 10) = sum(confVertices(i, 1 : numVertices(i)));% * d^2;
    elseif shape.cylinder == primitiveType(i)
        % cylinder cost: Normal consistency
        % (ax + by + cz)^2: (ax)^2 + (by)^2 + (cz)^2 + 2abxy + 2bcyz + 2caxz
        coefficients(i, 1) = dot(normalX(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 2) = dot(normalY(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 3) = dot(normalZ(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 4) = 2 * dot(normalX(i, 1 : numVertices(i)) .* normalY(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 5) = 2 * dot(normalY(i, 1 : numVertices(i)) .* normalZ(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 6) = 2 * dot(normalZ(i, 1 : numVertices(i)) .* normalX(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
    elseif shape.cone == primitiveType(i)
        % cone cost: Normal consistency
        % (ax + by + cz + d)^2: (ax)^2 + (by)^2 + (cz)^2 + 2abxy + 2bcyz + 2caxz + 2adx + 2bdy + 2cdz + d^2
        d = -cos(pi/2-inputParameters(i, 7));
        coefficients(i, 1) = dot(normalX(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 2) = dot(normalY(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 3) = dot(normalZ(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 4) = 2 * dot(normalX(i, 1 : numVertices(i)) .* normalY(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 5) = 2 * dot(normalY(i, 1 : numVertices(i)) .* normalZ(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 6) = 2 * dot(normalZ(i, 1 : numVertices(i)) .* normalX(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 7) = 2 * dot(normalX(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i))) * d;
        coefficients(i, 8) = 2 * dot(normalY(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i))) * d;
        coefficients(i, 9) = 2 * dot(normalZ(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i))) * d;
        coefficients(i, 10) = sum(confVertices(i, 1 : numVertices(i))) * d^2;
    end
end

initialFittingError = FittingError(inputParameters, numPrimitives, primitiveType, coefficients, collapseMap);

% optimization
options = optimset('Display', 'notify-detailed', 'Algorithm', 'interior-point', 'MaxFunEvals', Inf, 'MaxIter', maxIterNum, 'DerivativeCheck', 'off');
fobj = @(inputParameters)FittingError(inputParameters, numPrimitives, primitiveType, coefficients, collapseMap);
fcon = @(inputParameters)ConstrainNormal(inputParameters, numPrimitives, constraints, numConstraints, collapseMap);
tic;
[outputParameters, exitFittingError, exitFlag] = fmincon(fobj, inputParameters, [], [], [], [], [], [], fcon, options);
toc;

% copy parameters to collapsed primitives
for i = 1 : numPrimitives
    outputParameters(i, 1 : 3) = outputParameters(collapseMap(i), 1 : 3);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [energy] = FittingError(inputParameters, numPrimitives, primitiveType, coefficients, collapseMap)

LoadNameMap

% energy
energy = 0.0;

% fitting error
for i = 1 : numPrimitives
    if shape.plane == primitiveType(i)
        energy = energy + coefficients(i, 1) * inputParameters(collapseMap(i), 1) ^ 2 ...
            + coefficients(i, 2) * inputParameters(collapseMap(i), 2) ^ 2 ...
            + coefficients(i, 3) * inputParameters(collapseMap(i), 3) ^ 2 ...
            + coefficients(i, 4) * inputParameters(collapseMap(i), 1) * inputParameters(collapseMap(i), 2) ...
            + coefficients(i, 5) * inputParameters(collapseMap(i), 2) * inputParameters(collapseMap(i), 3) ...
            + coefficients(i, 6) * inputParameters(collapseMap(i), 3) * inputParameters(collapseMap(i), 1) ...
            + coefficients(i, 7) * inputParameters(collapseMap(i), 1) * inputParameters(i, 7) ...
            + coefficients(i, 8) * inputParameters(collapseMap(i), 2) * inputParameters(i, 7) ...
            + coefficients(i, 9) * inputParameters(collapseMap(i), 3) * inputParameters(i, 7) ...
            + coefficients(i, 10)* inputParameters(i, 7)^2;
    elseif shape.cylinder == primitiveType(i)
        energy = energy + coefficients(i, 1) * inputParameters(collapseMap(i), 1) ^ 2 ...
            + coefficients(i, 2) * inputParameters(collapseMap(i), 2) ^ 2 ...
            + coefficients(i, 3) * inputParameters(collapseMap(i), 3) ^ 2 ...
            + coefficients(i, 4) * inputParameters(collapseMap(i), 1) * inputParameters(collapseMap(i), 2) ...
            + coefficients(i, 5) * inputParameters(collapseMap(i), 2) * inputParameters(collapseMap(i), 3) ...
            + coefficients(i, 6) * inputParameters(collapseMap(i), 3) * inputParameters(collapseMap(i), 1);
    elseif shape.cone == primitiveType(i)
        energy = energy + coefficients(i, 1) * inputParameters(collapseMap(i), 1) ^ 2 ...
            + coefficients(i, 2) * inputParameters(collapseMap(i), 2) ^ 2 ...
            + coefficients(i, 3) * inputParameters(collapseMap(i), 3) ^ 2 ...
            + coefficients(i, 4) * inputParameters(collapseMap(i), 1) * inputParameters(collapseMap(i), 2) ...
            + coefficients(i, 5) * inputParameters(collapseMap(i), 2) * inputParameters(collapseMap(i), 3) ...
            + coefficients(i, 6) * inputParameters(collapseMap(i), 3) * inputParameters(collapseMap(i), 1) ...
            + coefficients(i, 7) * inputParameters(collapseMap(i), 1) ...
            + coefficients(i, 8) * inputParameters(collapseMap(i), 2) ...
            + coefficients(i, 9) * inputParameters(collapseMap(i), 3) ...
            + coefficients(i, 10);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c, ceq, gradc, gradceq] = ConstrainNormal(inputParameters, numPrimitives, constraints, numConstraints, collapseMap)

LoadNameMap

% Nonlinear inequality constraints
c = [];

% Nonlinear equality constraints
ceq = zeros(numPrimitives + numConstraints, 1);
for i = 1 : numPrimitives
    if i == collapseMap(i)
        ceq(i) = dot(inputParameters(i, 1 : 3), inputParameters(i, 1 : 3)) - 1.0;
    end
end

for i = 1 : numConstraints
    if relation.parallel == constraints(i, 1)
        ceq(numPrimitives+i) = 1.0 - dot(inputParameters(collapseMap(constraints(i, 2)), 1 : 3), inputParameters(collapseMap(constraints(i, 3)), 1 : 3)) ^ 2;
    elseif relation.orthogonal == constraints(i, 1)
        ceq(numPrimitives+i) = dot(inputParameters(collapseMap(constraints(i, 2)), 1 : 3), inputParameters(collapseMap(constraints(i, 3)), 1 : 3)) ^ 2;
    elseif relation.equalAngle == constraints(i, 1)
        ceq(numPrimitives+i) = dot(inputParameters(collapseMap(constraints(i, 2)), 1 : 3), inputParameters(collapseMap(constraints(i, 3)), 1 : 3)) ^ 2 ...
            - dot(inputParameters(collapseMap(constraints(i, 4)), 1 : 3), inputParameters(collapseMap(constraints(i, 5)), 1 : 3)) ^ 2;
    end
end

gradc = [];

gradceq = zeros(size(inputParameters, 2) * numPrimitives, numPrimitives + numConstraints);

end
