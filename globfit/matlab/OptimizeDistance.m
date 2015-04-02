function [outputParameters, initialFittingError, exitFittingError, exitflag] = OptimizeDistance(inputParameters, maxIterNum, numVertices, primitiveType, coordX, coordY, coordZ, normalX, normalY, normalZ, confVertices, constraints)

LoadNameMap

numPrimitives = size(numVertices, 1);
numConstraints = size(constraints, 1);
if 0 < numConstraints
    constraints(1 : numConstraints, 2 : 5) = constraints(1 : numConstraints, 2 : 5) + 1;
end

% compute parallel collapse map
collapseMap = 1 : 1 : numPrimitives;
for i = 1 : numConstraints
    if relation.coPlanar == constraints(i, 1)
        collapseMap(constraints(i, 3)) = constraints(i, 2);
    end
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
coefficients = zeros(numPrimitives, 3);
for i = 1 : numPrimitives
    if shape.plane == primitiveType(i)
        % plane cost: Euclidean distance
        % (ax + by + cz + d)^2: ((ax)^2 + (by)^2 + (cz)^2 + 2abxy + 2bcyz + 2caxz) + (2ax + 2by + 2cz)d + d^2
        a = fixedParameters(collapseMap(i), 1);
        b = fixedParameters(collapseMap(i), 2);
        c = fixedParameters(collapseMap(i), 3);
        coefficients(i, 1) = dot(coordX(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i))) * a ^ 2 ...
            + dot(coordY(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i))) * b ^ 2 ...
            + dot(coordZ(i, 1 : numVertices(i)) .^ 2, confVertices(i, 1 : numVertices(i))) * c ^ 2 ...
            + 2 * dot(coordX(i, 1 : numVertices(i)) .* coordY(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i))) * a * b ...
            + 2 * dot(coordY(i, 1 : numVertices(i)) .* coordZ(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i))) * b * c ...
            + 2 * dot(coordZ(i, 1 : numVertices(i)) .* coordX(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i))) * c * a;
        coefficients(i, 2) = 2 * dot(coordX(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i)))* a ...
            + 2 * dot(coordY(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i))) * b ...
            + 2 * dot(coordZ(i, 1 : numVertices(i)), confVertices(i, 1 : numVertices(i))) * c;
        coefficients(i, 3) = sum(confVertices(i, 1 : numVertices(i)));
    end
end

initialFittingError = FittingError(inputParameters, numPrimitives, primitiveType, coefficients, collapseMap);

% optimization
options = optimset('Display', 'iter-detailed', 'Algorithm', 'interior-point', 'MaxFunEvals', Inf, 'MaxIter', maxIterNum, 'DerivativeCheck', 'off');
fobj = @(inputParameters)FittingError(inputParameters, numPrimitives, primitiveType, coefficients, collapseMap);
fcon = @(inputParameters)ConstrainNormal(inputParameters, numPrimitives, constraints, numConstraints, collapseMap);
tic;
[outputParameters, exitFittingError, exitflag] = fmincon(fobj, inputParameters, [], [], [], [], [], [], fcon, options);
toc;

% copy parameters to collapsed primitives
for i = 1 : numPrimitives
    if shape.plane == primitiveType(i)
        outputParameters(i, 7) = outputParameters(collapseMap(i), 7);
        if dot(outputParameters(i, 1 : 3), outputParameters(collapseMap(i), 1 : 3)) < 0
            outputParameters(i, 7) = - outputParameters(i, 7);
        end
    end
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
        energy = energy + coefficients(i, 1)  ...
            + coefficients(i, 2) * inputParameters(collapseMap(i), 7) ...
            + coefficients(i, 3) * inputParameters(collapseMap(i), 7) ^ 2;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c, ceq, gradc, gradceq] = ConstrainNormal(inputParameters, numPrimitives, constraints, numConstraints, collapseMap)

LoadNameMap

% Nonlinear inequality constraints
c = [];

% Nonlinear equality constraints
ceq = zeros(numConstraints, 1);
for i = 1 : numConstraints
    if relation.equalLength == constraints(i, 1)
        d12 = inputParameters(collapseMap(constraints(i, 2)), 7) - inputParameters(collapseMap(constraints(i, 3)), 7);
        d34 = inputParameters(collapseMap(constraints(i, 4)), 7) - inputParameters(collapseMap(constraints(i, 5)), 7);
        ceq(i) = d12^2-d34^2;
    end
end


gradc = [];

gradceq = zeros(size(inputParameters, 2) * numPrimitives, numConstraints);

end