function [outputParameters, initialFittingError, exitFittingError, exitflag] = OptimizeRadius(inputParameters, maxIterNum, numVertices, primitiveType, coordX, coordY, coordZ, normalX, normalY, normalZ, confVertices, constraints)

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
    if shape.sphere == primitiveType(i) || shape.cylinder == primitiveType(i)
        if shape.sphere == primitiveType(i)
            vector = [coordX(i, 1:numVertices(i)) - fixedParameters(i, 4); coordY(i, 1:numVertices(i)) - fixedParameters(i, 5); coordZ(i, 1:numVertices(i)) - fixedParameters(i, 6)];
            vector = vector.^2;
            squaredDistance = sum(vector(:,:));
        else
            vector = [coordX(i, 1:numVertices(i)) - fixedParameters(i, 4); coordY(i, 1:numVertices(i)) - fixedParameters(i, 5); coordZ(i, 1:numVertices(i)) - fixedParameters(i, 6)];
            dotProduct = ((vector') * (fixedParameters(i, 1:3)'))';
            vector = vector.^2;
            squaredDistance = sum(vector(:,:)) - dotProduct.^2;
        end
        coefficients(i, 1) = dot(squaredDistance, confVertices(i, 1 : numVertices(i)));
        coefficients(i, 2) = -2*dot(sqrt(squaredDistance), confVertices(i, 1 : numVertices(i)));
        coefficients(i, 3) = sum(confVertices(i, 1 : numVertices(i)));
    end
end

initialFittingError = FittingError(inputParameters, numPrimitives, primitiveType, coefficients, collapseMap);

% optimization
options = optimset('Display', 'notify-detailed', 'LargeScale', 'off', 'MaxFunEvals', Inf, 'MaxIter', maxIterNum, 'DerivativeCheck', 'off');
fobj = @(inputParameters)FittingError(inputParameters, numPrimitives, primitiveType, coefficients, collapseMap);
tic;
[outputParameters, exitFittingError, exitflag] = fminunc(fobj, inputParameters, options);
toc;

% copy parameters to collapsed primitives
for i = 1 : numPrimitives
    if shape.sphere == primitiveType(i) || shape.cylinder == primitiveType(i) || shape.cone == primitiveType(i)
        outputParameters(i, 7) = outputParameters(collapseMap(i), 7);
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
    if shape.sphere == primitiveType(i) || shape.cylinder == primitiveType(i)
        radius = inputParameters(collapseMap(i), 7);
        energy = energy + coefficients(i, 1) + coefficients(i, 2) * radius + coefficients(i, 3) * radius^2;
    end
end

end