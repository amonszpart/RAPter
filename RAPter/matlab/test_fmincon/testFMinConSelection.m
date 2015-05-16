function x = testFMinConSelection(p)

USE_GA = 1;

%x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
if (~exist('p','var') )
    %p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/_1rect_input2_20140704_0852'
    %p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/_2rects'
    %p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/_rot_rects_noise2'
    p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/2rects_noise_0.003/';
end


%% objective function
%model       = gurobi_read( [ p filesep 'model.lp' ] );
global qo;
global Qo;
global A;
if ( ~exist('A','var') || size(A,1) == 0 )
    qo = readSparseMatrix( [ p filesep 'qo.cpp.txt' ])';
    Qo = readSparseMatrix( [ p filesep 'Qo.cpp.txt' ]); 
    A = readSparseMatrix( [ p filesep 'A.cpp.txt' ]);
    size(A)
    A = -unique(A,'rows');
    size(A)
end
sum(sum(A))
%qo = model.obj';
%Qo = model.Q;
%A = -model.A;

if ( USE_GA ) 
    %objectivef  = @(x) x * Qo * x' + qo * x';
    objectivef  = @(x) sum(x);
    %objectivef  = @(x) qo * x';
else
    objectivef  = @(x) x' * Qo * x + qo * x;
end
N = size(Qo,1);

%% x0
%x0 = randi([0 1],N,1);
%x0 = [ ones(8,1); zeros(N-8,1) ];
x0 = [ ones(1,1); zeros(N-1,1) ];
%[ p filesep 'starting_values.m']
%run( [ p filesep 'starting_values.m'] )
%rot_lines = perturb_rects( x0, Dim, [ 15, -25 ], [1,2,3,7; 4,5,6,8] ); %_2rects
b = -ones( size(A,1), 1 );

ub = ones(N,1) * 1;
lb = ones(N,1) * 0;
IntCon = [ 1 : N ];

fprintf('starting opt\n');
for it = 1 : 2
    %nlconf = @(x) my_nlcon2( x, model, bindings );
    nlconf = @(x) binary_nlcon( x );

    
    %% solve
    options = optimoptions('fmincon','MaxFunEvals',10*size(x0,1), 'TolFun', 1e-10, 'TolX', 1e-20 );
    %options = optimoptions();
    if USE_GA
        if ( size(x0,2) > 1 ) 
            x0 = x0'
        end
        if ( it == 1 )
            stall = 10;
        else
            stall = 500;
        end
        ga_opts = gaoptimset('StallGenLimit',stall,'TolFun',1e-10,'Generations',3000, ...
            'PlotFcns',@gaplotbestf,'InitialPopulation',x0');
        ga_opts = gaoptimset(ga_opts,'PlotFcns',@gaplotbestf,'Display','iter');
        ga_opts = gaoptimset( ga_opts, 'MutationFcn',@my_mutationadaptfeasible );
        fprintf('N: %d, size(x0^t): %d %d\n', N, size(x0') );
        if ( it == 1 )
            [x,fval,exitflag] = ga( objectivef,N,A,b,[],[],lb,ub,[],IntCon,ga_opts)
        else
            [x,fval,exitflag] = ga( objectivef,N,[],[],[],[],lb,ub,[],IntCon,ga_opts)
        end
        
        x0 = x;
        objectivef  = @(x) x * Qo * x' + qo * x';
    elseif it < 2
        x = fmincon( objectivef,  ... %fun
        x0,          ... %x0
        A,          ... %A
        b,          ... %b
        [], ...%Aeq,...%Aeq,         ... %Aeq
        [], ...%beq,         ... %beq
        lb,...%lb,          ... %lb
        ub,...%ub,          ... %ub
        nlconf,...%nlconf,      ... %optimoptions
        options      ...
        );
    else
        x = fmincon( objectivef,  ... %fun
        x0,          ... %x0
        [],          ... %A
        [],          ... %b
        Aeq, ...%Aeq,...%Aeq,         ... %Aeq
        beq, ...%beq,         ... %beq
        [],...%lb,          ... %lb
        [],...%ub,          ... %ub
        nlconf,      ... %optimoptions
        options      ...
        );
    end
    
    
    x0 = x;
    x
end
fprintf('sum(x) = %f\n', sum(x) );

%% save
% save lines

%write_files( fit_lines, [p filesep 'fit_lines.txt'] );
% save Q
%save( [p filesep 'x_solved.mat'], 'x' );
%Q = model.Q;
%Qf = full(Q);
%save( [p filesep 'Q.mat'], 'Q' );
         
end

function mutationChildren = my_mutationadaptfeasible(parents,options,GenomeLength, ...
    FitnessFcn,state,thisScore,thisPopulation,varargin)
    mutatonChildren = mutationadaptfeasible( parents, options, GenomeLength, fittnessFcn, state, thisScore, thisPopulation, varargin );
end

function [ c, ceq ] = binary_nlcon( x )
    c = 0;
    ceq= x .* (x-1);
end

% if normal substitution is used
function [ c, ceq ] = my_nlcon2( x, model, bindings )
    c = zeros( size(x,1),1);
    
%     ceq = zeros( numel(model.quadcon), 1 );
%     for i = 1 : size( ceq, 1 )
%         ceq(i) = x' * model.quadcon(i).Qc * x - model.quadcon(i).rhs;
%     end
    
end

function l = normal2line( n, d )
    dir = cross( [0,0,1], [n,1] );
    
    x0 = n * -d;
%     if ( x0(1) < 0 )
%         x0(1) = x0(1) * -1;
%     end
%      if ( x0(2) < 0 )
%          x0(2) = x0(2) * -1;
%      end
    
    % n . p + d = 0 ==> py = (d - nx * 0) / ny
    l = [ x0, 0, dir ];
    %fprintf('on_line == 0: %f\n', dot([x0,1], [n,d]) );
end
             