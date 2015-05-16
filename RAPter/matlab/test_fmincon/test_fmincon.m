function test_fmincon( p )

Dim          = 2;
stride       = Dim + 1;
USE_EQ_CNST  = 1; % 1, if Aeq * x = beq should be true
USE_DOT_CNST = 0; % 1, if Aeq * x = beq should be true
cosanglimit  = 0.01;

%x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
if (~exist('p','var') )
    %p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/_1rect_input2_20140704_0852'
    %p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/_2rects'
    %p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/_rot_rects_noise'
    %p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/_rot_rects_noise_0.015_20140711_1054';
    %p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/simple/';
    p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/2rects_noise_0.003/';
end

%% objective function
model       = gurobi_read( [ p filesep 'gf_model.lp' ] );
objectivef  = @(x) x' * model.Q * x;

%% x0
clear x0;
[ p filesep 'starting_values.m']
run( [ p filesep 'starting_values.m'] )
%rot_lines = perturb_rects( x0, Dim, [ 15, -25 ], [1,2,3,7; 4,5,6,8] ); %_2rects
rot_lines = perturb_rects( x0, Dim, [ 5 ], [1,2,3,4] );
rot_lines6d = [];
for l = 1 : size(x0,1) / stride
        n = rot_lines( (l-1)*stride+1 : (l-1)*stride + stride - 1 )';
        n = n / norm(n);
        
        d = rot_lines( (l-1)*stride + stride );
        rot_lines6d(l,:) = normal2line( n, d )
    end
write_files( rot_lines6d, [p filesep 'rot_lines.txt'] );
%x0 = rot_lines;

% match normals

perp_cnt = 0;
if USE_EQ_CNST
    bindings = zeros(0,3);
    Aeq = zeros( 0, size(x0,1) );
    beq = zeros( 0, 1 );
    for i = 1 : stride : size(x0,1)
        ni = [ x0(i), x0(i+1) ];
        for j = i+stride : stride : size(x0,1)
            nj = [ x0(j), x0(j+1) ];
            cosang = dot(ni,nj);
            if ( abs(cosang) >= 1-cosanglimit)
                % debug:
                %fprintf('n%d,n%d: %f\n', i, j, cosang );
                
                % log
                bindings( end+1, : ) = [i,j,1];
                
                % nx:
                Aeq( end+1, i   ) =  1; %Aeq( i, i ) + 1;
%                 if cosang > 0 
                    Aeq( end  , j   ) = -1;
%                 else
%                     Aeq( end  , j   ) = 1;
%                 end
                beq( end+1 )      = 0;
                
                % ny:
                Aeq( end+1, i+1 ) =  1;
%                 if cosang > 0 
                    Aeq( end  , j+1 ) = -1;
%                 else
%                     Aeq( end  , j+1 ) = 1;
%                 end
                beq( end+1 )      = 0;
            elseif ( abs(cosang) <= cosanglimit ) && perp_cnt < 4
                bindings( end+1, : ) = [i,j,0];
                perp_cnt = perp_cnt + 1;
            else
                fprintf('skipping %d-%d, since cosang %f\n', (i-1)/stride, (j-1)/stride, cosang );
            end
            
        end
    end
    
    for b = 1 : size(bindings,1)
        fprintf('bound %d-%d (%d)\n', (bindings(b,1)-1)/stride, (bindings(b,2)-1)/stride, bindings(b,3) );
    end
    
    Aeq
    beq
else
    Aeq = [];
    beq = [];
end


ub = ones(size(x0,1),1) * Inf;
lb = ones(size(x0,1),1) * -Inf;
for l = 1 : stride : size(x0,1)
    for d = 1 : Dim
        ub(l + d - 1) = 1.1;
        lb(l + d - 1) = -1.1;
    end
end
bindings

for it = 1 : 2
    nlconf = @(x) my_nlcon2( x, model, bindings );

    
    %% solve
    options = optimoptions('fmincon','MaxFunEvals',1000*size(x0,1), 'TolFun', 1e-20, 'TolX', 1e-20 );
    %options = optimoptions();
    if it < 2
        x = fmincon( objectivef,  ... %fun
            x0,          ... %x0
            [],          ... %A
            [],          ... %b
            [], ...%Aeq,...%Aeq,         ... %Aeq
            [], ...%beq,         ... %beq
            [],...%lb,          ... %lb
            [],...%ub,          ... %ub
            nlconf,      ... %optimoptions
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
    
    [c, ceq] = my_nlcon2( x, model, bindings );
    fprintf('constriants: %f\n', ceq);
    obj = x' * model.Q * x;
    fprintf('objective function: %f\n', obj );
    
    nLines    = numel( model.quadcon );
    fit_lines = zeros( nLines, 6 );
    l_width   = Dim + 1;
    for l = 1 : nLines
        n = x( (l-1)*l_width+1 : (l-1)*l_width + l_width - 1 )';
        n = n / norm(n);
        
        d = x( (l-1)*l_width + l_width );
        fit_lines(l,:) = normal2line( n, d );
        %fprintf('n . dir: %f\n', dot(fit_lines(l,4:5),[n]) );
    end
    
    [c, ceq] = my_nlcon2( x, model, bindings );
    x0 = x;
    x'
    %Aeq
    %beq
    % it
    %cosanglimit = 0.02;
    
end

%% save
% save lines

write_files( fit_lines, [p filesep 'fit_lines.txt'] );
% save Q
save( [p filesep 'x_solved.mat'], 'x' );
Q = model.Q;
Qf = full(Q);
save( [p filesep 'Q.mat'], 'Q' );
         
end

% if normal substitution is used
function [ c, ceq ] = my_nlcon2( x, model, bindings )
    c = [];
    
    ceq = zeros( numel(model.quadcon), 1 );
    for i = 1 : size( ceq, 1 )
        ceq(i) = x' * model.quadcon(i).Qc * x - model.quadcon(i).rhs;
    end
    
    for i = 1 : 1 %size( bindings, 1 )
        if ( bindings(i,3) == 0 )
            n1 = bindings(i,1);
            n2 = bindings(i,2);
            ceq(end+1) = abs(x(n1) * x(n2) + x(n1+1) * x(n2+1)) - bindings(i,3);
%         elseif ( bindings(i,3) == 1 )
%             n1 = bindings(i,1);
%             n2 = bindings(i,2);
%             %x(n1) * x(n2) + x(n1+1) * x(n2+1)
%             ceq(end+1) = abs(x(n1) * x(n2) + x(n1+1) * x(n2+1)) - bindings(i,3);
        end
     end
end

function l = normal2line( n, d )
    dir = cross( [0,0,1], [n,1] );
    %dir = dir / norm(dir);
%     if ( abs(n(1)) > abs(n(2)) )
%         x0 = [ -d / n(1), 0 ];
%     else
%         x0 = [ 0, -d / n(2) ];
%     end
%     x0 = [ -d / n(1) / 2, -d / n(2) / 2 ];
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
             