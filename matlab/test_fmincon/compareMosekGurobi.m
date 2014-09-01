%% this compares gurobi problem with mosek problem
p = '/home/bontius/workspace/SmartGeometry/ransacTest/build/2lines_rot_noise_0.003/';

% read mosek
Qo = readSparseMatrix( [ p filesep 'Qo.cpp.txt'] );
qo = readSparseMatrix( [ p filesep 'qo.cpp.txt'] );
A  = readSparseMatrix( [ p filesep 'A.cpp.txt'] );
% read gurobi
m = gurobi_read( [ p filesep 'cloud_20140809_1454_0/model.lp'] );

figure();
subplot(131);
imagesc( abs(full(m.obj) - full(qo)) )
title('qo');

subplot(132);
imagesc( abs(full(m.Q)' - full(Qo)) )
title('Qo');

subplot(133);
imagesc( abs(full(m.A) - full(A)) )
title('A');

