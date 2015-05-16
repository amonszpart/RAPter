p = '/home/bontius/ransacTest/build/stairs_noise_0.001';
qo = readSparseMatrix( [ p filesep 'problem/qo.csv'], 1 );
Qo = readSparseMatrix( [ p filesep 'problem/Qo.csv'], 1 );
A  = readSparseMatrix( [ p filesep 'problem/A.csv'], 1 );
size(unique(A,'rows'))

subplot(121)
hist(qo,100)
subplot(122)
hist(Qo(:)/100,100)



