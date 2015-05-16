function gurobiToSparseTxt( path, model_fname)
    if ( ~exist( 'path', 'var' ) )
        path = '/home/bontius/workspace/SmartGeometry/ransacTest/build/_1rect_input2_20140704_0852';
    end
    if ( ~exist( 'model_fname', 'var' ) )
        model_fname = 'model.lp';
    end
    m1 = gurobi_read( [ path  filesep model_fname ] );
    writeQ( full(m1.Q*2)', [path filesep 'Qo.txt'] );
    writeQ( m1.obj, [path filesep 'qo.txt'] );
    writeQ( [ full(m1.A), m1.rhs] , [path filesep 'Q.txt'] );
end