/*
   Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.

   File:      qcqo1.c

   Purpose:   To demonstrate how to solve a quadratic
              optimization problem using the MOSEK API.

              minimize  x_1^2 + 0.1 x_2^2 +  x_3^2 - x_1 x_3 - x_2
              s.t 1 <=  x_1 + x_2 + x_3 - x_1^2 - x_2^2 - 0.1 x_3^2 + 0.2 x_1 x_3
              x >= 0

 */

#include <stdio.h>
#include <iostream> // cout, cerr, endl
#include <vector>

//#include "globfit2/optimization/qp/mosekOpt.h"
#include "qcqpcpp/mosekOptProblem.h"

using qcqpcpp::MosekOpt;

inline MosekOpt<double>::BOUND toSGBound( MSKboundkeye bound )
{
    return MosekOpt<double>::BOUND( bound );
}

int main2( int argc, char** argv )
{
    typedef double      Scalar;
    typedef MSKrescodee ReturnType;

    const int numcon = 1; /* Number of constraints. */
    const int numvar = 3; /* Number of variables. */
    //const int numanz = 3; /* Number of non-zeros in A. */
    const int numqnz = 4; /* Number of non-zeros in Q. */

    MSKrescodee r;

    double q[] = {0.0,-1.0,0.0};

    MSKboundkeye bkc[] = {MSK_BK_LO};
    double blc[] = {1.0};
    double buc[] = {+MSK_INFINITY};

    MSKboundkeye bkx[] = {MSK_BK_LO,
                          MSK_BK_LO,
                          MSK_BK_LO};
    double blx[] = {0.0,
                    0.0,
                    0.0};
    double bux[] = {+MSK_INFINITY,
                    +MSK_INFINITY,
                    +MSK_INFINITY};

    MSKint32t aptrb[] = { 0 },
              aptre[] = { 3 },
              asub [] = { 0, 1, 2 };

    // x >= 0
    double aval[] = { 1.0, 1.0, 1.0 };
    MSKint32t qsubi[numqnz],
              qsubj[numqnz];
    double    qval[numqnz];

    qsubi[0] = 0; qsubj[0] = 0; qval[0] = 2.0;
    qsubi[1] = 1; qsubj[1] = 1; qval[1] = 0.2;
    qsubi[2] = 2; qsubj[2] = 0; qval[2] = -1.0;
    qsubi[3] = 2; qsubj[3] = 2; qval[3] = 2.0;

    MosekOpt<Scalar>::SparseMatrix QObj( numvar, numvar );
    QObj.insert(0,0) = 2.0;
    QObj.insert(1,1) = 0.2;
    QObj.insert(2,0) = -1.0;
    QObj.insert(2,2) = 2.0;

    MosekOpt<Scalar>::SparseMatrix QCon( numvar, numvar );
    QCon.insert(0,0) = -2.0;
    QCon.insert(1,1) = -2.0;
    QCon.insert(2,2) = -0.2;
    QCon.insert(2,0) =  0.2;

#if 1
    // test mosek wrapper
    std::vector<MosekOpt<Scalar>::Scalar> x_out, gt_out;
    {
        MosekOpt<Scalar> mosek;

        // set vars
        for ( int j = 0; j != numvar; ++j )
        {
            // add var
            mosek.addVariable( toSGBound(bkx[j]), blx[j], bux[j] );
            // obj function set
            mosek.setLinObjective( j, q[j] );
        }

        // set quadratic objective
        mosek.addQObjectives( QObj );

        // set constraints
        for ( int i = 0; i != numcon; ++i )
        {
            std::cout << "numcon: " << numcon << ", aptre[i] - aptrb[i]: " << aptre[i] - aptrb[i] << std::endl;
            // add var
            std::vector<double> coeffs( numvar );
            for ( int cid = 0; cid != aptre[i] - aptrb[i]; ++cid )
            {
                std::cout << "setting coeffs[" << asub[aptrb[i] + cid] << "] = aval[" <<  aptrb[i] + cid << "];" << std::endl;
                std::cout << "\t~==   coeffs[ asub[" << aptrb[i] << " + " << cid << "] ] = aval[ aptrb[" << i<< "] + " << cid << "];" << std::endl;
                coeffs[ asub[aptrb[i] + cid] ] = aval[ aptrb[i] + cid ];
            }
            std::cout<<"coeffs:";for(size_t vi=0;vi!=coeffs.size();++vi)std::cout<<coeffs[vi]<<" ";std::cout << "\n";
            mosek.addLinConstraint( toSGBound(bkc[i]), blc[i], buc[i], coeffs );
        }

        // set QConstraints
        {
            for ( int k = 0; k < QCon.outerSize(); ++k )
            {
                for ( MosekOpt<Scalar>::SparseMatrix::InnerIterator it(QCon, k); it; ++it )
                {
                    mosek.addQConstraint( 0, it.row(), it.col(), it.value() );
                    std::cout << "adding constraint " << it.row() << "," << it.col() << ", " << it.value() << std::endl;
                }
            }
        }

        int r = mosek.update();
        if ( r == MSK_RES_OK )
            std::cout << "mosek.finalize returns: " << r << ", MSK_RES_OK: " << (int)MSK_RES_OK << std::endl;
        else
            std::cerr << "mosek.finalize returns: " << r << ", MSK_RES_OK: " << (int)MSK_RES_OK << std::endl;

        mosek.printProblem();
        mosek.optimize( &x_out );
        std::cout<<"x_out:";for(size_t vi=0;vi!=x_out.size();++vi)std::cout<<x_out[vi]<<" ";std::cout << "\n";
        std::cout << std::endl << std::endl;
    }
#endif

    MSKint32t j,i;
    double xx[numvar];
    MSKenv_t env;
    MSKtask_t task;

    /* Create the mosek environment. */
    r = MSK_makeenv(&env,NULL);

    if ( r==MSK_RES_OK )
    {
        /* Create the optimization task. */
        r = MSK_maketask(env,numcon,numvar,&task);

        if ( r==MSK_RES_OK )
        {
            r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,qcqpcpp::mosekPrintStr);

            /* Append 'NUMCON' empty constraints.
   The constraints will initially have no bounds. */
            if ( r == MSK_RES_OK )
            {
                std::cout << "gt: MSK_appendcons(task,"<<numcon<<");" << std::endl;
                r = MSK_appendcons(task,numcon);
            }

            /* Append 'NUMVAR' variables.
   The variables will initially be fixed at zero (x=0). */
            if ( r == MSK_RES_OK )
            {
                std::cout << "gt: MSK_appendvars(task,"<<numvar<<");" << std::endl;
                r = MSK_appendvars(task,numvar);
            }

            /* Optionally add a constant term to the objective. */
            if ( r ==MSK_RES_OK )
            {
                std::cout << "gt: MSK_putcfix(task," << 0.0 << ");" << std::endl;
                r = MSK_putcfix(task,0.0);
            }
            for(j=0; j<numvar && r == MSK_RES_OK; ++j)
            {
                /* Set the linear term c_j in the objective.*/
                if(r == MSK_RES_OK)
                {
                    std::cout << "gt: putcj(task," << j << "," << q[j] << ")" << std::endl;
                    r = MSK_putcj(task,j,q[j]);
                }

                /* Set the bounds on variable j.
   blx[j] <= x_j <= bux[j] */
                if(r == MSK_RES_OK)
                {
                    r = MSK_putvarbound(task,
                                        j, /* Index of variable.*/
                                        bkx[j], /* Bound key.*/
                                        blx[j], /* Numerical value of lower bound.*/
                                        bux[j]); /* Numerical value of upper bound.*/
                    std::cout << "gt: MSK_putvarbound(task," << j << "," << bkx[j] << "," << blx[j] << "," << bux[j] << ");" << std::endl;
                }

                /* Input column j of A */
//                if(r == MSK_RES_OK)
//                {
//                    r = MSK_putacol(task,
//                                    j, /* Variable (column) index.*/
//                                    aptre[j]-aptrb[j], /* Number of non-zeros in column j.*/
//                                    asub+aptrb[j], /* Pointer to row indexes of column j.*/
//                                    aval+aptrb[j]); /* Pointer to Values of column j.*/
//                }
            }

            /* Set the bounds on constraints.
   for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
            for(i=0; i<numcon && r==MSK_RES_OK; ++i)
            {
//                r = MSK_putarow(task,
//                                j, /* Variable (column) index.*/
//                                aptre[j]-aptrb[j], /* Number of non-zeros in column j.*/
//                                asub+aptrb[j], /* Pointer to row indexes of column j.*/
//                                aval+aptrb[j]); /* Pointer to Values of column j.*/

                std::cout << "gt: MSK_putarow( task, " << i
                          << ", " << aptre[i] - aptrb[i]
                             << ", " << *(asub + aptrb[i])
                             << ", " << *(aval + aptrb[i])
                             << ");" << std::endl;
                r = MSK_putarow(task,
                                i,                 /* Row index.*/
                                aptre[i]-aptrb[i], /* Number of non-zeros in row i.*/
                                asub+aptrb[i],     /* Pointer to column indexes of row i.*/
                                aval+aptrb[i]);    /* Pointer to values of row i.*/

                std::cout << "my: MSK_putconbound( task, " << i
                          << ", " << bkc[i]
                          << ", " << blc[i]
                          << ", " << buc[i]
                          << ")" << std::endl;
                r = MSK_putconbound(task,
                                    i, /* Index of constraint.*/
                                    bkc[i], /* Bound key.*/
                                    blc[i], /* Numerical value of lower bound.*/
                                    buc[i]); /* Numerical value of upper bound.*/
            }

            if ( r==MSK_RES_OK )
            {
                /*
   * The lower triangular part of the Q^o
   * matrix in the objective is specified.
   */

//                qsubi[0] = 0; qsubj[0] = 0; qval[0] = 2.0;
//                qsubi[1] = 1; qsubj[1] = 1; qval[1] = 0.2;
//                qsubi[2] = 2; qsubj[2] = 0; qval[2] = -1.0;
//                qsubi[3] = 2; qsubj[3] = 2; qval[3] = 2.0;

                /* Input the Q^o for the objective. */
                std::cout<<"gt: MSK_putqobj(task," << numqnz << ",\n";
                for(size_t vi=0;vi!=numqnz;++vi)
                {
                    std::cout << qsubi[vi] << "," << qsubj[vi] << ", " << qval[vi] << std::endl;
                }
                std::cout << ");\n";
                r = MSK_putqobj(task,numqnz,qsubi,qsubj,qval);

            }

            if ( r==MSK_RES_OK )
            {
                /*
   * The lower triangular part of the Q^0
   * matrix in the first constraint is specified.
   This corresponds to adding the term
   - x_1^2 - x_2^2 - 0.1 x_3^2 + 0.2 x_1 x_3
   */

                qsubi[0] = 0; qsubj[0] = 0; qval[0] = -2.0;
                qsubi[1] = 1; qsubj[1] = 1; qval[1] = -2.0;
                qsubi[2] = 2; qsubj[2] = 2; qval[2] = -0.2;
                qsubi[3] = 2; qsubj[3] = 0; qval[3] = 0.2;

                /* Put Q^0 in constraint with index 0. */

                std::cout<<"gt: MSK_putqonk( task, " << 0
                        << ", " << 4
                        << ",\n";
                for(size_t vi=0;vi!=4;++vi)
                {
                    std::cout << qsubi[vi] << "," << qsubj[vi] << ", " << qval[vi] << std::endl;
                }
                std::cout << "); " << std::endl;
                r = MSK_putqconk(task,
                                 0,
                                 4,
                                 qsubi,
                                 qsubj,
                                 qval);
            }

            if ( r==MSK_RES_OK )
                r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);

            if ( r==MSK_RES_OK )
            {
                MSKrescodee trmcode;

                /* Run optimizer */
                r = MSK_optimizetrm(task,&trmcode);

                /* Print a summary containing information
   about the solution for debugging purposes*/
                MSK_solutionsummary (task,MSK_STREAM_LOG);

                if ( r==MSK_RES_OK )
                {
                    MSKsolstae solsta;
                    int j;

                    MSK_getsolsta(task,MSK_SOL_ITR,&solsta);

                    switch(solsta)
                    {
                        case MSK_SOL_STA_OPTIMAL:
                        case MSK_SOL_STA_NEAR_OPTIMAL:
                            MSK_getxx(task,
                                      MSK_SOL_ITR, /* Request the interior solution. */
                                      xx);

                            printf("Optimal primal solution\n");
                            gt_out.clear(); gt_out.reserve(numvar);
                            for(j=0; j<numvar; ++j)
                            {
                                gt_out.push_back( xx[j] );
                                printf("x[%d]: %e\n",j,xx[j]);
                            }

                            break;
                        case MSK_SOL_STA_DUAL_INFEAS_CER:
                        case MSK_SOL_STA_PRIM_INFEAS_CER:
                        case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
                        case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
                            printf("Primal or dual infeasibility certificate found.\n");
                            break;

                        case MSK_SOL_STA_UNKNOWN:
                            printf("The status of the solution could not be determined.\n");
                            break;
                        default:
                            printf("Other solution status.");
                            break;
                    }
                }
                else
                {
                    printf("Error while optimizing.\n");
                }
            }

            if (r != MSK_RES_OK)
            {
                /* In case of an error print error code and description. */
                char symname[MSK_MAX_STR_LEN];
                char desc[MSK_MAX_STR_LEN];

                printf("An error occurred while optimizing.\n");
                MSK_getcodedesc (r,
                                 symname,
                                 desc);
                printf("Error %s - '%s'\n",symname,desc);
            }
        }

        MSK_deletetask(&task);
    }
    MSK_deleteenv(&env);

    std::cout<<"x_out:\n";
    for ( size_t vi = 0; vi != x_out.size(); ++vi )
    {
        std::cout << "\tx(" << vi << ")=\t" << x_out[vi] << ((x_out[vi] == gt_out[vi]) ? " =OK= " : " =!NOT!= ") << gt_out[vi] << "(gt)" << std::endl;
    }
    std::cout << "\n";
    return ( r );
}

int main(int argc,char *argv[])
{
    return main2( argc, argv );

    const int numvar = 4,
            numcon = 3;

    double       c[]     = {3.0, 1.0, 5.0, 1.0};
    /* Below is the sparse representation of the A
         matrix stored by row. */
    MSKlidxt     aptrb[] = {0, 3, 7};
    MSKlidxt     aptre[] = {3, 7, 9};
    MSKidxt      asub[]  = { 0,1,2,
                             0,1,2,3,
                             1,3};
    double       aval[]  = { 3.0, 1.0, 2.0,
                             2.0, 1.0, 3.0, 1.0,
                             2.0, 3.0};

    /* Bounds on constraints. */
    MSKboundkeye bkc[]  = {MSK_BK_FX, MSK_BK_LO,     MSK_BK_UP    };
    double       blc[]  = {30.0,      15.0,          -MSK_INFINITY};
    double       buc[]  = {30.0,      +MSK_INFINITY, 25.0         };
    /* Bounds on variables. */
    MSKboundkeye bkx[]  = {MSK_BK_LO,     MSK_BK_RA, MSK_BK_LO,     MSK_BK_LO     };
    double       blx[]  = {0.0,           0.0,       0.0,           0.0           };
    double       bux[]  = {+MSK_INFINITY, 10.0,      +MSK_INFINITY, +MSK_INFINITY };
    MSKenv_t     env  = NULL;
    MSKtask_t    task = NULL;
    MSKrescodee  r;
    MSKidxt      i,j;

    // test mosek wrapper
    {
        MosekOpt<double> mosek;

        // set vars
        for ( int j = 0; j != numvar; ++j )
        {
            // add var
            mosek.addVariable( toSGBound(bkx[j]), blx[j], bux[j] );
            // obj function set
            mosek.setLinObjective( j, c[j] );
        }

        // set constraints
        for ( int i = 0; i != numcon; ++i )
        {
            // add var
            std::vector<double> coeffs( numvar );
            for ( int cid = 0; cid != aptre[i] - aptrb[i]; ++cid )
            {
                coeffs[ asub[aptrb[i] + cid] ] = aval[ aptrb[i] + cid ];
            }
            std::cout<<"coeffs:";for(size_t vi=0;vi!=coeffs.size();++vi)std::cout<<coeffs[vi]<<" ";std::cout << "\n";
            mosek.addLinConstraint( toSGBound(bkc[i]), blc[i], buc[i], coeffs );
        }

        mosek.update();
        mosek.printProblem();
        return mosek.optimize();
        std::cout << std::endl << std::endl;
    }

    /* Create the mosek environment. */
    r = MSK_makeenv(&env,NULL);

    if ( r==MSK_RES_OK )
    {
        /* Create the optimization task. */
        r = MSK_maketask(env,numcon,numvar,&task);

        /* Directs the log task stream to the 'printstr' function. */
        if ( r==MSK_RES_OK )
            r = MSK_linkfunctotaskstream( task,MSK_STREAM_LOG,NULL, qcqpcpp::mosekPrintStr );

        /* Append 'numcon' empty constraints.
         The constraints will initially have no bounds. */
        if ( r == MSK_RES_OK )
            r = MSK_appendcons(task,numcon);

        /* Append 'numvar' variables.
         The variables will initially be fixed at zero (x=0). */
        if ( r == MSK_RES_OK )
            r = MSK_appendvars(task,numvar);

        for(j=0; j<numvar && r == MSK_RES_OK; ++j)
        {
            /* Set the linear term c_j in the objective.*/
            if(r == MSK_RES_OK)
                r = MSK_putcj(task,j,c[j]);

            /* Set the bounds on variable j.
           blx[j] <= x_j <= bux[j] */
            if(r == MSK_RES_OK)
                r = MSK_putvarbound(task,
                                    j,           /* Index of variable.*/
                                    bkx[j],      /* Bound key.*/
                                    blx[j],      /* Numerical value of lower bound.*/
                                    bux[j]);     /* Numerical value of upper bound.*/
        }

        /* Set the bounds on constraints.
           for i=1, ...,numcon : blc[i] <= constraint i <= buc[i] */
        for(i=0; i<numcon && r==MSK_RES_OK; ++i)
        {
            r = MSK_putconbound(task,
                                i,           /* Index of constraint.*/
                                bkc[i],      /* Bound key.*/
                                blc[i],      /* Numerical value of lower bound.*/
                                buc[i]);     /* Numerical value of upper bound.*/

            /* Input row i of A */
            if(r == MSK_RES_OK)
                r = MSK_putarow(task,
                                i,                 /* Row index.*/
                                aptre[i]-aptrb[i], /* Number of non-zeros in row i.*/
                                asub+aptrb[i],     /* Pointer to column indexes of row i.*/
                                aval+aptrb[i]);    /* Pointer to values of row i.*/
        }

        /* Maximize objective function. */
        if (r == MSK_RES_OK)
            r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);

        if ( r==MSK_RES_OK )
        {
            MSKrescodee trmcode;

            /* Run optimizer */
            r = MSK_optimizetrm(task,&trmcode);

            /* Print a summary containing information
             about the solution for debugging purposes. */
            MSK_solutionsummary (task,MSK_STREAM_LOG);

            if ( r==MSK_RES_OK )
            {
                MSKsolstae solsta;

                if (r == MSK_RES_OK)
                    r = MSK_getsolsta (task,MSK_SOL_BAS,&solsta);
                switch(solsta)
                {
                    case MSK_SOL_STA_OPTIMAL:
                    case MSK_SOL_STA_NEAR_OPTIMAL:
                    {
                        double *xx = (double*) calloc(numvar,sizeof(double));
                        if ( xx )
                        {
                            MSK_getxx(task,
                                      MSK_SOL_BAS,    /* Request the basic solution. */
                                      xx);

                            printf("Optimal primal solution\n");
                            for(j=0; j<numvar; ++j)
                                printf("x[%d]: %e\n",j,xx[j]);
                        }
                        else
                        {
                            r = MSK_RES_ERR_SPACE;
                        }

                        free(xx);
                        break;
                    }
                    case MSK_SOL_STA_DUAL_INFEAS_CER:
                    case MSK_SOL_STA_PRIM_INFEAS_CER:
                    case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
                    case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
                        printf("Primal or dual infeasibility certificate found.\n");
                        break;
                    case MSK_SOL_STA_UNKNOWN:
                    {
                        char symname[MSK_MAX_STR_LEN];
                        char desc[MSK_MAX_STR_LEN];

                        /* If the solutions status is unknown, print the termination code
                   indicating why the optimizer terminated prematurely. */

                        MSK_getcodedesc(trmcode,
                                        symname,
                                        desc);

                        printf("The solutuion status is unknown.\n");
                        printf("The optimizer terminitated with code: %s\n",symname);
                        break;
                    }
                    default:
                        printf("Other solution status.\n");
                        break;
                }
            }
        }

        if (r != MSK_RES_OK)
        {
            /* In case of an error print error code and description. */
            char symname[MSK_MAX_STR_LEN];
            char desc[MSK_MAX_STR_LEN];

            printf("An error occurred while optimizing.\n");
            MSK_getcodedesc (r,
                             symname,
                             desc);
            printf("Error %s - '%s'\n",symname,desc);
        }

        /* Delete the task and the associated data. */
        MSK_deletetask(&task);
    }

    /* Delete the environment and the associated data. */
    MSK_deleteenv(&env);

    return r;
} // ... main()


