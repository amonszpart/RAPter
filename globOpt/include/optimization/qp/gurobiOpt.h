#include "gurobi_c++.h"

#include <random>
#include <limits>
#include <memory>  // shared_ptr
#include <fstream> // ofstream
#include <iomanip> // setfill, setw

#include "Eigen/Dense"

//#include "AMUtil2.h"
//#include "AMUtilPCL.h"

#include "globfit2/primitives/linePrimitive2.h"
#include "globfit2/visualization/visualizer.h"
#include "globfit2/io/io.h"
#include "globfit2/optimization/energyFunctors.h"
//#include "params.h"
#include "optimization/qp/MyGRBCallback.h"


//////////////////////////
/// GUROBIOPT
//////////////////////////

namespace GF2
{
    ///
    /// \sum_ij X_ij = 1, where j is any other dir_j that was copied to patch_i to create L_ij, selected, if X_ij == 1
    ///
    template <class PrimitiveContainerT, class PointContainerT>
    struct SimpleModelConstraintsStruct
    {
            static inline int
            eval ( GRBModel                                    & model
                   , std::vector<std::vector<GRBVar> >    const& vars
                   , PrimitiveContainerT                  const& lines
                   , PointContainerT                      const& /*points*/ )
            {
                // Constraints
                char cname[16];
                for ( size_t lid = 0; lid != lines.size(); ++lid )
                {
                    GRBLinExpr sum_constr = 0.0;
                    for ( size_t lid1 = 0; lid1 != lines[lid].size(); ++lid1 )
                    {
                        sum_constr += vars[lid][lid1];
                    }

                    sprintf( cname, "c_%d", static_cast<int>(lid) );
                    //model.addConstr( vars[lid][lid1] >= 1, cname );
                    model.addConstr( sum_constr, GRB_EQUAL, 1, cname );
                }

                return EXIT_SUCCESS;
            }
    };

    ///
    /// \sum_ij X_ij = 1, where i == patch_id, j \in line_ids that can explain patch_i inside scale
    ///
    template <class PrimitiveContainerT, class PointContainerT, class PointLineDistanceFunctor, typename Scalar>
    struct AdvancedModelConstraintsFunctor
    {
            static inline int
            eval ( GRBModel                                    & model
                   , std::vector<std::vector<GRBVar> >    const& vars
                   , PrimitiveContainerT                  const& lines
                   , PointContainerT                      const& points
                   , Scalar                               const  scale )
            {
                char cname[16];
                std::vector<GRBLinExpr> points_constraints( points.size() );

                for ( size_t lid = 0; lid != lines.size(); ++lid )
                    for ( size_t lid1 = 0; lid1 != lines[lid].size(); ++lid1 )
                    {
                        for ( size_t pid = 0; pid != points.size(); ++pid )
                        {
                            Scalar dist = PointLineDistanceFunctor::template eval<Scalar>( points[pid], lines[lid][lid1] );
                            if ( dist < scale )
                            {
                                points_constraints[pid] += vars[lid][lid1];
                            }
                        }
                    }
                for ( size_t pid = 0; pid != points_constraints.size(); ++pid )
                {
                    sprintf( cname, "c_%d", static_cast<int>(pid) );
                    model.addConstr( points_constraints[pid], GRB_GREATER_EQUAL, 1, cname );
                }

                return EXIT_SUCCESS;
            }
    };

    struct GurobiOptParams
    {
            float unary_weight      = 1.f;
            float pairwise_weight   = 1e2f;
            float complexity_weight = 1e2f;
            float time_chunk        = 20.f; // hot-start interval to dump intermediate results
            float time_limit        = 20.f; // overall run-time
            int   thread_count      = 0;
            std::string store_path  = "./";
            char  continuous        = GRB_CONTINUOUS;
    };

    template <typename Scalar, class PrimitiveT, class PointT>
    class GurobiOpt
    {
        public:
            typedef Eigen::Matrix                   <Scalar,3,1>                Vector;
            typedef typename std::vector            <std::vector<PrimitiveT> >  PrimitiveContainerT;
            typedef typename std::vector            <PointT                  >  PointContainerT;

            typedef MyPointPrimitiveDistanceFunctor                               PointPrimitiveDistanceFunctorT;
            typedef SqrtPrimitivePrimitiveEnergyFunctor<Scalar, PrimitiveT>         PrimitivePrimitiveDistanceFunctorT;
            typedef AdvancedModelConstraintsFunctor<
                PrimitiveContainerT, PointContainerT
                , PointPrimitiveDistanceFunctorT, Scalar >                      AddModelConstraintsFunctorT;

            static const std::string UNARY_WEIGHT;// ( "unary_weight" );
            static const std::string PAIRWISE_WEIGHT;//   = "pairwise_weight";

            //inline static int run();

            // template <class PrimitiveT, class PointT, class DistanceFunctor, typename Scalar>
            inline static int
            solve( PrimitiveContainerT            & out_lines
                   , PrimitiveContainerT      const& lines
                   , PointContainerT          const& points
                   , Scalar                   const  scale    = 0.05f
                   , std::vector<Scalar>      const& angles   = { Scalar(0), Scalar(M_PI_2), Scalar(M_PI) }
                   , bool                     const  verbose  = false
                   , GurobiOptParams          const* params   = NULL );

            inline static int
            globFit( PrimitiveContainerT            & out_lines
                     , PrimitiveContainerT      const& lines
                     , PointContainerT          const& points
                     , Scalar                   const  scale    = 0.05f
                     , std::vector<Scalar>      const& angles   = { Scalar(0), Scalar(M_PI_2), Scalar(M_PI) }
                     , bool                     const  verbose  = false
                     , GurobiOptParams          const  params   = GurobiOptParams() );


            inline static int
            storeResults( PrimitiveContainerT                        & out_lines
                          , PrimitiveContainerT                 const& lines
                          , PointContainerT                     const& points
                          , std::vector<std::vector<GRBVar> >   const& vars
                          , Scalar                              const  scale
                          , std::vector<Scalar>                 const& angles
                          , GurobiOptParams                     const& params
                          , Scalar                              const  time_stamp
                          , std::shared_ptr<MyGRBCallback>      const& grb_cb
                          , Scalar                              const  time_offset );

            inline static Scalar
            calculateEnergy( PrimitiveContainerT            const& lines_arg
                             , PointContainerT              const& points
                             , std::vector<Scalar>          const& angles
                             , Scalar                       const  scale
                             , Eigen::Matrix<Scalar,3,1>    const weights
                             , bool                         const  verbose = true);

            // unit test - gurobi example from the web
            inline static int test();
    };

} // ns gf2

#include <iostream>
namespace GF2
{
    template <typename Scalar, class PrimitiveContainerT, class PointContainerT> Scalar
    GurobiOpt<Scalar,PrimitiveContainerT,PointContainerT>::calculateEnergy( PrimitiveContainerT    const& lines_arg
                                                                            , PointContainerT      const& points
                                                                            , std::vector<Scalar>  const& angles
                                                                            , Scalar               const  scale
                                                                            , Eigen::Matrix<Scalar,3,1> const weights
                                                                            , bool                 const  verbose )
    {
        PrimitivePrimitiveDistanceFunctorT lineLineDistanceFunctor( angles );
        typedef typename PrimitiveContainerT::value_type::value_type PrimitiveT;

        // linearize
        std::vector<PrimitiveT> lines;
        {
            for ( size_t lid = 0; lid != lines_arg.size(); ++lid )
                for ( size_t lid1 = 0; lid1 != lines_arg[lid].size(); ++lid1 )
                    lines.push_back( lines_arg[lid][lid1] );
        }

        // costs
        Eigen::Matrix<Scalar,3,1> energies( Eigen::Matrix<Scalar,3,1>::Zero() );
        {
            // unary
            for ( size_t lid = 0; lid != lines.size(); ++lid )
            {
                Scalar unary_i = Scalar(0);
                for ( size_t pid = 0; pid != points.size(); ++pid )
                {
                    // if within scale, add unary cost
                    Scalar dist = PointPrimitiveDistanceFunctorT::eval<Scalar>( points[pid], lines[lid] );
                    if ( dist < scale )     unary_i += dist / scale;
                } // for points

                energies(0) += unary_i;
            }

            // pw
            for ( size_t lid = 0; lid != lines.size()-1; ++lid )
                for ( size_t olid = lid+1; olid != lines.size(); ++olid )
                {
                    energies(1) += lineLineDistanceFunctor.eval( lines[lid], lines[olid] );
                }

            // complx
            energies(2) = lines.size();
        }

        // print
        Scalar sum = weights.dot( energies );
        if ( verbose )
        {
            std::cout << "Energy: ";
            for ( int d = 0; d != energies.rows(); ++d )
            {
                std::cout << weights(d) << " * " << energies(d) << "  +  ";
            }

            std::cout << " = ";
            for ( int d = 0; d != energies.rows(); ++d )
                std::cout << weights(d) * energies(d) << "  +  ";

            std::cout << " = " << sum << std::endl;
        }

        return sum;
    } // ... calculateEnergy

    template <typename Scalar, class PrimitiveContainerT, class PointContainerT> int
    GurobiOpt<Scalar,PrimitiveContainerT,PointContainerT>::storeResults( PrimitiveContainerT                        & out_lines
                                                                         , PrimitiveContainerT                 const& lines
                                                                         , PointContainerT                     const& points
                                                                         , std::vector<std::vector<GRBVar> >   const& vars
                                                                         , Scalar                              const  scale
                                                                         , std::vector<Scalar>                 const& angles
                                                                         , GurobiOptParams                     const& params
                                                                         , Scalar                              const  time_stamp
                                                                         , std::shared_ptr<MyGRBCallback>      const& grb_cb
                                                                         , Scalar                              const  time_offset )
    {
        char str_time_stamp[64];
        sprintf( str_time_stamp, "%05d", (int)time_stamp );

        // out_lines
        out_lines.clear();
        {
            std::cout << "[" << __func__ << "]: " << "chosen primitives: ";
            for ( size_t lid = 0; lid != lines.size(); ++lid )
                for ( size_t lid1 = 0; lid1 != lines[lid].size(); ++lid1 )
                {
                    if ( vars[lid][lid1].get(GRB_DoubleAttr_X) > 0.5 )
                    {
                        out_lines.resize( lid+1 );
                        out_lines[lid].push_back( lines[lid][lid1] );

                        std::cout << vars[lid][lid1].get(GRB_StringAttr_VarName) << "(" << vars[lid][lid1].get(GRB_DoubleAttr_X) << "), ";
                    }
                }
            std::cout << std::endl;

            {
                std::stringstream ss;
                ss << params.store_path << "/" << "primitives_" << str_time_stamp << ".txt";
                io::savePrimitives( out_lines, ss.str() );
            }
        }

        if ( grb_cb )
        {
            std::stringstream ss;
            ss << params.store_path << "/" << "scores_" << str_time_stamp << ".txt";
            grb_cb->saveScores( ss.str(), /* gnuplot: */ false, time_offset );
        }

        return EXIT_SUCCESS;
    } // ... storeResults

    // template <class PrimitiveT, class PointT, class DistanceFunctor, typename Scalar>
    template <typename Scalar, class PrimitiveContainerT, class PointContainerT> int
    GurobiOpt<Scalar,PrimitiveContainerT,PointContainerT>::solve( PrimitiveContainerT         & out_lines
                                                                  , PrimitiveContainerT  const& lines
                                                                  , PointContainerT      const& points
                                                                  , Scalar               const  scale
                                                                  , std::vector<Scalar>  const& angles
                                                                  , bool                 const  verbose
                                                                  , GurobiOptParams      const* params  )
    {
        typedef typename PrimitiveContainerT::value_type::value_type PrimitiveT;

        if ( !params )
        {
            std::cerr << "[" << __func__ << "]: " << "params == NULL...exiting\n";
            return EXIT_FAILURE;
        }
        std::cout << "[" << __func__ << "]: "
                  << "\n\tUnary:\t" << params->unary_weight
                  << "\n\tPairwise:\t" << params->pairwise_weight
                  << "\n\tComplexity:\t" << params->complexity_weight
                  << "\n\tTime_chunk:" << params->time_chunk
                  << "\n\tTime_limit:" << params->time_limit
                  << "\n\tThreads: " << params->thread_count
                  << "\n\tstore: " << params->store_path
                  << std::endl;

        Eigen::Matrix<Scalar,3,1> weights; weights << params->unary_weight, params->pairwise_weight, params->complexity_weight;
        std::cout << "[" << __func__ << "]: " << "\tweights: " << weights.transpose() << std::endl;

        using std::cout;
        using std::endl;

        try
        {
            GRBEnv       env     = GRBEnv();
            GRBModel     model   = GRBModel(env);
            GRBQuadExpr  objectiveQExpr;           // objective function
            std::vector<std::vector<GRBVar> > vars;

            //MyGRBCallback cb;
            GRBVar *cb_vars;
            std::shared_ptr<MyGRBCallback>    grb_cb;

            PrimitivePrimitiveDistanceFunctorT lineLineDistanceFunctor( angles );

            // add costs
            {
                // variables
                char name[16];
                vars.resize( lines.size() );

                for ( size_t lid = 0; lid != lines.size(); ++lid )
                {
                    for ( size_t lid1 = 0; lid1 != lines[lid].size(); ++lid1 )
                    {
                        // add var
                        int gid     = lines[lid][lid1].getTag( PrimitiveT::GID );
                        int dir_gid = lines[lid][lid1].getTag( PrimitiveT::DIR_GID );
                        sprintf( name, "x_%d_%d", gid, dir_gid );
                        if ( GRB_CONTINUOUS == params->continuous )
                            vars[lid].emplace_back( model.addVar(0.0, 1.0, gid == dir_gid ? 1.0 : 0.0, GRB_CONTINUOUS, name) );
                        else
                            vars[lid].emplace_back( model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name) );
                    }
                }

                //model.addVars()
                model.update();

                // continuous lower bounds
                if ( GRB_CONTINUOUS == params->continuous )
                {
                    char cname[64];
                    for ( size_t lid = 0; lid != lines.size(); ++lid )
                    {
                        for ( size_t lid1 = 0; lid1 != lines[lid].size(); ++lid1 )
                        {
                            GRBLinExpr lb; lb += vars[lid][lid1];
                            sprintf( cname, "%lu_%lu_lb", lid, lid1 );
                            model.addConstr( lb, GRB_GREATER_EQUAL, 0, cname );
                        }
                    }
                }
                model.update();

                // contraints
                {
#               if 1
                    AddModelConstraintsFunctorT::eval( model, vars, lines, points, scale );
#               else
                    SimpleModelConstraintsStruct<PrimitiveContainerT,PointContainerT> addModelConstraintsFunctor;
                    addModelConstraintsFunctor.eval( model, vars, lines, points );
#               endif
                }

                // unary cost -> objective
                for ( size_t lid = 0; lid != lines.size(); ++lid )
                {
                    for ( size_t lid1 = 0; lid1 != lines[lid].size(); ++lid1 )
                    {
                        // unary cost
                        {
                            Scalar unary_i = Scalar(0);
                            for ( size_t pid = 0; pid != points.size(); ++pid )
                            {
                                // if within scale, add unary cost
                                Scalar dist = PointPrimitiveDistanceFunctorT::eval<Scalar>( points[pid], lines[lid][lid1] );
                                if ( dist < scale )
                                    unary_i += dist / scale;
                            } // for points

                            unary_i *= weights(0);
                            objectiveQExpr.addTerm( unary_i, vars[lid][lid1] );
                            objectiveQExpr.addTerm( weights(2), vars[lid][lid1] );
                            if ( verbose ) std::cout << "[" << __func__ << "]: " << "added " << unary_i << " * " << vars[lid][lid1].get(GRB_StringAttr_VarName) << std::endl;
                        }
                    }
                } // for lines
                model.update();

                // pairwise cost -> objective
                for ( size_t lid = 0; lid != lines.size(); ++lid )
                {
                    for ( size_t lid1 = 0; lid1 != lines[lid].size(); ++lid1 )
                    {
                        for ( size_t olid = 0; olid != lines.size(); ++olid )
                        {
                            for ( size_t olid1 = 0; olid1 != lines[olid].size(); ++olid1 )
                            {
#if 0
                                if ( (lid == olid) )
                                {
                                    if ( lid1 >= olid1 )    continue;
                                }
                                else if ( lid > olid )      continue;
#else
                                if ( (lid == olid) && (lid1 == olid1) ) continue;
#endif

                                Scalar dist = lineLineDistanceFunctor.eval( lines[lid][lid1], lines[olid][olid1] );
                                dist *= weights(1);
                                objectiveQExpr.addTerm( dist, vars[lid][lid1], vars[olid][olid1] );

                                if ( verbose ) std::cout << "[" << __func__ << "]: " << "added " << vars[lid][lid1].get(GRB_StringAttr_VarName)
                                                         << " * " << dist
                                                         << " * " << vars[olid][olid1].get(GRB_StringAttr_VarName) << std::endl;
                            } // ... olid1
                        } // ... olid
                    } // ... lid1
                } // ... lid
                std::cout << "qexpr.size(): " << objectiveQExpr.size() << std::endl;
                model.update();

                // set objective
                model.setObjective( objectiveQExpr, GRB_MINIMIZE );
                model.update();
            }

            // dump
            model.write( params->store_path + "/" + "model.lp" );
            std::cout << "wrote to " << params->store_path + "/" + "model.lp" << std::endl;

            // callback
            std::cout << "[" << __func__ << "]: " << "setting callback" << std::endl;
            {
                int numvars = model.get(GRB_IntAttr_NumVars);

                cb_vars = model.getVars();

                // Create a callback object and associate it with the model
                grb_cb.reset( new MyGRBCallback(cb_vars, numvars) );
                model.setCallback( grb_cb.get() );
            }

            // environment
            if ( params->thread_count > 0 ) model.getEnv().set( GRB_IntParam_Threads     , params->thread_count );
            if ( params->time_chunk   > 0 ) model.getEnv().set( GRB_DoubleParam_TimeLimit, params->time_chunk   );

            model.getEnv().set( GRB_IntParam_Presolve, GRB_PRESOLVE_CONSERVATIVE );

            // optimize
            float   curr_time       = 0.f;
            int     model_status    = GRB_INPROGRESS;
            do {
                std::cout << "[" << __func__ << "]: " << "optimizing " << curr_time << "s .." << curr_time + params->time_chunk << "s" << std::endl; fflush(stdout);
                model.optimize();

                storeResults( out_lines, lines, points, vars, scale, angles, *params, curr_time + params->time_chunk, grb_cb, curr_time );

                calculateEnergy( out_lines, points, angles, scale, weights );
                std::cout << "[" << __func__ << "]: " << "Minimized score: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;

                model_status = model.get( GRB_IntAttr_Status );
                std::cout << "[" << __func__ << "]: " << "model_status: " << model_status << std::endl;
            } while ( ((curr_time += params->time_chunk) < params->time_limit)
                      && (model_status != GRB_OPTIMAL)
                      && (model_status != GRB_SOLUTION_LIMIT) );

            // cleanup
            if ( cb_vars )
                delete [] cb_vars;

        } catch(GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch(...) {
            cout << "Exception during optimization" << endl;
        }

        return 0;
    } // ... GurobiOpt::solve()

    // template <class PrimitiveT, class PointT, class DistanceFunctor, typename Scalar>
    template <typename Scalar, class PrimitiveContainerT, class PointContainerT> int
    GurobiOpt<Scalar,PrimitiveContainerT,PointContainerT>::globFit( PrimitiveContainerT       & out_lines
                                                                  , PrimitiveContainerT  const& lines
                                                                  , PointContainerT      const& points
                                                                  , Scalar               const  scale
                                                                  , std::vector<Scalar>  const& angles
                                                                  , bool                 const  verbose
                                                                  , GurobiOptParams      const  params  )
    {
        const int Dim = 2; // 2: nx,ny, 3: +nz
        typedef typename PrimitiveContainerT::value_type::value_type PrimitiveT;
        typedef typename PointContainerT::value_type                 PointPrimitiveT;
        bool gurobi_log = false;

        try
        {
            GRBEnv       env     = GRBEnv();
            GRBModel     model   = GRBModel(env);
            GRBQuadExpr  objectiveQExpr;           // objective function
            std::vector<std::vector<std::vector<GRBVar> > > vars;
            std::vector<std::pair<std::string,Scalar> >     starting_values;
            std::vector< std::vector<std::pair<int,int> > > points_primitives;

            // add costs
            {
                // variables
                char name[16];
                vars.resize( lines.size() );

                Eigen::Matrix<Scalar,3,1> origin( Eigen::Matrix<Scalar,3,1>::Zero() );
                for ( size_t lid = 0; lid != lines.size(); ++lid )
                {
                    for ( size_t lid1 = 0; lid1 != lines[lid].size(); ++lid1 )
                    {
                        // add var
                        int gid     = lid;
                        int dir_gid = lid1;

                        Eigen::Matrix<Scalar,3,1> normal = lines[lid][lid1].template normal<Scalar>();
                        std::cout << "line_" << gid << "_" << dir_gid << ".n = " << normal.transpose() << std::endl;
                        vars[lid].emplace_back( std::vector<GRBVar>() );

                        sprintf( name, "nx_%d_%d", gid, dir_gid );
                        //vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, normal(0), GRB_CONTINUOUS, name) );
                        vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name) );
                        starting_values.push_back( std::pair<std::string,Scalar>(name,normal(0)) );


                        sprintf( name, "ny_%d_%d", gid, dir_gid );
                        vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name) );
                        starting_values.push_back( std::pair<std::string,Scalar>(name,normal(1)) );
                        //vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, normal(1), GRB_CONTINUOUS, name) );

                        if ( Dim > 2 )
                        {
                            sprintf( name, "nz_%d_%d", gid, dir_gid );
                            //vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, normal(2), GRB_CONTINUOUS, name) );
                            vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name) );
                            starting_values.push_back( std::pair<std::string,Scalar>(name,normal(2)) );
                        }

                        sprintf( name, "d_%d_%d", gid, dir_gid );
                        Scalar d = Scalar(-1) * lines[lid][lid1].getDistance( origin );
                        std::cout << "line_" << gid << "_" << dir_gid << ".d = " << d << std::endl;
                        // vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, lines[lid][lid1].point3Distance(origin), GRB_CONTINUOUS, name) );
                        vars[lid].back().emplace_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name) );
                        starting_values.push_back( std::pair<std::string,Scalar>(name,d) );
                    }
                }
                model.update();

                // |normal|^2 = 1
                {
                    char cname[64];
                    for ( size_t lid = 0; lid != lines.size(); ++lid )
                    {
                        for ( size_t lid1 = 0; lid1 != lines[lid].size(); ++lid1 )
                        {
                            GRBQuadExpr norm_constraint;
                            for ( int d = 0; d != Dim; ++d )
                                norm_constraint.addTerm( Scalar(1), vars[lid][lid1][d], vars[lid][lid1][d] );

                            sprintf( cname, "%lu_%lu_norm_ge", lid, lid1 );
                            model.addQConstr( norm_constraint, GRB_EQUAL, 1, cname );
                        }
                    }
                }
                model.update();

                /// cost -> objective
                // assignments

                if ( points_primitives.size() )
                    std::cerr << "[" << __func__ << "]: " << "warning, points_primitives not empty" << std::endl;

                // parse input associations instead of currently detecting point->primitive assocs
                {
                    int max_group_id = -1;
                    for ( size_t pid = 0; pid != points.size(); ++pid )
                        max_group_id = std::max( max_group_id, points[pid].getTag(PointPrimitiveT::GID) );
                    std::cout << "[" << __func__ << "]: " << "max_group_id: " << max_group_id << std::endl;
                    if ( max_group_id > 0 )
                    {
                        points_primitives.resize( points.size() );
                        for ( size_t pid = 0; pid != points.size(); ++pid )
                        {
                            // assume linear...TODO: 2d
                            int tag = points[pid].getTag(PointPrimitiveT::GID);
                            if ( tag >= 0 )
                                points_primitives[pid].push_back( std::pair<int,int>(0,tag) );
                        }
                    }
                }

                if ( !points_primitives.size() )
                {
                    points_primitives.resize( points.size() );
                    for ( size_t lid = 0; lid != lines.size(); ++lid )
                    {
                        for ( size_t lid1 = 0; lid1 != lines[lid].size(); ++lid1 )
                        {
                            Scalar dist = Scalar(0);
                            for ( size_t pid = 0; pid != points.size(); ++pid )
                            {
                                // if within scale, add point constraint to line
                                dist = PointPrimitiveDistanceFunctorT::eval<Scalar>( points[pid], lines[lid][lid1] );
                                if ( dist < scale )
                                {
                                    points_primitives[pid].push_back( std::pair<int,int>(lid,lid1) );
                                    if ( lid1 == 2 && points_primitives[pid].size() > 1 )
                                    {
                                        std::cout<<"[" << __func__ << "]: " << "oooo, points_primitives["<<pid<<"]:";
                                        for ( size_t vi = 0; vi!= points_primitives[pid].size(); ++vi )
                                            std::cout << "<" << points_primitives[pid][vi].first
                                                      << "," << points_primitives[pid][vi].second
                                                      << ">; ";
                                        std::cout << "\n";

                                    }
                                }
                            } // ... for points
                        } // ... for lines
                    } // ... for lines
                } // ... if points_primitives.empty()

                // add term
                Scalar coeff = Scalar(0);
                for ( size_t pid = 0; pid != points.size(); ++pid )
                {
                    // skip, if ambiguous assignment
                    if ( points_primitives[pid].size() != 1 )
                        continue;

                    const int lid  = points_primitives[pid][0].first;
                    const int lid1 = points_primitives[pid][0].second;

                    // debug
                    std::cout << "[" << __func__ << "]: " << "adding pid " << pid << " -> " << "lines[" << lid << "][" << lid1 << "]" << std::endl;

                    // x,y,z
                    for ( int dim = 0; dim != Dim; ++dim )
                    {
                        // p_x^2 . n_x^2
                        coeff = ((typename PointPrimitiveT::VectorType)points[pid])(dim);
                        coeff *= coeff;
                        objectiveQExpr.addTerm( coeff, vars[lid][lid1][dim], vars[lid][lid1][dim] );
                        if ( gurobi_log ) std::cout << "added qterm(pid: " << pid << "): "
                                                    << coeff << " * "
                                                    << vars[lid][lid1][dim].get(GRB_StringAttr_VarName) << " * "
                                                    << vars[lid][lid1][dim].get(GRB_StringAttr_VarName)
                                                    << std::endl;

                        // 2 . p_x . n_x . d
                        coeff = Scalar(2) * ((typename PointPrimitiveT::VectorType)points[pid])(dim);
                        objectiveQExpr.addTerm( coeff, vars[lid][lid1][dim], vars[lid][lid1].back() ); // back == d
                        if ( gurobi_log ) std::cout << "added qterm(pid: " << pid << "): "
                                                    << coeff << " * "
                                                    << vars[lid][lid1][dim].get(GRB_StringAttr_VarName) << " * "
                                                    << vars[lid][lid1].back().get(GRB_StringAttr_VarName)
                                                    << std::endl;
                    }

                    // d^2
                    coeff = Scalar(1);
                    objectiveQExpr.addTerm( coeff, vars[lid][lid1].back(), vars[lid][lid1].back() );
                    if ( gurobi_log ) std::cout << "added qterm(pid: " << pid << "): "
                                                << coeff << " * "
                                                << vars[lid][lid1].back().get(GRB_StringAttr_VarName) << " * "
                                                << vars[lid][lid1].back().get(GRB_StringAttr_VarName)
                                                << std::endl;

                    // 2 . px . py . nx . ny
                    coeff = Scalar(2) * ((typename PointPrimitiveT::VectorType)points[pid])(0) * ((typename PointPrimitiveT::VectorType)points[pid])(1);
                    objectiveQExpr.addTerm( coeff, vars[lid][lid1][0], vars[lid][lid1][1] );
                    if ( gurobi_log ) std::cout << "added qterm(pid: " << pid << "): "
                                                << coeff << " * "
                                                << vars[lid][lid1][0].get(GRB_StringAttr_VarName) << " * "
                                                << vars[lid][lid1][1].get(GRB_StringAttr_VarName)
                                                << std::endl;

                } // for points
                model.update();

                // set objective
                model.setObjective( objectiveQExpr, GRB_MINIMIZE );
                model.update();
            }

            // dump
            {
                // dump model
                model.write( params.store_path + "/" + "gf_model.lp" );
                std::cout << "wrote to " << params.store_path + "/" + "gf_model.lp" << std::endl;
                // dump x0 (starting values)
                {
                    std::string f_start_path = params.store_path + "/" + "starting_values.m";
                    ofstream f_start( f_start_path );
                    f_start << "x0 = [ ...\n";
                    for ( size_t i = 0; i != starting_values.size(); ++i )
                    {
                        f_start << "    " << starting_values[i].second << ",      % " << starting_values[i].first << "\n";
                    }
                    f_start << "];";
                    f_start.close();
                    std::cout << "[" << __func__ << "]: " << "wrote to " << f_start_path << std::endl;
                }
                // dump associations
                {
                    // reverse map for sorted output
                    std::vector< std::vector< std::vector<int> > > primitives_points( lines.size() );
                    for ( size_t pid = 0; pid != points_primitives.size(); ++pid )
                    {
                        if ( points_primitives[pid].size() == 1 )
                        {
                            const int lid  = points_primitives[pid][0].first;
                            const int lid1 = points_primitives[pid][0].second;
                            if ( primitives_points[lid].size() <= lid1 )
                                primitives_points[ lid ].resize( lid1 + 1 );
                            primitives_points[ lid ][ lid1 ].push_back( pid );
                            std::cout << "added " << "primitives_points[ "<<lid<<" ][ "<<lid1<<" ].back(): " << primitives_points[ lid ][ lid1 ].back() << std::endl;
                        }
                        else
                        {
                            std::cout << "skipping pid " << pid << ", since: ";
                            std::cout<<"points_primitives["<<pid<<"]:";
                            for(size_t vi=0;vi!=points_primitives[pid].size();++vi)std::cout<<points_primitives[pid][vi].first<<","
                                                                                           << ";" << points_primitives[pid][vi].second << ";  ";
                            std::cout << "\n";
                        }

                    }

                    std::string f_assoc_path = params.store_path + "/" + "points_primitives.txt";
                    if ( !boost::filesystem::exists(f_assoc_path) )
                    {
                        ofstream f_assoc( f_assoc_path );
                        f_assoc << "# point_id,primitive_pos_id,primitive_dir_id" << std::endl;
                        for ( size_t lid = 0; lid != primitives_points.size(); ++lid )
                            for ( size_t lid1 = 0; lid1 != primitives_points[lid].size(); ++lid1 )
                                for ( size_t pid_id = 0; pid_id != primitives_points[lid][lid1].size(); ++pid_id )
                                {
                                    f_assoc << primitives_points[lid][lid1][pid_id]
                                               << "," << lid
                                               << "," << lid1 << std::endl;
                                }
                        f_assoc.close();
                        std::cout << "[" << __func__ << "]: " << "wrote to " << f_assoc_path << std::endl;
                    }
                    else
                    {
                        std::cout << "[" << __func__ << "]: " << "did NOT write to " << f_assoc_path << ", since already exists" << std::endl;
                    }
                }
            }

            //model.getEnv().set( GRB_IntParam_Presolve, GRB_PRESOLVE_OFF );
            model.optimize();

        } catch(GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch(...) {
            cout << "Exception during optimization" << endl;
        }

        return 0;
    }

    template <typename Scalar, class PrimitiveT, class PointT> int
    GurobiOpt<Scalar,PrimitiveT,PointT>::test()
    {
        using namespace std;

        try
        {
            GRBEnv env = GRBEnv();

            GRBModel model = GRBModel(env);

            // Create variables

            GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "x");
            GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "y");
            GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "z");

            // Integrate new variables

            model.update();

            // Set objective

            GRBQuadExpr obj = x*x + x*y + y*y + y*z + z*z + 2*x;
            model.setObjective(obj);

            // Add constraint: x + 2 y + 3 z >= 4

            model.addConstr(x + 2 * y + 3 * z >= 4, "c0");

            // Add constraint: x + y >= 1

            model.addConstr(x + y >= 1, "c1");

            // Optimize model

            model.optimize();

            cout << x.get(GRB_StringAttr_VarName) << " "
                 << x.get(GRB_DoubleAttr_X) << endl;
            cout << y.get(GRB_StringAttr_VarName) << " "
                 << y.get(GRB_DoubleAttr_X) << endl;
            cout << z.get(GRB_StringAttr_VarName) << " "
                 << z.get(GRB_DoubleAttr_X) << endl;

            cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

            // Change variable types to integer

            x.set(GRB_CharAttr_VType, GRB_INTEGER);
            y.set(GRB_CharAttr_VType, GRB_INTEGER);
            z.set(GRB_CharAttr_VType, GRB_INTEGER);

            // Optimize model

            model.optimize();

            cout << x.get(GRB_StringAttr_VarName) << " "
                 << x.get(GRB_DoubleAttr_X) << endl;
            cout << y.get(GRB_StringAttr_VarName) << " "
                 << y.get(GRB_DoubleAttr_X) << endl;
            cout << z.get(GRB_StringAttr_VarName) << " "
                 << z.get(GRB_DoubleAttr_X) << endl;

            cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

        } catch(GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch(...) {
            cout << "Exception during optimization" << endl;
        }

        return EXIT_SUCCESS;
    }
} // ns gf2


