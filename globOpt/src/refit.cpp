#include <iostream>

#include "globfit2/globOpt_types.h"           // _2d, _3d namespaces
#include "globfit2/util/parse.h"              // GF2::console
#include "globfit2/io/io.h"
#include "globfit2/processing/angle_util.hpp" // AnglesT
#include "globfit2/io/inputParser.hpp"        // parseInput()
#include "globopt/util/impl/randUtil.hpp"     // randf()

#include "qcqpcpp/bonminOptProblem.h"
#include "globopt/primitives/impl/planePrimitive.hpp" // PlanePrimitive( pos ,normal )
#include "globfit2/processing/graph.hpp"

namespace globopt
{
    // ../refit -p primitives_it20.bonmin.csv -a points_primitives_it19.csv --scale 0.0175 --simple --3D
    template < class _PrimitiveMapT
             , class _PointContainerT
             >
    int refitNormals( _PrimitiveMapT             & outPrims
                   , _PrimitiveMapT        const& primitives
                   , _PointContainerT      const& points
                   , GF2::GidPidVectorMap  const& populations
                   , PidT                  const  targetPop
                   , bool                  const  verbose      = false )
    {

        typedef typename _PrimitiveMapT::PrimitiveT                         PrimitiveT;
    //        typedef typename _PrimitiveMapT::InnerContainerT                    InnerContainerT;
    //        typedef typename std::vector<InnerContainerT>                       PrimitiveVectorT;
        typedef typename PrimitiveT::Scalar                                 Scalar;
    //        typedef typename _PointContainerT::PrimitiveT                       PointPrimitiveT;
    //        typedef typename _PclCloudT::Ptr                                    PclCloudPtrT;
        typedef typename PrimitiveT::Position                               Position;
        typedef typename PrimitiveT::Direction                              Direction;
        typedef Eigen::Matrix<Scalar,4,1>                                   Vector4;

    //        typedef          Eigen::Map< Position >                             PositionMap;


        for ( typename _PrimitiveMapT::ConstIterator it(primitives); it.hasNext(); it.step() )
        {
            const GidT          gId     = it.getGid();
            const PrimitiveT&   prim    = *it;

            Position avgPos ( Position::Zero() );
            Vector4  normal ; normal  << prim.template normal(), Scalar(0.);
            Vector4  primPos; primPos << prim.template pos   (), Scalar(0.);
            for ( auto it = populations.at(gId).begin(); it != populations.at(gId).end(); ++it )
            {
                auto &pos = points[*it].template pos();

                Vector4 pt; pt << pos(0), pos(1), pos(2), Scalar(0.);
                // double k = (DOT_PROD_3D (points[i], p21) - dotA_B) / dotB_B;
                Scalar k =   (pt.dot(normal) - primPos.dot(normal))
                          /
                            normal.dot( normal );

                Vector4 pp = primPos + k * normal;
                // Calculate the projection of the point on the line (pointProj = A + k * B)
                avgPos += pp.template head<3>();
            } //...for inliers
            avgPos /= Scalar(populations.at(gId).size());

            PrimitiveT outPrim( avgPos, prim.template dir() );
            outPrim.copyTagsFrom( prim );
            GF2::containers::add( outPrims, gId, outPrim );
            std::cout << "moved [" << gId << "]: " << (outPrim.template pos() - prim.template pos()).norm() << std::endl;
        } //...for primitives

        return EXIT_SUCCESS;
    } //...refitSimple()


    // ../refit -p primitives_it20.bonmin.csv -a points_primitives_it19.csv --scale 0.0175 --simple --3D
    template < class _PrimitiveMapT
             , class _PointContainerT
             >
    int refitSimple( _PrimitiveMapT             & outPrims
                   , _PrimitiveMapT        const& primitives
                   , _PointContainerT      const& points
                   , GF2::GidPidVectorMap  const& populations
                   , PidT                  const  targetPop
                   , bool                  const  verbose      = false )
    {

        typedef typename _PrimitiveMapT::PrimitiveT                         PrimitiveT;
    //        typedef typename _PrimitiveMapT::InnerContainerT                    InnerContainerT;
    //        typedef typename std::vector<InnerContainerT>                       PrimitiveVectorT;
        typedef typename PrimitiveT::Scalar                                 Scalar;
    //        typedef typename _PointContainerT::PrimitiveT                       PointPrimitiveT;
    //        typedef typename _PclCloudT::Ptr                                    PclCloudPtrT;
        typedef typename PrimitiveT::Position                               Position;
        typedef typename PrimitiveT::Direction                              Direction;
        typedef Eigen::Matrix<Scalar,4,1>                                   Vector4;

    //        typedef          Eigen::Map< Position >                             PositionMap;

        for ( typename _PrimitiveMapT::ConstIterator it(primitives); it.hasNext(); it.step() )
        {
            const GidT          gId     = it.getGid();
            const PrimitiveT&   prim    = *it;

            Position avgPos ( Position::Zero() );
            Vector4  normal ; normal  << prim.template normal(), Scalar(0.);
            Vector4  primPos; primPos << prim.template pos   (), Scalar(0.);
            for ( auto it = populations.at(gId).begin(); it != populations.at(gId).end(); ++it )
            {
                auto &pos = points[*it].template pos();

                Vector4 pt; pt << pos(0), pos(1), pos(2), Scalar(0.);
                // double k = (DOT_PROD_3D (points[i], p21) - dotA_B) / dotB_B;
                Scalar k =   (pt.dot(normal) - primPos.dot(normal))
                          /
                            normal.dot( normal );

                Vector4 pp = primPos + k * normal;
                // Calculate the projection of the point on the line (pointProj = A + k * B)
                avgPos += pp.template head<3>();
            } //...for inliers
            avgPos /= Scalar(populations.at(gId).size());

            PrimitiveT outPrim( avgPos, prim.template dir() );
            outPrim.copyTagsFrom( prim );
            GF2::containers::add( outPrims, gId, outPrim );
            std::cout << "moved [" << gId << "]: " << (outPrim.template pos() - prim.template pos()).norm() << std::endl;
        } //...for primitives

        return EXIT_SUCCESS;
    } //...refitSimple()

    template < class _PrimitiveMapT
             , class _PointContainerT
             >
    int refitNonLin( _PrimitiveMapT             & outPrims
                   , _PrimitiveMapT        const& primitives
                   , _PointContainerT      const& points
                   , GF2::GidPidVectorMap  const& populations
                   , PidT                  const  targetPop
                   , bool                  const  verbose      = false )
    {

        typedef typename _PrimitiveMapT::PrimitiveT                         PrimitiveT;
//        typedef typename _PrimitiveMapT::InnerContainerT                    InnerContainerT;
//        typedef typename std::vector<InnerContainerT>                       PrimitiveVectorT;
        typedef typename PrimitiveT::Scalar                                 Scalar;
//        typedef typename _PointContainerT::PrimitiveT                       PointPrimitiveT;
//        typedef typename _PclCloudT::Ptr                                    PclCloudPtrT;
        typedef typename PrimitiveT::Position                               Position;
//        typedef          Eigen::Map< Position >                             PositionMap;

        typedef double                        OptScalar;
        typedef qcqpcpp::BonminOpt<OptScalar> OptProblemT;
        OptProblemT problem;

        const int Dims = 2; // 2: nx,ny, 3: +nz
        typedef          std::pair<LidT,LidT>                        LIdPair;
        typedef          std::map<LIdPair, std::vector<LidT> >       PrimsVarsT;

        std::vector<std::pair<std::string,Scalar> >     starting_values;    // [var_id] = < var_name, x0 >
        PrimsVarsT                                      prims_vars;         // < <lid,lid1>, var_id >, associates primitives to variables

        char                                            name[255];          // variable name
        const Position                                  origin( Position::Zero() ); // to calculate d

        // add variables
        {
            for ( typename _PrimitiveMapT::ConstIterator it(primitives); it.hasNext(); it.step() )
            {
                const LidT          lId0    = it.getLid0();
                const LidT          lId1    = it.getLid1();
                const GidT          gId     = it.getGid();
                const DidT          dId     = it.getDid();
                const PrimitiveT&   prim    = *it;
                const LIdPair       varKey  = LIdPair( lId0, lId1 );


                Position normal = prim.template normal();
                std::cout << "line_" << gId << "_" << dId << ".n = " << normal.transpose() << std::endl;

                sprintf( name, "nx_%ld_%ld", gId, dId );
                prims_vars[ varKey ].push_back( problem.addVariable(OptProblemT::BOUND::RANGE, -problem.getINF(), problem.getINF(), OptProblemT::VAR_TYPE::CONTINUOUS, OptProblemT::LINEAR, name) );
                starting_values.push_back( std::pair<std::string,Scalar>(name,normal(0)) );

                sprintf( name, "ny_%ld_%ld", gId, dId );
                prims_vars[ varKey ].push_back( problem.addVariable(OptProblemT::BOUND::RANGE, -problem.getINF(), problem.getINF(), OptProblemT::VAR_TYPE::CONTINUOUS, OptProblemT::LINEAR, name) );
                starting_values.push_back( std::pair<std::string,Scalar>(name,normal(1)) );

                if ( Dims > 2 )
                {
                    sprintf( name, "nz_%ld_%ld", gId, dId );
                    prims_vars[ varKey ].push_back( problem.addVariable(OptProblemT::BOUND::RANGE, -problem.getINF(), problem.getINF(), OptProblemT::VAR_TYPE::CONTINUOUS, OptProblemT::LINEAR, name) );
                    starting_values.push_back( std::pair<std::string,Scalar>(name,normal(2)) );
                }

                sprintf( name, "d_%ld_%ld", gId, dId );
                Scalar d = Scalar(-1.) * prim.getDistance( origin );
                std::cout << "line_" << gId << "_" << dId << ".d = " << d << std::endl;
                prims_vars[ varKey ].push_back( problem.addVariable(OptProblemT::BOUND::RANGE, -problem.getINF(), problem.getINF(), OptProblemT::VAR_TYPE::CONTINUOUS, OptProblemT::LINEAR, name) );
                starting_values.push_back( std::pair<std::string,Scalar>(name,d) );
            } //...for primitives
        } //...add variables

        // add constraints: |normal|^2 = 1
        {
            for ( typename _PrimitiveMapT::ConstIterator it(primitives); it.hasNext(); it.step() )
            {
                const LidT          lId0    = it.getLid0();
                const LidT          lId1    = it.getLid1();
        //                const GidT          gId     = it.getGid();
        //                const DidT          dId     = it.getDid();
        //                const PrimitiveT&   prim    = *it;
                const LIdPair       varKey  = LIdPair( lId0, lId1 );

                // sparse quadratic matrix
                OptProblemT::SparseMatrix norm_constraint( problem.getVarCount(), problem.getVarCount() );

                // 1 * nx * nx + 1 * ny * ny + 1 *  nz * nz
                for ( int dim = 0; dim != Dims; ++dim )
                    norm_constraint.insert( prims_vars[varKey][dim]
                                          , prims_vars[varKey][dim] ) = Scalar( 1. );

                // add constraint instance
                problem.addConstraint  ( OptProblemT::BOUND::EQUAL, /* >= 1 */ Scalar(1.), /* <= 1 */ Scalar(1.), /* linear constraint coeffs: */ NULL );
                // add quadratic coefficients
                problem.addQConstraints( norm_constraint );
            } //...for primitives
        } //...constraints

        /// cost -> objective: minimize \sum_n \sum_p ((n_j . p_i) + d)^2  where point p_i is assigned to line with normal n_j
        //if ( Dims == 3 ) throw new std::runtime_error("datafit unimplemented for 3D");
        // n0^2 p0^2 + n1^2 p1^2 + n2^2 p2^2 +
        // 2 n0 p0 d + 2 n1 p1 d + 2 n2 p2 d +
        // 2 n0 n1 p0 p1 + d^2 +
        // 2 n0 n2 p0 p2 + 2 n1 n2 p1 p2
        {
            // for each line
            Scalar coeff = Scalar( 0. );
            for ( typename _PrimitiveMapT::ConstIterator it(primitives); it.hasNext(); it.step() )
            {
                const LidT          lId0    = it.getLid0();
                const LidT          lId1    = it.getLid1();
                const GidT          gId     = it.getGid();
        //                const DidT          dId     = it.getDid();
        //                const PrimitiveT&   prim    = *it;
                const LIdPair       varKey  = LIdPair( lId0, lId1 );

                // for each assigned point
                Scalar rat = std::min( Scalar(1.), targetPop / Scalar(populations.at(gId).size()) );
                for ( size_t pIdId = 0; pIdId != populations.at(gId).size(); ++pIdId )
                {
                    if ( randf<Scalar>() > rat ) continue;
                    const PidT pid = populations.at(gId)[pIdId];

                    // debug
                    if ( verbose ) std::cout << "[" << __func__ << "]: " << "adding pid " << pid << " -> " << "lines[" << lId0 << "][" << lId1 << "]" << std::endl;

                    // for each dimension: x,y,z
                    for ( int dim = 0; dim != Dims; ++dim )
                    {
                        // (p_x)^2 . (n_x)^2
                        coeff = points[pid].pos()( dim );                                                   // (p_x)
                        coeff *= coeff;                                                                     // (p_x)^2
                        problem.addQObjective( prims_vars[varKey][dim], prims_vars[varKey][dim], coeff );   // (p_x)^2 . (n_x)^2

                        // debug
                        if ( verbose ) std::cout << "added qterm(pid: " << pid << "): "
                                                    << coeff << " * "
                                                    << problem.getVarName( prims_vars[varKey][dim] ) << " * "
                                                    << problem.getVarName( prims_vars[varKey][dim] )
                                                    << std::endl;

                        // 2 . p_x . n_x . d
                        coeff = Scalar(2.) * points[pid].pos()( dim );              // 2 . p_x
                        problem.addQObjective( /* n_x: */ prims_vars[varKey][dim]
                                             , /*   d: */ prims_vars[varKey][Dims]
                                             ,            coeff );                  // 2 . p_x . n_x . d

                        // debbug
                        if ( verbose ) std::cout << "added qterm(pid: " << pid << "): "
                                                 << coeff << " * "
                                                 << problem.getVarName( prims_vars[varKey][dim] ) << " * "
                                                 << problem.getVarName( prims_vars[varKey][Dims] )
                                                 << std::endl;
                    }

                    // d^2
                    coeff = Scalar(1);
                    problem.addQObjective( /* d: */ prims_vars[varKey][Dims], /* d: */ prims_vars[varKey][Dims], coeff );

                    // debug
                    if ( verbose ) std::cout << "added qterm(pid: " << pid << "): "
                                             << coeff << " * "
                                             << problem.getVarName( prims_vars[varKey][Dims] ) << " * "
                                             << problem.getVarName( prims_vars[varKey][Dims] )
                                             << std::endl;

                    // 2 . px . py . nx . ny
                    coeff = Scalar(2.) * points[pid].pos()(0) * points[pid].pos()(1);
                    problem.addQObjective( prims_vars[varKey][0], prims_vars[varKey][1], coeff );

                    // debug
                    if ( verbose ) std::cout << "added qterm(pid: " << pid << "): "
                                                << coeff << " * "
                                                << problem.getVarName( prims_vars[varKey][0] ) << " * "
                                                << problem.getVarName( prims_vars[varKey][1] )
                                                << std::endl;

                    // 2 n0 n2 p0 p2 + 2 n1 n2 p1 p2
                    if ( Dims > 2 )
                    {
                        // 2 n0 n2 p0 p2
                        coeff = Scalar(2.) * points[pid].pos()(0) * points[pid].pos()(2);
                        problem.addQObjective( prims_vars[varKey][0], prims_vars[varKey][2], coeff );
                        // 2 n1 n2 p1 p2
                        coeff = Scalar(2.) * points[pid].pos()(1) * points[pid].pos()(2);
                        problem.addQObjective( prims_vars[varKey][1], prims_vars[varKey][2], coeff );
                    }
                } // for points
            } //...for primitives
        } //...add objective

        // add perpendicular constraints
        struct Triplet { LIdPair prim0, prim1; Scalar rhs;
                         Triplet( LIdPair const& prim0, LIdPair const& prim1, Scalar rhs )
                            : prim0(prim0), prim1(prim1), rhs(rhs) {}
                       };

        typedef GF2::graph::EdgeT<Scalar>                                       EdgeT;
        typedef GF2::graph::EdgeListT<Scalar>                                   EdgeListT;
        typedef typename GF2::Graph<Scalar, typename GF2::MyGraphConfig<Scalar>::UndirectedGraph> GraphT;

        EdgeListT edgeList;


        const Scalar angTolRad = 0.01 / 180. * M_PI;
        std::vector< Triplet > constraints;
        int perpCount = 0, paralCount = 0;
        for ( typename _PrimitiveMapT::ConstIterator it0(primitives); it0.hasNext(); it0.step() )
            for ( typename _PrimitiveMapT::ConstIterator it1(primitives); it1.hasNext(); it1.step() )
            {
                // skip not same groupID ones
                if ( it0.getDid() != it1.getDid() ) continue;
                if ( it0.getGid() > it1.getGid() ) continue; // increasing order
                if ( it0.getGid() == it1.getGid() &&
                     it0.getLid1() == it1.getLid1() ) continue;

                Scalar angle = std::abs( GF2::angleInRad(it0->template dir(), it1->template dir()) );
                while ( angle > M_PI ) angle -= M_PI;
                if ( !perpCount && std::abs( M_PI_2 - angle ) < angTolRad ) // 90degs
                {
                    constraints.push_back( Triplet( LIdPair(it0.getLid0(),it0.getLid1()),
                                                    LIdPair(it1.getLid0(),it1.getLid1()),
                                                    0.) );
                    //++perpCount;
                    edgeList.insert( EdgeT(it0.getGid(), it1.getGid(), 20.)  );
                    std::cout << "perpconstr: <"
                              << constraints.back().prim0.first << ","
                              << constraints.back().prim0.second << ","
                              << prims_vars.at( constraints.back().prim0 ).at(0) << ","
                              << prims_vars.at( constraints.back().prim0 ).at(1)
                              << "> . "
                              << constraints.back().prim1.first << ","
                              << constraints.back().prim1.second << ","
                              << prims_vars.at( constraints.back().prim1 ).at(0) << ","
                              << prims_vars.at( constraints.back().prim1 ).at(1)
                              << "> ==  "
                              << constraints.back().rhs
                              << " (angle: " << angle * 180.0 / M_PI << ")"
                              << std::endl;
                }
                else if ( paralCount < 1 &&
                          ( (angle < angTolRad) || (std::abs(M_PI-angle) < angTolRad) )
                        ) // parallel
                {
                    constraints.push_back( Triplet( LIdPair(it0.getLid0(),it0.getLid1()),
                                                    LIdPair(it1.getLid0(),it1.getLid1()),
                                                    1.) );
                    //++paralCount;
                    edgeList.insert( EdgeT(it0.getGid(), it1.getGid(), 1.)  );
                    std::cout << "paralconstr: <"
                              << constraints.back().prim0.first << ","
                              << constraints.back().prim0.second << ","
                              << prims_vars.at( constraints.back().prim0 ).at(0) << ","
                              << prims_vars.at( constraints.back().prim0 ).at(1)
                              << "> . "
                              << constraints.back().prim1.first << ","
                              << constraints.back().prim1.second << ","
                              << prims_vars.at( constraints.back().prim1 ).at(0) << ","
                              << prims_vars.at( constraints.back().prim1 ).at(1)
                              << "> ==  "
                              << constraints.back().rhs
                              << " (angle: " << angle * 180.0 / M_PI << ")"
                              << std::endl;
                }
            }

        GraphT g(edgeList);
        typedef std::pair<GidT,GidT> GidPair;
        std::set< GidPair > mstEdges;
        g.spanningTree( mstEdges );

        for ( int i = 0; i != constraints.size(); ++i ) // perps.size()
        {
            OptProblemT::SparseMatrix perp_constraint( problem.getVarCount(), problem.getVarCount() );
            const std::vector<LidT>& prim0VarIds = prims_vars.at( constraints[i].prim0 );
            const std::vector<LidT>& prim1VarIds = prims_vars.at( constraints[i].prim1 );
            auto getGid = [&primitives]( LidT const lid0, LidT const lid1 ) {
                auto it = primitives.begin();
                std::advance( it, lid0 );
                return it->second.at( lid1 );
            };
            const GidT gid0 = getGid( constraints[i].prim0.first, constraints[i].prim0.second ).getTag( PrimitiveT::TAGS::GID );
            const GidT gid1 = getGid( constraints[i].prim1.first, constraints[i].prim1.second ).getTag( PrimitiveT::TAGS::GID );
            if ( mstEdges.find( GidPair(std::min(gid0,gid1), std::max(gid0,gid1) ) ) == mstEdges.end() )
                continue;
            std::cout << "constraint: " << gid0 << " - " << gid1 << " = " << constraints[i].rhs << std::endl;

            // 1 * nx0 * nx1 + 1 * ny0 * ny1 = 1/0
            for ( int dim = 0; dim != Dims; ++dim )
                perp_constraint.insert( prim1VarIds.at(dim), prim0VarIds.at(dim) ) = Scalar( 1. ); // reverse order for lower triangular

            // add constraint instance
            problem.addConstraint  ( OptProblemT::BOUND::EQUAL, constraints[i].rhs, constraints[i].rhs, /* linear constraint coeffs: */ NULL );
            // add quadratic coefficients
            problem.addQConstraints( perp_constraint );
        }

        // starting point
        {
            if ( starting_values.size() != problem.getVarCount() )
            {
                std::cerr << "[" << __func__ << "]: " << "starting_values.size() " << starting_values.size() << " != " << problem.getVarCount() << " problem.getVarCount()" << std::endl;
                return EXIT_FAILURE;
            }

            OptProblemT::SparseMatrix x0( problem.getVarCount(), 1 );
            for ( LidT i = 0; i != starting_values.size(); ++i )
            {
                x0.insert( i, 0 ) = starting_values[i].second;
            }
            problem.setStartingPoint( x0 );
        } //...starting values

        // save
        {
            problem.write( "./datafit_problem" );
        }

        int err = EXIT_SUCCESS;
        // solve
        {
            OptProblemT::ReturnType r = problem.getOkCode();

            if ( problem.getOkCode() == r )
            {
                if ( verbose ) { std::cout << "[" << __func__ << "]: " << "calling problem optimize...\n"; fflush(stdout); }
                r = problem.update();
                if ( verbose ) { std::cout << "[" << __func__ << "]: " << "finished problem optimize...\n"; fflush(stdout); }
            }

            // optimize
            std::vector<OptScalar> x_out;
            if ( r == problem.getOkCode() )
            {
                // log
                if ( verbose ) { std::cout << "[" << __func__ << "]: " << "calling problem optimize...\n"; fflush(stdout); }

                // work
                r = problem.optimize( &x_out, OptProblemT::OBJ_SENSE::MINIMIZE );

                // check result
                if ( r != problem.getOkCode() )
                {
                    std::cerr << "[" << __func__ << "]: " << "ooo...optimize didn't work with code " << r << std::endl; fflush(stderr);
                    err = r;
                }

                // output result
                //for ( LidT i = 0; i < x_out.size(); i += 4 )
                LidT i = 0;
                for ( typename _PrimitiveMapT::ConstIterator it(primitives); it.hasNext() && i < x_out.size(); it.step(), i += (Dims+1) )
                {
                    PrimitiveT outPrim;
                    Position normal( Position::Zero() );
                    for ( int d = 0; d != Dims; ++d )
                        normal(d) = x_out[ i+d ];
                    PrimitiveT::generateFrom( outPrim, (normal).eval(), Scalar(x_out[i+Dims]) );
                    outPrim.copyTagsFrom( *it );
                    GF2::containers::add( outPrims, it.getGid(), outPrim );
                }
            } //...optimize
        } //...solve

        return err;
    } //...refitNonLin()

    //! \brief Unfinished function. Supposed to do GlobFit.
    template < class _PrimitiveMapT
             , class _PointContainerT
             , class _PclCloudT
             >
    int
    refit( int    argc
         , char** argv )
    {
        typedef typename _PrimitiveMapT::PrimitiveT                         PrimitiveT;
        typedef typename _PrimitiveMapT::InnerContainerT                    InnerContainerT;
        typedef typename std::vector<InnerContainerT>                       PrimitiveVectorT;
        typedef typename PrimitiveT::Scalar                                 Scalar;
        typedef typename _PointContainerT::PrimitiveT                       PointPrimitiveT;
        typedef typename _PclCloudT::Ptr                                    PclCloudPtrT;
        typedef          GF2::PidT                                          PidT;
        typedef          GF2::LidT                                          LidT;
        typedef          GF2::GidT                                          GidT;
        typedef          GF2::DidT                                          DidT;
        typedef typename PrimitiveT::Position                               Position;
        typedef          Eigen::Map< Position >                             PositionMap;

        int                     err             = EXIT_SUCCESS;
        Scalar                  scale           = 0.05f;
        std::string             cloud_path      = "cloud.ply",
                                primitives_path,
                                associations_path = "";
        GF2::AnglesT            angle_gens( {GF2::AnglesT::Scalar(90.)} );
        bool                    verbose         = false;
        PidT                    targetPop       = 100;

        std::cout << "hello refit" << std::endl;

        _PointContainerT        points;
        PclCloudPtrT            pclCloud;
        PrimitiveVectorT        primitivesVector;
        _PrimitiveMapT          primitives;
        bool                    valid_input = true, simple = false;
        struct Params { Scalar scale; } params;

        // parse
        std::string primsPath(""), assocPath("");
        {
            valid_input = (EXIT_SUCCESS ==
                           GF2::parseInput<InnerContainerT,_PclCloudT>(points, pclCloud, primitivesVector, primitives, params, argc, argv) );
            if ( !valid_input )
            {
                std::cerr << "[" << __func__ << "]: " << "could not parse input..." << std::endl;
                return EXIT_FAILURE;
            }

            GF2::console::parse_argument( argc, argv, "--target-pop", targetPop );
            std::cout << "[" << __func__ << "]: " << "counting " << targetPop << " points from each primitive, change with --target-pop\n";

            simple = GF2::console::find_switch( argc, argv, "--simple" );
            std::cout << "[" << __func__ << "]: " << "performint simple fitting: " << ( simple ? "YES" : "NO" ) << std::endl;

            primsPath = GF2::parsePrimitivesPath( argc, argv );
            assocPath = GF2::parseAssocPath     ( argc, argv );
        }

        // assignments
        GF2::GidPidVectorMap populations;
        GF2::processing::getPopulations( populations, points );

        // WORK
        _PrimitiveMapT outPrims;
        if ( !simple )
            refitNonLin( outPrims, primitives, points, populations, targetPop, verbose );
        else
            refitSimple( outPrims, primitives, points, populations, targetPop, verbose );

        // save
        std::string outName = primsPath.substr( 0, primsPath.rfind(".csv") );
        std::stringstream ssPrims; ssPrims << outName << ".refit.csv";
        GF2::io::savePrimitives<PrimitiveT,typename InnerContainerT::const_iterator>( outPrims, ssPrims.str() );
        std::cout << "results written to " << ssPrims.str() << "\n";
        std::cout << "../show.py -s " << scale << " -p " << ssPrims.str() << " -a " << assocPath << " --save-poly" << std::endl;

        // visualize
        pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
        vptr->setBackgroundColor( .1, .1, .1 );
        vptr->addPointCloud<GF2::PclPointT>( pclCloud, "cloud" );


        char name[255];
        for ( typename _PrimitiveMapT::ConstIterator it(outPrims); it.hasNext(); it.step() )
        {
            sprintf( name, "line_%ld", it.getUniqueId() );
            PrimitiveT::template draw<PointPrimitiveT>( /*   primitive: */ *it
                                                      , /*      points: */ points
                                                      , /*   threshold: */ scale
                                                      , /*     indices: */ &populations.at( it.getGid() )
                                                      , /*      viewer: */ vptr
                                                      , /*   unique_id: */ name
                                                      , /*      colour: */ 1., 0., 0.
                                                      , /* viewport_id: */ 0
                                                      );
        }

        vptr->spin();
        return err;
    } // ...refit()

} //...ns globopt

int main( int argc, char** argv )
{
//    typedef float                                                           Scalar;
//    typedef GF2::graph::EdgeT<Scalar>                                       EdgeT;
//    typedef GF2::graph::EdgeListT<Scalar>                                   EdgeListT;
//    typedef GF2::Graph<Scalar, GF2::MyGraphConfig<Scalar>::UndirectedGraph> GraphT;

//    EdgeListT edgeList;
//    edgeList.insert( EdgeT(0,1,1) );
//    edgeList.insert( EdgeT(0,2,1) );
//    edgeList.insert( EdgeT(1,2,1) );
//    edgeList.insert( EdgeT(1,3,1) );
//    edgeList.insert( EdgeT(3,4,1) );
//    edgeList.insert( EdgeT(4,0,1) );
//    GraphT g(edgeList);
//    g.spanningTree();
//    return EXIT_SUCCESS;

    if ( GF2::console::find_switch(argc,argv,"--3D") )
    {
        return globopt::refit< GF2::_3d::PrimitiveMapT
                             , GF2::PointContainerT
                             , GF2::PclCloudT
                             >( argc, argv );
    }
    else
    {
        return globopt::refit< GF2::_2d::PrimitiveMapT
                             , GF2::PointContainerT
                             , GF2::PclCloudT
                             >( argc, argv );
    } //...if find_switch
} //...main

