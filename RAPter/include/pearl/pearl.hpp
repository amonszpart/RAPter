#ifndef __RAPTER_PEARL_HPP__
#define __RAPTER_PEARL_HPP__

#include <limits>

#include "omp.h"

#include "Eigen/Dense"
#include "gco/GCoptimization.h"

#include "pcl/visualization/pcl_visualizer.h" //DEBUG

#include "pearl/pearl.h"

#include "rapter/typedefs.h"
#include "rapter/simpleTypes.h"
#include "rapter/processing/util.hpp"   // GidPidVectorMap
#include "rapter/util/containers.hpp"   // containers::add
#include "rapter/util/impl/pclUtil.hpp" // smartgeometry::
#include "rapter/processing/util.hpp"   // processing::getNeighbourhoodIndices

namespace am
{
    template <class _PclCloudT, class _PointsContainerT, class _PrimitiveT>
    inline int
    Pearl::run( std::vector<int>                       & labels
              , std::vector<_PrimitiveT>               & lines
              , _PclCloudT                        const& cloud
              , _PointsContainerT                 const& points
              , std::vector<int>                  const* indices
              , Params                            const& params
              , std::vector<std::vector<int> >         * label_history
              , std::vector<std::vector<_PrimitiveT> > * line_history
              , int                               const  nPropose
              )
    {
        int err = EXIT_SUCCESS;

        // prepare
        //lines.clear();
        labels.clear();

        // propose
        if ( !lines.size() )
        {
            err += Pearl::propose( lines, cloud, indices, params, nPropose );
            if ( err != EXIT_SUCCESS ) return err;
            std::cout << "lines.size(): " << lines.size() << std::endl;
        }
        else
        {
            std::cout << "[" << __func__ << "]: " << "received " << lines.size() << "input primitives, skipping propose step!" << std::endl;
        }

        int iteration_id = 0;
        std::vector<int> prev_labels;
        do
        {
            std::cout << "\n[" << __func__ << "]: " << "iteration " << iteration_id << std::endl;
            prev_labels = labels;

            // expand
            err += Pearl::expand( /* [in/out]      labels: */ labels
                                , /* [in]        pclCloud: */ cloud
                                , /* [in] pointPrimitives: */ points
                                , /* [in]         indices: */ NULL
                                , /* [in]      primitives: */ lines
                                , /* [in]        gammasqr: */ params.gammasqr // 50*50
                                , /* [in]            beta: */ params.beta     // params.scale*100
                                , /* [in]      parameters: */ params );
            if ( err != EXIT_SUCCESS ) return err;

            if ( label_history )
                label_history->emplace_back( labels );
            if ( line_history )
                line_history->emplace_back( lines );

            // refit
            err += Pearl::refit( lines, labels, *cloud, params ); // 25, 1000.f, 5e6f
            if ( err != EXIT_SUCCESS ) return err;

            if ( prev_labels.size() )
            {
                int label_diff = 0;
                for ( size_t i = 0; i != labels.size(); ++i )
                    label_diff += (labels[i] != prev_labels[i]);
                std::cout << "it[" << iteration_id << "] labeldiff: " << label_diff << std::endl;
            }

        } while ( /**/ (iteration_id++ != params.max_pearl_iterations)
                  &&   (!prev_labels.size() || !std::equal( labels.begin(), labels.end(), prev_labels.begin() )) );

        return err;
    }

    template <typename PointT, class TLine>
    inline int
    Pearl::propose( std::vector<TLine>                            & lines
                    , boost::shared_ptr<pcl::PointCloud<PointT> >   cloud
                    , std::vector<int>                       const* indices
                    , Params                                 const& params
                    , int const nPropose
                    )
    {
        typedef typename TLine::Scalar Scalar;

        if ( indices ) { std::cerr << __PRETTY_FUNCTION__ << "]: indices must be NULL, not implemented yet..." << std::endl; return EXIT_FAILURE; }

        // get neighbourhoods
        std::vector<std::vector<int   > > neighs;
        std::vector<std::vector<Scalar> > sqr_dists;
        rapter::processing::getNeighbourhoodIndices( neighs
                                                , cloud
                                                , indices
                                                , &sqr_dists
                                                , params.max_neighbourhood_size   // 5
                                                , params.scale // 0.02f
                                                , /* soft_radius: */ true
                                                );

        // every point proposes primitive[s] using its neighbourhood
        float chance = std::min( float(nPropose) / float(neighs.size()), 1.f );
        for ( size_t pid = 0; pid != neighs.size(); ++pid )
        {
            if ( nPropose && (rand() / float(RAND_MAX) > chance) ) continue;
#if 1
            // can't fit a line to 0 or 1 points
            if ( neighs[pid].size() < 2 ) continue;

            Eigen::Matrix<Scalar,TLine::Dim,1> line;
            if ( TLine::EmbedSpaceDim == 2 )
            {
                int err = smartgeometry::geometry::fitLinearPrimitive( /*           output: */ line
                                                                       , /*         points: */ *cloud
                                                                       , /*          scale: */ params.scale
                                                                       , /*        indices: */ &(neighs[pid])
                                                                       , /*    refit times: */ 2
                                                                       , /* use input line: */ false  );
                if ( err == EXIT_SUCCESS )
                    lines.emplace_back( TLine(line) );
            }
            else
            {
                Eigen::Matrix<Scalar,4,1> plane;
                int err = smartgeometry::geometry::fitLinearPrimitive<pcl::PointCloud<PointT>, float, 4>
                          ( /*           output: */ plane
                                                                       , /*         points: */ *cloud
                                                                       , /*          scale: */ params.scale
                                                                       , /*        indices: */ &(neighs[pid])
                                                                       , /*    refit times: */ 2
                                                                       , /* use input line: */ false
                                                                       );
                if ( err == EXIT_SUCCESS )
                    lines.emplace_back( TLine( -plane.template head<3>() * plane(3), plane.template head<3>()) );
            }


#else
            // skip, if no neighbours found this point won't contribute a primitive for now
            if ( neighs[pid].size() < 2 )  continue;

            // compute neighbourhood cov matrix
            Eigen::Matrix<Scalar,3,3> cov;
            smartgeometry::computeCovarianceMatrix<pcl::PointCloud<PointT>,float>( cov, *cloud, &(neighs[pid]), NULL, NULL );
            // solve for neighbourhood biggest eigen value
            Eigen::SelfAdjointEigenSolver< Eigen::Matrix<Scalar, 3, 3> > es;
            es.compute( cov );
            // get eigen vector for biggest eigen value
            const int max_eig_val_id = std::distance( es.eigenvalues().data(), std::max_element( es.eigenvalues().data(), es.eigenvalues().data()+3 ) );
            Eigen::Matrix<Scalar,3,1> eig2 = es.eigenvectors().col( max_eig_val_id ).normalized();
            // compute line direction perpendicular to eigen vector
            Eigen::Matrix<Scalar,3,1> p0 = cloud->at(pid).getVector3fMap();                                      // point on line

            lines.emplace_back( TLine(p0, eig2) );
#endif
        }

        return EXIT_SUCCESS;
    }

    template <class _PointsContainerT, typename _PclPointT, class _PrimitiveT, typename Scalar >
    inline int
    Pearl::expand(
            std::vector<int>                 & labels
            , boost::shared_ptr<pcl::PointCloud<_PclPointT> > cloud
            , _PointsContainerT         const& points
            , std::vector<int>          const* indices
            , std::vector<_PrimitiveT>  const& lines
            , Scalar                    const  gammasqr
            , Scalar                    const  beta
            , Params                    const& params )
    {
        using rapter::PidT;
        using rapter::LidT;
        using rapter::ULidT;

        typedef typename _PointsContainerT::PrimitiveT PointPrimitiveT;

        if ( indices ) { std::cerr << __PRETTY_FUNCTION__ << "]: indices must be NULL, not implemented yet..." << std::endl; return EXIT_FAILURE; }

        const size_t num_pixels = indices ? indices->size() : cloud->size();
        const size_t num_labels = lines.size();
        bool noAssignmentsYet = false;

        rapter::GidPidVectorMap populations;
        if ( !labels.size() )
        {
            std::cout << "first" << std::endl;
            noAssignmentsYet = true;
        }
        else
        {
            assert( points.size() == labels.size() );
            for ( size_t pid = 0; pid != points.size(); ++pid )
            {
                rapter::containers::add( populations, labels.at(pid), static_cast<PidT>(pid) );
            }
        }

        labels.resize( num_pixels );

        // first set up the array for data costs
        typename _PrimitiveT::ExtremaT extrema;

        const float intMultSqr = params.int_mult * params.int_mult;
        std::cout << "intMult: " << params.int_mult << std::endl;
        Scalar *data = new Scalar[ num_pixels * num_labels ];
        int data_zero_cnt = 0, data_max_cnt = 0;
        for ( size_t pid = 0; pid != num_pixels; ++pid )
            for ( size_t line_id = 0; line_id != num_labels; ++line_id )
            {
                //float dist = std::abs( lines[line_id].getDistance( /*     point: */ cloud->at( indices ? (*indices)[pid] : pid ).getVector3fMap() ) ); // abs() added by Aron 18/1/2015
                float dist = std::numeric_limits<Scalar>::max();
                if ( noAssignmentsYet || populations.find(line_id) == populations.end() || populations[line_id].size() == 0 )
                {
                    dist = std::abs( lines[line_id].getDistance( /* point: */ cloud->at( indices ? (*indices)[pid] : pid ).getVector3fMap() ) ); // abs() added by Aron 18/1/2015
                }
                else
                {
                    lines[line_id].template getExtent<PointPrimitiveT>( extrema, points, params.scale, &populations.at(line_id), /* force_axis_aligned: */ true );
                    dist = std::abs( lines[line_id].getFiniteDistance( extrema,
                                                                       cloud->at( indices ? (*indices)[pid] : pid ).getVector3fMap() ) ); // abs() added by Aron 18/1/2015
                }

                //std::cout << "dist: " << dist;
                dist *= dist * intMultSqr;
                // std::cout << ", data: " << dist << std::endl;

                if ( (dist != dist) || (dist < 0) || (dist > params.data_max) )
                {
                    dist = params.data_max;
                    ++data_max_cnt;
                }
                else if (dist == 0)
                    ++data_zero_cnt;

                data[ pid * num_labels + line_id ] = dist;
            }
        std::cout << "Zero datacosts: " << data_zero_cnt << "/" << num_pixels << "(" << Scalar(data_zero_cnt)/num_pixels*Scalar(100.) << "%)"
                  << ", capped datacost: " << data_max_cnt << "/" << num_pixels << "(" << Scalar(data_max_cnt)/num_pixels*Scalar(100.) << "%)\n";

        // next set up the array for smooth costs
        //const int smooth_2 = params.lambdas(2)/2;
        Scalar *smooth = new Scalar[ num_labels * num_labels ];
        for ( ULidT l1 = 0; l1 < num_labels; ++l1 )
            for ( ULidT l2 = 0; l2 < num_labels; ++l2 )
                smooth[l1+l2*num_labels] = /*smooth_2 * */ (l1 != l2); // dirac/potts

        std::vector<std::vector<int  > > neighs;
        std::vector<std::vector<float> > sqr_dists;
        rapter::processing::getNeighbourhoodIndices( neighs
                                                   , cloud
                                                   , indices
                                                   , &sqr_dists
                                                   , params.max_neighbourhood_size   // 10
                                                   , params.scale // 0.01f
                                                   );
        try
        {
            gco::GCoptimizationGeneralGraph *gc = new gco::GCoptimizationGeneralGraph(num_pixels,num_labels);
            gc->setDataCost  ( data   ); // unary
            gc->setSmoothCost( smooth ); // pairwise labelwise
            gc->setLabelCost ( beta   ); // complexity ( number of labels)

            // set neighbourhoods
            std::vector<float> neighvals; neighvals.reserve( neighs.size() * 15 );
            for ( size_t pid = 0; pid != neighs.size(); ++pid )
                for ( size_t pid2 = 1; pid2 != neighs[pid].size(); ++pid2 ) // don't count own
                {
                    float distsqr  = sqr_dists[pid][pid2] * params.int_mult * params.int_mult;
                    float neighval = std::abs( params.lambdas(2) * exp( -1.f * distsqr / gammasqr) );

                    if ( neighval > (long)INT_MAX )
                        std::cerr << "neighbours exceeds limits: " << neighval << std::endl;

                    if ( neighval < 0 )
                        std::cerr << "neighval: " << neighval << std::endl;

                    gc->setNeighbors( pid, neighs[pid][pid2], neighval ); // pairwise - point wise
                    neighvals.push_back( neighval );
                }

            // debug
            {
                float sum = std::accumulate( neighvals.begin(), neighvals.end(), 0 );
                std::cout << "neighval avg: " << sum / (float)neighvals.size();
                std::sort( neighvals.begin(), neighvals.end() );
                std::cout << ", median: " << neighvals[ neighvals.size() / 2 ];
                float n_min = neighvals[0];
                std::cout << ", neighval_min: " << n_min << ", neighval_max: " << neighvals.back() << std::endl;
            }

            if ( neighvals[ neighvals.size() / 2 ] > 0 )
            {
                std::cout << params.beta << ", " << params.lambdas(2) << ", " << params.pottsweight;
                if ( (params.lambdas(2) * params.pottsweight < 1e10f) )
                {

                    printf("\tBefore optimization energy is %f", gc->compute_energy() ); fflush(stdout);
                    gc->expansion( 10 );// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
                    printf("\tAfter optimization energy is %f\n",gc->compute_energy());

                    // copy output
                    for ( size_t pid = 0; pid != num_pixels; ++pid )
                    {
                        labels[ pid ] = gc->whatLabel( pid );
                    }
                }
                else
                {
                    std::cout << "skipping " << std::endl;
                    for ( size_t pid = 0; pid != num_pixels; ++pid )
                    {
                        labels[ pid ] = 0;
                    }
                }
            }
            else
                std::cerr << "[" << __func__ << "]: " << "effective pairwise median !> 0, so not running" << std::endl;

            // cleanup
            if ( gc     ) { delete gc; gc = NULL; }
        }
        catch ( gco::GCException e )
        {
            e.Report();
        }
        catch ( std::exception &e )
        {
            std::cerr << "[" << __func__ << "]: " << e.what() << std::endl;
        }

        // cleanup
        if ( data   ) { delete [] data  ; data   = NULL; }
        if ( smooth ) { delete [] smooth; smooth = NULL; }

        return EXIT_SUCCESS;
    }

    template <typename _PointContainerT, class _PrimitiveT> inline int
    Pearl::refit( std::vector<_PrimitiveT>        &lines
                  , std::vector<int> const& labels
                  , _PointContainerT          const& cloud
                  , Params           const& params
                  //, std::vector<int> const* indices // TODO: this is not used...
                  )
    {
        using std::vector;
        typedef typename _PrimitiveT::Scalar Scalar;
        typedef typename _PointContainerT::PointType PclPointT;
        //typedef typename _PointContainerT::value_type PointPrimitiveT;

//        typedef pcl::PointXYZRGB _PclPointT;
//        typedef pcl::PointCloud<_PclPointT> _PclCloudT;
//        typedef typename _PclCloudT::Ptr _PclCloudPtrT;


//        _PclCloudPtrT pclCloud( new _PclCloudT() );
//        PointPrimitiveT::template toCloud<_PclCloudPtrT, _PointContainerT, GF2::PCLPointAllocator<PointPrimitiveT::Dim> >
//                ( pclCloud, points );

        // assign points to lines
        vector<vector<int> > lines_points( lines.size() );
        for ( size_t pid = 0; pid != labels.size(); ++pid )
            lines_points[ labels[pid] ].push_back( pid );

        if ( _PrimitiveT::EmbedSpaceDim == 2 )
        {
            for ( size_t line_id = 0; line_id != lines.size(); ++line_id )
            {
                if ( lines_points[line_id].size() > 1 )
                    smartgeometry::geometry::fitLinearPrimitive<_PointContainerT,Scalar,6>( /*           output: */ lines[line_id].coeffs()
                                                                                 , /*         points: */ cloud
                                                                                 , /*          scale: */ params.scale
                                                                                 , /*        indices: */ &(lines_points[line_id])
                                                                                 , /*    refit times: */ 2
                                                                                 , /* use input line: */ true );
            }//...for
        }//...embedspace=2
        else //PlANES, embedspace=3
        {
#if 0
            pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
            vptr->addPointCloud<PclPointT>( cloud.makeShared(), "cloud" );
            pcl::PointCloud<pcl::PointXYZRGB>::Ptr colourCloud( new pcl::PointCloud<pcl::PointXYZRGB>() );
            vptr->setBackgroundColor( .7, .7, .7 );
            pcl::ModelCoefficients modelCoeffs; modelCoeffs.values.resize(4);
#endif
            #pragma omp parallel for
            for ( size_t planeId = 0; planeId < lines.size(); ++planeId )
            {
#if 0
                vptr->removePointCloud( "colourCloud" );
                vptr->removeAllShapes();
#endif
                if ( lines_points[planeId].size() > 1 )
                {
                    //Eigen::Map< Eigen::Matrix<Scalar,4,1> > map( lines[planeId].coeffs().data(), 4 );
                    Eigen::Matrix<Scalar,4,1> plane4; plane4 << lines[planeId].normal(), lines[planeId].getDistance( Eigen::Matrix<Scalar,3,1>::Zero() );
#if 0
                    std::copy( plane4.data(), plane4.data()+4, modelCoeffs.values.data() );
                    vptr->addPlane( modelCoeffs, lines[planeId].template pos()(0)
                                    , lines[planeId].template pos()(1)
                                    , lines[planeId].template pos()(2)
                                    , "plane4" );
                    //vptr->spin();
#endif

                    smartgeometry::geometry::fitLinearPrimitive<_PointContainerT,Scalar,4>( /*           output: */ plane4 //lines[planeId].coeffs().template head<4>()
                                                                                   , /*         points: */ cloud
                                                                                   , /*          scale: */ params.scale
                                                                                   , /*        indices: */ &(lines_points[planeId])
                                                                                   , /*    refit times: */ 2
                                                                                   , /* use input line: */ true );
#if 0
                    std::cout << "refit from " << lines[planeId].toString() << ", to " << plane4.transpose() ;
                    {
                        pcl::copyPointCloud( cloud, lines_points[planeId], *colourCloud );
                        vptr->addPointCloud( colourCloud, "colourCloud" );
                         //lines_points[planeId]
                    }
#endif

                    _PrimitiveT plane6( /*     x0: */ -plane4.template head<3>() * plane4(3)
                                      , /* normal: */ plane4.template head<3>() );
                    plane6.copyTagsFrom( lines[planeId] );
                    lines[planeId] = plane6;
                    //std::cout << " and then \n\t" << plane6.toString() << std::endl;
#if 0
                    std::copy( plane4.data(), plane4.data()+4, modelCoeffs.values.data() );
                    vptr->addPlane( modelCoeffs, plane6.template pos()(0)
                                    , plane6.template pos()(1)
                                    , plane6.template pos()(2)
                                    ,"refit4" );
#endif

                    //vptr->spin();
                }
            } //...for
        } //...planes

        return EXIT_SUCCESS;
    } //...refit()


    inline int
    Pearl::getActiveLabelsCount( std::vector<int>         const& labels
                                 , std::vector<unsigned>       * active_labels_arg )
    {
        std::vector<unsigned> active_labels;
        for ( size_t pid = 0; pid != labels.size(); ++pid )
        {
            if ( static_cast<int>(active_labels.size()) <= labels[pid] )
                active_labels.resize( labels[pid]+1, 0 );
            active_labels[ labels[pid] ] = 1;
        }

        if ( active_labels_arg )
            *active_labels_arg = active_labels;

        return std::accumulate( active_labels.begin(), active_labels.end(), 0 );
    }
} // nsam

#endif // __RAPTER_PEARL_HPP__
