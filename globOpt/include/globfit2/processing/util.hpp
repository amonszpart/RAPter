#ifndef GF2_PROC_UTIL_HPP
#define GF2_PROC_UTIL_HPP

#include <map>
#include <set>
#include <vector>
#include "pcl/search/kdtree.h"
#include "globfit2/util/containers.hpp" // add()

/*!   \brief GlobOpt main namespace. */
namespace GF2 {

    typedef std::map   < int, int >         GidIntMap;
    typedef std::set   < int >              PidSet;
    typedef std::vector< int >              PidVector;
    typedef std::map   < int, PidSet >      GidPidSetMap;
    typedef std::map   < int, PidVector >   GidPidVectorMap;

    /*! \brief Various 3D processing snippets that don't fit elsewhere. */
    namespace processing
    {

        /*! \brief Calculate the number of points that are assigned to a GID (group id).
         *
         *  \tparam                 Concept: std::map< int, int >. Key: gid, value: population (\#points assigned to a group id (GID).
         *  \param[out] populations Contains an int value for each GID that occurred in points. The value is the count, how many points had this GID.
         *  \param[in]  points      The points containing gids retrievable by \code points[pid].getTag( _PointPrimitiveT::GID ) \endcode.
         *  \return                 EXIT_SUCCESS
        */
        template <class _GidIntMap, class _PointContainerT > inline int
        calcPopulations( _GidIntMap & populations, _PointContainerT const& points )
        {
            typedef typename _PointContainerT::value_type _PointPrimitiveT;

            for ( size_t pid = 0; pid != points.size(); ++pid )
            {
                ++populations[ points[pid].getTag(_PointPrimitiveT::GID) ];
            }

            return EXIT_SUCCESS;
        } //...calcPopulation

        /*! \brief Calculate the number of points that are assigned to a GID (group id).
        *   \tparam                 Concept: std::map< int, std::set<int> >. Key: GID, value: list of point ids that have that GID. (\#points assigned to a group id (GID).
        *   \param[out] populations Contains an int value for each GID that occurred in points. The value is the count, how many points had this GID.
        *   \param[in]  points      The points containing gids retrievable by \code points[pid].getTag( _PointPrimitiveT::GID ) \endcode.
        *   \return                 EXIT_SUCCESS
        */
        template <class _GidIntSetMap, class _PointContainerT > inline int
        getPopulations( _GidIntSetMap & populations, _PointContainerT const& points )
        {
            typedef typename _PointContainerT::value_type _PointPrimitiveT;

            for ( size_t pid = 0; pid != points.size(); ++pid )
            {
                const int gid = points[pid].getTag( _PointPrimitiveT::GID );
                containers::add( populations, gid, static_cast<int>(pid) );
            }

            return EXIT_SUCCESS;
        } //...getPopulations

        /*! \brief Adds angles from generator to in/out vector \p angles.
         *
         *  \tparam _AnglesContainerT Concept: std::vector<_Scalar>.
         *  \tparam _Scalar           Floating point precision type.
         *  \param[in,out] angles     Container (vector) to append to.
         *  \param[in] angle_gen      Generator element. 0, \p angle_gen, 2 * \p angle_gen, ..., M_PI will be appended.
         *  \param[in] verbose        Enable logging.
         */
        template <class _AnglesContainerT, class _Scalar> inline int
        appendAnglesFromGenerators( _AnglesContainerT &angles, std::vector<_Scalar> &angle_gens, char verbose = true )
        {
            std::cout << "[" << __func__ << "]: " << "ASSUMING DEGREES" << std::endl;
            std::set<_Scalar> angles_set;
            // copy
            angles_set.insert( angles.begin(), angles.end() );

            // insert 0 element
            angles_set.insert( _Scalar(0) );

            for ( int i = 0; i != angle_gens.size(); ++i )
            {
                _Scalar angle_gen_rad = angle_gens[i] * M_PI/_Scalar(180.);
                // generate
                for ( _Scalar angle = angle_gen_rad; angle < M_PI; angle+= angle_gen_rad )
                    angles_set.insert( angle );
            }

            // insert 0 element
            angles_set.insert( _Scalar(M_PI) );

            angles.resize( angles_set.size() );
            std::copy( angles_set.begin(), angles_set.end(), angles.begin() );

            // print
            if ( verbose )
            {
                std::cout << "Desired angles: {";
                for ( size_t vi=0;vi!=angles.size();++vi)
                    std::cout << angles[vi] << "(" << angles[vi] * _Scalar(180.) / M_PI << ")" << ((vi==angles.size()-1) ? "" : ", ");
                std::cout << "}\n";
            }

            return EXIT_SUCCESS;
        } //...appendAnglesFromGenerator()
        
        /*! \brief Applies a functor to each primitive.
         *  \tparam _PrimitiveT          Primitive wrapper class. Concept: LinePrimitive2.
         *  \tparam _inner_iterator      Iterates over the mapped_value of _PrimitiveContainerT. Concept: std::vector<int>::iterator.
         *  \tparam _PrimitiveContainerT Stores all lines of a patch under it's GID key. Concept: std::map< int, std::vector<_PrimitiveT> >.
         *  \tparam _FunctorT            Has an eval(_PrimitiveT &) function to be called for each primitive.
         *  \param[in,out] prims         The primitives to transform.
         *  \param[in] functor           An instance of the functor that should be called for each primitive.
         */
        template < class _PrimitiveT
                 , typename _inner_iterator
                 , class _FunctorT
                 , class _PrimitiveContainerT
                 >
        int transformPrimitivesMap( _PrimitiveContainerT       & primitives
                                  , _FunctorT             const& functor = NULL
                                  )
        {
            typedef typename _PrimitiveContainerT::iterator outer_iterator;

            int ret = 0;

            // for all patches
            for ( outer_iterator outer_it  = primitives.begin();
                                 outer_it != primitives.end();
                               ++outer_it )
            {
                for ( _inner_iterator inner_it  = containers::valueOf<_PrimitiveT>(outer_it).begin();
                                      inner_it != containers::valueOf<_PrimitiveT>(outer_it).end();
                                    ++inner_it )
                {
                    // apply functor to primitive
                    ret += functor.eval( *inner_it );
                } //...inner loop (for each direction in patch)
            } //...outer loop (for each patch)

            return ret;
        } //...TransformPrimitivesMap

        /*! \brief Applies a functor to each primitive.
         *  \tparam _PrimitiveT          Primitive wrapper class. Concept: LinePrimitive2.
         *  \tparam _inner_iterator      Iterates over the mapped_value of _PrimitiveContainerT. Concept: std::vector<int>::iterator.
         *  \tparam _PrimitiveContainerT Stores all lines of a patch under it's GID key. Concept: std::map< int, std::vector<_PrimitiveT> >.
         *  \tparam _FunctorT            Has an eval(_PrimitiveT &) function to be called for each primitive.
         *  \param[in,out] prims         The primitives to transform.
         *  \param[in] functor           An instance of the functor that should be called for each primitive.
         */
        template < class _PrimitiveT
                 , typename _inner_const_iterator
                 , class _FunctorT
                 , class _PrimitiveContainerT
                 >
        int filterPrimitives( _PrimitiveContainerT  const& primitives
                            , _FunctorT                  & functor
                            )
        {
            typedef typename _PrimitiveContainerT::const_iterator outer_const_iterator;

            int ret = 0;

            // for all patches
            for ( outer_const_iterator outer_it  = primitives.begin();
                                       outer_it != primitives.end();
                                     ++outer_it )
            {
                int lid = 0;
                for ( _inner_const_iterator inner_it  = containers::valueOf<_PrimitiveT>(outer_it).begin();
                                            inner_it != containers::valueOf<_PrimitiveT>(outer_it).end();
                                          ++inner_it, ++lid )
                {
                    // apply functor to primitive
                    ret += functor.eval( *inner_it, lid );
                } //...inner loop (for each direction in patch)
            } //...outer loop (for each patch)

            return ret;
        } //...TransformPrimitivesMap

        /*! \brief Evaluate a functor to each primitive to test if it have to be destroyed
         *  \tparam _PrimitiveT          Primitive wrapper class. Concept: LinePrimitive2.
         *  \tparam _inner_iterator      Iterates over the mapped_value of _PrimitiveContainerT. Concept: std::vector<int>::iterator.
         *  \tparam _PrimitiveContainerT Stores all lines of a patch under it's GID key. Concept: std::map< int, std::vector<_PrimitiveT> >.
         *  \tparam _FunctorT            Has a bool eval(_PrimitiveT &) function to be called for each primitive.
         *  \param[in,out] prims         The primitives to transform.
         *  \param[in] functor           An instance of the functor that should be called for each primitive.
         *  \return                      Number of erased primitives
         */
        template < class _PrimitiveT
                 , typename _inner_iterator
                 , class _FunctorT
                 , class _PrimitiveContainerT
                 >
        int erasePrimitives( _PrimitiveContainerT  & primitives
                            , _FunctorT            & functor
                            )
        {
            typedef typename _PrimitiveContainerT::iterator outer_iterator;

            int ret = 0;

            // for all patches
            for ( outer_iterator outer_it  = primitives.begin();
                                 outer_it != primitives.end();
                                       /* do nothing */ )
            {
                for ( _inner_iterator inner_it  = containers::valueOf<_PrimitiveT>(outer_it).begin();
                                            inner_it != containers::valueOf<_PrimitiveT>(outer_it).end();
                                          /* Do nothing */)
                {
                    // apply functor to primitive
                    if( functor.eval( *inner_it ) ){
                        ++ret;
                        inner_it = containers::valueOf<_PrimitiveT>(outer_it).erase(inner_it);
                    }else{
                        ++inner_it ;
                    }
                } //...inner loop (for each direction in patch)

                if(containers::valueOf<_PrimitiveT>(outer_it).size() == 0){
                    outer_it = primitives.erase(outer_it);
                }else
                    ++outer_it;

            } //...outer loop (for each patch)

            return ret;
        } //...TransformPrimitivesMap

        /*!
         * \brief Discard primitives not assigned to any point
         * Internally, use #erasePrimitives with nested functor
         *  \tparam _PrimitiveT          Primitive wrapper class. Concept: LinePrimitive2.
         *  \tparam _inner_iterator      Iterates over the mapped_value of _PrimitiveContainerT. Concept: std::vector<int>::iterator.
         *  \tparam _PointContainerT     Stores all points. Concept: GF2::PointPrimitive
         *  \tparam _PrimitiveContainerT Stores all lines of a patch under it's GID key. Concept: std::map< int, std::vector<_PrimitiveT> >.
         *  \param[in,out] prims         The primitives to transform.
         *  \param[in] functor           An instance of the functor that should be called for each primitive.
         *  \return                      Number of erased primitives
         */
        template < class _PrimitiveT
                 , typename _inner_iterator
                 , class _PointContainerT
                 , class _PrimitiveContainerT
                 >
        int eraseNonAssignedPrimitives( _PrimitiveContainerT   & primitives,
                                         _PointContainerT const & points) {
            // In a first run over all points, collect all used GID

            typedef typename _PointContainerT::value_type PointT;

            struct {
                std::set<int> gids;

                //! Return true when the id cannot be found in the set -> induce destruction
                inline
                bool eval(const _PrimitiveT& p) const{
                    return gids.find( p.getTag(_PrimitiveT::GID) ) == gids.end();
                }
            } nestedFunctor;

            for(typename _PointContainerT::const_iterator it  = points.begin();
                                                          it != points.end();
                                                        ++it)
                nestedFunctor.gids.insert((*it).getTag(PointT::GID));

            return erasePrimitives<_PrimitiveT, _inner_iterator>(primitives, nestedFunctor);
        }
        
        
        /**
         * @brief fitLine               [Re]Fits 3D line to a [part of a] pointcloud.
         * @param[out] line             Output line, and possibly input line to refit, if \param start_from_input_line is true.
         * @param cloud                 Points to fit to. Must have methods operator[] and pos()->Eigen::Vector3f. If \param p_indices!=NULL, must have at least max(*p_indices) points.
         * @param scale                 Distance where point get's zero weight
         * @param p_indices             Indices to use from cloud. Can be NULL, in which case the whole cloud is used.
         * @param refit                 How many refit iterations. 0 means fit once, and refit 0 times, obviously (TODO to fix...).
         * @param initial_line          Use this line to calculate weights on the 0th iteration already.
         * @param fit_pos_only          Use this to refit only the position of the line or both position and orientation. Requires initial_line!=NULL
         */
        template <int rows, class PrimitiveT, class PointsT, typename Scalar> inline int
        fitLinearPrimitive( PrimitiveT                      & primitive
                            , PointsT                  const& cloud
                            , Scalar                          scale
                            , std::vector<int>         const* p_indices             = NULL
                            , int                             refit                 = 0
                            , PrimitiveT               const* initial_line          = NULL
                            , bool                            fit_pos_only          = false
                            , bool                            debug                 = false )
        {
            //SG_STATIC_ASSERT( (rows == 4) || (rows == 6), smartgeometry_fit_linear_model_rows_not_4_or_6 );
            typedef Eigen::Matrix<Scalar,3,1> Position;

            // number of points to take into account
            const int N = p_indices ? p_indices->size() : cloud.size();

            // skip, if not enought points found to fit to
            if ( N < 2 ) { std::cerr << "[" << __func__ << "]: " << "can't fit line to less then 2 points..." << std::endl; return EXIT_FAILURE; }

            if ( fit_pos_only && initial_line == NULL ){ std::cerr << "[" << __func__ << "]: " << "can't fit position only without initial line..." << std::endl; return EXIT_FAILURE; }

            // copy input, if exists
            if ( initial_line )
                primitive = *initial_line;

            int iteration = 0; // track refit iterations
            do
            {
                // LeastSquares weights
                std::vector<Scalar> weights( N, 1.f );

                // calculate weights, if value in "line" already meaningful
                if ( initial_line || (iteration > 0) )
                {
                    // calculate distance from all points
                    for ( size_t point_id = 0; point_id != N; ++point_id )
                    {
                        weights[point_id] = primitive.getDistance( cloud[ p_indices ? (*p_indices)[point_id] : point_id ].pos() );
                    }

                    // the farther away, the smaller weight -->
                    // w_i = f( dist_i / scale ), dist_i < scale; f(x) = (x^2-1)^2
                    for ( size_t wi = 0; wi != weights.size(); ++wi )
                    {
                        if ( weights[wi] < scale )
                        {
                            weights[wi] /= scale;                                               // x = dist_i / scale
                            weights[wi] = (weights[wi] * weights[wi] - static_cast<Scalar>(1)); // x^2-1
                            weights[wi] *= weights[wi];                                         // (x^2-1)^2
                        }
                        else
                            weights[wi] = static_cast<Scalar>(0);                               // outside scale, truncated to 0

                    }
                }

                // compute centroid of cloud or selected points
                Position centroid( Position::Zero() );
                Scalar sumW = 0;
                for ( size_t point_id = 0; point_id != N; ++point_id )
                {
                    const unsigned int id = p_indices ? (*p_indices)[point_id] : point_id;
                    centroid += cloud[ id ].pos() * weights[ point_id ];
                    sumW     +=                     weights[ point_id ];
                }
                if ( sumW > Scalar(0.) )
                    centroid /= sumW;

                if ( fit_pos_only )
                {
                    primitive = PrimitiveT(centroid.template head<3>(),
                                           initial_line->dir());
                    continue; // we can stop now and go to next iteration
                }

                // compute neighbourhood covariance matrix
                Eigen::Matrix<Scalar,3,3> cov( Eigen::Matrix<Scalar,3,3>::Zero() );
                for ( size_t point_id = 0; point_id != N; ++point_id )
                {
                    const unsigned int id = p_indices ? (*p_indices)[point_id] : point_id;
                    Position pos = cloud[ id ].template pos() - centroid; // eigen expression template
                    cov   += pos * pos.transpose() * weights[ point_id ];
                }                
                cov /= sumW;

                // solve for neighbourhood biggest eigen value
                Eigen::SelfAdjointEigenSolver< Eigen::Matrix<Scalar, 3, 3> > es;
                es.compute( cov );

                if ( rows == 6 ) // line -> dir ==
                {
                    // get eigen vector for biggest eigen value
                    const int max_eig_val_id = std::distance( es.eigenvalues().data(), std::max_element( es.eigenvalues().data(), es.eigenvalues().data()+3 ) );

                    // output line
                    primitive = PrimitiveT( centroid.template head<3>(),
                                            es.eigenvectors().col(max_eig_val_id).normalized() );
                }
                else if ( rows == 4 ) // plane
                {
                    // get eigen vector for biggest eigen value
                    const int min_eig_val_id = std::distance( es.eigenvalues().data(), std::min_element( es.eigenvalues().data(), es.eigenvalues().data()+3 ) );

                    // output line
                    primitive = PrimitiveT( centroid.template head<3>(),
                                            es.eigenvectors().col(min_eig_val_id).normalized() );
                }
                else
                    std::cerr << "[" << __func__ << "]: " << "lines(rows==6) or planes(rows==4), not rows == " << rows << std::endl;

#if DEBUG_FITLINE
                for ( int d = 0; d < 6; ++d )
                {
                    if ( primitive(d) != primitive(d) )
                    {
                        std::cerr << "centroid: \n" << centroid.transpose() << std::endl
                                  << "cov\n: " << cov << std::endl
                                  << "es.eigenvectors():\n" << es.eigenvectors() << std::endl;
                    }
                }

                if ( debug )
                {
                    std::cout << "[" << __func__ << "]: " << line.transpose() << std::endl;
                    std::cout << "[" << __func__ << "]: " << es.eigenvectors().col(0).transpose() << std::endl;
                    std::cout << "[" << __func__ << "]: " << es.eigenvectors().col(1).transpose() << std::endl;
                    std::cout << "[" << __func__ << "]: " << es.eigenvectors().col(2).transpose() << std::endl;

                }
#endif
            }
            while ( iteration++ < refit );

            return EXIT_SUCCESS;
        } // ... fitline
        
        /**
     * @brief getNeighbourhoodIndices   Get's a list of neighbours for each point in pointcloud/indices_arg, indices untested
     * @param[out] neighbour_indices    List of list of neighbour indices. One list for each point in cloud.
     * @param[in ] cloud                3D point cloud.
     * @param[in ] indices_arg          Optional, if given, selects point from cloud (untested).
     * @param[out] p_distances          Optional, neighbourhood squared distances arranged as neighbour_indices.
     * @param[in ] K                    Maximum number of neighbours
     * @param[in ] radius               Optional, maximum radius to look for K neighbours in
     * @param[in ] soft_radius          Return K neighbours even if some outside radius
     */
        template <typename MyPointT>
        inline int
        getNeighbourhoodIndices( std::vector<std::vector<int> >                            & neighbour_indices
                                , boost::shared_ptr<pcl::PointCloud<MyPointT> >              cloud
                                , std::vector<int>                                    const* indices_arg        = NULL
                                , std::vector<std::vector<float> >                         * p_distances        = NULL
                                , int                                                        K                  = 15
                                , float                                                      radius             = -1.f
                                , bool                                                       soft_radius        = false
                                )
        {
            // prepare output
            const int   N              = indices_arg ? indices_arg->size() : cloud->size();
            const bool  doRadiusSearch = radius > 0.f;

            neighbour_indices.resize( N );
            if ( p_distances )    p_distances->resize( N );

            // create KdTree
            typename pcl::search::KdTree<MyPointT>::Ptr tree( new pcl::search::KdTree<MyPointT> );
            if ( indices_arg )
            {
                pcl::IndicesPtr indices_ptr( new std::vector<int>() );
                *indices_ptr = *indices_arg; // copy indices
                tree->setInputCloud( cloud, indices_ptr );
            }
            else
                tree->setInputCloud( cloud );

            MyPointT            searchPoint;
            std::vector<float>  sqr_dists;
            int                 found_points_count = 0;
            for ( size_t pid = 0; pid != N; ++pid )
            {
                // copy search point
                searchPoint = indices_arg ? cloud->at( (*indices_arg)[pid] )
                                          : cloud->at( pid                 );

                // calculate neighbourhood indices
                if ( doRadiusSearch ) found_points_count = tree->radiusSearch  ( searchPoint, radius, neighbour_indices[pid], sqr_dists, /* all: */ 0 );
                else                  found_points_count = tree->nearestKSearch( searchPoint,      K, neighbour_indices[pid], sqr_dists    );

                if ( (found_points_count <2 ) && (soft_radius) )
                    found_points_count = tree->nearestKSearch( searchPoint,      K, neighbour_indices[pid], sqr_dists    );

                // output distances
                if ( found_points_count > 0 )
                {
                    if ( p_distances )
                    {
                        p_distances->at(pid) = sqr_dists;
                    }
                }
                else
                {
                    // clear results for this point
                    neighbour_indices[pid].resize(0);
                    if ( p_distances )
                        p_distances->at(pid).resize(0);
                    // report error
                    std::cerr << __func__ << ": no neighs found for point " << pid << std::endl;
                }
            }

            return EXIT_SUCCESS;
        }

        template < class _Scalar
                 , class _PointContainerT
                 , class _IndicesContainerT >
        inline Eigen::Matrix<_Scalar,3,1> getCentroid( _PointContainerT const& points, _IndicesContainerT const& indices )
        {
            Eigen::Matrix<_Scalar,3,1> centroid( Eigen::Matrix<_Scalar,3,1>::Zero() );
            for ( int pid_id = 0; pid_id != indices.size(); ++pid_id )
            {
                centroid += points[ indices[pid_id] ].template pos();
            }

            if ( indices.size() )
                centroid /= _Scalar( indices.size() );
            else
                std::cout << "empty..." << std::endl;
            return centroid;
        }

    } //...ns processing
} //...ns GF2

#endif // GF2_PROC_UTIL_HPP
