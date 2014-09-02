#ifndef __GF2_CANDIDATEGENERATOR_H__
#define __GF2_CANDIDATEGENERATOR_H__

#include <opencv2/core/core.hpp>        // Mat
#include <opencv2/highgui/highgui.hpp>  // imread
#include <opencv2/imgproc/imgproc.hpp>  // cvtColor
#include "pcl/common/intersections.h"   // lineWithLineIntersection
#include "pcl/kdtree/kdtree.h"          // nearestneighboursearch

#include "pcltools/util.hpp"            // addGaussianNoise, fitLinearPrimitive
#include "globfit2/optimization/patchDistanceFunctors.h" // FullLinkagePointPatchDistanceFunctor

namespace GF2
{
    struct CandidateGeneratorParams
    {
            enum RefitMode { SPATIAL, AVG_DIR };

            float     angle_limit                 = 0.08f;   //!< \brief angle threshold for similar lines
            float     angle_limit_div             = 10.f;    //!< \brief Determines, how much the angle_limit is divided by, to get the patch-similarity threshold
            float     patch_dist_limit            = 1.f;     //!< \brief Patchify takes "patch_dist_limit * scale" as maximum spatial distance
            int       nn_K                        = 15;      //!< \brief Number of points used for primitive fitting
            int       patch_population_limit      = 10;      //!< \brief Threshold, on when a patch can distribute its directions
            RefitMode refit_mode                  = AVG_DIR; //!< \brief Determines, how a patch gets it's final direction. A NL-LSQ to the points, or the average of local orientations.
            bool      show_fit_lines              = false;
            bool      show_candidates             = true;

            inline std::string printRefitMode() const
            {
                switch ( refit_mode )
                {
                    case AVG_DIR: return "avg_dir"; break;
                    case SPATIAL: return "spatial"; break;
                    default:      return "UNKNOWN"; break;
                }
            } // ...printRefitMode()

            inline int parseRefitMode( std::string const& refit_mode_string )
            {
                int err = EXIT_SUCCESS;

                if ( !refit_mode_string.compare("avg_dir") )
                {
                    this->refit_mode = AVG_DIR;
                }
                else if ( !refit_mode_string.compare("spatial") )
                {
                    this->refit_mode = SPATIAL;
                }
                else
                {
                    err = EXIT_FAILURE;
                    this->refit_mode = SPATIAL;
                    std::cerr << "[" << __func__ << "]: " << "Could NOT parse " << refit_mode_string << ", assuming SPATIAL" << std::endl;
                }

                return err;
            } // ...parseRefitMode()
    }; // ...struct CandidateGeneratorParams

    // generate
    /// (1) patchify
    //// (1.1) propose
    //// (1.2) agglomerative
    /// (2)
    class CandidateGenerator
    {
            typedef std::pair<int,int>      PidLid;

            template <typename _Scalar, typename _PrimitiveT>
            struct Patch : public std::vector<PidLid>
            {
                public:
                    typedef _Scalar     Scalar;
                    typedef _PrimitiveT PrimitiveT;

                    //using std::vector<PidLid>::vector;
                    Patch()
                        : _representative( Eigen::Matrix<_Scalar,3,1>::Zero(), Eigen::Matrix<_Scalar,3,1>::Ones() )
                        ,_n(0) {}
                    Patch( PidLid const& elem )
                        : _representative( Eigen::Matrix<_Scalar,3,1>::Zero(), Eigen::Matrix<_Scalar,3,1>::Ones() )
                        , _n(0) { this->push_back( elem ); }

                    inline _PrimitiveT      & getRepresentative()       { return _representative; }
                    inline _PrimitiveT const& getRepresentative() const { return _representative; }

                    template <class _PointContainerT>
                    inline void update( _PointContainerT const& points )
                    {
                        // take average of all points
                        Eigen::Matrix<_Scalar,3,1> pos( _representative.pos() * _n );
                        Eigen::Matrix<_Scalar,3,1> dir( _representative.dir() * _n );
                        for ( int pid_id = _n; pid_id < this->size(); ++pid_id, ++_n )
                        {
                            const int pid = this->operator []( pid_id ).first;
                            pos += points[ pid ].pos();
                            dir += points[ pid ].dir();
                        } // ... for all new points

                        _representative = _PrimitiveT( (pos / _n), (dir / _n).normalized() );
                    }

                    inline void update( _PrimitiveT const& other, _Scalar other_n )
                    {
                        _Scalar scalar_n( _n );
                        // take average of all points
                        Eigen::Matrix<_Scalar,3,1> pos( _representative.pos() * scalar_n );
                        Eigen::Matrix<_Scalar,3,1> dir( _representative.dir() * scalar_n );
                        pos += other.pos() * other_n;
                        dir += other.dir() * other_n;
                        scalar_n += other_n;

                        _representative = _PrimitiveT( pos / scalar_n, (dir / scalar_n).normalized() );
                        _n = scalar_n;
                    }

                    template <class _PointT>
                    inline void updateWithPoint( _PointT const& pnt )
                    {
                        _Scalar scalar_n( _n );
                        // take average of all points
                        Eigen::Matrix<_Scalar,3,1> pos( _representative.pos() * scalar_n );
                        Eigen::Matrix<_Scalar,3,1> dir( _representative.dir() * scalar_n );
                        pos += pnt.pos();
                        dir += pnt.dir();
                        scalar_n += Scalar(1);

                        _representative = _PrimitiveT( pos / scalar_n, (dir / scalar_n).normalized() );
                        _n = scalar_n;
                    }

                    inline typename Eigen::Matrix<Scalar,3,1> pos() const { return _representative.pos(); }
                    inline typename Eigen::Matrix<Scalar,3,1> dir() const { return _representative.dir(); }

                protected:
                    _PrimitiveT _representative;
                    int         _n;              //!< \brief How many points are averaged in _representative
            }; // ... struct Patch

        public:
            //! \brief Main functionality to generate lines from points.
            template <  class       PrimitivePrimitiveAngleFunctorT // concept: energyFunctors.h::PrimitivePrimitiveAngleFunctor
                      , class       PointPatchDistanceFunctorT      // concept: FullLinkagePointPatchDistanceFunctorT
                      , class       PrimitiveContainerT             // concept: std::vector<std::vector<LinePrimitive2>>
                      , class       PointContainerT                 // concept: std::vector<PointPrimitive>
                      , class       PrimitiveT                      = typename PrimitiveContainerT::value_type::value_type  // concept: LinePrimitive2
                      , typename    Scalar                          = typename PrimitiveT::Scalar                           // concept: float
                      >
            static inline int
            generate( PrimitiveContainerT             &  lines
                      , PointContainerT               &  points // non-const to be able to add group tags
                      , Scalar                   const   scale
                      , std::vector<Scalar>      const&  angles
                      , CandidateGeneratorParams const&  params );

            //! \brief patchify     Groups points into patches based on their local fit lines and positions.
            template < /*class        _PrimitivePrimitiveAngleFunctorT
                      ,*/ class       _PointPatchDistanceFunctorT
                      , class       _PointContainerT
                      , class       _PrimitiveT
                      , typename    _Scalar>
            static inline int patchify( std::vector<_PrimitiveT>           & patch_lines
                                        , _PointContainerT                 & points
                                        , _Scalar                     const  scale
                                        , std::vector<_Scalar>        const& angles
                                        , _PointPatchDistanceFunctorT const& pointPatchDistanceFunctor
                                        , CandidateGeneratorParams    const& params
                                        );

            //! \brief propose      Create local fits for local neighbourhoods.
            template <  class       TLines
                      , class       PointsPtrT>
            static inline int
            propose(  TLines                 & lines
                    , PointsPtrT        const  cloud
                    , std::vector<int>  const* indices
                    , int               const  K
                    , float             const  radius
                    , bool              const  soft_radius
                    , std::vector<int>       * mapping
                    );

            //! \brief                              Greedy region growing
            //! \tparam _PointPatchDistanceFunctorT Concept: FullLinkagePointPatchDistanceFunctor
            //! \tparam _PatchesT
            //! \tparam _PrimitiveContainerT        Concept: std::vector<LinePrimitive2>
            //! \tparam _PointContainerT            Concept: std::vector<PointPrimitive>
            //! \tparam _PrimitiveT                 Concept: LinePrimitive2
            //! \tparam _Scalar                     Concept: float
            //! \tparam _PointT                     Concept: PointPrimitive
            template < class       _PatchPatchDistanceFunctorT
                     , class       _PatchesT
                     , class       _PrimitiveContainerT
                     , class       _PointContainerT
                     , class       _PrimitiveT          = typename _PrimitiveContainerT::value_type
                     , typename    _Scalar              = typename _PrimitiveT::Scalar
                     , class       _PointT              = typename _PointContainerT::value_type
                     >
            static inline int
            regionGrow( _PointContainerT                 & points
                      , _PrimitiveContainerT        const& lines
                      , _Scalar                     const  scale
                      , std::vector<_Scalar>        const& /*angles*/
                      , CandidateGeneratorParams    const  /*params*/
                      , _PatchPatchDistanceFunctorT const& patchPatchDistanceFunctor
                      , int                         const  gid_tag_name              = _PointT::GID
                      , std::vector<int>            const* point_ids_arg             = NULL
                      , _PatchesT                        * groups_arg                = NULL );

            //! \brief image_2_2DCloud
            template <  class       PointAllocatorFunctorT
                      , class       PointContainerT>
            static inline int
            image_2_2DCloud( PointContainerT         & cloud
                             , std::string             img_path
                             , int                     N_samples
                             , float const             Z
                             , float const             scale );

            //! \brief image_2_2DCloud
            template <  class PointAllocatorFunctorT
                      , class PointContainerT>
            static inline int
            image_2_2DCloud( PointContainerT         & cloud
                           , cv::Mat         const & img
                           , int   const             N_samples
                           , float const             Z
                           , float const             scale );
        protected:
            template < typename _Scalar
                     , class _PointPatchDistanceFunctorT
                     , class _PatchT
                     , class _PointContainerT
                     >  static inline int
            _tagPointsFromGroups( _PointContainerT                 & points
                                , _PatchT                     const& groups
                                , _PointPatchDistanceFunctorT const& pointPatchDistanceFunctor
                                , int                         const  gid_tag_name );
    }; //...class CandidateGenerator
} // ...ns::GF2

//_____________________________________________________________________________________________________________________
// HPP
#include "globfit2/my_types.h"          //  PCLPointAllocator
//#include "globfit2/visualization/visualization.h"   // debug

namespace GF2
{
    template <  class       _PrimitivePrimitiveAngleFunctorT
              , class       _PatchPatchDistanceFunctorT                                                  // concept: FullLinkagePointPatchDistanceFunctor
              , class       _PrimitiveContainerT
              , class       _PointContainerT
              , class       _PrimitiveT
              , typename    _Scalar> int
    CandidateGenerator::generate( _PrimitiveContainerT             & lines
                                , _PointContainerT                 & points // non-const to be able to add group tags
                                , _Scalar                     const  scale
                                , std::vector<_Scalar>        const& angles
                                , CandidateGeneratorParams    const& params )
    {
        typedef typename _PointContainerT::value_type _PointPrimitiveT;

        _PatchPatchDistanceFunctorT patchPatchDistanceFunctor( params.patch_dist_limit * scale, params.angle_limit, scale );
        std::cout << "[" << __func__ << "]: " << "PatchPatchDistance by " << patchPatchDistanceFunctor.toString() << std::endl;

        // (1) Create patches
        std::vector<_PrimitiveT> patch_lines;
        patchify<_PatchPatchDistanceFunctorT>( patch_lines // tagged lines at GID with patch_id
                                               , points    // filled points with directions and tagged at GID with patch_id
                                               , scale
                                               , angles
                                               , patchPatchDistanceFunctor
                                               , params );
        if ( params.show_candidates )
        {
            _PrimitiveContainerT patches( patch_lines.size() );
            for ( int i = 0; i != patches.size(); ++i )
                patches[i].push_back( patch_lines[i] );
            //pcl::visualization::PCLVisualizer::Ptr vptr =
//                    Visualizer<_PrimitiveContainerT,_PointContainerT>::template show<_Scalar>( patches
//                                                                                               , points
//                                                                                               , /*     scale: */ scale
//                                                                                               , /*    colour: */ { 0.f, 0.f, 1.f}
//                                                                                               , /*      spin: */ true
//                                                                                               , /* show cons: */ NULL
//                                                                                               , /*       ids: */ false
//                                                                                               , /* tags_only: */ true );
            //vptr->spin();
        }
        // (2) Mix and Filter

        // count patch populations
        std::vector<int> populations( patch_lines.size() ); // populations[patch_id] = all points with GID==patch_id
        for ( size_t pid = 0; pid != points.size(); ++pid )
        {
            ++populations[ points[pid].getTag(_PointPrimitiveT::GID) ];

            if ( patch_lines[ points[pid].getTag(_PointPrimitiveT::GID) ].getTag(_PrimitiveT::GID) != points[pid].getTag(_PointPrimitiveT::GID)  )
                std::cerr << "[" << __func__ << "]: " << "wtf" << std::endl;
        }

        // convert local fits
        const _Scalar angle_limit( params.angle_limit / params.angle_limit_div );
        int nlines = 0;
        lines.resize( patch_lines.size() );
        for ( size_t patch_id0 = 0; patch_id0 != patch_lines.size(); ++patch_id0 ) // patch_id
        {
            lines[patch_id0].reserve( lines.size() );
            for ( size_t patch_id1 = patch_id0; patch_id1 != patch_lines.size(); ++patch_id1 )
            {
                bool add0, add1;
                add0 = add1 = (patch_id0 == patch_id1);

                _Scalar closest_angle = 0.f;
                {
                    _Scalar angdiff = _PrimitivePrimitiveAngleFunctorT::template eval<_Scalar>(patch_lines[patch_id0], patch_lines[patch_id1], angles, &closest_angle );
                    add0 = add1 = ( angdiff < angle_limit );
                }

                // filter by population
                if ( params.patch_population_limit > 0 )
                {
                    // don't copy direction to patch0 from too small patch1
                    if ( (patch_id0 != patch_id1) && (populations[patch_id1] < params.patch_population_limit) )
                    {
                        //std::cout << "skipping " << patch_id0 << "<--" << patch_id1 << "since pop is " << populations[patch_id1] << std::endl;
                        add0 = false;
                    }

                    // don't copy direction to patch1 from too small patch0
                    if ( (patch_id0 != patch_id1) && populations[patch_id0] < params.patch_population_limit )
                    {
                        //std::cout << "skipping " << patch_id0 << "-->" << patch_id1 << "since pop is " << populations[patch_id0] << std::endl;
                        add1 = false;
                    }

                    // debug check
                    if ( !populations[patch_id0] )
                        std::cerr << "[" << __func__ << "]: " << "pop[" << patch_id0 << "] is " << populations[patch_id0] << std::endl;
                    // debug check
                    if ( !populations[patch_id1] )
                        std::cerr << "[" << __func__ << "]: " << "pop[" << patch_id1 << "] is " << populations[patch_id1] << std::endl;
                }

                if ( !add0 && !add1 )
                    continue;

                Eigen::Matrix<_Scalar,3,1> dir0 = patch_lines[patch_id1].dir(), dir1 = patch_lines[patch_id0].dir();
                if ( (closest_angle > 0.f) && (closest_angle < M_PI) )
                {
                    dir0 = Eigen::AngleAxisf(-closest_angle, Eigen::Matrix<_Scalar,3,1>::UnitZ() ) * dir0;
                    dir1 = Eigen::AngleAxisf( closest_angle, Eigen::Matrix<_Scalar,3,1>::UnitZ() ) * dir1;
                }

                _PrimitiveT cand0 = _PrimitiveT( patch_lines[patch_id0].pos(), dir0 );
                _PrimitiveT cand1 = _PrimitiveT( patch_lines[patch_id1].pos(), dir1 );

                // copy line from pid to pid1
                if ( add0 )
                {
                    lines[patch_id0].emplace_back( cand0 );
                    lines[patch_id0].back().setTag( _PrimitiveT::GID    , patch_lines[patch_id0].getTag(_PrimitiveT::GID) );
                    lines[patch_id0].back().setTag( _PrimitiveT::DIR_GID, patch_lines[patch_id1].getTag(_PrimitiveT::GID) );
                    ++nlines;
                }

                if ( (patch_id0 != patch_id1) && add1 )
                {
                    // copy line from pid1 to pid
                    lines[patch_id1].emplace_back( cand1 );
                    lines[patch_id1].back().setTag( _PrimitiveT::GID    , patch_lines[patch_id1].getTag(_PrimitiveT::GID) );
                    lines[patch_id1].back().setTag( _PrimitiveT::DIR_GID, patch_lines[patch_id0].getTag(_PrimitiveT::GID) );
                    ++nlines;
                }
            }
        }
        //vptr->spin();

        std::cout << "[" << __func__ << "]: " << "finished generating, we now have " << nlines << " candidates" << std::endl;
        if ( params.show_candidates )
        {

//            Visualizer<_PrimitiveContainerT,_PointContainerT>::template show<_Scalar>( lines, points
//                                                                                     , /*     scale: */ scale
//                                                                                     , /*    colour: */ { 0.f, 0.f, 1.f}
//                                                                                     , /*      spin: */ true
//                                                                                     , /*    angles: */ NULL
//                                                                                     , /*  show_ids: */ false
//                                                                                     , /* tags_only: */ false );
        }

        return EXIT_SUCCESS;
    } // ...CandidateGenerator::generate()

    //! \brief patchify Groups unoriented points into oriented patches represented by a single primitive
    //!                 (1) local fits
    //!                 (2) group to patches
    //!                 (3) refit lines to patches
    //! \param patch_lines Patch representing lines output.
    //! \param points      Input points that get assigned to the patches.
    //! \param scale       Spatial scale to use for fits.
    //! \param angles      Desired angles to use for groupings
    template < class       _PointPatchDistanceFunctorT
             , class       _PointContainerT
             , class       _PrimitiveT
             , typename    _Scalar> int
    CandidateGenerator::patchify( std::vector<_PrimitiveT>         & patch_lines
                                , _PointContainerT                 & points // non-const, because group tags the points with primitive ids
                                , _Scalar                     const  scale
                                , std::vector<_Scalar>        const& angles
                                , _PointPatchDistanceFunctorT const& pointPatchDistanceFunctor
                                , CandidateGeneratorParams    const& params
                                )
    {
        typedef Patch<_Scalar,_PrimitiveT> PatchT;
        typedef std::vector< PatchT >      PatchesT;

        typedef typename _PointContainerT::value_type PointPrimitiveT;
        typedef pcl::PointCloud<pcl::PointXYZ>        CloudXYZ;

        // to pcl cloud
        CloudXYZ::Ptr cloud( new CloudXYZ() );
        PointPrimitiveT::template toCloud<CloudXYZ::Ptr, _PointContainerT, PCLPointAllocator<PointPrimitiveT::Dim> >( cloud, points );

        // (1) local fit lines from pcl cloud
        std::vector<_PrimitiveT> fit_lines;
        std::vector<int       > point_ids;
        {
            propose( /* [out]       lines: */ fit_lines
                     , /*            points: */ cloud
                     , /*           indices: */ NULL
                     , /*              nn_K: */ params.nn_K
                     , /*         nn_radius: */ scale
                     , /*       soft_radius: */ true
                     , /* [out]     mapping: */ &point_ids ); // contains point id for fit_line

            // copy line direction into point
            for ( int pid_id = 0; pid_id != point_ids.size(); ++pid_id )
            {
                const int pid = point_ids[pid_id];
                points[pid].coeffs().template segment<3>(3) = fit_lines[pid_id].dir();
            }
        } // ... (1) local fit

        // (2) group
        PatchesT groups;
        {
            std::cout << "[" << __func__ << "]: " << "grouping started..." << std::endl; fflush(stdout);
            regionGrow                                  ( /* [in/out]  points/pointsWGIDTag: */ points
                                                           , /* [in]                     lines: */ fit_lines
                                                           , /* [in]                     scale: */ scale
                                                           , /* [in]            desired_angles: */ angles
                                                           , /* [in]  CandidateGeneratorParams: */ params
                                                           , /* [in] PointPatchDistanceFunctor: */ pointPatchDistanceFunctor
                                                           , /* [in]              gid_tag_name: */ PointPrimitiveT::GID
                                                           , /* [in]            lines_pointids: */ &point_ids
                                                           , /* [out]          groups_pointids: */ &groups );
            std::cout << "[" << __func__ << "]: " << "group finished..." << std::endl; fflush(stdout);
        } // ... (2) group

        // (3) refit lines to patches
        if ( params.refit_mode == CandidateGeneratorParams::AVG_DIR )
        {
            patch_lines.reserve( groups.size() );
            for ( size_t gid = 0; gid != groups.size(); ++gid )
            {
                patch_lines.emplace_back( groups[gid].getRepresentative() );
                patch_lines.back().setTag( _PrimitiveT::GID, gid );
            }
            std::cout << "[" << __func__ << "]: " << "tagging of " << groups.size() << "patches finished" << std::endl;
        }
        else // params.refit_mode == CandidateGeneratorParams::SPATIAL )
        {
            patch_lines.reserve( groups.size() );
            for ( size_t gid = 0; gid != groups.size(); ++gid )
            {
                // gather indices from group
                std::vector<int> indices( groups[gid].size() );
                {
                    for ( size_t pid_id = 0; pid_id != groups[gid].size(); ++pid_id )
                        indices[ pid_id ] = groups[gid][pid_id].first;
                }

                // refit line to points
                if ( indices.size() > 1 ) // if more points assigned to group
                {

                    Eigen::Matrix<_Scalar,_PrimitiveT::Dim,1> line;
                    int err = smartgeometry::geometry::fitLinearPrimitive( /*           output: */ line
                                                                           , /*         points: */ *cloud
                                                                           , /*          scale: */ scale
                                                                           , /*        indices: */ &(indices)
                                                                           , /*    refit times: */ 3
                                                                           , /* use input line: */ false
                                                                           );
                    if ( EXIT_SUCCESS == err )
                    {
                        patch_lines.emplace_back( _PrimitiveT(line) );
                        patch_lines.back().setTag( _PrimitiveT::GID, gid );
                    }
                    else
                    {
                        std::cerr << "[" << __func__ << "]: " << "this shouldn't happen..." << std::endl; fflush(stderr);
                    }
                }
                else if ( indices.size() ) // if at least one point is assigned to group
                {
                    // gets here, if indices.size() == 1, just use the fit_line in this case
                    int point_id = indices[0];
                    int lid      = std::distance( point_ids.begin(), std::find(point_ids.begin(), point_ids.end(), point_id) );
                    if ( lid < 0 || lid >= fit_lines.size() )
                        std::cerr <<"[" << __func__ << "]: " << "could not find point_id " << point_id << " in fit_lines.point_ids...not copying line...\n";
                    else
                    {
                        patch_lines.push_back( fit_lines[lid] );
                        patch_lines.back().setTag( _PrimitiveT::GID, gid );
                    }
                } // ... if at least one point is assigned to group
            } // ... for groups
        } // ... (3) refit

        return EXIT_SUCCESS;
    } // ...CandidateGenerator::patchify()

    //! \brief Greedy region growing
    template < class       _PatchPatchDistanceFunctorT
             , class       _PatchesT
             , class       _PrimitiveContainerT
             , class       _PointContainerT
             , class       _PrimitiveT
             , typename    _Scalar
             , class       _PointT> int
    CandidateGenerator::regionGrow( _PointContainerT                 & points
                                  , _PrimitiveContainerT        const& lines
                                  , _Scalar                     const  scale
                                  , std::vector<_Scalar>        const& /*angles*/
                                  , CandidateGeneratorParams    const  /*params*/
                                  , _PatchPatchDistanceFunctorT const& patchPatchDistanceFunctor
                                  , int                         const  gid_tag_name
                                  , std::vector<int>            const* point_ids_arg
                                  , _PatchesT                        * groups_arg               )
    {
        std::cout << "[" << __func__ << "]: " << "running at " << patchPatchDistanceFunctor.getSpatialThreshold() << " radius threshold" << std::endl;
        if ( point_ids_arg && (lines.size() != point_ids_arg->size()) )
        {
            std::cerr << "[" << __func__ << "]: " << "lines.size() != point_ids_arg->size()" << "...exiting" << std::endl;
            return EXIT_FAILURE;
        }

        typedef typename _PatchesT::value_type   PatchT;
        typedef          std::vector<PatchT>     Patches;

        // create patches with a single point in them
        std::deque<int> starting_cands( lines.size() );
        for ( int lid = 0; lid != lines.size(); ++lid )
        {
            const int pid = point_ids_arg ? (*point_ids_arg)[ lid ] : lid;
            starting_cands.push_back( pid );
        }

        // debug show
        {
            std::vector<_PrimitiveContainerT> llines = { lines };
//            Visualizer<std::vector<_PrimitiveContainerT>,_PointContainerT>::template show<_Scalar>( llines
//                                                                                                  , points
//                                                                                                  , /*    scale: */ 0.05f
//                                                                                                  , /*   colour: */ { 0.f, 0.f, 1.f}
//                                                                                                  , /*     spin: */ true
//                                                                                                  , /*   angles: */ NULL
//                                                                                                  , /* show_ids: */ false );
        }

        // prebulid ann cloud
        pcl::PointCloud<pcl::PointXYZ>::Ptr ann_cloud( new pcl::PointCloud<pcl::PointXYZ>() );
        {
            ann_cloud->reserve( points.size() );
            for ( size_t pid = 0; pid != points.size(); ++pid )
            {
                pcl::PointXYZ pnt;
                pnt.x = points[pid].template pos()(0); // convert PointPrimitive to Eigen::Matrix, and get (0)
                pnt.y = points[pid].template pos()(1);
                pnt.z = points[pid].template pos()(2);
                ann_cloud->push_back( pnt );
            }
        }
        typename pcl::search::KdTree<pcl::PointXYZ>::Ptr tree( new pcl::search::KdTree<pcl::PointXYZ> );
        tree->setInputCloud( ann_cloud );

        Patches patches; patches.reserve( std::max(1.5*sqrt(points.size()),10.) );

        // get unassigned point
        std::vector<bool> assigned( points.size(), false );
        std::vector<bool> visited( points.size(), false );

        const int K = 20;
        std::vector<float>  sqr_dists( K );
        std::vector< int >  neighs( K );
        int                 found_points_count  = 0;
        pcl::PointXYZ       searchPoint;
        const _Scalar       spatial_thresh      = patchPatchDistanceFunctor.getSpatialThreshold();
//        const _Scalar       ang_thresh          = patchPatchDistanceFunctor.getAngularThreshold();
        const _Scalar       dist_weight         = 0.5*0.5; //sqrt(spatial_thresh); //(patchPatchDistanceFunctor.getAngularThreshold() * patchPatchDistanceFunctor.getAngularThreshold()) / (_Scalar(1.5) * spatial_thresh * spatial_thresh );
        const _Scalar       sqr_spatial_thresh  = spatial_thresh * spatial_thresh;
        const _Scalar       sqr_ang_thresh      = patchPatchDistanceFunctor.getAngularThreshold() * patchPatchDistanceFunctor.getAngularThreshold();
        const _Scalar       max_dist            = patchPatchDistanceFunctor.getSpatialThreshold() * _Scalar(3.5); // longest axis of ellipse)

        // look for neighbours, merge most similar
        while ( starting_cands.size() )
        {
            // remove point from unassigned
            int pid = starting_cands.front();
            starting_cands.pop_front();

            if ( visited[pid] && !assigned[pid] )
                std::cerr << "[" << __func__ << "]: " << " visited but not assigned...:JSDL:FJSD:LKFSDK:SFDJ" << std::endl;

            if ( visited[pid] )
                continue;

            // add to new cluster
            if ( !assigned[pid] )
            {
                patches.push_back( (PatchT){ PidLid(pid,-1) } );
                patches.back().update( points );
                assigned[ pid ] = true;
            }

            // look for unassigned neighbours
            searchPoint.x = points[ pid ].template pos()(0);
            searchPoint.y = points[ pid ].template pos()(1);
            searchPoint.z = points[ pid ].template pos()(2);
            found_points_count = tree->radiusSearch( searchPoint, max_dist, neighs, sqr_dists, 0);

            if ( neighs.size() != found_points_count )
                std::cerr << "as;dlkjfa;lkdsjfl;asdjkf neighs.size() " << neighs.size() << " != " << found_points_count << " found_points_count" << std::endl;
            for ( size_t pid_id = 0; pid_id != neighs.size(); ++pid_id )
            {
                const int pid2 = neighs[ pid_id ];
                if ( !assigned[pid2] )
                {
                    _Scalar diff = GF2::angleInRad( patches.back().template dir(), points[pid2].template dir() );
                    //if ( (sqrt(sqr_dists[pid_id]) < spatial_thresh) && (diff < ang_thresh) )                                      // original condition
                    if ( (dist_weight * sqr_dists[pid_id] / sqr_spatial_thresh + diff * diff / sqr_ang_thresh ) <= _Scalar(1) )     // ellipse longer along line
                    //if ( ((2*ang_thresh-diff)*(2*ang_thresh-diff) / sqr_ang_thresh) - (sqr_dists[pid_id] / sqr_spatial_thresh ) ) // hyperbole?
                    {
                        patches.back().push_back( PidLid(pid2,-1) );
                        patches.back().updateWithPoint( points[pid2] );
                        assigned[pid2] = true;
                        // enqueue for visit
                        starting_cands.push_front( pid2 );
                    }
                }
            }

            visited[pid] = true;
        }

        // copy patches to groups
        _PatchesT tmp_groups;                                       // local var, if output not needed
        _PatchesT *groups = groups_arg ? groups_arg : &tmp_groups;  // relay, if output needed
        (*groups).insert( (*groups).end(), patches.begin(), patches.end() );

        _tagPointsFromGroups<_Scalar>( points
                                     , *groups, patchPatchDistanceFunctor, gid_tag_name );

        return EXIT_SUCCESS;
    } // ...CandidateGenerator::regionGrow()

    template < typename _Scalar
             , class    _PointPatchDistanceFunctorT
             , class    _PatchT
             , class    _PointContainerT
             >  int
    CandidateGenerator::_tagPointsFromGroups( _PointContainerT                 & points
                                            , _PatchT                     const& groups
                                            , _PointPatchDistanceFunctorT const& pointPatchDistanceFunctor
                                            , int                         const  gid_tag_name )
    {
        // set group tag for grouped points
        std::vector<bool> visited( points.size(), false );
        for ( int gid = 0; gid != groups.size(); ++gid )
            for ( size_t pid_id = 0; pid_id != groups[gid].size(); ++pid_id )
            {
                const int pid = groups[gid][pid_id].first;

                visited.at(pid) = true; // TODO: change to []
                points[pid].setTag( gid_tag_name, gid );
            }

        // add left out points to closest patch
        for ( int pid = 0; pid != points.size(); ++pid )
        {
            if ( visited[pid] )
                continue;

            _Scalar     min_dist = std::numeric_limits<_Scalar>::max(), dist;
            int         min_gid  = 0;
            for ( int gid = 1; gid != groups.size(); ++gid )
            {
                if ( (dist = pointPatchDistanceFunctor.evalSpatial(pid, groups[gid], points)) < min_dist )
                {
                    min_dist = dist;
                    min_gid  = gid;
                }
            }

            points[pid].setTag( gid_tag_name, min_gid );
        }

        return EXIT_SUCCESS;
    } // ...CandidateGenerator::tagPointsFromGroups()

    template <  class TLines
              , class PointsPtrT> int
    CandidateGenerator:: propose(
            TLines                  & lines
            , PointsPtrT       const  cloud
            , std::vector<int> const* indices
            , int              const  K
            , float            const  radius
            , bool             const  soft_radius
            , std::vector<int>      * point_ids
            )
    {
        using std::vector;
        typedef typename TLines::value_type         TLine;
        typedef typename TLine::Scalar              Scalar;
        typedef typename PointsPtrT::element_type   PointsT;

        if ( indices ) { std::cerr << __PRETTY_FUNCTION__ << "]: indices must be NULL, not implemented yet..." << std::endl; return EXIT_FAILURE; }

        // get neighbourhoods
        std::vector< std::vector<int   > > neighs;
        std::vector< std::vector<Scalar> > sqr_dists;
        smartgeometry::getNeighbourhoodIndices( /*   [out] neighbours: */ neighs
                                                , /* [in]  pointCloud: */ cloud
                                                , /* [in]     indices: */ indices
                                                , /* [out]  sqr_dists: */ &sqr_dists
                                                , /* [in]        nn_K: */ K              // 5
                                                , /* [in]      radius: */ radius         // 0.02f
                                                , /* [in] soft_radius: */ soft_radius    // true
                                                );

        // only use, if more then 2 data-points
        if ( std::count_if( neighs.begin(), neighs.end(), [] (vector<int> const& n1) { return n1.size() > 2; } ) < 2 )
        {
            std::cerr << "[" << __func__ << "]: " << "not enough to work with (<2)...change scale " << radius << std::endl;
            return EXIT_SUCCESS;
        }

        // every point proposes primitive[s] using its neighbourhood
        int skipped = 0;
        for ( size_t pid = 0; pid != neighs.size(); ++pid )
        {
            // can't fit a line to 0 or 1 points
            if ( neighs[pid].size() < 2 )
            {
                ++skipped;
                continue;
            }

            Eigen::Matrix<Scalar,TLine::Dim,1> line;
            int err = smartgeometry::geometry::fitLinearPrimitive<PointsT,Scalar,TLine::Dim>( /*           output: */ line
                                                                                              , /*         points: */ *cloud
                                                                                              , /*          scale: */ radius
                                                                                              , /*        indices: */ &(neighs[pid])
                                                                                              , /*    refit times: */ 2
                                                                                              , /* use input line: */ false
                                                                                              );
            // debug
            //            {
            //                // debug
            //                pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
            //                {
            //                    vptr->setBackgroundColor( .5, .5, .6 );
            //                    PointsPtrT tmp_cloud( new PointsT() );
            //                    for ( int nid = 0; nid != neighs[pid].size(); ++nid )
            //                        tmp_cloud->push_back( cloud->at(neighs[pid][nid]) );
            //                    vptr->addPointCloud( tmp_cloud );

            //                }
            //                char plane_name[255]; sprintf( plane_name, "plane%d", pid );
            //                TLine::draw( line, cloud, radius, &(neighs[pid]), vptr, plane_name, .7, .7, .2, 0 );
            //                pcl::ModelCoefficients plane_coeffs;
            //                plane_coeffs.values.resize(4); std::copy( line.data(), line.data() + 4, plane_coeffs.values.begin() );

            //                sprintf( plane_name, "plane2_%d", pid );
            //                vptr->addPlane( plane_coeffs, cloud->at(pid).x, cloud->at(pid).y, cloud->at(pid).z, plane_name, 0 );

            //                vptr->spin();
            //            }

            if ( err == EXIT_SUCCESS )      lines.emplace_back( TLine(line) );
            if ( point_ids )
            {
                point_ids->emplace_back( pid );
            }
        }
        std::cout << "[" << __func__ << "]: "
                  << skipped << "/" << neighs.size() << ": " << skipped / static_cast<float>(neighs.size()) * 100.f << "% of points did not produce primitives, so the primitive count is:"
                  << lines.size() << " = " << lines.size() / static_cast<float>(neighs.size()) *100.f << "%" << std::endl;

        return EXIT_SUCCESS;
    } // ...CandidateGenerator::propose()

    /**
     * @brief GlobFit2::image_2_2DCloud Randomly samples a 2D gray image with 3D points at black pixels to create a 3D cloud.
     * @param cloud
     * @param img_path
     * @param N_samples
     * @param Z
     * @return EXIT_SUCCESS/EXIT_FAILURE
    */
    template <  class PointAllocatorFunctorT
              , class PointContainerT> int
    CandidateGenerator::image_2_2DCloud( PointContainerT              &cloud
                                         , std::string                 img_path
                                         , int                         N_samples
                                         , const float                 Z
                                         , const float scale )
    {
        return image_2_2DCloud( cloud, cv::imread(img_path,cv::IMREAD_UNCHANGED), N_samples, Z, scale );
    } // ...CandidateGenerator::image_2_2DCloud()

    template <  class PointAllocatorFunctorT
              , class PointContainerT> int
    CandidateGenerator::image_2_2DCloud( PointContainerT                &cloud
                                         , cv::Mat                const & img
                                         , const int                      N_samples
                                         , const float                    Z
                                         , const float                    scale )
    {
        typedef typename PointContainerT::PointType PointT;
        //   const float scale = 1.f;

//        if ( !cloud )
//            cloud = pcl::PointCloud<MyPoint>::Ptr( new pcl::PointCloud<MyPoint>() );
        cloud.reserve( N_samples );

        cv::Mat gray;
        if ( img.channels() > 1 )   cv::cvtColor( img, gray, cv::COLOR_RGB2GRAY );
        else                        gray = img;

        cv::threshold( gray, gray, 200., 255., cv::THRESH_BINARY );
        cv::Mat sampled( gray.clone() ); sampled.setTo(0);

        int   pix_count = gray.cols * gray.rows - cv::countNonZero( gray );
        float prob      = N_samples  / (float)pix_count;
        while ( cloud.size() < static_cast<size_t>(N_samples) )
        {
            pix_count = gray.cols * gray.rows - cv::countNonZero( gray );
            prob      = (N_samples - cloud.size())  / (float)pix_count;

            for ( int y = 0; (y < gray.rows) && (cloud.size() != (size_t)N_samples); ++y  )
            {
                for ( int x = 0; (x < gray.cols) && (cloud.size() != (size_t)N_samples); ++x )
                {
//                    float depth = (Eigen::Vector3f() << x, (gray.rows-y), 0).finished().norm();
//                    if ( RANDF() > depth / scale / sqrt(2.f) ) continue;

                    if ( (gray.at<uchar>(y,x) < 255) && (rand()/static_cast<double>(RAND_MAX) < prob) && (sampled.at<uchar>(y,x) == 0) )
                    {
                        //std::cout << "cv::sum( gray(cv::Range(y-1,y+1), cv::Range(x-1, x+1))): " << cv::sum( sampled(cv::Range(std::max(0,y-1),std::min(gray.rows,y+1)),
                        //                                                                                             cv::Range(std::max(0,x-1),std::min(gray.cols,x+1))) )[0] << std::endl;
                        if ( cv::sum( sampled(cv::Range(std::max(0,y-1),std::min(gray.rows,y+1)),
                                              cv::Range(std::max(0,x-1),std::min(gray.cols,x+1))) )[0] >= 3 ) continue;

//                        MyPoint pnt;
//                        pnt.x = x             / (float)gray.cols * scale;
//                        pnt.y = (gray.rows-y) / (float)gray.rows * scale;
//                        pnt.z = Z;

//                        am::util::pcl::setPointColor( pnt, 1, 0, 0 );
//                        cloud.push_back( pnt );


                        cloud.push_back( PointAllocatorFunctorT::template create<PointT>(
                                             (Eigen::Matrix<float,3,1>() <<
                                              x             / (float)gray.cols * scale,
                                              (gray.rows-y) / (float)gray.rows * scale,
                                              Z
                                              ).finished()
                                             )
                                         );

                        sampled.at<uchar>(y,x) = 1;
                    }
                }
            }
        }

        cv::imshow( "img", gray );
        cv::waitKey( 20 );

        return EXIT_SUCCESS;
    } // ...CandidateGenerator::image_2_2DCloud()
} // ... ns GF2


#endif // __GF2_CANDIDATEGENERATOR_H__

