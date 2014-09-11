#ifndef GF2_SEGMENTATION_H
#define GF2_SEGMENTATION_H

#include <utility> // pair
#include <vector>
#include "Eigen/Dense"

namespace GF2 {

namespace segmentation {
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
}

class Segmentation
{
    public:
        /*! \brief                  Step 1. Generates primitives from a cloud. Reads "cloud.ply" and saves "candidates.txt".
         *  \param argc             Contains --cloud cloud.ply, and --scale scale.
         *  \param argv             Contains --cloud cloud.ply, and --scale scale.
         *  \return                 EXIT_SUCCESS.
         *  \post                   "patches.txt" and "points_primitives.txt" on disk in "cloud.ply"'s parent path.
         */
        template < class _PrimitiveT
                 , class _PrimitiveContainerT
                 , class _PointPrimitiveT
                 , class _PointContainerT
                 , typename _Scalar
                 >
        static inline int
        segmentCli( int argc, char** argv );

        //! \param[in/out] points
        //! \param[in]     scale    Fit radius
        //! \param[in]     nn_K     Nearest neighbour count to fit primitive to.
        template < class     _PointPrimitiveT
                 , class     _PrimitiveT
                 , typename  _Scalar
                 , class     _PointContainerT
                 >
        static inline int
        orientPoints( _PointContainerT       &points
                    , _Scalar          const  scale
                    , int              const  nn_K   );

        /*!
         * \brief patchify Groups unoriented points into oriented patches represented by a single primitive
         *                 (2) group to patches
         *                 (3) refit lines to patches
         * \tparam _PrimitiveContainerT    Concept: std::map< int, std::vector<_PrimitiveT> >. Groups primitives by their GID.
         * \tparam _PatchPatchDistanceFunctorT Has an eval( point, primitive ) function for all points and primitives. Concept: \ref MyPointPatchDistanceFunctor
         * \param[out] patches                   Ready to process patches, each of them with one direction only. ( { 0: [Primitive00] }, { 1: [Primitive11] }, ... )
         * \param[in,out] points                 Input points that get assigned to the patches by their GID tag set.
         * \param[in]  scale                     Spatial scale to use for fits.
         * \param[in]  angles                    Desired angles to use for groupings.
         * \param[in]  patchPatchDistanceFunctor #regionGrow() uses the thresholds encoded to group points. The evalSpatial() function is used to assign orphan points.
         * \param[in]  nn_K                      Number of nearest neighbour points looked for in #regionGrow().
         */
        template <
                 class       _PrimitiveT
                 , typename    _Scalar
                 , class       _PrimitiveContainerT
                 , class       _PointContainerT
                 , class       _PatchPatchDistanceFunctorT
                 > static inline int
        patchify( _PrimitiveContainerT                   & patches
                , _PointContainerT                       & points
                , _Scalar                           const  scale
                , std::vector<_Scalar>              const& angles
                , _PatchPatchDistanceFunctorT       const& patchPatchDistanceFunctor
                , int                               const  nn_K
                );

        /*! \brief                               Greedy region growing
         *  \tparam _PrimitiveContainerT         Concept: std::vector<\ref GF2::LinePrimitive2>
         *  \tparam _PointContainerT             Concept: std::vector<\ref GF2::PointPrimitive>
         *  \tparam _PointPatchDistanceFunctorT  Concept: \ref RepresentativeSqrPatchPatchDistanceFunctorT.
         *  \tparam _PatchesT                    Concept: vector< \ref segmentation::Patch <_Scalar,_PrimitiveT> >
         *  \tparam _PrimitiveT                  Concept: \ref GF2::LinePrimitive2
         *  \tparam _Scalar                      Concept: float
         *  \tparam _PointT                      Concept: \ref GF2::PointPrimitive
         *  \param[in,out] points                Input points to create patches from. The \p gid_tag_name field of the points will be set according to their patch assignment.
         *  \param[out] groups_arg               Holds the point-id-groups, that can then be refit to to get a patch location and direction. Concept: vector<\ref segmentation::Patch>.
         *  \param[in] scale                     Spatial extent of the input. In practice unused, since the \p patchPatchDistanceFunctor was constructed with it.
         *  \param[in] patchPatchDistanceFunctor Takes two patches, and decides, whether they are similar enough to be merged.
         *                                       In practice, takes a patch that is currently grown, and a temporary patch
         *                                       that only contains a neighbouring point, and decides. See in \ref RepresentativeSqrPatchPatchDistanceFunctorT.
         *  \param[in] gid_tag_name              The key value of GID in _PointT. Suggested to be: _PointT::GID.
         *  \param[in] nn_K                      Number of nearest neighbour points looked for.
         */
        template < class       _PrimitiveT
                 , class       _PointContainerT
                 , class       _PatchPatchDistanceFunctorT
                 , class       _PatchesT
                 , typename    _Scalar              = typename _PrimitiveT::Scalar
                 , class       _PointT              = typename _PointContainerT::value_type
                 >
        static inline int
        regionGrow( _PointContainerT                 & points
                  , _PatchesT                        & groups_arg
                  , _Scalar                     const  /*scale*/
                  , _PatchPatchDistanceFunctorT const& patchPatchDistanceFunctor
                  , int                         const  gid_tag_name              //= _PointT::GID
                  , int                         const  nn_K );

        //! \brief propose      Create local fits to local neighbourhoods, these will be the point orientations.
        template <  class _PrimitiveContainerT
                  , class _PointContainerPtrT>
        static inline int
        fitLocal( _PrimitiveContainerT      & lines
                , _PointContainerPtrT  const  cloud
                , std::vector<int>     const* indices
                , int                  const  K
                , float                const  radius
                , bool                 const  soft_radius
                , std::vector<int>          * mapping
                );
    protected:
        template < class    _PointPrimitiveT
                 , typename _Scalar
                 , class    _PointPatchDistanceFunctorT
                 , class    _PatchT
                 , class    _PointContainerT
                 >  static inline int
        _tagPointsFromGroups( _PointContainerT                 & points
                            , _PatchT                     const& groups
                            , _PointPatchDistanceFunctorT const& pointPatchDistanceFunctor
                            , int                         const  gid_tag_name );
}; //...Segmentation

} //...ns GF2

#ifndef GF2_INC_SEGMENTATION_HPP
#   define GF2_INC_SEGMENTATION_HPP
#   include "globfit2/optimization/impl/segmentation.hpp"
#endif

#endif // GF2_SEGMENTATION_H
