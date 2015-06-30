#include "rapter/typedefs.h"
#include "rapter/optimization/segmentation.h"
#include "rapter/optimization/impl/segmentation.hpp"
#include "rapter/primitives/impl/planePrimitive.hpp"

namespace rapter
{
    template int
    Segmentation::segmentCli< rapter::_3d::PrimitiveT
                            , rapter::_3d::PrimitiveVectorT
                            , rapter::PointPrimitiveT
                            , rapter::PointContainerT
                            , rapter::Scalar
                            >
                            ( int argc, char** argv );

    template int
    Segmentation::segmentCli< rapter::_2d::PrimitiveT
                            , rapter::_2d::PrimitiveVectorT
                            , rapter::PointPrimitiveT
                            , rapter::PointContainerT
                            , rapter::Scalar
                            >
                            ( int argc, char** argv );

    template int
    Segmentation::orientPoints< rapter::PointPrimitiveT
                              , rapter::_3d::PrimitiveT
                              , rapter::Scalar
                              , rapter::PointContainerT
                              >
                              ( rapter::PointContainerT       &points
                              , rapter::Scalar          const  scale
                              , int                  const  nn_K
                              , int                  const  verbose );

    template int
    Segmentation::orientPoints< rapter::PointPrimitiveT
                              , rapter::_2d::PrimitiveT
                              , rapter::Scalar
                              , rapter::PointContainerT
                              >
                              ( rapter::PointContainerT       &points
                              , rapter::Scalar          const  scale
                              , int                  const  nn_K
                              , int                  const  verbose );


    template int
    Segmentation::patchify< rapter::_2d::PrimitiveT
            , rapter::Scalar
            , rapter::_2d::PrimitiveVectorT
            , rapter::PointContainerT
            , typename rapter::_2d::PatchPatchDistanceFunctorT
            >
            ( rapter::_2d::PrimitiveVectorT                 & patches
            , rapter::PointContainerT                       & points
            , rapter::Scalar                           const  scale
            , std::vector<rapter::Scalar>              const& angles
            , rapter::_2d::PatchPatchDistanceFunctorT  const& patchPatchDistanceFunctor
            , int                                      const  nn_K
            , int                                      const  verbose
            , size_t                                   const  patchPopLimit
            );

    template int
    Segmentation::patchify< rapter::_3d::PrimitiveT
                          , rapter::Scalar
                          , rapter::_3d::PrimitiveVectorT
                          , rapter::PointContainerT
                          , typename rapter::_3d::PatchPatchDistanceFunctorT
                          >
                          ( rapter::_3d::PrimitiveVectorT                 & patches
                          , rapter::PointContainerT                       & points
                          , rapter::Scalar                           const  scale
                          , std::vector<rapter::Scalar>              const& angles
                          , rapter::_3d::PatchPatchDistanceFunctorT  const& patchPatchDistanceFunctor
                          , int                                      const  nn_K
                          , int                                      const  verbose
                          , size_t                                   const  patchPopLimit
                          );

    namespace segm_templinst
    {
        namespace _2d
        {
            typedef segmentation::Patch<rapter::Scalar, rapter::_2d::PrimitiveT> PatchT;
            typedef std::vector< PatchT > PatchesT;
        }
        namespace _3d
        {
            typedef segmentation::Patch<rapter::Scalar, rapter::_3d::PrimitiveT> PatchT;
            typedef std::vector< PatchT > PatchesT;
        }
    }
    template int
    Segmentation::regionGrow  < rapter::_2d::PrimitiveT
                              , rapter::PointContainerT
                              , rapter::_2d::PatchPatchDistanceFunctorT
                              , segm_templinst::_2d::PatchesT
                              , rapter::Scalar
                              , rapter::PointPrimitiveT
                              >
                              ( rapter::PointContainerT             & points
                              , segm_templinst::_2d::PatchesT    & groups_arg
                              , rapter::Scalar                 const  /*scale*/
                              , rapter::_2d::PatchPatchDistanceFunctorT const& patchPatchDistanceFunctor
                              , GidT                        const  gid_tag_name              //= _PointT::GID
                              , int                         const  nn_K
                              , bool                        const  verbose
                              );
    template int
    Segmentation::regionGrow  < rapter::_3d::PrimitiveT
                              , rapter::PointContainerT
                              , rapter::_3d::PatchPatchDistanceFunctorT
                              , segm_templinst::_3d::PatchesT
                              , rapter::Scalar
                              , rapter::PointPrimitiveT
                              >
                              ( rapter::PointContainerT             & points
                              , segm_templinst::_3d::PatchesT    & groups_arg
                              , rapter::Scalar                 const  /*scale*/
                              , rapter::_3d::PatchPatchDistanceFunctorT const& patchPatchDistanceFunctor
                              , GidT                        const  gid_tag_name              //= _PointT::GID
                              , int                         const  nn_K
                              , bool                        const  verbose
                              );

    template int
    Segmentation::fitLocal  ( rapter::_2d::InnerPrimitiveContainerT & lines
                            , pcl::PointCloud<pcl::PointXYZ>::Ptr const  cloud
                            , std::vector<int>     const* indices
                            , int                  const  K
                            , float                const  radius
                            , bool                 const  soft_radius
                            , std::vector<PidT>         * mapping
                            , int                  const  verbose
                            );

    template int
    Segmentation::fitLocal  ( rapter::_3d::InnerPrimitiveContainerT & lines
                            , pcl::PointCloud<pcl::PointXYZ>::Ptr const  cloud
                            , std::vector<int>     const* indices
                            , int                  const  K
                            , float                const  radius
                            , bool                 const  soft_radius
                            , std::vector<PidT>         * mapping
                            , int                  const  verbose
                            );

    template int
    Segmentation::_tagPointsFromGroups< rapter::PointPrimitiveT
                                      , rapter::Scalar
                                      , rapter::_2d::PatchPatchDistanceFunctorT
                                      , segm_templinst::_2d::PatchesT
                                      , rapter::PointContainerT
                                      >
                                      ( rapter::PointContainerT                      & points
                                      , segm_templinst::_2d::PatchesT        const& groups
                                      , rapter::_2d::PatchPatchDistanceFunctorT const& pointPatchDistanceFunctor
                                      , GidT                                 const  gid_tag_name );

    template int
    Segmentation::_tagPointsFromGroups< rapter::PointPrimitiveT
                                      , rapter::Scalar
                                      , rapter::_3d::PatchPatchDistanceFunctorT
                                      , segm_templinst::_3d::PatchesT
                                      , rapter::PointContainerT
                                      >
                                      ( rapter::PointContainerT                      & points
                                      , segm_templinst::_3d::PatchesT        const& groups
                                      , rapter::_3d::PatchPatchDistanceFunctorT const& pointPatchDistanceFunctor
                                      , GidT                                 const  gid_tag_name );
}
