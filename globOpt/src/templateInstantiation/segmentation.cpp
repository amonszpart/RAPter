#include "globfit2/optimization/impl/segmentation.hpp"
#include "globfit2/globOpt_types.h"

namespace GF2
{
    template int
    Segmentation::segmentCli< GF2::_3d::PrimitiveT
                            , GF2::_3d::PrimitiveVectorT
                            , GF2::PointPrimitiveT
                            , GF2::PointContainerT
                            , GF2::Scalar
                            >
                            ( int argc, char** argv );

    template int
    Segmentation::segmentCli< GF2::_2d::PrimitiveT
                            , GF2::_2d::PrimitiveVectorT
                            , GF2::PointPrimitiveT
                            , GF2::PointContainerT
                            , GF2::Scalar
                            >
                            ( int argc, char** argv );

    template int
    Segmentation::orientPoints< GF2::PointPrimitiveT
                              , GF2::_3d::PrimitiveT
                              , GF2::Scalar
                              , GF2::PointContainerT
                              >
                              ( GF2::PointContainerT       &points
                              , GF2::Scalar          const  scale
                              , int                  const  nn_K
                              , int                  const  verbose );

    template int
    Segmentation::orientPoints< GF2::PointPrimitiveT
                              , GF2::_2d::PrimitiveT
                              , GF2::Scalar
                              , GF2::PointContainerT
                              >
                              ( GF2::PointContainerT       &points
                              , GF2::Scalar          const  scale
                              , int                  const  nn_K
                              , int                  const  verbose );


    template int
    Segmentation::patchify< GF2::_2d::PrimitiveT
            , GF2::Scalar
            , GF2::_2d::PrimitiveVectorT
            , GF2::PointContainerT
            , typename GF2::_2d::PatchPatchDistanceFunctorT
            >
            ( GF2::_2d::PrimitiveVectorT                 & patches
            , GF2::PointContainerT                       & points
            , GF2::Scalar                           const  scale
            , std::vector<GF2::Scalar>              const& angles
            , GF2::_2d::PatchPatchDistanceFunctorT  const& patchPatchDistanceFunctor
            , int                                   const  nn_K
            , int                                   const  verbose
            , int                                   const  patchPopLimit
            );

    template int
    Segmentation::patchify< GF2::_3d::PrimitiveT
                          , GF2::Scalar
                          , GF2::_3d::PrimitiveVectorT
                          , GF2::PointContainerT
                          , typename GF2::_3d::PatchPatchDistanceFunctorT
                          >
                          ( GF2::_3d::PrimitiveVectorT                 & patches
                          , GF2::PointContainerT                       & points
                          , GF2::Scalar                           const  scale
                          , std::vector<GF2::Scalar>              const& angles
                          , GF2::_3d::PatchPatchDistanceFunctorT  const& patchPatchDistanceFunctor
                          , int                                   const  nn_K
                          , int                                   const  verbose
                          , int                                   const  patchPopLimit
                          );

    namespace segm_templinst
    {
        namespace _2d
        {
            typedef segmentation::Patch<GF2::Scalar,GF2::_2d::PrimitiveT> PatchT;
            typedef std::vector< PatchT > PatchesT;
        }
        namespace _3d
        {
            typedef segmentation::Patch<GF2::Scalar,GF2::_3d::PrimitiveT> PatchT;
            typedef std::vector< PatchT > PatchesT;
        }
    }
    template int
    Segmentation::regionGrow  < GF2::_2d::PrimitiveT
                              , GF2::PointContainerT
                              , GF2::_2d::PatchPatchDistanceFunctorT
                              , segm_templinst::_2d::PatchesT
                              , GF2::Scalar
                              , GF2::PointPrimitiveT
                              >
                              ( GF2::PointContainerT             & points
                              , segm_templinst::_2d::PatchesT    & groups_arg
                              , GF2::Scalar                 const  /*scale*/
                              , GF2::_2d::PatchPatchDistanceFunctorT const& patchPatchDistanceFunctor
                              , GidT                        const  gid_tag_name              //= _PointT::GID
                              , int                         const  nn_K
                              , bool                        const  verbose
                              );
    template int
    Segmentation::regionGrow  < GF2::_3d::PrimitiveT
                              , GF2::PointContainerT
                              , GF2::_3d::PatchPatchDistanceFunctorT
                              , segm_templinst::_3d::PatchesT
                              , GF2::Scalar
                              , GF2::PointPrimitiveT
                              >
                              ( GF2::PointContainerT             & points
                              , segm_templinst::_3d::PatchesT    & groups_arg
                              , GF2::Scalar                 const  /*scale*/
                              , GF2::_3d::PatchPatchDistanceFunctorT const& patchPatchDistanceFunctor
                              , GidT                        const  gid_tag_name              //= _PointT::GID
                              , int                         const  nn_K
                              , bool                        const  verbose
                              );

    template int
    Segmentation::fitLocal  ( GF2::_2d::InnerPrimitiveContainerT & lines
                            , pcl::PointCloud<pcl::PointXYZ>::Ptr const  cloud
                            , std::vector<int>     const* indices
                            , int                  const  K
                            , float                const  radius
                            , bool                 const  soft_radius
                            , std::vector<PidT>         * mapping
                            , int                  const  verbose
                            );

    template int
    Segmentation::fitLocal  ( GF2::_3d::InnerPrimitiveContainerT & lines
                            , pcl::PointCloud<pcl::PointXYZ>::Ptr const  cloud
                            , std::vector<int>     const* indices
                            , int                  const  K
                            , float                const  radius
                            , bool                 const  soft_radius
                            , std::vector<PidT>         * mapping
                            , int                  const  verbose
                            );

    template int
    Segmentation::_tagPointsFromGroups< GF2::PointPrimitiveT
                                      , GF2::Scalar
                                      , GF2::_2d::PatchPatchDistanceFunctorT
                                      , segm_templinst::_2d::PatchesT
                                      , GF2::PointContainerT
                                      >
                                      ( GF2::PointContainerT                      & points
                                      , segm_templinst::_2d::PatchesT        const& groups
                                      , GF2::_2d::PatchPatchDistanceFunctorT const& pointPatchDistanceFunctor
                                      , GidT                                 const  gid_tag_name );

    template int
    Segmentation::_tagPointsFromGroups< GF2::PointPrimitiveT
                                      , GF2::Scalar
                                      , GF2::_3d::PatchPatchDistanceFunctorT
                                      , segm_templinst::_3d::PatchesT
                                      , GF2::PointContainerT
                                      >
                                      ( GF2::PointContainerT                      & points
                                      , segm_templinst::_3d::PatchesT        const& groups
                                      , GF2::_3d::PatchPatchDistanceFunctorT const& pointPatchDistanceFunctor
                                      , GidT                                 const  gid_tag_name );
}
