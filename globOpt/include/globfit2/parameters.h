#ifndef GF2_PARAMETERS_H
#define GF2_PARAMETERS_H

#include <vector>
#include "Eigen/Dense"

namespace GF2 {

    //! \brief Collection of parameters for \ref CandidateGenerator::generate.
    template <typename _Scalar>
    struct CandidateGeneratorParams
    {
            enum RefitMode { SPATIAL    //!< \deprecated Spatial refit to the member points does not make sense for oriented points.
                           , AVG_DIR    //!< \brief Patch representative direction will be the average position and average direction of member oriented points.
                           };
            enum PatchPatchDistMode {
#if GF2_WITH_FULL_LINKAGE
                                      FULL_MIN      //!< \brief \ref FullLinkagePatchPatchDistanceFunctorT, \ref SpatialPatchPatchMinDistanceFunctorT
                                    , FULL_MAX      //!< \brief \ref FullLinkagePatchPatchDistanceFunctorT, \ref SpatialPatchPatchMaxDistanceFunctorT
                                    , SQR_MIN       //!< \brief \ref SquaredPatchPatchDistanceFunctorT
                                    , REPR_MIN,      //!< \brief \ref RepresentativePatchPatchDistanceFunctorT
#endif // GF2_WITH_FULL_LINKAGE
                                      REPR_SQR      //!< \brief \ref RepresentativeSqrPatchPatchDistanceFunctorT
                                    };

            _Scalar   scale                       = 0.05;     //!< \brief Scale parameter of input.
            std::vector<_Scalar> angles           = {0, M_PI_2, M_PI }; //!< \brief Desired angles.
            _Scalar   angle_limit                 = 0.08f;    //!< \brief angle threshold for similar lines
            _Scalar   parallel_limit              = _Scalar(1e-6); //!< \brief Two lines are parallel, if their angle is smaller than this. Used in \ref Merging.
            _Scalar   angle_limit_div             = 10.f;     //!< \brief Determines, how much the angle_limit is divided by, to get the patch-similarity threshold
            _Scalar   patch_dist_limit_mult       = 1.f;      //!< \brief Patchify takes "patch_dist_limit * scale" as maximum spatial distance
            _Scalar   patch_spatial_weight        = _Scalar(0.5); //!< \brief Weight of spatial term in \f$ patch\_spatial\_weight^2 \cdot \frac{spat\_dist^2}{spat\_thresh} + \frac{ang\_diff^2}{ang\_thresh} < 1 \f$ regionGrowing distance functor. \sa \ref GF2::RepresentativeSqrPatchPatchDistanceFunctorT.
            int       nn_K                        = 15;       //!< \brief Number of points used for primitive fitting
            int       patch_population_limit      = 10;       //!< \brief Threshold, on when a patch can distribute its directions
            RefitMode refit_mode                  = AVG_DIR;  //!< \brief Determines, how a patch gets it's final direction. A NL-LSQ to the points, or the average of local orientations.

            //! \brief Regiongrow uses these options to determine which functor to use for a patch-patch distance.
            //! \sa \ref FullLinkagePatchPatchDistanceFunctorT, \ref SquaredPatchPatchDistanceFunctorT, \ref RepresentativePatchPatchDistanceFunctorT,
            //! \sa \ref RepresentativeSqrPatchPatchDistanceFunctorT, \ref SpatialPatchPatchMinDistanceFunctorT, \ref SpatialPatchPatchMaxDistanceFunctorT
            PatchPatchDistMode patch_dist_mode    = REPR_SQR;

            //___Parsers___

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

            inline int parsePatchDistMode( std::string const& patch_dist_mode_string )
            {
                int err = EXIT_SUCCESS;
#if GF2_WITH_FULL_LINKAGE
                if ( !patch_dist_mode_string.compare("full_min") )
                {
                    this->patch_dist_mode = FULL_MIN;
                }
                else if ( !patch_dist_mode_string.compare("full_max") )
                {
                    this->patch_dist_mode = FULL_MAX;
                }
                else if ( !patch_dist_mode_string.compare("squared_min") )
                {
                    this->patch_dist_mode = SQR_MIN;
                }
                else if ( !patch_dist_mode_string.compare("representative_min") )
                {
                    this->patch_dist_mode = REPR_MIN;
                }
                else
#endif // GF2_WITH_FULL_LINKAGE
                    if ( !patch_dist_mode_string.compare("representative_sqr") )
                {
                    this->patch_dist_mode = REPR_SQR;
                }
                else
                {
                    err = EXIT_FAILURE;
                    std::cerr << "[" << __func__ << "]: " << "Could NOT parse " << patch_dist_mode_string << ", assuming " << printPatchDistMode() << std::endl;
                }

                return err;
            } // ...parsePatchDistMode()

            inline std::string printPatchDistMode() const
            {
                switch ( patch_dist_mode )
                {
#if GF2_WITH_FULL_LINKAGE
                    case FULL_MIN: return std::string("full_min"); break;
                    case FULL_MAX: return std::string("full_max"); break;
                    case SQR_MIN:  return std::string("squared_min"); break;
                    case REPR_MIN: return std::string("representative_min"); break;
#endif
                    case REPR_SQR: return std::string("representative_sqr"); break;
                    default:       return "UNKNOWN"; break;
                }
            } // ...printPatchDistMode()

    }; // ...struct CandidateGeneratorParams

    //! \brief Collection of parameters the \ref ProblemSetup::formulate needs.
    template <typename _Scalar>
    struct ProblemSetupParams
    {
        //! \brief Options to calculate data cost of primitive by.
        enum DATA_COST_MODE { ASSOC_BASED    = 0 //!< \brief \sa \ref problemSetup::associationBasedDataCost
                            , BAND_BASED     = 1 //!< \brief \sa \ref problemSetup::bandBasedDataCost
                            , INSTANCE_BASED = 2 //!< \brief \sa \ref problemSetup::instanceBasedDataCost
                            };
        //! \brief Options to calculate constraints for patches in optimization.
        enum CONSTR_MODE    { PATCH_WISE     = 0 //!< \brief \sa \ref problemSetup::everyPatchNeedsDirectionConstraint
                            , POINT_WISE     = 1 //!< \brief \sa \ref problemSetup::everyPointNeedsPatchConstraint
                            , HYBRID         = 2 //!< \brief \sa \ref problemSetup::largePatchesNeedDirectionConstraint
                            };

        //! \brief Scale parameter of input.
        _Scalar                      scale           = 0.05;
        //! \brief Desired angles.
        std::vector<_Scalar>         angles          = {0, M_PI_2, M_PI };
        //! \brief Optimization weights. [ 0: unary, 1: pairwise, 2: complexity ]
        Eigen::Matrix<_Scalar,-1,1>  weights         = (Eigen::Matrix<_Scalar,-1,1>(3,1) << 10000, 10, 1).finished();
        //! \brief \copydoc CONSTR_MODE
        CONSTR_MODE    constr_mode     = PATCH_WISE;
        //! \brief \copydoc DATA_COST_MODE
        DATA_COST_MODE data_cost_mode  = ASSOC_BASED;
        //! \brief Point-count threshold, that decides if a patch is big or small.
        int                          patch_pop_limit = 5;
    };

    template <typename _Scalar>
    struct MergeParams : public CandidateGeneratorParams<_Scalar>
    {
        char do_adopt = 2; //!< \brief Adopt orphaned points greedily, 1: unambiguous only, 2: all.
//        //! \brief Scale parameter of input.
//        _Scalar                      scale           = 0.05;
//        //! \brief Desired angles.
//        std::vector<_Scalar>         angles          = {0, M_PI_2, M_PI };
    };

}

#endif // GF2_PARAMETERS_H
