#ifndef GF2_PARAMETERS_H
#define GF2_PARAMETERS_H

#include <vector>
#include "Eigen/Dense"

namespace GF2 {

    //! \brief Collection of parameters for \ref CandidateGenerator::generate.
    template <typename _Scalar>
    struct CandidateGeneratorParams
    {
            enum RefitMode { SPATIAL, AVG_DIR };
            enum PatchPatchDistMode { FULL_MIN, FULL_MAX, SQR_MIN, REPR_MIN };

            _Scalar   scale                       = 0.05;     //!< \brief Scale parameter of input.
            std::vector<_Scalar> angles           = {0, M_PI_2, M_PI }; //!< \brief Desired angles.
            _Scalar   angle_limit                 = 0.08f;    //!< \brief angle threshold for similar lines
            _Scalar   angle_limit_div             = 10.f;     //!< \brief Determines, how much the angle_limit is divided by, to get the patch-similarity threshold
            _Scalar   patch_dist_limit            = 1.f;      //!< \brief Patchify takes "patch_dist_limit * scale" as maximum spatial distance
            int       nn_K                        = 15;       //!< \brief Number of points used for primitive fitting
            int       patch_population_limit      = 10;       //!< \brief Threshold, on when a patch can distribute its directions
            RefitMode refit_mode                  = AVG_DIR;  //!< \brief Determines, how a patch gets it's final direction. A NL-LSQ to the points, or the average of local orientations.

            //! \brief Regiongrow uses these options to determine which functor to use for a patch-patch distance.
            //! \sa \ref FullLinkagePatchPatchDistanceFunctorT, \ref SquaredPatchPatchDistanceFunctorT, \ref RepresentativePatchPatchDistanceFunctorT, \ref SpatialPatchPatchMinDistanceFunctorT, \ref SpatialPatchPatchMaxDistanceFunctorT
            PatchPatchDistMode patch_dist_mode    = REPR_MIN;

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
                    case FULL_MIN: return std::string("full_min"); break;
                    case FULL_MAX: return std::string("full_max"); break;
                    case SQR_MIN:  return std::string("squared_min"); break;
                    case REPR_MIN: return std::string("representative_min"); break;
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
    struct MergeParams
    {
        //! \brief Scale parameter of input.
        _Scalar                      scale           = 0.05;
        //! \brief Desired angles.
        std::vector<_Scalar>         angles          = {0, M_PI_2, M_PI };
    };

}

#endif // GF2_PARAMETERS_H
