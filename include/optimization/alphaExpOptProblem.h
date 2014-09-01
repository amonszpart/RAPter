#ifndef __GF2_ALPHAEXPOPTPROBLEM_H__
#define __GF2_ALPHAEXPOPTPROBLEM_H__

#include <vector>
#include <deque>
#include "Eigen/Dense"

#include "GCoptimization.h"
#include "pcltools/util.hpp"

#include "my_types.h"
#include "optimization/optProblem.h"
#include "primitives/linePrimitive.h"

namespace am
{
    class AlphaExpOptProblem : public OptProblem
    {
        public:
            typedef float Scalar;

            AlphaExpOptProblem( std::vector<LinePrimitive>         const& lines
                              , pcl::PointCloud<MyPoint>::ConstPtr        cloud
                              , double                             const  threshold
                              , std::vector<int>                   const* indices
                              , Eigen::VectorXf                    const& lambdas
                              , int                                       max_step_count = 65536
                              , float                              const  trunc_pw_at_angle = 0.4
                              );

            virtual Scalar              getEnergy();
            virtual int                 next();
            virtual bool                hasNext();
            virtual int                 stepBack();
            virtual int                 init( std::vector<int> const& config );
            virtual std::vector<int>    getConfig();
            virtual ~AlphaExpOptProblem();

            template <class PrimitivesT, class PointsPtrT>
            static inline int solve( MaskType                            & opt_mask
                                     , PrimitivesT                  const& primitives
                                     , PointsPtrT                          cloud
                                     , Scalar                              working_scale
                                     , Eigen::Matrix<Scalar,-1,1>          lambdas
                                     , Scalar                       const  trunc_at
                                     , std::vector<Scalar>          const& desired_angles
                                     , std::vector<int>                  * labels         = NULL
                                     , Scalar                       const  angle_thresh   = 0.2  );

        protected:
            pcl::PointCloud<MyPoint>::ConstPtr  _cloud;
            std::vector<int>                    _indices;
            std::vector<LinePrimitive>          _lines;
            Eigen::VectorXf                     _lambdas;
//            std::deque<std::pair<std::vector<int>,Scalar> >  _configs_energies;
//            const size_t                        _ids_energies_capacity;

//            const double                        _threshold;

            bool                                _energy_up_to_date;
//            int                                 _step_count, _max_step_count;
            GCoptimizationGeneralGraph          *_gc;


    };
} // ns am

#ifndef __GF2_ALPHAEXPOPTPROBLEM_INC_HPP__
#   define __GF2_ALPHAEXPOPTPROBLEM_INC_HPP__
#   include "alphaExpOptProblem.hpp"
#endif

#endif // __GF2_ALPHAEXPOPTPROBLEM_H__
