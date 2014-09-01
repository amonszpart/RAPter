#ifndef __GF2_OPTPROBLEM_H__
#define __GF2_OPTPROBLEM_H__

#include <vector>
#include "my_types.h"

namespace am
{

#ifndef __GF2_POWERSET__
#   define __GF2_POWERSET__
    extern std::vector<MaskType> powerSet( int );
#endif

    class OptProblem
    {
        public:
            typedef float Scalar;

            virtual Scalar getEnergy( void const* p_desired_angles = NULL ) = 0;
            virtual int    next     ( void const* p_desired_angles = NULL ) = 0;
            virtual bool   hasNext()   = 0;
            virtual int    stepBack()  = 0;
            virtual int    init( std::vector<int> const& config ) = 0;
            virtual std::vector<int> getConfig() = 0;
    };

} // ns am

#endif // __GF2_OPTPROBLEM_H__
