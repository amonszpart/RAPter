#include "my_types.h"
#include <vector>
#include <stddef.h>
#include "optimization/simAnnOptProblem.h"

namespace am
{
#if 0
    std::vector<MaskType>
    simAnnOptProblem<powerSet( int K )
    {
        if ( K == 1 )
            return {{0},{1}};

        std::vector< std::vector<int> > prev = powerSet( K-1 );

        std::vector< std::vector<int> > out;
        out.resize( prev.size() << 1 );

        int offset = 0;
        for ( int prefix = 0; prefix != 2; ++prefix )
        {
            for ( size_t i = 0; i != prev.size(); ++i )
            {
                out[offset+i].resize( prev[i].size() + 1 );
                out[offset+i][0] = prefix;
                std::copy( prev[i].begin(), prev[i].end(), out[offset+i].begin()+1 );
            }
            offset += prev.size();
        }

        return out;
    } // func powerSet
#endif
} // ns am
