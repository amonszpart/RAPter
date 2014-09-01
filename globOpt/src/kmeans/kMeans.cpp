#include <stdlib.h>
namespace am
{
    namespace kmeans
    {

        double randf( double m )
        {
            return m * rand() / (RAND_MAX - 1.);
        }

    } // ns kmeans
} //ns am
