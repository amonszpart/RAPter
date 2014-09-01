#define NANOSVG_IMPLEMENTATION
#include "/home/bontius/workspace/3rdparty/nanosvg/src/nanosvg.h"
#include "stdio.h"
#include <iostream>

int main(int argc, char *argv[])
{
    if ( argc < 2 )
    {
        std::cout << "usage:\t" << "x.svg" << std::endl;
        return 1;
    }

    // Load
    struct NSVGimage* image = NULL;
    std::cout << "parsing " << argv[1] << std::endl;
    image = nsvgParseFromFile( argv[1], "px", 96);
    if ( !image )
    {
        std::cout << "image == NULL...exiting" << std::endl;
        return 1;
    }
    printf( "size: %f x %f\n", image->width, image->height );
    // Use...
    for ( NSVGshape* shape = image->shapes; shape != NULL; shape = shape->next)
    {
        std::cout << "shape\n";
        for ( NSVGpath* path = shape->paths; path != NULL; path = path->next)
        {
            for ( int i = 0; i < path->npts-1; i += 3)
            {
                std::cout << "path->pts[2*i]: " << path->pts[2*i] << std::endl;
                std::cout << "path->pts[2*i+1]: " << path->pts[2*i+1] << std::endl;
//                float* p = &path->pts[i*2];
//                //drawCubicBez(p[0],p[1], p[2],p[3], p[4],p[5], p[6],p[7]);
//                std::cout << "p: ";
//                for ( int j = 0; j != 8; ++j )
//                    std::cout << "p[" << j << "]: " << p[j] << "\t";
//                std::cout << std::endl;
            }
        }
    }
    // Delete
    if ( image ) nsvgDelete(image);
    std::cout << "finished" << std::endl;

    return 0;
}
