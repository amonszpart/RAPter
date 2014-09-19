#ifndef GF2_DEPTHIO_HPP
#define GF2_DEPTHIO_HPP

namespace GF2 {
namespace io {

/*! \brief Load a depth image stored as a .dat file.
 *          Assumes, that first two int entries are height and width (in this order).
 *
// Tout is the format of the data file, Timg is the format of cv::Mat
*/
template <typename Timg, typename Tout>
int
loadDepthAsBinary( ::cv::Mat &img, std::string path, float alpha = 1.f )
{
    FILE *fp = fopen( path.c_str(), "rb" );
    if (!fp)
    {
        printf( "loadDepthAsBinary failed, when opening %s\n", path.c_str() );
        return false;
    }

    int w,h, read_size = 0;
    read_size += fread( &h,sizeof(int),1,fp );
    read_size += fread( &w,sizeof(int),1,fp );
    if ( read_size != 2 )
    {
        std::cerr << "[" << __func__ << "]: " << "w or h not read properly: " << w << ", " << h << std::endl;
    }

    Tout *p_depth = new Tout[ h * w ];

    read_size = fread( p_depth, sizeof(Tout), w*h, fp );
    if ( read_size != w*h )
        std::cerr << "[" << __func__ << "]: " << "read size != " << w*h << std::endl;

    img.create( h, w, CV_16UC1 );
    int p_depth_index = 0;
    for ( unsigned y = 0; y < img.rows; ++y  )
    {
        for ( unsigned x = 0; x < img.cols; ++x, ++p_depth_index )
        {
            img.at<Timg>(y,x) = p_depth[ p_depth_index ] * alpha;
        }
    }

    fclose( fp );
    if ( p_depth )
        delete [] p_depth;

    return EXIT_SUCCESS;
}

} //...ns io
} //...ns GF2

#endif // GF2_DEPTHIO_HPP
