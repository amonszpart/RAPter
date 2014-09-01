#ifndef GF2_MERGING_HPP
#define GF2_MERGING_HPP

#include <iostream>
#include "globfit2/visualization/visualizer.h"

namespace GF2 {

class Merging
{
    public:

        template < class    _PrimitiveContainerT
                 , class    _PointContainerT
                 , typename _Scalar              = float
                 , class    _PointPrimitiveT = typename _PointContainerT::value_type >
        static inline int merge( int argc, char** argv );

};

} //... namespace GF2

//____________________________________________________

namespace GF2 {

template < class    _PrimitiveContainerT
         , class    _PointContainerT
         , typename _Scalar
         , class    _PointPrimitiveT
         >
int
Merging::merge( int argc, char** argv )
{
    std::string cloud_path = "cloud.ply",
                prims_path = "primitives.bonmin.txt";
    _Scalar     scale      = 0.05;
    _Scalar     angle_gen  = M_PI_2;
    // parse params
    {
        bool valid_input = true;
        valid_input &= pcl::console::parse_argument( argc, argv, "--scale", scale          ) >= 0;
        valid_input &= pcl::console::parse_argument( argc, argv, "--cloud", cloud_path     ) >= 0;
        valid_input &= pcl::console::parse_argument( argc, argv, "--prims", prims_path     ) >= 0;
        pcl::console::parse_argument( argc, argv, "--angle-gen", angle_gen );

        std::cerr << "[" << __func__ << "]: " << "Usage:\t gurobi_opt --formulate\n"
                  << " --scale " << scale << "\n"
                  << " --cloud " << cloud_path << "\n"
                  << " --prims " << prims_path << "\n"
                  << std::endl;

        if ( !valid_input || pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
        {
            std::cerr << "[" << __func__ << "]: " << "--scale, --prims, --cloud are compulsory" << std::endl;
            return EXIT_FAILURE;
        }
    } // ... parse params

    // read primitives
    _PrimitiveContainerT prims;
    {
        io::readPrimitives( prims, prims_path );
    } //...read primitives

    // Read points
    _PointContainerT points;
    {
        io::readPoints( points, cloud_path );
    }

    // Read desired angles
    std::vector<_Scalar> angles = { _Scalar(0) };
    {
        for ( _Scalar angle = angle_gen; angle < M_PI; angle+= angle_gen )
            angles.push_back( angle );
        angles.push_back( M_PI );

        std::cout << "Desired angles: {";
        for ( size_t vi=0;vi!=angles.size();++vi)
            std::cout << angles[vi] << ((vi==angles.size()-1) ? "" : ", ");
        std::cout << "}\n";
    } // ... read angles

    std::vector<std::pair<int,int> > points_primitives;
    io::readAssociations( points_primitives, "points_primitives.txt", NULL );
    for ( size_t i = 0; i != points.size(); ++i )
    {
        if ( i > points_primitives.size() )
        {
            std::cerr << "more points than associations..." << std::endl;
            return EXIT_FAILURE;
        }
        points[i].setTag( _PointPrimitiveT::GID, points_primitives[i].first );

        if ( (points[i].getTag( _PointPrimitiveT::GID ) >= static_cast<int>(prims.size())) || (points[i].getTag( _PointPrimitiveT::GID ) < 0) )
            std::cerr << "points[" << i << "].getTag(GID) >= prims.size() || < 0: " << points[i].getTag( _PointPrimitiveT::GID ) << " >= " << prims.size() << std::endl;
    }

    auto vptr = GF2::Visualizer<_PrimitiveContainerT,_PointContainerT>::template show<_Scalar>( /* primitives: */ prims
                                                                                  , /*     points: */  points
                                                                                  , /*      scale: */ scale
                                                                                  , /*     colour: */ {0,0,1}
                                                                                  , /*       spin: */ false
                                                                                  , /*     angles: */ &angles
                                                                                  , /*   show_ids: */ false
                                                                                  , /*   use_tags: */ true );
    // todo: if close, add edge

    return EXIT_SUCCESS;
}//...Merging::merge()

} //...namespace GF2

#endif // GF2_MERGING_HPP
