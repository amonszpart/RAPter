#ifndef GF2_MERGING_HPP
#define GF2_MERGING_HPP

#include "globfit2/optimization/merging.h"

#include "globfit2/parameters.h"
#include "globfit2/visualization/visualization.h"
#include "globfit2/io/io.h"
//_______________________HPP_____________________________

namespace GF2 {

template < class    _PrimitiveContainerT
         , class    _PointContainerT
         , typename _Scalar
         , class    _PointPrimitiveT
         , class    _PrimitiveT
         >
int
Merging::mergeCli( int argc, char** argv )
{
    MergeParams<_Scalar> params;

    std::string cloud_path = "cloud.ply",
                prims_path = "primitives.bonmin.txt";
    _Scalar     angle_gen  = M_PI_2;
    // parse params
    {
        bool valid_input = true;

        valid_input &= pcl::console::parse_argument( argc, argv, "--scale", params.scale ) >= 0;
        valid_input &= pcl::console::parse_argument( argc, argv, "--prims", prims_path   ) >= 0;

        // cloud
        pcl::console::parse_argument( argc, argv, "--cloud", cloud_path );
        valid_input &= boost::filesystem::exists( cloud_path );

        pcl::console::parse_argument( argc, argv, "--angle-gen", angle_gen );

        std::cerr << "[" << __func__ << "]: " << "Usage:\t gurobi_opt --formulate\n"
                  << "\t--scale " << params.scale << "\n"
                  << "\t--prims " << prims_path << "\n"
                  << "\t--cloud " << cloud_path << "\n"
                  << "\t[--angle-gen " << angle_gen << "]\n"
                  << std::endl;

        if ( !valid_input || pcl::console::find_switch(argc,argv,"--help") || pcl::console::find_switch(argc,argv,"-h") )
        {
            std::cerr << "[" << __func__ << "]: " << "--scale, --prims are compulsory, --cloud needs to exist" << std::endl;
            return EXIT_FAILURE;
        }
    } // ... parse params

    // read primitives
    _PrimitiveContainerT prims;
    typedef std::vector<_PrimitiveT> PatchT;
    std::map< int, std::vector<_PrimitiveT> > prims_map;
    {
        io::readPrimitives<_PrimitiveT, PatchT>( prims, prims_path, &prims_map );
    } //...read primitives

    // Read points
    _PointContainerT points;
    {
        io::readPoints<_PointPrimitiveT>( points, cloud_path );
    }

    // Read desired angles
    params.angles = { _Scalar(0) };
    {
        for ( _Scalar angle = angle_gen; angle < M_PI; angle+= angle_gen )
            params.angles.push_back( angle );
        params.angles.push_back( M_PI );

        // log
        std::cout << "[" << __func__ << "]: " << "Desired angles: {";
        for ( size_t vi=0;vi!=params.angles.size();++vi)
            std::cout << params.angles[vi] << ((vi==params.angles.size()-1) ? "" : ", ");
        std::cout << "}\n";
    } // ... read angles

    // associations
    std::vector<std::pair<int,int> > points_primitives;
    io::readAssociations( points_primitives, "points_primitives.txt", NULL );
    for ( size_t i = 0; i != points.size(); ++i )
    {
        // error check
        if ( i > points_primitives.size() )
        {
            std::cerr << "more points than associations..." << std::endl;
            return EXIT_FAILURE;
        }

        // store association in point
        points[i].setTag( _PointPrimitiveT::GID, points_primitives[i].first );

        // error check 2
        if ( points[i].getTag(_PointPrimitiveT::GID) == -1 )
            std::cerr << "[" << __func__ << "]: " << "point assigned to patch with id -1" << std::endl;
    }

    // select unassigned points
    for ( size_t i = 0; i != points.size(); ++i )
    {
        int gid = points[i].getTag( _PointPrimitiveT::GID );
        if ( !prims_map[gid].size() )
        {
            points[i].setTag( _PointPrimitiveT::GID, -2 );
        }
    }


    io::writeAssociations<_PointPrimitiveT>( points, "test_assoc.csv" );

    return EXIT_SUCCESS;
}//...Merging::merge()

} //...namespace GF2

#endif // GF2_MERGING_HPP
