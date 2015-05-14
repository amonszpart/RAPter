/**
  * inputParser.hpp
  * : Class name (if applicable)
  * Author: Aron Monszpart <a.monszpart@cs.ucl.ac.uk>
  * Created: 10/11/2014
  */

#ifndef GF2_INPUTPARSER_HPP
#define GF2_INPUTPARSER_HPP

#include <string>
#include "boost/filesystem.hpp"
#include "globfit2/simple_types.h"
#include "globfit2/util/parse.h"
#include "globfit2/io/io.h"
#include "globfit2/processing/util.hpp"
#include "globfit2/processing/angle_util.hpp" // appendAngles

namespace GF2
{

inline std::string parseAssocPath( int argc, char** argv )
{
    std::string input_prims_path;
    if ( GF2::console::parse_argument( argc, argv, "-a"     , input_prims_path) < 0 )
        GF2::console::parse_argument( argc, argv, "--assoc", input_prims_path);

    return input_prims_path;
} //...parsePrimitivesPath()

inline std::string parsePrimitivesPath( int argc, char** argv )
{
    std::string input_prims_path;
    if ( GF2::console::parse_argument( argc, argv, "-p"     , input_prims_path) < 0 )
        GF2::console::parse_argument( argc, argv, "--prims", input_prims_path);

    return input_prims_path;
} //...parsePrimitivesPath()

template < class _InnerPrimitiveContainerT
         , class _PclCloudT
         , class _PointContainerT
         , class _PrimitiveContainerT
         , class _PrimitiveMapT
         , typename _ParamsT >
inline int parseInput( _PointContainerT         &points
                     , typename _PclCloudT::Ptr &pcl_cloud
                     , _PrimitiveContainerT     &initial_primitives
                     , _PrimitiveMapT           &patches
                     , _ParamsT                 &params
                     , int                       argc, char **argv
                     , bool                      read_assoc = true )
{
    typedef typename _PointContainerT::value_type           PointPrimitiveT;
    typedef typename _InnerPrimitiveContainerT::value_type  PrimitiveT;
    typedef typename PrimitiveT::Scalar                     Scalar;

    bool valid_input = true;
    std::string cloud_path("cloud.ply"), input_prims_path, associations_path;

    if (    (GF2::console::parse_argument( argc, argv, "--cloud", cloud_path) < 0)
         && (!boost::filesystem::exists(cloud_path)) )
    {
        std::cerr << "[" << __func__ << "]: " << "--cloud does not exist: " << cloud_path << std::endl;
        valid_input = false;
    }

    // primitives
    input_prims_path = parsePrimitivesPath(argc, argv);
    if (    /*(GF2::console::parse_argument( argc, argv, "-p"     , input_prims_path) < 0)
         && (GF2::console::parse_argument( argc, argv, "--prims", input_prims_path) < 0)*/
         !input_prims_path.size()
         || (!boost::filesystem::exists(input_prims_path)) )
    {
        std::cerr << "[" << __func__ << "]: " << "-p or --prims is compulsory" << std::endl;
        valid_input = false;
    }

    associations_path = parseAssocPath(argc,argv);
    if (    /*(pcl::console::parse_argument( argc, argv, "-a", associations_path) < 0)
         && (pcl::console::parse_argument( argc, argv, "--assoc", associations_path) < 0)*/
            (!boost::filesystem::exists(associations_path)) )
    {
        if ( read_assoc )
        {
            std::cerr << "[" << __func__ << "]: " << "-a or --assoc is compulsory" << std::endl;
            valid_input = false;
        }
    }

    // scale
    if (    (GF2::console::parse_argument( argc, argv, "--scale", params.scale) < 0)
         && (GF2::console::parse_argument( argc, argv, "-sc"    , params.scale) < 0) )
    {
        std::cerr << "[" << __func__ << "]: " << "--scale is compulsory" << std::endl;
        valid_input = false;
    }

    // READ
    int err     = EXIT_SUCCESS;
    pcl_cloud   = typename _PclCloudT::Ptr( new _PclCloudT() );
    if ( EXIT_SUCCESS == err )
    {
        err = GF2::io::readPoints<PointPrimitiveT>( points, cloud_path, &pcl_cloud );
        if ( err != EXIT_SUCCESS )
        {
            std::cerr << "[" << __func__ << "]: " << "readPoints returned error " << err << std::endl;
            valid_input = false;
        }
    } //...read points

    if ( EXIT_SUCCESS == err )
    {
        std::cout << "[" << __func__ << "]: " << "reading primitives from " << input_prims_path << "...";
        err = GF2::io::readPrimitives<PrimitiveT, _InnerPrimitiveContainerT>( initial_primitives, input_prims_path, &patches );
        if ( EXIT_SUCCESS == err )
            std::cout << "[" << __func__ << "]: " << "reading primitives ok (#: " << initial_primitives.size() << ")\n";
        else
            std::cerr << "[" << __func__ << "]: " << "reading primitives failed" << std::endl;
    } //...read primitives

    // read assoc
    if ( read_assoc )
    {
        std::vector< std::pair<PidT,LidT> > points_primitives;
        GF2::io::readAssociations( points_primitives, associations_path, NULL );
        for ( size_t i = 0; i != points.size(); ++i )
        {
            // store association in point
            points[i].setTag( PointPrimitiveT::TAGS::GID, points_primitives[i].first );
        }
    }

    return err + !valid_input;
}

inline int parseAngles( AnglesT &angles, int argc, char** argv, AnglesT* angle_gens_out = NULL, bool no_paral = false, bool verbose = false, bool inRad = false )
{
    AnglesT angle_gens;
    if ( pcl::console::parse_x_arguments( argc, argv, "--angle-gens", angle_gens ) < 0 )
        return EXIT_FAILURE;

    angles::appendAnglesFromGenerators( angles, angle_gens, no_paral, verbose, inRad );

    if ( angle_gens_out ) *angle_gens_out = angle_gens;

    return EXIT_SUCCESS;
}

} //...ns GF2


#endif // INPUTPARSER_HPP
