#include <iostream>
#include <vector>
#include <string>

#include <my_types.h>
#include "primitives/linePrimitive.h"
#include "primitives/planePrimitive.h"
#include "pcl/io/pcd_io.h"
#include "pcl/io/ply_io.h"
#include "pcl/console/parse.h"

#include "post_opt/postAligner.h"
#include "localAnalysis.h"
#include "globFit2.h"
#include "optimization/alphaExpOptProblem.h"


Eigen::Matrix<float,3,1> getLine2DNormal( am::LinePrimitive const& line )
{
    return line().template segment<3>(3).cross( Eigen::Vector3f::UnitZ() ).normalized();
}

template <class PrimitivesT, class CloudPtrT> inline int
work( CloudPtrT cloud
      , std::string const& prims_path
      , int argc
      , char** argv )
{
    using std::vector;
    using std::string;
    using std::cout;
    using std::endl;
    using am::LinePrimitive;
    using am::MyCloud;
    using am::AlphaExpOptProblem;
    using Eigen::Matrix;

    typedef typename PrimitivesT::value_type    PrimitiveT;
    typedef typename PrimitiveT::Scalar         Scalar;

    float gap_length = 10.f;
    pcl::console::parse_argument( argc, argv, "--gap", gap_length );

    const Scalar scale      = 0.01f;
    const Scalar trunc_at   = 0.2f;
    Matrix<Scalar,4,1> lambdas; lambdas << /* label_cost: */ 1.f, /* align: */ 100000.f, /* smooth: */ .1f, 0.f;
    std::vector<Scalar> desired_angles = { 0.f, M_PI_2, M_PI};
    Scalar angle_threshold = 0.01f;

    vector<int> lines_clusters;
    PrimitivesT prims( am::Primitive<PrimitiveT::Dim>::template read<PrimitiveT>(prims_path,&lines_clusters) );
//    for ( size_t i = 0; i != prims.size(); ++i )
//        cout << "[" << lines_clusters[i] << "] "
//                << prims[i]()(0) << ", " << prims[i]()(1) << ", " << prims[i]()(2) << endl;


    // solve GCO
    am::MaskType     opt_mask( prims.size(), 0 );
    std::vector<int> labels;
    {
        AlphaExpOptProblem::solve( opt_mask
                                   , prims
                                   , cloud
                                   , scale
                                   , lambdas
                                   , trunc_at
                                   , desired_angles
                                   , &labels
                                   , angle_threshold );
        opt_mask.resize( prims.size(), 0 );
        for ( size_t pid = 0; pid != labels.size(); ++pid )
            opt_mask[ labels[pid] ] = 1;
    }

    am::GlobFit2::visualize( prims, opt_mask, cloud, scale * 2, {0,M_PI_2,M_PI},NULL, NULL, NULL, true, true, gap_length );

//    vector<int> points_lines;
//    GF2::PostAligner::points2Prims( points_lines
//                                    , cloud
//                                    , prims
//                                    , &(LinePrimitive::point3Distance) );

//    int err = GF2::PostAligner::run< vector<LinePrimitive>, am::MyCloud::Ptr, LinePrimitive::Scalar, 2>
//            ( prims, cloud, points_lines, lines_clusters, &getLine2DNormal );
    return 0;
}

int main( int argc, char **argv )
{
    using std::vector;
    using std::string;
    using std::cout;
    using std::endl;
    using am::LinePrimitive;
    using am::MyCloud;
    using am::AlphaExpOptProblem;
    using Eigen::Matrix;
    using namespace GF2;

    cout << "Hello post-opt!" << endl;

    string prims_path = "lines.txt";
    if ( argc > 1 ) prims_path = string( argv[1] );
    string cloud_path = "cloud.pcd";
    if ( argc > 2 ) cloud_path = string( argv[2] );

    MyCloud::Ptr cloud( new MyCloud() );
    pcl::io::loadPCDFile( cloud_path, *cloud );
//    for ( size_t pid = 0; pid != cloud->size(); ++pid )
//    {
//        cout << cloud->at( pid ).x << ", " << cloud->at(pid).y << ", " << cloud->at(pid).z << endl;
//    }

    if ( prims_path.find("planes.txt") != std::string::npos )
        return work<vector<am::PlanePrimitive>, am::MyCloud::Ptr>( cloud, prims_path, argc, argv );
    else
        return work<vector<am::LinePrimitive>, am::MyCloud::Ptr>( cloud, prims_path, argc, argv );

    //return EXIT_SUCCESS;
}
