#include <iostream>

#include "pcl/io/pcd_io.h"
#include "pcl/io/ply_io.h"
#include "pcl/io/obj_io.h"
#include "pcl/PolygonMesh.h"

#include "boost/program_options.hpp"

#include "amDefines.h"
#include "pcltools/util.hpp"
#include "CVTools.h"
#include "primitives/planePrimitive.h"
#include "primitives/houghLine.h"

#include "kmeans/kMeans.h"
#include "cuboidRansac.h"
#include "optimization.h"
#include "globFit2.h"
#include "localAnalysis.h"
#include "ground_truth/ground_truth.h"
#include "globfit2/io/io.h"

#include "smartgeometry/util.h"

#ifdef USE_PEARL
#   include "pearl/pearl.hpp"
#endif // use_pearl

using namespace std;
using namespace am;
using namespace GF2;

enum METHOD_TYPE { PEARL_METHOD, GF2_METHOD } method = GF2_METHOD;

// predecl
int autodiff_main();
int testCuda();

template <class PrimitivesT, class PointsPtrT, typename Scalar = typename PrimitivesT::value_type::Scalar> int
sanityCheck( PrimitivesT    const& gt_lines
             , PrimitivesT  const& lines
             , am::MaskType const& opt_mask
             , PointsPtrT          cloud
             , Eigen::Matrix<Scalar,-1,1> const& lambdas
             , Scalar       const  gf2_scale
             , std::vector<Scalar> const& desired_angles
             , Scalar       const  trunc_pw_angle_at)
{
    using std::vector;
    typedef typename PrimitivesT::value_type PrimitiveT;

    vector< vector<PrimitiveT> const* > final_lines = { &lines /*, &gt_lines*/ };
    vector<int> gt_mask( gt_lines.size(), 1 );
    vector< am::MaskType const* > final_masks;
    final_masks.push_back( &opt_mask );
    //final_masks.push_back( &gt_mask );
    vector< string > captions = { "OPTIMIZED_ENERGY: ",
                                  "       GT_ENERGY: " };

    Eigen::Matrix<Scalar,-1,1> energies(4,1);
    for ( size_t i = 0; i != final_lines.size(); ++i )
    {
        // get energy
        std::cout << captions[i] << GlobFit2::queryEnergy( energies, *(final_masks[i]), *(final_lines[i]), cloud, lambdas, gf2_scale, desired_angles, 0.01, trunc_pw_angle_at ) << " = ";
        // print energy
        for ( int row = 0; row != energies.rows(); ++row )
            std::cout << energies(row) * lambdas(row) << ((row == energies.rows()-1) ? "" : " + ");
        std::cout << " = ";
        for ( int row = 0; row != energies.rows(); ++row )
            std::cout << energies(row) << " * " << lambdas(row) << ((row == energies.rows()-1) ? "" : " + ");
        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}

/// SOLVE ///
template <class PrimitiveT> std::vector<PrimitiveT>
gt_primitives_from_lines( std::vector<LinePrimitive> const& gt_lines );
template <> std::vector<LinePrimitive>
gt_primitives_from_lines<LinePrimitive>( std::vector<LinePrimitive> const& gt_lines )
{
    return gt_lines;
};
template <> std::vector<PlanePrimitive>
gt_primitives_from_lines<PlanePrimitive>( std::vector<LinePrimitive> const& gt_lines )
{
    typedef typename PlanePrimitive::Scalar Scalar;

    std::vector<PlanePrimitive> gt_primitives; gt_primitives.reserve( gt_lines.size() );
    for ( size_t lid = 0; lid != gt_lines.size(); ++lid )
    {
        LinePrimitive               line = gt_lines.at( lid );
        Eigen::Matrix<Scalar,3,1>   dir = line.dir().cross( Eigen::Matrix<Scalar,3,1>::UnitZ() );
        gt_primitives.push_back( PlanePrimitive(line.pos(),dir) );
    }

    return gt_primitives;
};

template <class PrimitiveT, class PointsPtrT,typename Scalar = typename PrimitiveT::value_type::Scalar >
int solve( PointsPtrT cloud
           , Scalar                     const working_scale
           , std::vector<Scalar>        const& gf2_desired_angles
           , std::vector<LinePrimitive> const& gt_lines
           , bool                       const show_clusters
           , int                        const desiredK
           , Eigen::Matrix<float,-1,1>  const lambdas
           , int                        const gf2_nThreads
           , int                        const gf2_thread_depth
           , Scalar                     const gf2_trunc_pw_angle_at
           , METHOD_TYPE                const method
           , float                      const gf2_k_explore_ratio
           , float                      const gf2_filter_coeff
           , float                      const gf2_angle_similarity
           , int                        const gf2_add_n_friends
           , int                        const argc
           , char**                     const argv
           #if USE_PEARL
           , Pearl::Params              const* p_pearl_params = NULL
           #endif
           )
{
    typedef vector<PrimitiveT> PrimitivesT;

    PrimitivesT gt_primitives = gt_primitives_from_lines<PrimitiveT>( gt_lines );
    const int         DIM       = (PrimitiveT::Dim == 6) ? 2 : 3;
    const std::string timestamp = smartgeometry::util::timestamp2Str();
    const std::string out_dir = timestamp.substr(1) + "/";
    boost::filesystem::create_directory( out_dir );
    auto lambdas2 = lambdas;
    //lambdas2(0) *= 0.1f;

    vector<PrimitiveClustering<PrimitivesT> > clusterings, clusterings2;
    PrimitivesT candidates, candidates2;

    // output
    MaskType opt_mask( candidates.size(), 0 ), opt_mask2;
    vector<int> labels;

    struct
    {
        int         k       = -1;
        float       energy  = FLT_MAX;
        MaskType    opt_mask;
    } best_mask;

    if ( method == PEARL_METHOD )
    {
#if USE_PEARL
        std::cout << "candidates.size(): " << candidates.size() << std::endl;
        std::cout << "running pearl with " << "beta: " << p_pearl_params->beta
                  << "smooth: " << p_pearl_params->lambdas(2) << std::endl;

        // run Pearl
        Pearl::run( labels
                    , candidates
                    , cloud
                    , NULL
                    , *p_pearl_params
                    , NULL /*beta_range.size() > 1 ? NULL : &label_history
                                                      , beta_range.size() > 1 ? NULL : &line_history */);
        //target_mask.resize( candidates.size(), 0 );

        // parse output labels to opt_mask
        std::cout << "returned candidates.size(): " << candidates.size() << std::endl;
        opt_mask.resize( candidates.size(), 0 );
        for ( size_t pid = 0; pid != labels.size(); ++pid )     opt_mask[ labels[pid] ] = 1;
#else
        std::cerr << "[" << __func__ << "]: " << "define USE_PEARL 1 to run this" << std::endl;
#endif
    }
    else // GF2
    {
        const float filter_coeff    = gf2_filter_coeff      < 0.f ? ((DIM==2)   ?      .1f :   .75f) : gf2_filter_coeff;
        const float angle_sim       = gf2_angle_similarity  < 0.f ? ((DIM == 2) ?  0.0001f : 0.005f) : gf2_angle_similarity;
        const float add_n_friends   = gf2_add_n_friends     < 0.f ? ((DIM==2)   ?       40 :     2 ) : gf2_add_n_friends;

        {
            LocalAnalysis::runSimpler( /*       out_clusterings: */ clusterings
                                       , /*              points: */ cloud
                                       , /*       working_scale: */ working_scale
                                       , /* angle_sim_threshold: */ angle_sim
                                       , /*      max_neigh_size: */ ((DIM == 2) ? 10   : 20  )
                                       , /*      desired_angles: */ gf2_desired_angles
                                       , /*      is_friend_func: */ &(LocalAnalysis::matching<Scalar>)
                                       , /*         addNFriends: */ add_n_friends
                                       , /*        filter_coeff: */ filter_coeff
                                       , /*    input_primitives: */ static_cast<PrimitivesT*>( NULL )
                                       , /* clust_friends_trunc: */ gf2_trunc_pw_angle_at
                                       );

            if ( show_clusters )
                LocalAnalysis::display( clusterings[0].lines, cloud, 0.1f, &(clusterings[0]) ); // display new lines temporarily
        }

        // store candidate generation output
        candidates = clusterings[0].lines;

        std::cout << "desiredK: " << desiredK << std::endl;

        std::ofstream ks_energies( out_dir + "ks_energies.txt" );

        int start_k = std::min( desiredK  , static_cast<int>(floor((2.f - gf2_k_explore_ratio) * desiredK)) );
        int stop_k  = std::max( desiredK+1, static_cast<int>(ceil (       gf2_k_explore_ratio  * desiredK)) );
        std::cout << "start_k: " << start_k << ", stop_k: " << stop_k << std::endl;

        // put desiredK in front, calculate neighbours afterwards
        std::vector<int> ks;
        {
            ks.push_back( desiredK );
            for ( int k = start_k; k != stop_k; ++k )
                if ( k != desiredK )
                    ks.push_back( k );
        }

        for ( size_t k_id = 0; k_id != ks.size(); ++k_id )
        {
            const int k = ks[ k_id ];
            std::cout << "running opt with " << k << std::endl;
            float tmp_energy = FLT_MAX;

            // optimize
            GlobFit2::optimize( /*   [out]             mask: */ opt_mask
                                , /* [in ]       candidates: */ candidates
                                , /* [in ]       clustering: */ &(clusterings[0])
                                , /* [in ]           points: */ cloud
                                , /* [in ]          lambdas: */ lambdas
                                , /* [in ]            scale: */ working_scale
                                , /* [in ]   desired_angles: */ gf2_desired_angles
                                , /* [in ]    restart_count: */ gf2_nThreads
                                , /* [in ]    sa run_length: */ gf2_thread_depth
                                , /* [in ]    inlier thresh: */ 0.01
                                , /* [in ]   trunc at angle: */ gf2_trunc_pw_angle_at
                                , /* [in ]      pnt indices: */ NULL
                                , /* [in ]    target #lines: */ k
                                , /* [out]      assignments: */ /*beta_range.size() > 1 ? NULL : &labels */ NULL
                                , /* [out]  solution_energy: */ &tmp_energy
                                , /* [out] output_directory: */ &out_dir );
            

            //            // split
//            PrimitivesT split_candidates;
//            GlobFit2::split( split_candidates
//                             , opt_mask
//                             , candidates
//                             , cloud );

            // second
            {
                //GlobFit2::visualize( candidates, opt_mask, cloud, working_scale, gf2_desired_angles, NULL, NULL, labels.size() ? &(labels) : NULL, false, true );
                for ( size_t prim_id = 0; prim_id != opt_mask.size(); ++prim_id )
                {
                    if ( !opt_mask[prim_id] )   continue;
                    candidates2.push_back( candidates[prim_id] );
                }

                LocalAnalysis::runSimpler( /*       out_clusterings: */ clusterings2
                                           , /*              points: */ cloud
                                           , /*       working_scale: */ working_scale
                                           , /* angle_sim_threshold: */ -1.f
                                           , /*      max_neigh_size: */ ((DIM == 2) ? 10 : 20  )
                                           , /*      desired_angles: */ gf2_desired_angles
                                           , /*      is_friend_func: */ &(LocalAnalysis::matching<Scalar>)
                                           , /*         addNFriends: */ -1
                                           , /*        filter_coeff: */ 0.f // filter_coeff
                                           , /*    input_primitives: */ &candidates2
                                           , /* clust_friends_trunc: */ gf2_trunc_pw_angle_at
                                           );
                LocalAnalysis::display( clusterings2[0].lines, cloud, 0.1f, &(clusterings2[0]) ); // display new lines temporarily

                // optimize
                GlobFit2::optimize( /*   [out]             mask: */ opt_mask2
                                    , /* [in ]       candidates: */ clusterings2[0].lines
                                    , /* [in ]       clustering: */ &(clusterings2[0])
                                    , /* [in ]           points: */ cloud
                                    , /* [in ]          lambdas: */ lambdas2
                                    , /* [in ]            scale: */ working_scale
                                    , /* [in ]   desired_angles: */ gf2_desired_angles
                                    , /* [in ]    restart_count: */ gf2_nThreads
                                    , /* [in ]    sa run_length: */ gf2_thread_depth
                                    , /* [in ]    inlier thresh: */ 0.01
                                    , /* [in ]   trunc at angle: */ gf2_trunc_pw_angle_at
                                    , /* [in ]      pnt indices: */ NULL
                                    , /* [in ]    target #lines: */ k
                                    , /* [out]      assignments: */ /*beta_range.size() > 1 ? NULL : &labels */ NULL
                                    , /* [out]  solution_energy: */ &tmp_energy
                                    , /* [out] output_directory: */ &out_dir );
                std::cout << "opt_mask2: ";
                std::cout<<"opt_mask2:";
                for(size_t vi=0;vi!=opt_mask2.size();++vi)
                    if ( opt_mask2[vi] )
                    {
                        std::cout<<vi<<" ";
                    }
                std::cout << "\n";

                for(size_t vi=0;vi!=opt_mask2.size();++vi)
                    if ( opt_mask2[vi] )
                        std::cout << "clusterings2[0].lines[ " << vi << "]: " << clusterings2[0].lines[vi]().transpose() << std::endl;

            }


            if ( ks.size() > 1 ) GlobFit2::visualize( candidates, opt_mask, cloud, working_scale, gf2_desired_angles, NULL, NULL, labels.size() ? &(labels) : NULL, false, true );

            ks_energies << k << " " << tmp_energy << std::endl;
            if ( tmp_energy < best_mask.energy )
            {
                best_mask.energy    = tmp_energy;
                best_mask.k         = std::accumulate( opt_mask.begin(), opt_mask.end(), 0 );
                best_mask.opt_mask  = opt_mask;

                std::cout << "storing k " << k << " as best, with mask: ";
                std::cout <<"opt_mask:";for(size_t vi=0;vi!=best_mask.opt_mask.size();++vi) if(best_mask.opt_mask[vi]) std::cout << vi << " ";std::cout << "\n";
            }
        }

        {
            ks_energies.close();
            if ( stop_k - start_k > 1 )
            {
                int err = 0;
                err = system( ("gnuplot -e \"plot '" + out_dir + "ks_energies.txt' u 1:2 w points\" -p").c_str() );
                if ( err != EXIT_SUCCESS ) std::cerr << "system returned code " << err << std::endl;
            }
        }

        opt_mask =  best_mask.opt_mask;
        std::cout << "min_k: " << best_mask.k<< " with energy " << best_mask.energy << std::endl;
        std::cout << "opt_mask:";for(size_t vi=0;vi!=opt_mask.size();++vi) if(opt_mask[vi]) std::cout << vi << " ";std::cout << "\n";
    } // if method == gf2...

    sanityCheck<PrimitivesT, MyCloud::Ptr, float> ( gt_primitives, candidates, opt_mask, cloud, lambdas, working_scale, gf2_desired_angles, gf2_trunc_pw_angle_at );
    if ( clusterings2.size() && clusterings2[0].lines.size() ) sanityCheck<PrimitivesT, MyCloud::Ptr, float> ( gt_primitives, clusterings2[0].lines, opt_mask2, cloud, lambdas2, working_scale, gf2_desired_angles, gf2_trunc_pw_angle_at );

    // dump
    {
        // primitives
        std::string primitives_out_path = out_dir + ( (DIM == 3) ? "planes" : "lines") + ".txt";
        Primitive<PrimitiveT::Dim>::dump( primitives_out_path
                                          , candidates
                                          , clusterings.size() ? &(clusterings[0].lines_clusters)
                                                               : (labels.size() ? &labels
                                                                                : NULL)
                                          , NULL );

        if ( method == GF2_METHOD )
        {
            primitives_out_path = out_dir + ( (DIM == 3) ? "planes2" : "lines2") + ".txt";
            Primitive<PrimitiveT::Dim>::dump( primitives_out_path
                                              , clusterings2[0].lines
                    , clusterings2.size() ? &(clusterings2[0].lines_clusters)
                : (labels.size() ? &labels
                : NULL)
                , NULL );
        }

        // cloud
        pcl::io::savePCDFile( out_dir + "cloud.pcd", *cloud, am::util::pcl::allIndicesOf(cloud)->indices );

        // solution
        {
            ofstream opt_f( out_dir + "solution.txt" );
            for ( size_t k = 0; k != candidates.size(); ++k )
            {
                if ( !opt_mask[k] ) continue;
                opt_f << k << "\t";
                opt_f << candidates[k]().transpose() << std::endl;
            }

            // scale
            opt_f << "scale\t" << working_scale << "\n";

            // desired angles
            opt_f << "desired_angles\t";
            for ( size_t angi = 0; angi != gf2_desired_angles.size(); ++angi )
                opt_f << gf2_desired_angles[angi] << "\t";
            opt_f << std::endl;

            // args
            opt_f << "#";
            for ( int argi = 0; argi != argc; ++argi )
                opt_f << argv[argi] << " ";
            opt_f << "\n";

            opt_f.close();
        }
    }

    // visualize output
    {
        GlobFit2::visualize( candidates, opt_mask, cloud, working_scale, gf2_desired_angles, NULL, NULL, labels.size() ? &(labels) : NULL, true, true );
        if ( method == GF2_METHOD )
            GlobFit2::visualize( clusterings2[0].lines, opt_mask2, cloud, working_scale, gf2_desired_angles, NULL, NULL, NULL, true, true );
        //GlobFit2::visualize( gt_lines  , gt_mask , cloud, 1., NULL, NULL, labels.size() ? &(labels) : NULL       );
    }

    return EXIT_SUCCESS;
}


int main( int argc, char** argv )
{
    std::cout << "entering " << std::endl; fflush(stdout);
    return 0;

    using std::vector;
    typedef float Scalar;

    /// ---------------------------------------------------------------------------------------------------------------
    /// PARAMETERS
    /// for documentation, see "CLI ARGS" section or run with "--help"

    // METHOD
    std::string     method_name         = "gf2";
    // GT
    std::string     gt_name             = "pearl0";
    float           gt_spacing          = 50.f;
    int             gt_nLines           = 5;
    float           gt_decay            = .1f;
    int             gt_nRects           = 2;
    int             gt_nAngles          = 1;
    float           gt_noise            = 0.003f;
    int             gt_nPoints          = 100;
    float           gt_scene_size       = 1.f;
    int             gt_clutter_lines    = 0;
    bool            gt_3D               = false;
    //std::vector<float> gf2_desired_angles = {0,60.f*M_PI/180.f, M_PI_2, 120.f*M_PI/180.f, M_PI}; // {0, M_PI_2, M_PI };
    std::vector<float> gf2_desired_angles = { 0, M_PI_2, M_PI }; // {0, M_PI_2, M_PI };
    //std::vector<float> gf2_desired_angles = { 0, M_PI_4, M_PI_2, M_PI_4 * 3.0, M_PI }; // {0, M_PI_2, M_PI };
    float const     Z                   = 0.f;
    int             seed                = 123456;
    // GF2
    float           gf2_scale           = .03f; // main algorithm input parameter
    float           gf2_trunc_pw_angle_at = 0.25f;
    Eigen::Matrix<float,-1,1> lambdas(4,1); lambdas << 0.f, 1.f, 10.f, 0.f;
    int             gf2_nThreads        = 16;
    int             gf2_thread_depth    = 300;
    float           gf2_k_explore_ratio = 1.5;
    float           gf2_filter_coeff    = -1.f;
    float           gf2_angle_similarity = -1.f;
    int             gf2_add_n_friends   = 40;
    bool            show_clusters       = false;

#if USE_PEARL
    Pearl::Params   pearl_params;
#endif
    /// ---------------------------------------------------------------------------------------------------------------

    // --method gf2 --datacost 0.01 --pwcost 1000 --nthreads 128 --thread_depth 500 --sample_k 1.1 --sample_scale 1 --scale 0.05 --seed 123456 --3D --input-cloud /home/bontius/workspace/SmartGeometry/ransacTest/data/ --K 4
    // --method gf2 --gt stairs --gt-nlines 6 --gt-spacing 4 --gt-decay 0.2 --nrects 8 --nangles 3 --noise 0.001 --datacost 0.01 --pwcost 1000 --npoints 2000 --nthreads 128 --thread_depth 500 --sample_k 1.1 --sample_scale 1 --nodebug --scale 0.05 --seed 123456 --3D --input-cloud /home/bontius/Downloads/house_050_noisy00000.pcd --K 4
    // --method gf2 --gt stairs --gt-nlines 6 --gt-spacing 4 --gt-decay 0.15 --nrects 2 --nangles 2 --noise 0.005 --datacost 0.5 --pwcost 50 --npoints 300 --nthreads 128 --thread_depth 300 --sample_k 1.15 --sample_scale 1.5 --nodebug --scale 0.05 --show-clusters --trunc_at 0.4
    // --method gf2 --gt pearl0 --gt-nlines 6 --gt-spacing 60 --gt-decay 0.15 --nrects 2 --nangles 2 --noise 0.005 --datacost 0.5 --pwcost 50 --npoints 300 --nthreads 128 --thread_depth 300 --sample_k 1.15 --sample_scale 1.5 --nodebug --scale 0.05 --show-clusters --trunc_at 0.4
    // --method gf2 --gt stratified --gt-nlines 6 --gt-spacing 60 --gt-decay 0.15 --nrects 2 --nangles 2 --noise 0.005 --datacost 0.5 --pwcost 50 --npoints 300 --nthreads 128 --thread_depth 300 --sample_k 1.15 --sample_scale 1.5 --nodebug --scale 0.05 --show-clusters
    // --method pearl --beta 1000 --smooth_weight 30000 --gt pearl0 --gt-nlines 6 --gt-spacing 60 --gt-decay 0.15 --nrects 8 --nangles 3 --noise 0.005 --datacost 0.5 --pwcost 50 --npoints 400 --nthreads 128 --thread_depth 300 --sample_k 1.0 --sample_scale 2 --nodebug --scale 0.05 --seed 123456

    // Sunday: --method gf2 --gt stairs --gt-nlines 6 --gt-spacing 4 --gt-decay 0.2 --nrects 8 --nangles 3 --noise 0.001 --datacost 0.5 --pwcost 1000 --npoints 2000 --nthreads  32 --thread_depth 500 --sample_k 1.2 --sample_scale 1 --nodebug --scale 0.05 --seed 123456 --3D --input-cloud /home/bontius/Downloads/house_050_noisy00000.pcd --K 3
    // CLI ARGS
    namespace po = boost::program_options;
    po::variables_map vm;
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options() // --nodebug --gt stratified --nrects 4 --nangles 2 --noise 0.002 --method gf2 --datacost 2500 --pwcost 10 --npoints 500 --nthreads 8 --sample_k 1.2 --thread_depth 500 --sample_scale 1.5
                ("method"       , po::value<std::string>(&method_name     )->default_value("gf2"   ), "[method]    : pearl, __GF2__" )
                ("seed"         , po::value<int        >(&seed            )->default_value(123456  ), "[GT]        : random seed")
                ("gt"           , po::value<std::string>(&gt_name         )->default_value("pearl0"), "[GT]        : __pearl0__, stairs, stratified, (old: rects, 2rects)" )
                ("gt-spacing"   , po::value<float      >(&gt_spacing      )->default_value(50.f    ), "[GT][pearl] : how much space between \"fence\" lines; [stairs]: how much longer a line is than a column (i.e. 1.5)" )
                ("gt-nlines"    , po::value<int        >(&gt_nLines       )->default_value(5       ), "[GT][pearl] : how many \"fence\" lines              ; [stairs]: how many lines" )
                ("gt-decay"     , po::value<float      >(&gt_decay        )->default_value(.1f     ), "[GT]                                                  [stairs]: how much line lenght decays" )
                ("nrects"       , po::value<int        >(&gt_nRects       )->default_value(2       ), "[GT][rects] : how many autorectangles" )
                ("nangles"      , po::value<int        >(&gt_nAngles      )->default_value(1       ), "[GT][rects] : how many angles for the rectangles" )
                ("noise"        , po::value<float      >(&gt_noise        )->default_value(0.003f  ), "[GT]        : how much noise after sampling points" )
                ("npoints"      , po::value<int        >(&gt_nPoints      )->default_value(100     ), "[GT]        : how many points to use to sample the scene")
                ("sample_scale" , po::value<float      >(&gt_scene_size   )->default_value(1.f     ), "[GT]        : scene size (0.f..1.f)" )
                ("clutter-lines", po::value<int        >(&gt_clutter_lines)->default_value(0       ), "[GT]        : add n randomly oriented lines crossing the center" )
                ("input-cloud"  , po::value<std::string>(                 )                         , "[GT]        : load cloud instead of subsampling image - assumes that the generator image is the same!")
                ("input-lines"  , po::value<std::string>(                 )                         , "[GT]        : use input lines instead of fitting locally")
                ("2rect_angle"  , po::value<float      >(                 )                         , "[GT]        : old, dunno" )
                ("scale"        , po::value<float      >(&gf2_scale       )->default_value(0.03f   ), "[GF2]       : operation scale" )
                ("datacost"     , po::value<float      >(&lambdas(0)      )->default_value(-1.f    ), "[GF2]       : lambda weight for GF2 data cost" )
                ("pwcost"       , po::value<float      >(&lambdas(2)      )->default_value(-1.f    ), "[GF2]       : lambda weight for GF2 pairwise cost"  )
                ("trunc-pw"     , po::value<float      >(&gf2_trunc_pw_angle_at)->default_value(0.25f), "[GF2]     : Truncate pw cost at this angle (in radians)" )
                ("nthreads"     , po::value<int        >(&gf2_nThreads    )->default_value(16      ), "[GF2][opt]  : how many initialisations" )
                ("thread_depth" , po::value<int        >(&gf2_thread_depth)->default_value(300     ), "[GF2][opt]  : how long to run one search thread" )
                ("filter-candidates", po::value<float  >(&gf2_filter_coeff)->default_value(-1.f    ), "[GF2]       : candidate filtering strictness: 0.f - no filtering, 0.1f - allowing filtering, 1.f - filtering at scale, 100.f strict filtering")
                ("sample_k"     , po::value<float      >(&gf2_k_explore_ratio)->default_value(1.2  ), "[GF2][opt]  : how big range to sampl around GT #lines" )
                ("angle-thresh" , po::value<float      >(&gf2_angle_similarity)->default_value(-1.f), "[GF2]       : how strict is the clustering" )
                ("n-friends"    , po::value<int        >(&gf2_add_n_friends)->default_value(-1     ), "[GF2]       : how many 'perpendicular' clusters" )
                ("K"            , po::value<int        >(                 )                        )
                ("nodebug"                                                                          , "[PEARL]     : don't allow debug messages" )
                ("show-clusters"                                                                    , "[2D]        : show proposed lines")
                ("3D"                                                                               , "[3D]        : planes instead")
                ("gt-only"                                                                          , "return after showing gt")
                ("help"                                                                             , "produce help message" )
        #if USE_PEARL
                ("maxit"        , po::value<int        >(&pearl_params.max_pearl_iterations), "[pearl]: how many expand-refit iterations are allowed")
                ("beta"         , po::value<float      >(&pearl_params.beta      )->default_value(-1.f) )
                ("smooth_weight", po::value<float      >(&pearl_params.lambdas(2))->default_value(-1.f) )
        #endif
                ;

        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if ( vm.count("help" ) ) { std::cout << desc << "\n"; return EXIT_FAILURE; }
        if ( vm.count("gt"   ) )   std::cout << "gt: "    << gt_name    << std::endl;
        if ( vm.count("noise") )   std::cout << "noise: " << gt_noise << std::endl;
        show_clusters = vm.count("show-clusters");
        gt_3D         = vm.count("3D");
#if USE_PEARL
        if ( vm.count("nodebug") )
            pearl_params.debug = false;
#endif
    }

    srand( seed );

    // METHOD
    if ( vm.count("method") )
    {
        std::transform( method_name.begin(), method_name.end(), method_name.begin(), ::tolower );
        if      ( !(method_name.compare("pearl")) )    method = PEARL_METHOD;
        else if ( !(method_name.compare("gf2"  )) )    method = GF2_METHOD;
        else                                           { std::cerr << "unknown method..." << method_name << "...exiting" << std::endl; return EXIT_FAILURE; }
    }

    // GT
    cv::Mat                       img;
    pcl::PointCloud<MyPoint>::Ptr cloud;
    vector<LinePrimitive>         *p_lines = NULL;
    int                           desiredK;
    vector<float>                 angles;
    if ( vm.count("gt") )
    {
        p_lines = new std::vector<LinePrimitive>();

        if      (   !gt_name.compare("pearl0")     ) { gtPearl( img, *p_lines, gt_spacing, gt_nLines, gt_scene_size ); }
        else if (   !gt_name.compare("stairs")     ) { gtStairs( img, *p_lines, gt_nLines, gt_spacing, gt_decay, gt_scene_size ); } // --method gf2 --gt stairs --gt-nlines 6 --gt-spacing 4 --gt-decay 0.15 --nrects 2 --nangles 2 --noise 0.01 --datacost 0.1 --pwcost 0.6 --npoints 400 --nthreads 8 --sample_k 1.2 --thread_depth 500 --sample_scale 1.5 --nodebug
        else if (   !gt_name.compare("stratified") ) { gtImg_stratified( img, *p_lines, angles, gt_nRects, gt_nAngles, gt_scene_size ); }
        // OLD
        else if (   !gt_name.compare("2rects")     ) { gtImg( img, *p_lines, NULL ); }
        else if (   !gt_name.compare("2rects2")
                 || !gt_name.compare("2rects3")    )
        {
            vector<float> tmp_angles = {0.f,
                                        vm.count("2rect_angle") ? vm["2rect_angle"].as<float>()
                                                                : 40.f };
            vector<cv::Point2i> points =
            {
                cv::Point2i(20 ,20), cv::Point2i(150,220),
                cv::Point2i(190,40), cv::Point2i(240,240),
            };

            gtImg( img, *p_lines, &tmp_angles, gt_name.compare("2rects3") ? NULL : &points );
        }
        else if ( !gt_name.compare("2rects")     ) { gtImg2( img, *p_lines, gt_nRects ); }
        else                                       { std::cerr << "unrecognized GT type..." << gt_name << "...exiting\n"; return EXIT_FAILURE; }

        if ( gt_clutter_lines )
        {
            gtRandomLines( img, *p_lines, gt_clutter_lines, gt_scene_size );
            cv::imshow( "img", img );
            cv::waitKey( 100 );
        }

        if ( vm.count("gt-only") )
        {
            cv::waitKey();
            return EXIT_SUCCESS;
        }

        // count GT k
        desiredK = p_lines->size();
        cv::imwrite( "input.png", img );
    }

    // IMAGE -> CLOUD
    if ( !cloud )
    {
        if ( vm.count("input-cloud") )
        {
            cloud.reset( new MyCloud() );
            std::string path( vm["input-cloud"].as<string>() );
            if ( !boost::filesystem::exists( path ) )
            {
                std::cerr << "input-cloud: " << path << " doest not exist!!! exiting..." << std::endl;
                return EXIT_FAILURE;
            }

            std::cout << "loading " << path << std::endl;
            if ( path.find("pcd") != std::string::npos )
                pcl::io::loadPCDFile( path, *cloud );
            else if ( path.find("ply") != std::string::npos )
            {
                bool success = true;
                try
                {
                    pcl::io::loadPLYFile( path, *cloud );
                    success = true;
                }
                catch (std::exception const& e)
                {
                    std::cerr << e.what() << std::endl;
                    success = false;
                }

                if ( !success )
                {
                    std::cout << "attempting mesh..." << std::endl;

                    pcl::PolygonMesh testMesh;
                    pcl::io::loadPLYFile( path, testMesh );
                    std::cout << "mesh loaded, " << testMesh.cloud.width << " points" << std::endl;
                    return EXIT_FAILURE;
                }
            }
            else if ( path.find("obj") != std::string::npos )
            {
                pcl::io::loadOBJFile( path, *cloud );
            }
            else
            {
                std::cerr << "unrecognized cloud extension: " << path << std::endl;
                return EXIT_FAILURE;
            }

            std::cout << "loaded: " << cloud->size() << " points" << std::endl;
            if ( !cloud->size() ) return EXIT_FAILURE;

            pcl::visualization::PCLVisualizer::Ptr vptr( new pcl::visualization::PCLVisualizer() );
            {
                vptr->setBackgroundColor( .5, .5, .6 );
                vptr->addPointCloud( cloud );
                vptr->spinOnce();
            }

            if ( cloud->size() != static_cast<size_t>(gt_nPoints) )
            {
                std::cerr << "[" << __func__ << "]: " << "cloud->size() != gt_nPoints...\n";

                // subsample cloud
                if ( cloud->size() > static_cast<size_t>(gt_nPoints) )
                {
                    MyCloud::Ptr tmp_cloud( new MyCloud() ); tmp_cloud->reserve( gt_nPoints + 10 );
                    float prob = gt_nPoints / static_cast<float>( cloud->size() );
                    for ( size_t pid = 0; pid != cloud->size(); ++pid )
                    {
                        if ( (rand()/static_cast<float>(RAND_MAX)) > prob ) continue;
                        tmp_cloud->push_back( cloud->at(pid) );
                    }
                    cloud = tmp_cloud;
                }

                // resize cloud
                {
                    MyPoint min_pt, max_pt;
                    pcl::getMinMax3D( *cloud, min_pt, max_pt );
                    const float rescale = (max_pt.getVector3fMap() - min_pt.getVector3fMap()).norm() / gt_scene_size;
                    for ( size_t pid = 0; pid != cloud->size(); ++pid )
                    {
                        cloud->at( pid ).x = (cloud->at(pid).x - min_pt.x) / rescale;
                        cloud->at( pid ).y = (cloud->at(pid).y - min_pt.y) / rescale;
                        cloud->at( pid ).z = (cloud->at(pid).z - min_pt.z) / rescale;
                    }
                }

                {
                    pcl::visualization::PCLVisualizer::Ptr vptr1( new pcl::visualization::PCLVisualizer() );
                    vptr1->setBackgroundColor( .5, .5, .6 );
                    vptr1->addPointCloud( cloud );
                    vptr1->addSphere( cloud->at(cloud->size()/2), gf2_scale, "scale_sphere", 0 );
                    vptr1->setShapeRenderingProperties( pcl::visualization::PCL_VISUALIZER_OPACITY, .5, "scale_sphere" );
                    vptr1->spinOnce(5000);
                }

                //return EXIT_FAILURE;
            }

            if ( vm.count("K") )
                desiredK = vm["K"].as<int>();
        }
        else
        {
            if ( img.empty() )
            {
                std::cerr << "image input disables...I need gt lines provided..." << std::endl;
                return EXIT_FAILURE;
            }

            const int n_points = gt_3D ? pow(gt_nPoints,1/3.) : gt_nPoints;
            GlobFit2::image_2_2DCloud( /*  out_cloud: */ cloud
                                       , /*     path: */ img
                                       , /* N_points: */ (gt_3D ? n_points * n_points : n_points)
                                       , /*        Z: */ Z
                                       , /* scene scale: */ gt_scene_size );

            if ( gt_3D )
            {
                pcl::PointCloud<MyPoint>::Ptr cloud2( new pcl::PointCloud<MyPoint>() );
                for ( size_t pid = 0; pid != cloud->size(); ++pid )
                {
                    int points_new = 0;
                    while ( points_new < n_points )
                    {
                        cloud2->push_back( cloud->at(pid) );
                        float uniform = (rand()/static_cast<float>(RAND_MAX));
                        // inv_cdf(x^2): x = (3*uniform)^(1/3); a^(1/3) = exp(log(a)/3);
                        cloud2->back().z = Z + exp(log(3.f*uniform)/3.f) * gt_scene_size;
                        ++points_new;
                    }
                }
                pcl::copyPointCloud( *cloud2, *cloud );
            }

            Eigen::Vector3f sensor_origin( .0f, .0f, .0f );
            smartgeometry::addGaussianNoise<MyPoint>( cloud
                                                      , NULL
                                                      , { gt_noise, gt_noise, gt_3D ? gt_noise : 0.f }
                                                      , Eigen::Vector3f::Zero()
                                                      , &sensor_origin );
        }
    }

    // load lines?
    if ( vm.count("input-lines") )
    {

    }

#if USE_PEARL
    pearl_params.max_neighbourhood_radius = gf2_scale;
#endif
    std::cout << "entering " << std::endl; fflush(stdout);
    return 0;

    if ( gt_3D )
        solve<PlanePrimitive>( cloud
                               , gf2_scale
                               , gf2_desired_angles
                               , *p_lines
                               , show_clusters
                               , desiredK
                               , lambdas
                               , gf2_nThreads
                               , gf2_thread_depth
                               , gf2_trunc_pw_angle_at
                               , method
                               , gf2_k_explore_ratio
                               , gf2_filter_coeff
                               , gf2_angle_similarity
                               , gf2_add_n_friends
                               , argc
                               , argv
                               , &pearl_params );
    else
        solve<LinePrimitive>( cloud
                              , gf2_scale
                              , gf2_desired_angles
                              , *p_lines
                              , show_clusters
                              , desiredK
                              , lambdas
                              , gf2_nThreads
                              , gf2_thread_depth
                              , gf2_trunc_pw_angle_at
                              , method
                              , gf2_k_explore_ratio
                              , gf2_filter_coeff
                              , gf2_angle_similarity
                              , gf2_add_n_friends
                              , argc
                              , argv
                              , &pearl_params );

    SAFE_DELETE( p_lines );

    return 0;
}

#if USE_CUDA
#include "gf2cuda/gf2Cuda.h"
#include "my_cuda_util/myCudaUtil.h"

int testCuda()
{
    std::vector<float> test_series;
    int limit = 1;
    for ( int s = 0; s != 13; ++s )
    {
        test_series.clear();
        for ( int i = 0; i < limit; ++i )
            test_series.push_back( rand() % 2 );
        std::cout << "testing " << test_series.size() << ":\n";
        smartgeometry::cuda_util::testReduce( test_series );

        if ( limit == 1 ) limit = 2;
        else limit *= 10;
    }

    return 0;
}
#endif // USE_CUDA
