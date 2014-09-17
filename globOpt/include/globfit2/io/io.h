#ifndef __GF2_IO_H__
#define __GF2_IO_H__

#include <string>
#include <iomanip>
#include <fstream> // ifstream
#include "Eigen/Sparse"

#if GF2_USE_PCL
#   include <pcl/io/ply_io.h>
#   include <pcl/io/pcd_io.h>
#endif // GF2_USE_PCL

#include "globfit2/util/pcl_util.hpp"
#include "globfit2/util/containers.hpp"


namespace GF2
{
    namespace io
    {
        //! \brief Dumps primitives with GID and DIR_GID to disk.
        //! \tparam PrimitiveT Concept: PrimitiveContainerT::value_type::value_type aka GF2::LinePrimitive2.
        //! \tparam PrimitiveContainerT Concept: vector< vector< GF2::LinePrimitive2 > >.
        template <class PrimitiveT, class _inner_const_iterator, class PrimitiveContainerT> inline int
        savePrimitives( PrimitiveContainerT const& primitives
                      , std::string out_file_name
                      , bool verbose = false
                      )
        {
            typedef typename PrimitiveContainerT::const_iterator outer_const_iterator;

            //const int Dim = PrimitiveT::Dim;
            typedef typename PrimitiveT::VectorType VectorType;

            // out_lines
            std::string parent_path = boost::filesystem::path(out_file_name).parent_path().string();
            if ( !parent_path.empty() )
                if ( !boost::filesystem::exists(parent_path) )
                    boost::filesystem::create_directory( boost::filesystem::path(parent_path) );

            std::ofstream out_file( out_file_name );
            if ( !out_file.is_open() ) { std::cerr << "could not open file..." << out_file_name << std::endl; return EXIT_FAILURE; }

            int lid = 0;
            outer_const_iterator gid_end_it = primitives.end();
            for ( outer_const_iterator gid_it = primitives.begin(); gid_it != gid_end_it; ++gid_it, ++lid )
            //for ( size_t lid = 0; lid != primitives.size(); ++lid )
            {
                int lid1 = 0;
                _inner_const_iterator lid_end_it = containers::valueOf<PrimitiveT>(gid_it).end();
                for ( _inner_const_iterator lid_it = containers::valueOf<PrimitiveT>(gid_it).begin(); lid_it != lid_end_it; ++lid_it, ++lid1 )
                //for ( size_t lid1 = 0; lid1 != primitives[lid].size(); ++lid1 )
                {
                    //for ( int d = 0; d != Dim; ++d )
                        //out_file << std::setprecision(9) << ((VectorType)*lid_it)(d) << ",";
                    //out_file << std::setprecision(9) << ((VectorType)primitives.at(lid).at(lid1))(d) << ",";

                    // saves primitive part
                    out_file << (*lid_it).toFileEntry();
                    // save taggable part (todo: merge these)
                    out_file << lid_it->getTag( PrimitiveT::GID ) << ",";
                    out_file << lid_it->getTag( PrimitiveT::DIR_GID ) << ",";
                    out_file << lid_it->getTag( PrimitiveT::CHOSEN )  << "\n";
                }
            }
            out_file.close();
            if ( verbose ) std::cout << "[" << __func__ << "]: " << "saved " << out_file_name << std::endl;

            return EXIT_SUCCESS;
        }

        //! \brief Reads primitives with their GIDs and dir_GIDs from file.
        //! \tparam PatchT Concept: vector< \ref GF2::LinePrimitive2 >.
        template <
                   class       PrimitiveT          /*= typename PrimitiveContainerT::value_type::value_type*/
                 , class       PatchT
                 , class       PrimitiveContainerT
                 >
        inline int readPrimitives( PrimitiveContainerT &lines
                                 , std::string const& path
                                 //, void(*add)(PrimitiveContainerT & prims, PrimitiveT const& prim, int gid, int dir_gid ) = NULL
                                 , std::map<int, typename PrimitiveContainerT::value_type> *patches = NULL
                                 )
        {
            typedef typename PrimitiveT::Scalar Scalar;
            //enum {                              Dim    = PrimitiveT::Dim }; // Important! This decides how many floats to read
            const int Dim = PrimitiveT::getFileEntryLength();
            //typedef typename PrimitiveContainerT::value_type PatchT;
            typedef std::map<int, PatchT>                    PatchMap; // <GID, vector<primitives> >

            // open file
            std::ifstream file( path.c_str() );
            if ( !file.is_open() )
            {
                std::cerr << "[" << __func__ << "] couldn't open file" << std::endl;
                return EXIT_FAILURE;
            }

            std::map<int, PatchT> tmp_lines;
            int lid        = 0; // deprecated, tracks linear id
            int line_count = 0;
            std::string line;
            while ( getline(file, line) )
            {
                if ( line[0] == '#' )          continue;   // skip comments
                if ( line.empty()   ) { ++lid; continue; } // deprecated, groups used to be separated by empty lines

                std::vector<Scalar> floats;
                std::istringstream iss( line );
                std::string        tmp_str;
                while ( /**/ (floats.size() < Dim)
                        &&   (std::getline(iss, tmp_str, ',') || std::getline(iss, tmp_str)) )
                {
                    floats.push_back( atof(tmp_str.c_str()) );
                }

                // error check
                if ( floats.size() < Dim )
                    std::cerr << "[" << __func__ << "]: " << "not good, floats.size() < Dim..." << std::endl;

                // rest
                int gid = -1, dir_gid = -1, chosen = -1;
                if ( !iss.eof() )
                {
                    if ( std::getline(iss, tmp_str, ',') )  gid     = atoi( tmp_str.c_str() );
                    if ( std::getline(iss, tmp_str, ',') )  dir_gid = atoi( tmp_str.c_str() );
                    if ( std::getline(iss, tmp_str, ',') )  chosen  = atoi( tmp_str.c_str() );
                } // if patch information

                // insert into proper patch, if gid specified
                if ( gid > -1 )
                {
                    //tmp_lines[ gid ].push_back( PrimitiveT(floats) );
                    tmp_lines[ gid ].push_back( PrimitiveT::fromFileEntry(floats) );
                    tmp_lines[ gid ].back().setTag( PrimitiveT::GID    , gid     );
                    tmp_lines[ gid ].back().setTag( PrimitiveT::DIR_GID, dir_gid );
                    tmp_lines[ gid ].back().setTag( PrimitiveT::CHOSEN , chosen  );
                }
                else // just make a new patch for it
                {
                    throw new std::runtime_error("[io::readPrims] code not up to date to handle gid==-1 cases, please add proper GID to primitives");
//                    lines.push_back( PatchT() );
//                    lines[lid].push_back( PrimitiveT(floats) );
//                    lines[lid].back().setTag( PrimitiveT::GID, line_count ); // 1D indexing only
                }

                ++line_count;
            } // while getline

            // copy all patches from map to vector (so that there are no empty patches in the vector)
            //lines.reserve( lines.size() + tmp_lines.size() );
            typename PatchMap::const_iterator end_it = tmp_lines.end();
            for ( typename PatchMap::const_iterator it = tmp_lines.begin(); it != end_it; ++it )
            {
//                if ( add )
//                {
//                    typename PatchT::const_iterator end_it2 = it->second.end();
//                    for ( typename PatchT::const_iterator it2 = it->second.begin(); it2 != end_it2; ++it2 )
//                    {
//                        int tmp_gid     = it2->getTag(PrimitiveT::GID);
//                        int tmp_dir_gid = it2->getTag(PrimitiveT::DIR_GID);
//                        add( lines, *it2, tmp_gid, tmp_dir_gid );

//                        if ( tmp_gid != it->first )
//                            std::cerr << "[" << __func__ << "]: " << "it->first " << it->first << " != " << tmp_gid << " it2->second.gid" << std::endl;
//                    }
//                }
//                else
//                {
                    // assume vector<vector< PrimitiveT > >
                    lines.push_back( PatchT() );
                    lines.back().insert( lines.back().end(), it->second.begin(), it->second.end() );
//                }
            }
            if ( patches )
            {
                *patches = tmp_lines;
            }

            file.close();

            return EXIT_SUCCESS;
        } // ... readPrimitives()

        //! \brief                      Reads point-primitive associations from file
        //! \param points_primitives    [Out]    points_primitives[pid] = pair<lid,lid1>
        //! \param path                 [In]     Path of file to read
        //! \param linear_indices       [In/Out] If not null, will get filled like this: linear_indices[pid] = line_id
        //! \return                     EXIT_SUCCESS if could open file
        inline int readAssociations( std::vector<std::pair<int,int> >      & points_primitives
                                     , std::string                    const& path
                                     , std::map<int,int>                   * linear_indices )
        {
            std::ifstream f( path.c_str() );
            if ( !f.is_open() )
            {
                std::cerr << "[" << __func__ << "] couldn't open file" << std::endl;
                return EXIT_FAILURE;
            }

            std::set< std::pair<int,int> > lines;
            std::string line;
            while ( getline(f, line) )
            {
                if ( line[0] == '#' )
                    continue;

                std::vector<int> ints;
                std::istringstream iss( line );
                std::string        tmp_str;
                while ( /**/ (ints.size() < 3)
                        &&   (std::getline(iss, tmp_str, ',') || std::getline(iss, tmp_str)) )
                {
                    ints.push_back( atoi(tmp_str.c_str()) );
                }
                if ( ints.size() != 3 )
                {
                    std::cerr << "[" << __func__ << "]: " << "skipping " << line << ", since could not read 3 ints" << std::endl;
                    continue;
                }

                // 2d indices (lid,lid1)
                if ( static_cast<int>(points_primitives.size()) <= ints[0] )
                    points_primitives.resize( ints[0]+1, std::pair<int,int>(-1,-1) );
                points_primitives[ ints[0] ] = std::pair<int,int>( ints[1], ints[2] );
                lines.insert( std::pair<int,int>(ints[1], ints[2]) );

                // 1d indices (line_id)
                if ( linear_indices )
                {
                    (*linear_indices)[ ints[0] ] = lines.size() - 1; // assumes sorted set
                }

            } // while getline

            return EXIT_SUCCESS;
        } // ... readAssociations

        //! \brief                       Write points' associations to GID and DIR_GID.
        //! \tparam     _PointPrimitiveT Concept: \ref GF2::PointPrimitive.
        //! \tparam     _PointContainerT Concept: vector< \ref GF2::PointPrimitive >
        //! \param[in]  points           Output point vector
        //! \param[in]  path             PLY source path
        //! \return                      EXIT_SUCCESS
        template < class _PointPrimitiveT
                 , class _PointContainerT
                 >
        inline int writeAssociations( _PointContainerT const& points
                                    , std::string const& f_assoc_path )
        {
            // open
            std::ofstream f_assoc( f_assoc_path );
            if ( !f_assoc.is_open() ) { std::cerr << "[" << __func__ << "]: " << "could not open " << f_assoc_path << " for writing..." << std::endl; return EXIT_FAILURE; }

            // preamble
            f_assoc << "# point_id,primitive_gid,primitive_dir_gid" << std::endl;
            // write points
            for ( size_t pid = 0; pid != points.size(); ++pid )
            {
                f_assoc << pid
                           << "," << points[pid].getTag( _PointPrimitiveT::GID )
                           << "," << -1  // assigned to patch, but no direction
                           << std::endl;
            }
            // finish
            f_assoc.close();

            std::cout << "[" << __func__ << "]: " << "wrote to " << f_assoc_path << std::endl;

            return EXIT_SUCCESS;
        } //...writeAssociations

        //! \brief                    Read stored points, and convert them to non-PCL format.
        //! \param[out] points        Output point vector
        //! \param[in]  path          PLY source path
        //! \return                   EXIT_SUCCESS
        template < class _PointT /*= typename _PointContainerT::value_type */
                   , class _PointContainerT
                 >
        inline int
        readPoints( _PointContainerT &points
                  , std::string        path
                  , pcl::PointCloud<pcl::PointNormal>::Ptr *cloud_arg = NULL )
        {
            pcl::PointCloud<pcl::PointNormal>::Ptr cloud;

            // sample image
            if ( !cloud ) cloud.reset( new pcl::PointCloud<pcl::PointNormal>() );
            if ( path.find("ply") != std::string::npos )
                pcl::io::loadPLYFile( path, *cloud );
            else
                pcl::io::loadPCDFile( path, *cloud );

            // convert to raw vector
            std::vector<typename _PointT::VectorType> raw_points;
            pclutil::cloudToVector<typename _PointT::RawAllocator>( raw_points, cloud );

            // convert to tagged vector
            points.reserve( raw_points.size() );
            for ( size_t pid = 0; pid != raw_points.size(); ++pid )
            {
                raw_points[pid](3) = cloud->at(pid).normal_x;
                raw_points[pid](4) = cloud->at(pid).normal_y;
                raw_points[pid](5) = cloud->at(pid).normal_z;

                points.emplace_back( _PointT(raw_points[pid]) );
                points.back().setTag( _PointT::PID, pid );
                points.back().setTag( _PointT::GID, pid );
            }
            if ( cloud_arg ) *cloud_arg = cloud;

            return EXIT_SUCCESS;
        } // ...Solver::readPoints()

        //! \brief                   Write points to almost PLY
        //! \param[in] points        Points
        //! \param[in] path          PLY destination path
        //! \return EXIT_SUCCESS
        template < class _PointT /*= typename _PointContainerT::value_type */
                   , class _PointContainerT
                 >
        inline int
        writePoints( _PointContainerT &points
                   , std::string       path )
        {
            std::ofstream file( path.c_str() );
            if ( ! file.is_open() ) { std::cerr << "[" << __func__ << "]: " << "could not open " << path << std::endl; return EXIT_FAILURE; }

            file << "ply\n"
                 << "format ascii 1.0\n"
                 << "comment Aron generated\n"
                 << "element vertex " << points.size() << "\n"
                 << "property float x\n"
                 << "property float y\n"
                 << "property float z\n"
                 << "property float nx\n"
                 << "property float ny\n"
                 << "property float nz\n"
                 << "end_header\n";

            for ( typename _PointContainerT::const_iterator it = points.begin(); it != points.end(); ++it )
            {
                Eigen::Matrix< typename _PointT::Scalar, 3, 1 > const& pos = it->pos();
                Eigen::Matrix< typename _PointT::Scalar, 3, 1 > const& dir = it->dir();
                file << pos(0) << " " << pos(1) << " " << pos(2) << " " << dir(0) << " " << dir(1) << " " << dir(2) << std::endl;
            }

            file.close();

            return EXIT_SUCCESS;
        } //...writePoints()

    } // ... ns io
} // ... ns GF2

#ifndef __GF2_IO_INC_HPP__
#   define __GF2_IO_INC_HPP__
#   include "globfit2/io/io.hpp"
#endif

#endif // __GF2_IO_H__

