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


namespace GF2
{
    namespace io
    {
        //! \brief Dumps primitives with GID and DIR_GID to disk.
        //! \tparam PrimitiveT Concept: PrimitiveContainerT::value_type::value_type aka GF2::LinePrimitive2.
        //! \tparam PrimitiveContainerT Concept: vector< vector< GF2::LinePrimitive2 > >.
        template <class PrimitiveT /* = typename PrimitiveContainerT::value_type::value_type*/, class PrimitiveContainerT> inline int
        savePrimitives( PrimitiveContainerT const& primitives, std::string out_file_name, bool verbose = false )
        {
            const int Dim = PrimitiveT::Dim;
            typedef typename PrimitiveT::VectorType VectorType;

            // out_lines
            std::ofstream out_file( out_file_name );
            for ( size_t lid = 0; lid != primitives.size(); ++lid )
            {
                for ( size_t lid1 = 0; lid1 != primitives[lid].size(); ++lid1 )
                {
                    for ( int d = 0; d != Dim; ++d )
                        //out_file << ((VectorType)primitives[lid][lid1])(d) << ((d!=Dim-1)?",":"\n");
                        out_file << std::setprecision(9) << ((VectorType)primitives[lid][lid1])(d) << ",";
                    out_file << primitives[lid][lid1].getTag( PrimitiveT::GID ) << ",";
                    out_file << primitives[lid][lid1].getTag( PrimitiveT::DIR_GID ) << "\n";
                }
            }
            out_file.close();
            if ( verbose ) std::cout << "[" << __func__ << "]: " << "saved " << out_file_name << std::endl;

            return EXIT_SUCCESS;
        }

        template <
                   class       PrimitiveT          /*= typename PrimitiveContainerT::value_type::value_type*/
                 , class       PrimitiveContainerT
                 >
        inline int readPrimitives( PrimitiveContainerT &lines, std::string const& path )
        {
            typedef typename PrimitiveT::Scalar Scalar;
            enum {                              Dim    = PrimitiveT::Dim };
            typedef typename PrimitiveContainerT::value_type PatchT;
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
                if ( floats.size() < Dim )
                    std::cerr << "[" << __func__ << "]: " << "not good, floats.size() < Dim..." << std::endl;

                // rest
                int gid = -1, dir_gid = -1;
                if ( !iss.eof() )
                {
                    if ( std::getline(iss, tmp_str, ',') )
                    {
                        gid = atoi( tmp_str.c_str() );
                        //lines[lid].back().setTag( PrimitiveT::GID, gid );
                    }
                    if ( std::getline(iss, tmp_str, ',') )
                    {
                        dir_gid = atoi( tmp_str.c_str() );
                        //lines[lid].back().setTag( PrimitiveT::DIR_GID, dir_gid );
                    }
                } // if patch information

                // insert into proper patch, if gid specified
                if ( gid > -1 )
                {
//                    if ( gid >= lines.size() )
//                        lines.resize( gid+1 );

//                    lines[ gid ].emplace_back( PrimitiveT(floats) );
//                    lines[ gid ].back().setTag( PrimitiveT::GID, gid );
//                    lines[ gid ].back().setTag( PrimitiveT::DIR_GID, dir_gid );
                    tmp_lines[ gid ].push_back( PrimitiveT(floats) );
                    tmp_lines[ gid ].back().setTag( PrimitiveT::GID    , gid     );
                    tmp_lines[ gid ].back().setTag( PrimitiveT::DIR_GID, dir_gid );
                }
                else // just make a new patch for it
                {
                    lines.push_back( PatchT() );
                    lines[lid].push_back( PrimitiveT(floats) );
                    lines[lid].back().setTag( PrimitiveT::GID, line_count ); // 1D indexing only
                }

                ++line_count;
            } // while getline

            // copy all patches from map to vector (so that there are no empty patches in the vector)
            lines.reserve( lines.size() + tmp_lines.size() );
            typename PatchMap::const_iterator end_it = tmp_lines.end();
            for ( typename PatchMap::const_iterator it = tmp_lines.begin(); it != end_it; ++it )
            {
                lines.push_back( PatchT() );
                lines.back().insert( lines.back().end(), it->second.begin(), it->second.end() );
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

        //! \brief Solver::readPoints Read stored points, and convert them to non-PCL format.
        //! \param [out]points        Output point vector
        //! \param [in] path          PLY source path
        //! \return EXIT_SUCCESS
        template < class _PointT /*= typename _PointContainerT::value_type */
                   , class _PointContainerT
                 >
        inline int
        readPoints( _PointContainerT &points
                  , std::string        path
                  , pcl::PointCloud<pcl::PointXYZRGB>::Ptr *cloud_arg = NULL )
        {
            pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud;

            // sample image
            if ( !cloud ) cloud.reset( new pcl::PointCloud<pcl::PointXYZRGB>() );
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
                points.emplace_back( _PointT(raw_points[pid]) );
                points.back().setTag( _PointT::PID, pid );
                points.back().setTag( _PointT::GID, pid );
            }
            if ( cloud_arg ) *cloud_arg = cloud;

            return EXIT_SUCCESS;
        } // ...Solver::readPoints()

#if 0 // deprecated
        extern int
        readSolution( std::string const path );

        template <class PrimitivesT, typename Scalar = typename PrimitivesT::value_type::Scalar> int
        saveSolution( std::string path
                      , MaskType const& opt_mask
                      , PrimitivesT const& candidates
                      , Scalar const working_scale
                      , std::vector<Scalar> const& gf2_desired_angles
                      , int  argc
                      , char **argv);
#endif
    } // ... ns io
} // ... ns GF2

#ifndef __GF2_IO_INC_HPP__
#   define __GF2_IO_INC_HPP__
#   include "globfit2/io/io.hpp"
#endif

#endif // __GF2_IO_H__
