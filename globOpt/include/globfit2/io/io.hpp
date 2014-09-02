#include "io.h"
#include <fstream>

#ifdef WITH_QCQPCPP
#   include "qcqpcpp/io/io.h"
#endif

namespace GF2
{
    namespace io
    {
#if 0 // deprecated
        template <class PrimitivesT, typename Scalar = typename PrimitivesT::value_type::Scalar> int
        saveSolution( std::string path
                      , MaskType const& opt_mask
                      , PrimitivesT const& candidates
                      , Scalar const working_scale
                      , std::vector<Scalar> const& gf2_desired_angles
                      , int  argc
                      , char **argv)
        {
            std::ofstream opt_f( path );
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

            return EXIT_SUCCESS;
        } // ... saveSolution
        template <class PrimitivesT> int
        readSolution( MaskType           & mask
                      , std::string const& path
                      , PrimitivesT      * primitives )
        {
#if 0

            MaskType opt_mask(planes.size(),0);

            std::vector<float> desired_angles;
            float scale = -1.f;
            {
                std::ifstream f_solution( solution_path );
                if ( !f_solution.is_open() )
                {
                    std::cerr << "[" << __func__ << "]: " << "could not open " << solution_path << "...exiting\n";
                    return EXIT_FAILURE;
                }

                std::string line;
                enum READ_STATE { PRIMS, SCALE, DESIRED_ANGLES } read_state = PRIMS;
                while ( getline(f_solution, line) )
                {
                    // skip comments
                    if ( line[0] == '#') continue;

                    // read word by word
                    std::istringstream  iss( line );
                    std::string         word;
                    int                 word_id = 0;
                    while ( std::getline(iss, word, '\t') )
                    {
                        // set state
                        if ( word_id == 0 )
                        {
                            if ( line.find("scale") != std::string::npos )
                                read_state = READ_STATE::SCALE;
                            else if ( line.find("desired_angles") != std::string::npos )
                                read_state = READ_STATE::DESIRED_ANGLES;
                            else
                                read_state = READ_STATE::PRIMS;
                        }

                        switch ( read_state )
                        {
                            case READ_STATE::PRIMS:
                                if ( !word_id ) // first number is the primitive id
                                    opt_mask[ atoi(word.c_str()) ] = 1;
                                break;
                            case READ_STATE::SCALE:
                                if ( word_id == 1 )
                                    scale = atof(word.c_str());
                                break;
                            case READ_STATE::DESIRED_ANGLES:
                                if ( word_id > 0)
                                {
                                    float angle = atof(word.c_str());
                                    std::cout << "adding angle " << angle << std::endl; fflush(stdout);
                                    desired_angles.push_back( angle );
                                    std::cout << "angles is noow ";
                                    for(size_t vi=0;vi!=desired_angles.size();++vi)std::cout<<desired_angles[vi]<<" ";std::cout << "\\n"; fflush(stdout);
                                }
                                break;
                        }
                        ++word_id;
                    }
                }

                f_solution.close();

                std::cout << "[" << __func__ << "]: " << "read opt_mask: ";
                for ( size_t vi=0; vi != opt_mask.size(); ++vi)
                    if ( opt_mask[vi] )
                        std::cout << vi << " ";
                std::cout << "\n";
                std::cout << "[" << __func__ << "]: " << "read scale: " << scale << std::endl;
                std::cout << "[" << __func__ << "]: " << "read desired_angles: ";
                for ( size_t vi=0; vi != desired_angles.size(); ++vi)
                    std::cout << desired_angles[vi] << " ";
                std::cout << "\n";
            }

            if ( (scale < 0.f) || !desired_angles.size() )
            {
                std::cerr << "[" << __func__ << "]: " << "scale < 0.f or no angles...exiting\\n";
                return EXIT_FAILURE;
            }
#endif

            return EXIT_SUCCESS;
        } // ... readSolution
#endif // deprecated
    } // ... ns io
} // ... ns GF2
