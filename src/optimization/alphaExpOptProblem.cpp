#include "optimization/alphaExpOptProblem.h"

namespace am
{
    AlphaExpOptProblem::~AlphaExpOptProblem()
    {
        if ( _gc )
        {
            try
            {
                delete _gc;
                _gc = NULL;
            }
            catch (GCException e)
            {
                e.Report();
            }
        }
    }

    AlphaExpOptProblem::AlphaExpOptProblem( std::vector<LinePrimitive>           const& lines
                                            , pcl::PointCloud<MyPoint>::ConstPtr        cloud
                                            , double                             const  threshold
                                            , std::vector<int>                   const* indices
                                            , Eigen::VectorXf                    const& lambdas
                                            , int                                       max_step_count
                                            , float                              const  trunc_pw_at_angle
                                            )
        : _lines                ( lines )
        //        , _ids_energies_capacity( 10 )
        , _lambdas              ( lambdas )
        , _cloud                ( cloud )
        //        , _threshold            ( threshold )
        , _energy_up_to_date    ( false )
        //        , _step_count           ( -1 )
        //        , _max_step_count       ( max_step_count )
        //        , _trunc_pw_at_angle    ( trunc_pw_at_angle )
    {
        if ( indices )
            _indices = *indices;
    }

    int
    AlphaExpOptProblem::init( std::vector<int> const& config )
    {
        //_configs_energies.push_front( std::pair<std::vector<int>,Scalar>(config,-1) );
#if 0
        const int num_pixels = _cloud->size();
        const int num_labels = _lines.size();

        // first set up the array for data costs
        int *data = new int[ num_pixels * num_labels ];
        for ( int pid = 0; pid != num_pixels; ++pid )
            for ( int line_id = 0; line_id != num_labels; ++line_id )
                data[pid*num_labels+line_id] = _lines[line_id].distanceToPoint( /*     point: */ _cloud->at(pid).getVector4fMap()
                                                                                , /* squared: */ true                             );

        // next set up the array for smooth costs
        int *smooth = new int[num_labels*num_labels];
        for ( int l1 = 0; l1 < num_labels; l1++ )
            for (int l2 = 0; l2 < num_labels; l2++ )
                smooth[l1+l2*num_labels] = (l1-l2)*(l1-l2) <= 4  ? (l1-l2)*(l1-l2):4;

        try
        {
            _gc = new GCoptimizationGeneralGraph( num_pixels,num_labels );
            _gc->setDataCost(data);
            _gc->setSmoothCost(smooth);

            // now set up a grid neighborhood system
            // first set up horizontal neighbors
            for (int y = 0; y < height; y++ )
                for (int  x = 1; x < width; x++ )
                    gc->setNeighbors(x+y*width,x-1+y*width);

            // next set up vertical neighbors
            for (int y = 1; y < height; y++ )
                for (int  x = 0; x < width; x++ )
                    gc->setNeighbors(x+y*width,x+(y-1)*width);
        }
        catch (GCException e){
            e.Report();
        }

        delete [] smooth;
        delete [] data;

#endif
        return EXIT_SUCCESS;
    }

#if 0
    AlphaExpOptProblem::Scalar
    AlphaExpOptProblem::calcEnergy( std::vector<int>             const& line_mask
                                    , Eigen::Matrix<Scalar,-1,1>      * energies )
    {
        const int N = _cloud->size();
        const int K = _lines.size();

        Eigen::Matrix<Scalar,-1,1> curr_energies( Eigen::Matrix<Scalar,3,1>::Zero() );

        // E_complexity
        {
            curr_energies( 1 ) = std::accumulate( line_mask.begin(), line_mask.end(), 0 );
            //std::cout << "complexity term is " << e_complexity[comb_id] << std::endl;
            if ( curr_energies(1) <= 0.f )   return curr_energies(0);
        }

        // E_data
        {
            // get line inliers
            std::vector<int>    points_lines ( N, -1      );
            std::vector<double> min_distances( N, DBL_MAX );
            {
                // init sac model
                pcl::SampleConsensusModelLine<MyPoint> sacline( _cloud );
                if ( _indices.size() )
                    sacline.setIndices( _indices );
                std::vector<double> tmp_distances;

                // assignments - get closest line for each point
                for ( int line_id = 0; line_id != K; ++line_id )
                {
                    if ( !line_mask[line_id] ) continue; // skip unchosen lines
                    LinePrimitive const& line = _lines[ line_id ];

                    sacline.getDistancesToModel( line.coeffs(), tmp_distances );
                    for ( int pnt_id = 0; pnt_id != N; ++pnt_id )
                    {
                        if ( tmp_distances[pnt_id] < min_distances[pnt_id] )
                        {
                            min_distances[pnt_id] = tmp_distances[pnt_id];
                            points_lines [pnt_id] = line_id;
                        }
                    } // for pnt_id
                } // for line_id
            } // get line inliers

            curr_energies(0) = std::inner_product( min_distances.begin(), min_distances.end(), min_distances.begin(), 0.f );
            //std::cout << "data term is " << e_data[comb_id] << std::endl;

        } // E_data

        // E_pw
        {
            // for every i-j pair
            for ( int i = 0; i != K-1; ++i )
            {
                // skip unchosen lines
                if ( !line_mask[i] ) continue;
                // cache
                LinePrimitive const& line0 = _lines[ i ];

                for ( int j = i+1; j != K; ++j )
                {
                    if ( !line_mask[j] ) continue;
                    LinePrimitive const& line1 = _lines[ j ];

                    float cos_angle    = line0.dir().dot( line1.dir() );
                    float min_ang_diff = std::min( fabs(1.f - cos_angle), fabs(cos_angle) ); // TODO: generalize
                    curr_energies(2) += exp( min_ang_diff ) - 1.f;
                }
            }
            //std::cout << "pw term is " << e_pw[comb_id] << std::endl;
        } // E_pw

        // final E
        if ( energies )
        {
            *energies = curr_energies;
        }

        return _lambdas.dot( curr_energies );
    }
#endif

    AlphaExpOptProblem::Scalar
    AlphaExpOptProblem::getEnergy()
    {
#if 0
        //        if ( _energy_up_to_date )
        //            return _configs_energies.front().second;
        Scalar ret_energy;

        try
        {
            printf("\nBefore optimization energy is %Ld",gc->compute_energy());
            gc->expansion(2);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
            ret_energy = gc->compute_energy();
            printf("\nAfter optimization energy is %Ld\n", ret_energy );

            for ( int  i = 0; i < num_pixels; i++ )
                result[i] = gc->whatLabel(i);

            for ( int y = 0; y < height; ++y )
            {
                for ( int x = 0; x < width; ++x )
                {
                    std::cout << result[x+y*width] << ", ";
                }

                std::cout << std::endl;
            }

            // E_complexity
            {
                curr_energies( 1 ) = std::accumulate( line_mask.begin(), line_mask.end(), 0 );
                //std::cout << "complexity term is " << e_complexity[comb_id] << std::endl;
                if ( curr_energies(1) <= 0.f )   return curr_energies(0);
            }

        }
        catch (GCException e)
        {
            e.Report();
        }

        //        // selection
        //        std::vector<int> const& line_mask = _configs_energies.front().first;

        //        // final E
        //        _configs_energies.front().second = calcEnergy( line_mask );
        //        _energy_up_to_date               = true;

        //        return _configs_energies.front().second;

        _energy_up_to_date = true;;
        return ret_energy;
#endif
        return -1.f;
    }

    bool
    AlphaExpOptProblem::hasNext()
    {
        return !_energy_up_to_date;
    }

    int
    AlphaExpOptProblem::next()
    {
        std::cerr << "[" << __func__ << "] " << "not implemented" << std::endl;
        return EXIT_SUCCESS;
    }

    int
    AlphaExpOptProblem::stepBack()
    {
        std::cerr << "[" << __func__ << "] " << "not implemented" << std::endl;
        return EXIT_SUCCESS;
    }

    std::vector<int>
    AlphaExpOptProblem::getConfig()
    {
        // TODO
        std::cerr << "[" << __func__ << "] " << "todo" << std::endl;
        //return _configs_energies.front().first;

        return std::vector<int>();
    }

} // ns am
