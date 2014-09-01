#ifndef __GF2_EXHAUSTIVEOPTPROBLEM_H__
#define __GF2_EXHAUSTIVEOPTPROBLEM_H__

#include "optimization/optProblem.h"

#include <vector>
#include <deque>
#include "Eigen/Dense"

#include "my_types.h"
#include "primitives/houghLine.h"

namespace am
{
    class ExhaustiveOptProblem : public OptProblem
    {
        public:
            typedef float Scalar;

            ExhaustiveOptProblem( std::vector<LinePrimitive>           const& lines
                                  , pcl::PointCloud<MyPoint>::ConstPtr        cloud
                                  , double                             const  threshold
                                  , std::vector<int>                   const* indices
                                  , Eigen::VectorXf                    const& lambdas );

            virtual Scalar  getEnergy();
            virtual int     next();
            virtual bool    hasNext();
            virtual int     stepBack();
            virtual int     init( std::vector<int> const& config );
            virtual std::vector<int> getConfig();

        protected:
            std::vector< LinePrimitive    >     _lines;
            std::vector< std::vector<int> >     _combs;
            std::deque<std::pair<unsigned,Scalar> >  _ids_energies;
            const size_t                        _ids_energies_capacity;
            Eigen::VectorXf                     _lambdas;

            pcl::PointCloud<MyPoint>::ConstPtr  _cloud;
            std::vector<int>                    _indices;
            const double                        _threshold;

            bool                                _energy_up_to_date;



    };

    std::vector<std::vector<int> >
    powerSet( int K )
    {
        if ( K == 1 )
            return {{0},{1}};

        std::vector< std::vector<int> > prev = powerSet( K-1 );

        std::vector< std::vector<int> > out;
        out.resize( prev.size() << 1 );

        int offset = 0;
        for ( int prefix = 0; prefix != 2; ++prefix )
        {
            for ( size_t i = 0; i != prev.size(); ++i )
            {
                out[offset+i].resize( prev[i].size() + 1 );
                out[offset+i][0] = prefix;
                std::copy( prev[i].begin(), prev[i].end(), out[offset+i].begin()+1 );
            }
            offset += prev.size();
        }

        return out;
    }

    ExhaustiveOptProblem::ExhaustiveOptProblem( std::vector<LinePrimitive>           const& lines
                                                , pcl::PointCloud<MyPoint>::ConstPtr        cloud
                                                , double                             const  threshold
                                                , std::vector<int>                   const* indices
                                                , Eigen::VectorXf                    const& lambdas )
        : _lines( lines )
        , _combs( powerSet(lines.size()) )
        , _ids_energies_capacity( 10 )
        , _lambdas( lambdas )
        , _cloud( cloud )
        , _threshold        ( threshold )
        , _energy_up_to_date( false )
    {
        _ids_energies.push_front( std::pair<unsigned,Scalar>({0,-1}) );
        if ( indices )
            _indices = *indices;
    }

    ExhaustiveOptProblem::Scalar
    ExhaustiveOptProblem::getEnergy()
    {
        const int N = _cloud->size();
        const int K = _lines.size();

        Eigen::Matrix<Scalar,-1,1> curr_energies( Eigen::Matrix<Scalar,3,1>::Zero() );

        // selection
        std::vector<int> const& line_mask = _combs[ _ids_energies.front().first ];

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
        _ids_energies.front().second = _lambdas.dot( curr_energies );
        _energy_up_to_date           = true;

        return _ids_energies.front().second;
    }

    bool
    ExhaustiveOptProblem::hasNext()
    {
        return _ids_energies.front().first+1 < _combs.size();
    }

    int
    ExhaustiveOptProblem::next()
    {
        // update
        _ids_energies.push_front( std::make_pair<int,Scalar>( _ids_energies.front().first+1,
                                                              static_cast<Scalar>(-1)       ) );
        _energy_up_to_date = false;

        // limit deque
        if ( _ids_energies.size() > _ids_energies_capacity )
            _ids_energies.pop_back();

        // debug
        if ( _ids_energies.front().first % 100 )
            std::cout << _ids_energies.front().first << "/" << _combs.size() << std::endl;

        return EXIT_SUCCESS;
    }

    int
    ExhaustiveOptProblem::stepBack()
    {
        _ids_energies.pop_front();
        _energy_up_to_date = true;

        return EXIT_SUCCESS;
    }

    std::vector<int>
    ExhaustiveOptProblem::getConfig()
    {
        return _combs[ _ids_energies.front().first ];
    }

    int ExhaustiveOptProblem::init( std::vector<int> const& config )
    {
        _ids_energies.push_front( std::pair<int,Scalar>(0,-1) );

        return EXIT_SUCCESS;
    }

//private void evaluateNextStepPosition() {

//		// current value of the cost function
//		double currentCost = getGlobalCost(costFactor);
//		double bestCost = currentCost;
//		status.reset();
//		status.setTemperature(getCurrentTemperature());
//		status.incrementNumber();

//		for (long index = 1; index <= iterations; index++) {

//			final double currentCostOLD = currentCost;
//			currentCost = evaluateNodePositionCandidate(currentCost);
//			if (currentCost > bestCost) {
//				final double temp = (bestCost - currentCost) / status.getTemperature();
//				if (temp > -15.0 && Math.exp(temp) > Math.random()) {
//					bestCost = currentCost;
//					status.incrementWorse();
//				} else {
//					swapNodes();
//					status.incrementRejected();
//					currentCost = currentCostOLD;
//				}
//			} else if (currentCost < bestCost) {
//				bestCost = currentCost;
//				status.incrementBetter();
//			} else {
//				status.incrementConst();
//			}
//		}
//	}

} // ns am

#endif // __GF2_EXHAUSTIVEOPTPROBLEM_H__
