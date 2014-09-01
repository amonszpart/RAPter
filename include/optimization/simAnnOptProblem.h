#ifndef __GF2_SIMANNOPTPROBLEM_H__
#define __GF2_SIMANNOPTPROBLEM_H__

#include <vector>
#include <deque>
#include "Eigen/Dense"

#include "my_types.h"
#include "optimization/optProblem.h"
#include "primitives/linePrimitive.h"

#if GF2_USE_CUDA
#include "gf2cuda/gf2Cuda.h"
#endif
#include "lineClustering.hpp"

namespace am
{
    struct SimAnnStats
    {
            SimAnnStats( int max_steps ) { hamming_hits.reserve(max_steps); rev_jumps.reserve(max_steps); }

            std::vector< std::pair<int,float> > hamming_hits;
            std::vector< std::pair<int,float> > rev_jumps;

            void addHammingHit( std::pair<int,float> hit ) { hamming_hits.push_back(hit); }
            void addReverseJump( std::pair<int,float> jump ) { rev_jumps.push_back(jump); }
    };

    template <typename Scalar>
    struct PointPrimitiveDistanceCache
    {
            PointPrimitiveDistanceCache( int n_lines, int max_size = 50 )
                //: distances ( n_lines, std::vector<Scalar>() )
                : max_size  ( max_size )
                , ready     ( n_lines, false ) {}

            std::map<int,std::vector< Scalar > >    distances; // stores distances from point point_id to line line_id at distances[line_id][point_id]
            int                                  max_size;  // how many entries are allowed in distances
            std::deque< int >                    history;   // front was added most recently to distances, back is the oldest entry
            std::vector<bool>                    ready;     // is it ok to read distances at this index?

            bool contains( int id ) const
            {
                return ( ready[id] ); //&& (distances[id].size() > 0) );
            }

            std::vector<Scalar> const& get( int id ) const
            {
                if ( !ready[id] /*|| (id > distances.size())*/ )
                    std::cerr << "[" << __func__ << "]: " << " cache does not have this many lines or not ready..." << id << " vs. " << distances.size() << std::endl;
//                else
//                    std::cout << "cache hit " << id << std::endl;

                return distances.at( id );
            }

            void push( int id, std::vector<Scalar> const& dists )
            {
#               pragma omp critical
                {
                    if ( distances.find(id) != distances.end() )
                        std::cerr << "[" << __func__ << "]: " << "something is wrong, distances[" << id << "]: " << distances.at(id).size() << "\n";

                    //distances.insert( std::pair<id, dists> );

                    distances.insert( std::pair<int,std::vector<Scalar> >(id, dists) );
                    ready[ id ] = true;
                    history.push_front( id );

                    if ( history.size() > max_size )
                    {
                        //std::cout << "dropping " << history.back() << std::endl;
                        if ( !distances[ history.back() ].size() ) std::cerr << "[" << __func__ << "]: " << " dropping empty " << history.back() << std::endl;
                        distances.erase( distances.find(history.back()) );
                        ready[ history.back() ] = false;
                        history.pop_back();
                    }
                }
            }
    };

    template <class PrimitivesT>
    class SimAnnOptProblem : public OptProblem
    {
        public:
            typedef typename PrimitivesT::value_type PrimitiveT;
            typedef typename PrimitiveT::Scalar Scalar;

            SimAnnOptProblem( PrimitivesT                          const& lines
                              , pcl::PointCloud<MyPoint>::ConstPtr        cloud
                              , double                             const  threshold
                              , std::vector<int>                   const* indices
                              , Eigen::VectorXf                    const& lambdas
                              , float                                     scale
                              , int                                       max_step_count = 65536
                              , float                              const  trunc_pw_at_angle = .25f
                              , PointPrimitiveDistanceCache<Scalar>     * p_points_prims_dist_cache = NULL );

            virtual Scalar              getEnergy( void const* p_desired_angles );
            virtual int                 next     ( void const* p_desired_angles );
            virtual bool                hasNext  ();
            virtual int                 stepBack ();
            virtual int                 init( MaskType const& config );
            virtual MaskType    getConfig();

            int
            initWithClusters( PrimitiveClustering<PrimitivesT> const& clustering
                              , MaskType                       const& config
                              , int                            const  K         );

            static Scalar
            calcDataTerms(pcl::PointCloud<MyPoint>::ConstPtr        cloud
                          , PrimitivesT                      const& lines
                          , MaskType                         const& line_mask
                          , std::vector<int>                 const* indices = NULL
                          , float                            const* scale   = NULL
                          , std::vector<int>                      * labels  = NULL
                          , PointPrimitiveDistanceCache<Scalar>   * cache   = NULL
                          , Eigen::Matrix<float, 3, 1>       const* sensor  = NULL );

            static Scalar inline
            calcPairwiseTerm( PrimitiveT    const& l1
                              , PrimitiveT  const& l2
                              , std::vector<Scalar> const &desired_angles
                              , Scalar              const trunc_at  = 0.4f
                              , bool                const verbose   = false
                              );

            SimAnnOptProblem::Scalar
            calcEnergy( std::vector<int>               const& line_mask
                        , std::vector<Scalar>          const& desired_angles
                        , Eigen::Matrix<Scalar, -1, 1>      * energies          = NULL
                        , bool                         const  verbose           = false ) const;

            int getStepCount() { return _step_count; }
            int setSubsampleLinesRatio( Scalar rat ) { _subsample_lines_ratio = rat; return EXIT_SUCCESS; }

            SimAnnStats report() const;

            inline std::vector<MaskType> powerSet(int);
            inline Scalar
            getExhaustiveBest( MaskType &opt_mask, std::vector<Scalar> const& desired_angles ) const;

        protected:
            typedef std::deque<std::pair<MaskType,Scalar> > HistoryT;

            PrimitivesT                                 _lines;
            HistoryT                                    _configs_energies;
            const size_t                                _configs_energies_capacity;
            Eigen::VectorXf                             _lambdas;

            pcl::PointCloud<MyPoint>::ConstPtr          _cloud;
            std::vector<int>                            _indices;
            const double                                _threshold; // probably deprected
            const float                                 _scale;     // the real thing

            bool                                        _energy_up_to_date;
            int                                         _step_count, _max_step_count;
            Eigen::Matrix<Scalar,-1,-1>                 _smooth_scores;
            Scalar                                      _subsample_lines_ratio; // used in "next()", we don't want to check every line on greedy discard

#           if GF2_USE_CUDA
            GF2::Gf2Cuda                                _gf2_cuda;
#           endif

            PrimitiveClustering<PrimitivesT>            _clustering;
            //mutable std::vector< std::vector<float> >   _lines_points_dist_cache;
            mutable PointPrimitiveDistanceCache<Scalar>* _points_prims_dist_cache;
            Scalar                                      _trunc_pw_at_angle;
        public:
            int                                         _omp_thread_id; // used to debug outputs
            Eigen::Matrix<Scalar,3,1>                   _sensor;
            SimAnnStats                                 _statistics;
    };
} // ns am

#ifndef __GF2_INC_HOUGHLINE_HPP__
#   define __GF2_INC_HOUGHLINE_HPP__
#   include "simAnnOptProblem.hpp"
#endif // __GF2_INC_HOUGHLINE_HPP__

#endif // __GF2_SIMANNOPTPROBLEM_H__
