#ifndef GF2_ANGLE_UTIL_HPP
#define GF2_ANGLE_UTIL_HPP

namespace GF2
{

    template <typename _Scalar, class _AnglesT>
    inline int deduceGenerators( _AnglesT &angle_gens, _AnglesT const& angles )
    {

        std::cout<<"[" << __func__ << "]: " << "angles:";
        for(size_t vi=0;vi!=angles.size();++vi)std::cout<<angles[vi]<<" ";std::cout << "\n";

        // first element is always the parallel generator
        angle_gens.push_back( _Scalar(0.) );

        // check every input angle, if it's already a multiple of an output angle generator
        // if not, add it to the generators

        // for every input angle
        for ( int i = 0; i != angles.size(); ++i )
        {
            // skip parallel
            if ( (angles[i] == _Scalar(0.)) || (angles[i] == _Scalar(M_PI)) )
                continue;

            bool found_divisor = false;

            // for every generator
            // start from 1, we don't want division by 0
            for ( int j = 0; j < angle_gens.size(); ++j )
            {
                if ( (angle_gens[j] == _Scalar(0.)) || (angle_gens[j] == _Scalar(M_PI)) )
                    continue;

                double intpart = 0.;
                double rest = modf( angles[i] / angle_gens[j], &intpart );
                std::cout << "[" << __func__ << "]: "
                          << angles[i] * 180. / M_PI
                          << " / " << angle_gens[j] * 180. / M_PI
                          << " = " << intpart
                          << " * " << angle_gens[j] * 180. / M_PI
                          << " + " << rest
                          << std::endl;

                if ( (intpart != 0.) && (rest == 0.) )
                {
                    found_divisor = true;
                    break;
                }

            } //...angle_gens

            if ( !found_divisor )
                angle_gens.push_back( angles[i] );
        } //...angles
    } //...deduceGenerators

    typedef std::map< int, std::map<int,int> > DirAngleMapT; // <did, <angle_id, count> >


    template <typename _PrimitiveT, typename _PrimitiveContainerT, typename _AnglesT>
    struct DirAnglesFunctorInner
    {
        typedef typename _PrimitiveT::Scalar Scalar;

        _PrimitiveT const& _prim0;
        _AnglesT    const& _angles; // input all angles

        DirAngleMapT &_dirAngles; // accumulator out

        inline DirAnglesFunctorInner( _PrimitiveT const& prim0, _AnglesT const& angles, DirAngleMapT &dirAngles ) : _prim0( prim0 ), _angles( angles ), _dirAngles(dirAngles) {}

        inline int eval( _PrimitiveT const& prim1, int /*lid*/ )
        {
            if (
                    ( _prim0.getTag(_PrimitiveT::TAGS::DIR_GID) != prim1.getTag(_PrimitiveT::TAGS::DIR_GID) ) // have same ID
                 || (    (_prim0.getTag( _PrimitiveT::STATUS)!=_PrimitiveT::STATUS_VALUES::ACTIVE)
                      || ( prim1.getTag( _PrimitiveT::STATUS)!=_PrimitiveT::STATUS_VALUES::ACTIVE) )
               )
                return 0;

            int closest_angle_id = 0;
            Scalar angle = std::abs( MyPrimitivePrimitiveAngleFunctor::eval(_prim0, prim1, _angles, &closest_angle_id) );
            if ( angle < Scalar(1e-6) )
            {
                _dirAngles[ _prim0.getTag(_PrimitiveT::TAGS::DIR_GID) ][ closest_angle_id ]++;
            }

            return 1;
        }
    }; //...DirAnglesFunctorInner

    template <typename _PrimitiveT, typename _PrimitiveContainerT, typename _AnglesT>
    struct DirAnglesFunctorOuter
    {
        DirAnglesFunctorOuter( _PrimitiveContainerT const& primitives, _AnglesT const& angles, DirAngleMapT &dirAngles )
            : _primitives( primitives ), _angles( angles ), _dirAngles( dirAngles ) {}

        _PrimitiveContainerT const& _primitives;
        _AnglesT             const& _angles;
        DirAngleMapT              & _dirAngles; // output

        inline int eval( _PrimitiveT const& prim, int /*lid*/ )
        {
            DirAnglesFunctorInner<_PrimitiveT,_PrimitiveContainerT,_AnglesT> _innerFunctor( prim, _angles, _dirAngles );
            return processing::filterPrimitives<_PrimitiveT,typename _PrimitiveContainerT::mapped_type>( _primitives, _innerFunctor );
        }
    }; //... DirAnglesFunctorOuter

    template <typename _Scalar>
    inline void genAngles( AnglesT &single_gen, _Scalar const& angle, AnglesT angleGensInRad )
    {
        int i = 0;
        while ( i != angleGensInRad.size() )
        {
            if ( angleGensInRad[i] == _Scalar(0.) || angleGensInRad[i] == _Scalar(M_PI) )
            {
                ++i;
                continue;
            }

            if ( fmod(angle, angleGensInRad[i]) == _Scalar(0.) )
            {
                single_gen.push_back( angleGensInRad[i] );
                break;
            }
            ++i;
        }

        if ( single_gen.size() )
            std::cout << "[" << __func__ << "]: " << "deduced " << single_gen[0] * 180. / M_PI << " from " << angle << std::endl;
    }

    /*! \brief Decides, which perfect angles to allow for each direction id
     * \return allowedAngles filled with custom "angles" for each direction id
     */
    inline void selectAngles( std::map< int, AnglesT >            & allowedAngles
                            , DirAngleMapT                   const& dirAngles
                            , AnglesT const& angles
                            , AnglesT const& angle_gens )
    {
        typedef typename AnglesT::value_type Scalar;
        typedef DirAngleMapT::mapped_type AngleHistogramT; // vector<angle_id>

        std::map<Scalar,int> votes; // <angle_gen, votecount>

        DirAngleMapT::const_iterator end_it = dirAngles.end();
        for ( DirAngleMapT::const_iterator it = dirAngles.begin(); it != end_it; ++it )
        {
            AngleHistogramT const& hist = it->second;
            int j = 0;
            for ( AngleHistogramT::const_iterator hist_it = hist.begin(); hist_it != hist.end(); ++hist_it, ++j )
            {
                const int angle_id = hist_it->first;

                // skip votes on parallel
                if ( (angles[angle_id] == Scalar(0.)) || (angles[angle_id] == Scalar(M_PI)) )
                    continue;
                else if ( j == 0 || j == angles.size() )
                    std::cout << "not skipping " << angles[angle_id] << std::endl;

                // for each generator, check, if the angle is an integer divisor of it
                for ( int i = 0; i != angle_gens.size(); ++i )
                {
                    // skip parallel
                    if ( (angle_gens[i] == Scalar(0.)) || (angle_gens[i] == Scalar(M_PI)) )
                        continue;
                    else if ( i == 0 || i == angle_gens.size() )
                            std::cout << "not skipping " << angle_gens[i] << std::endl;

                    double intpart = 0;
                    Scalar fractpart = modf( angles[angle_id] / angle_gens[i], &intpart );

                    if ( (fractpart < 1e-6) && intpart )
                    {
                        votes[ i ] += hist_it->second;
                        std::cout << "added vote " << hist_it->second << " to angle_gen " << angle_gens[i] << " due to angle " << angles[angle_id] << std::endl;
                    }
                }
            } // for voted angles

            int max_vote = 0, max_vote_id = 0;
            std::cout << "votes:";
            for ( auto it = votes.begin(); it != votes.end(); ++it )
            {
                std::cout << "<" << it->first << "," << it->second << ">, ";
                if ( it->second > max_vote )
                {
                    max_vote = it->second;
                    max_vote_id = it->first;
                }
            }
            std::cout << std::endl;

            if ( max_vote )
            {
                AnglesT single_gen;
                genAngles( single_gen, angles[max_vote_id], angle_gens );
                processing::appendAnglesFromGenerators( allowedAngles[it->first], single_gen, /* no_paral: */ false, true, /*inRad:*/ true );
            }

            votes.clear();
        }
    }
} //...GF2

#endif // GF2_ANGLE_UTIL_HPP
