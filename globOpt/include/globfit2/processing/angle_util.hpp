#ifndef GF2_ANGLE_UTIL_HPP
#define GF2_ANGLE_UTIL_HPP

#include "globfit2/optimization/energyFunctors.h"

namespace GF2
{
    namespace angles
    {
        /*! \brief Adds angles from generator to in/out vector \p angles.
         *
         *  \tparam _AnglesContainerT Concept: std::vector<_Scalar>.
         *  \tparam _Scalar           Floating point precision type.
         *  \param[in,out] angles     Container (vector) to append to.
         *  \param[in] angle_gen      Generator element. 0, \p angle_gen, 2 * \p angle_gen, ..., M_PI will be appended.
         *  \param[in] verbose        Enable logging.
         */
        template <class _AnglesContainerT> inline int
        appendAnglesFromGenerators( _AnglesContainerT &angles, _AnglesContainerT angle_gens, bool no_parallel, char verbose, bool inRad = false )
        {
            typedef typename _AnglesContainerT::Scalar Scalar;

            if ( !inRad )
                std::cout << "[" << __func__ << "]: " << "ASSUMING DEGREES" << std::endl;

            std::set<Scalar> angles_set;
            // copy
            angles_set.insert( angles.begin(), angles.end() );

            // insert 0 element
            if ( !no_parallel )
                angles_set.insert( Scalar(0) );

            for ( int i = 0; i != angle_gens.size(); ++i )
            {
                Scalar angle_gen_rad = angle_gens[i] * (inRad ? Scalar(1.) : M_PI/Scalar(180.));

                if ( angle_gen_rad == Scalar(0.) )
                {
                    std::cerr << "[" << __func__ << "]: " << "skipping 0 as generator" << std::endl;
                    continue;
                }

                // generate
                for ( Scalar angle = angle_gen_rad; angle < M_PI; angle+= angle_gen_rad )
                    angles_set.insert( angle );
            }

            // insert 0 element
            if ( !no_parallel )
                angles_set.insert( Scalar(M_PI) );

            angles.resize( angles_set.size() );
            std::copy( angles_set.begin(), angles_set.end(), angles.begin() );

            // print
            if ( verbose )
            {
                std::cout << "[" << __func__ << "]: " << "Desired angles: {";
                for ( size_t vi=0;vi!=angles.size();++vi)
                    std::cout << angles[vi] << "(" << angles[vi] * Scalar(180.) / M_PI << ")" << ((vi==angles.size()-1) ? "" : ", ");
                std::cout << "}\n";
            }

            return EXIT_SUCCESS;
        } //...appendAnglesFromGenerator()
    } //...ns angles

    template <typename _Scalar, class _AnglesT>
    inline int deduceGenerators( _AnglesT &angle_gens, _AnglesT const& angles, const int verbose = 0 )
    {
        if ( verbose )
        {
            std::cout<<"[" << __func__ << "]: " << "angles:";
            for(size_t vi=0;vi!=angles.size();++vi)std::cout<<angles[vi]*180./M_PI<<" ";std::cout << "\n";
        }

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
                if ( verbose )
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

    typedef std::map< DidT, std::map<int,int> > DirAngleMapT; // <did, <angle_id, count> >

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
                 || (    (_prim0.getTag( _PrimitiveT::TAGS::STATUS)!=_PrimitiveT::STATUS_VALUES::ACTIVE)
                      || ( prim1.getTag( _PrimitiveT::TAGS::STATUS)!=_PrimitiveT::STATUS_VALUES::ACTIVE) )
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
    inline void genAngles( AnglesT &single_gen, _Scalar const& angle, AnglesT angleGensInRad, bool const verbose = false )
    {
        if ( angle == _Scalar(0.) || angle == _Scalar(M_PI) )
            return;

        int i = 0;
        while ( i != angleGensInRad.size() )
        {
            if ( angleGensInRad[i] == _Scalar(0.) || angleGensInRad[i] == _Scalar(M_PI) )
            {
                ++i;
                continue;
            }

            if ( (fmod(angle, angleGensInRad[i]) == _Scalar(0.)) && (angleGensInRad[i] != _Scalar(0.)) )
            {
                single_gen.push_back( angleGensInRad[i] );
                break;
            }
            ++i;
        }

        if ( single_gen.size() )
            if ( verbose ) std::cout << "[" << __func__ << "]: " << "deduced " << single_gen[0] * 180. / M_PI << " from " << angle * 180. / M_PI << std::endl;
    }

    /*! \brief Decides, which perfect angles to allow for each direction id
     * \return allowedAngles filled with custom "angles" for each direction id
     */
    inline void selectAngles( std::map< DidT, AnglesT >            & allowedAngles
                            , DirAngleMapT                   const& dirAngles
                            , AnglesT                        const& angles
                            , AnglesT                        const& angle_gens
                            , bool                           const  verbose = false )
    {
        /*! \brief { <angle_id, count>, ... }
         *         i.e. { <0, 10>, <2, 6> } means:
         *         -  10 pairs are parallel and 6 pairs are 180 (also parallel), if angles == { 0, 90, 180 }
         *         -  but means 6 pairs with 90, if angle_gens == { 60, 90 }, angles == { 0, 60, 90, 120, 180 }
         */
        typedef          DirAngleMapT::mapped_type  AngleHistogramT;
        typedef typename AnglesT::value_type        Scalar;

        std::map<Scalar,int> votes; // { <angle_gen, votecount>, ... }

        // for all did-s
        DirAngleMapT::const_iterator end_it = dirAngles.end();
        for ( DirAngleMapT::const_iterator it = dirAngles.begin(); it != end_it; ++it )
        {
            // first: did
            // second: map< angle_id, count >
            DidT            const  did  = it->first;
            AngleHistogramT const& hist = it->second;

            // Skip, if already set in earlier iteration (and added from GEN_ANGLE tag)
            if ( allowedAngles.find(did) != allowedAngles.end() ) continue;

            // deduce generators \todo Aron: replace with deduceGenerators? maybe not doing the same...
            int j = 0;
            for ( AngleHistogramT::const_iterator hist_it = hist.begin(); hist_it != hist.end(); ++hist_it, ++j )
            {
                // first: angle_id
                // second: ocurrence count
                const int angle_id = hist_it->first;

                // skip votes on parallel
                if ( (angles[angle_id] == Scalar(0.)) || (angles[angle_id] == Scalar(M_PI)) )
                    continue;
                // debug
                else if ( j == 0 || j == angles.size() )
                    std::cout << "[" << __func__ << "]: " << "not skipping " << angles[angle_id] << std::endl;

                // for each generator, check, if the angle is an integer divisor of it
                for ( int i = 0; i != angle_gens.size(); ++i )
                {
                    // skip parallel
                    if ( (angle_gens[i] == Scalar(0.)) || (angle_gens[i] == Scalar(M_PI)) )
                        continue;
                    // debug
                    else if ( i == 0 || i == angle_gens.size() )
                            std::cout << "[" << __func__ << "]: " << "not skipping " << angle_gens[i] << std::endl;

                    // check common denominator
                    double intpart = 0;
                    Scalar fractpart = modf( angles[angle_id] / angle_gens[i], &intpart );

                    /* if k x gen + ~0, then this is a multiple of the generator angle_gens[i]
                     * i.e. angles[angle_id] == 60, then 1 x 60 + 0 --> vote for 60, if angle_gens = { 60, 90 }
                     * , but if angles[angle_id] == 120, angle_gens = { 90 }, then 1 x 90 + 0.5 will not vote for 90
                     */
                    if ( (fractpart < 1e-6) && intpart )
                    {
                        votes[ i ] += hist_it->second;
                        if ( verbose ) std::cout << "[" << __func__ << "]: " << "added vote " << hist_it->second
                                                 << " to angle_gen " << angle_gens[i]
                                                 << " due to angle " << angles[angle_id] << std::endl;
                    }
                } //...for generators
            } //...for voted angle_ids

            // select max voted angle \todo change to probabilistic
            int max_vote = 0, max_vote_id = 0;
            if ( verbose ) std::cout << "[" << __func__ << "]: " << "votes:";
            for ( auto vote_it = votes.begin(); vote_it != votes.end(); ++vote_it )
            {
                // first: angle_gen
                // second: vote count
                if ( verbose ) std::cout << "<" << vote_it->first << "," << vote_it->second << ">, ";
                if ( vote_it->second > max_vote )
                {
                    max_vote    = vote_it->second;
                    max_vote_id = vote_it->first;
                }
            }
            if ( verbose ) std::cout << std::endl;

            // Add, if has at least 1 non-parallel vote (parallel votes are not put into votes)
            if ( max_vote )
            {
                /* Generate multiple from single generator
                 * (i.e. 60 -> 0, 60, 120, 180, or 90 -> 0, 90, 180 )
                 */
                AnglesT single_gen;
                genAngles( single_gen, angles[max_vote_id], angle_gens );

                // store angles from single generator for did
                angles::appendAnglesFromGenerators( /*      out: */ allowedAngles[ did ]
                                                  , /*       in: */ single_gen
                                                  , /* no_paral: */ false
                                                  , /*  verbose: */ true
                                                  , /*    inRad: */ true );
            }

            // clear votes to prepare for next did
            votes.clear();
        } //...for all dids
    } //...selectAngles()
} //...GF2

#endif // GF2_ANGLE_UTIL_HPP
