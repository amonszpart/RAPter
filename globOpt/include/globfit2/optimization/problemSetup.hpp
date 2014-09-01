#ifndef GF2_PROBLEMSETUP_HPP
#define GF2_PROBLEMSETUP_HPP

#include <vector>

namespace GF2
{
    namespace ProblemSetup
    {
        //! \brief Adds unary costs to problem based on point to primitive associations.
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
                 , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        associationBasedDataCost( _OptProblemT              & problem
                                , _PrimitiveContainerT const& prims
                                , _PointContainerT     const& points
                                , _AssocT              const& lids_varids
                                , _WeightsT            const& weights
                                , _Scalar              const  /*scale*/ )
        {
            typedef typename _AssocT::key_type IntPair;

            int err = EXIT_SUCCESS;

            for ( size_t lid = 0; lid != prims.size(); ++lid )
            {
                // check, if any directions for patch
                if ( !prims[lid].size() )
                {
                    std::cerr << "[" << __func__ << "]: " << "no directions for patch[" << lid << "]! This shouldn't happen. Skipping patch..." << std::endl;
                    continue;
                }

                // cache patch group id to match with point group ids
                const int gid = prims[lid][0].getTag( _PrimitiveT::GID );

                // for each direction
                for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
                {
                    // point count for normalization
                    unsigned cnt = 0;
                    // data-cost coefficient (output)
                    _Scalar unary_i = _Scalar(0);
                    // for each point, check if assigned to main patch (TODO: move assignment test to earlier, it's indep of lid1)
                    for ( size_t pid = 0; pid != points.size(); ++pid )
                    {
                        //if ( points[pid].getTag( PointPrimitiveT::GID ) == static_cast<int>(lid) )
                        if ( points[pid].getTag( _PointPrimitiveT::GID ) == gid )
                        {
                            // if within scale, add unary cost
                            _Scalar dist = _PointPrimitiveDistanceFunctor::template eval<_Scalar>( points[pid], prims[lid][lid1] );
                            unary_i += dist;
                            ++cnt;              // normalizer
                        }
                    } // for points

                    _Scalar coeff = cnt ? /* complx: */ weights(2) + /* unary: */ weights(0) * unary_i / _Scalar(cnt)
                                        : /* complx: */ weights(2) + /* unary: */ weights(0) * _Scalar(2);            // add large weight, if no points assigned

                    // add to problem
                    problem.addLinObjective( /* var_id: */ lids_varids.at( IntPair(lid,lid1) )
                                           , /*  value: */ coeff );
                } //...for each direction
            } //...for each patch

            return err;
        } //...associationBasedDataCost

        //! \brief
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
                 , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        everyPatchNeedsDirectionConstraint( _OptProblemT              & problem
                                          , _PrimitiveContainerT const& prims
                                          , _PointContainerT     const& /*points*/
                                          , _AssocT              const& lids_varids
                                          , _WeightsT            const& /*weights*/
                                          , _Scalar              const  /*scale*/ )
        {
            typedef typename _AssocT::key_type IntPair;

            int err = EXIT_SUCCESS;

            // one direction / patch needs to be choosen
            std::vector< double              > coeffs ( problem.getVarCount(), 0 ); // constraint coefficient line in A
            std::set   < std::vector<double> > uniqueA;                             // to ensure unique constraints
            // for all patches
            for ( size_t lid = 0; lid != prims.size(); ++lid )
            {
                // init row to 0
                std::fill( coeffs.begin(), coeffs.end(), 0 );

                // add 1 for each direction patch -> at least one direction has to be chosen for this patch
                for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
                {
                    coeffs[ /* varid: */ lids_varids.at(IntPair(lid,lid1)) ] = 1.0;
                }

                // unique insertion
                unsigned prev_size = uniqueA.size(); // store prev set size to check for unique insert
                uniqueA.insert( coeffs );            // insert into set
                if ( uniqueA.size() > prev_size )    // if line was unique, add to problem
                    problem.addLinConstraint( _OptProblemT::BOUND::GREATER_EQ, 1, problem.getINF(), coeffs ); // 1 <= A( lid, : ) * X <= INF
            } // ... for each patch

            return err;
        } // ...everyPatchNeedsDirectionConstraint

        //! \brief Hybrid mode, large patches need one direction, points in small patches need to be assigned to one.
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
                 , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        largePatchesNeedDirectionConstraint( _OptProblemT              & problem
                                          , _PrimitiveContainerT const& prims
                                          , _PointContainerT     const& /*points*/
                                          , _AssocT              const& lids_varids
                                          , _WeightsT            const& /*weights*/
                                          , _Scalar              const  /*scale*/
                                          , int                  const  pop_limit )
        {
            typedef typename _AssocT::key_type IntPair;

            int err = EXIT_SUCCESS;

            // one direction / patch needs to be choosen
            std::vector< double              > coeffs ( problem.getVarCount(), 0 ); // constraint coefficient line in A
            std::set   < std::vector<double> > uniqueA;                             // to ensure unique constraints
            // for all patches
            for ( size_t lid = 0; lid != prims.size(); ++lid )
            {
                // init row to 0
                std::fill( coeffs.begin(), coeffs.end(), 0 );

                // add 1 for each direction patch -> at least one direction has to be chosen for this patch
                for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
                {
                    coeffs[ /* varid: */ lids_varids.at(IntPair(lid,lid1)) ] = 1.0;
                }

                // unique insertion
                unsigned prev_size = uniqueA.size(); // store prev set size to check for unique insert
                uniqueA.insert( coeffs );            // insert into set
                if ( uniqueA.size() > prev_size )    // if line was unique, add to problem
                    problem.addLinConstraint( _OptProblemT::BOUND::GREATER_EQ, 1, problem.getINF(), coeffs ); // 1 <= A( lid, : ) * X <= INF
            } // ... for each patch

            return err;
        } // ...everyPatchNeedsDirectionConstraint

        //! \brief Nic's version
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
                 , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        instanceBasedDataCost( _OptProblemT              & problem
                                , _PrimitiveContainerT const& prims
                                , _PointContainerT     const& points
                                , _AssocT              const& lids_varids
                                , _WeightsT            const& weights
                                , _Scalar              const  /*scale*/ )
        {
                int err = EXIT_SUCCESS;
                std::cerr << "[" << __func__ << "]: " << "UNIMPLEMENTED" << std::endl;
    //        std::map<int, std::pair<Scalar,int> > repr_dirs; // [lid] = <cost,npoints>
    //        // unary cost -> lin objective
    //        for ( size_t lid = 0; lid != prims.size(); ++lid )
    //        {
    //            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
    //            {
    //                unsigned cnt = 0;
    //                Scalar unary_i = Scalar(0);
    //                for ( size_t pid = 0; pid != points.size(); ++pid )
    //                {
    //                    if ( points[pid].getTag( PointPrimitiveT::GID ) == static_cast<int>(lid) )
    //                    {
    //                        // if within scale, add unary cost
    //                        Scalar dist = PointPrimitiveDistanceFunctor::eval<Scalar>( points[pid], prims[lid][lid1] );
    //                        unary_i += dist;
    //                        ++cnt;
    //                    }
    //                } // for points

    //                if ( prims[lid].getTag(PrimitiveT::GID) != lid )
    //                    std::cerr << "[" << __func__ << "]: " << "a;sldjfa;lkjdf, " << prims[lid].getTag(PrimitiveT::GID) << "lid != tag" << prims[lid].getTag(PrimitiveT::GID) << std::endl;

    //                if ( prims[lid].getTag(PrimitiveT::GID) == prims[lid1].getTag(PointPrimitiveT::GID) )
    //                    repr_dirs[ lid ] = std::pair<Scalar,int>( unary_i, cnt );
    //            }
    //        } // for lines

    //        for ( size_t lid = 0; lid != prims.size(); ++lid )
    //        {
    //            for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
    //            {
    //                Scalar coeff = /* complx: */ weights(2) + /* unary: */ weights(0) * unary_i / Scalar(cnt);
    //                if ( !cnt )
    //                    coeff = weights(2) + weights(0) * 2;

    //                problem.addLinObjective( /* var_id: */ lids_varids[Key(lid,lid1)]
    //                        , coeff );
    //            }
    //        }
            return err;
        } //...instanceBasedDataCost()

        //! \brief Adds unary costs to problem based on scale wide band assocation.
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
                 , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        bandBasedDataCost( _OptProblemT              & problem
                         , _PrimitiveContainerT const& prims
                         , _PointContainerT     const& points
                         , _AssocT              const& lids_varids
                         , _WeightsT            const& weights
                         , _Scalar              const  scale )
        {
            typedef typename _AssocT::key_type IntPair;

            int err = EXIT_SUCCESS;

            for ( size_t lid = 0; lid != prims.size(); ++lid )
            {
                for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
                {
                    // unary cost
                    {
                        unsigned cnt = 0;
                        _Scalar unary_i = _Scalar(0);
                        for ( size_t pid = 0; pid != points.size(); ++pid )
                        {
                            // if within scale, add unary cost
                            _Scalar dist = _PointPrimitiveDistanceFunctor::template eval<_Scalar>( points[pid], prims[lid][lid1] );
                            if ( dist < scale ) //TODO: patch membership
                            {
                                unary_i += dist;
                                ++cnt;
                            }
                                //unary_i += dist / scale;
                        } // for points

                        _Scalar coeff = /* complx: */ weights(2) + /* unary: */ weights(0) * unary_i / _Scalar(cnt);
                        if ( !cnt )
                            coeff = weights(2) + weights(0) * 2;

                        problem.addLinObjective( /* var_id: */ lids_varids.at( IntPair(lid,lid1) )
                                               , coeff );
                    }
                }
            } //...for patches

            return err;
        } //...bandBasedDataCost

        //! \brief
        //! \tparam _AssocT     Associates a primitive identified by <lid,lid1> with a variable id in the problem. Default: std::map< std::pair<int,int>, int >
        template < class _PointPrimitiveDistanceFunctor
                 , class _PrimitiveT        /* = typename _PrimitiveContainerT::value_type::value_type */
                 , class _PointPrimitiveT   /* = typename _PointContainerT::value_type */
                 , typename _Scalar
                 , class _OptProblemT
                 , class _PrimitiveContainerT
                 , class _PointContainerT
                 , class _AssocT
                 , class _WeightsT
                 >
        static inline int
        everyPointNeedsPatchConstraint( _OptProblemT              & problem
                                      , _PrimitiveContainerT const& prims
                                      , _PointContainerT     const& points
                                      , _AssocT              const& lids_varids
                                      , _WeightsT            const& /*weights*/
                                      , _Scalar              const  scale )
        {
            typedef typename _AssocT::key_type IntPair;

            int err = EXIT_SUCCESS;

            std::set< std::vector<double> > uniqueA;
            // each point needs at least one line
            for ( size_t pid = 0; pid != points.size(); ++pid )
            {
                // work
                std::vector<double> coeffs( problem.getVarCount(), 0 );
                bool non_zero = false;
                for ( size_t lid = 0; lid != prims.size(); ++lid )
                    for ( size_t lid1 = 0; lid1 != prims[lid].size(); ++lid1 )
                    {
                        _Scalar dist = _PointPrimitiveDistanceFunctor::template eval<_Scalar>( points[pid], prims[lid][lid1] );
                        if ( dist < scale )
                        {
                            coeffs[ /* mosek.varid: */ lids_varids.at(IntPair(lid,lid1)) ] = 1.0;
                            non_zero = true;
                        }
                    }

                if ( non_zero )
                {
                    unsigned prev_size = uniqueA.size();
                    uniqueA.insert( coeffs );
                    if ( uniqueA.size() > prev_size )
                        problem.addLinConstraint( _OptProblemT::BOUND::GREATER_EQ, 1, problem.getINF(), coeffs );
                }
            } //... for each point

            return err;
        } //...everyPointNeedsPatchConstraint

    } //... namespace ProblemSetup
} //...namespace GF2

#endif // GF2_PROBLEMSETUP_HPP
