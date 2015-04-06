#ifndef __AM__PRIMITIVE_HPP__
#define __AM__PRIMITIVE_HPP__

#include <vector>
#include "Eigen/Dense"

#include <iostream> // cout, endl
#include <fstream>

#include "globopt/primitives/taggable.h"

namespace GF2
{
    template <class StoredT>
    struct CachedField
    {
        class CachedFieldException : public std::runtime_error
        {
            public:
                explicit CachedFieldException( const std::string& __arg )
                    : std::runtime_error( __arg )
                {
                    std::cerr << "[CachedFieldException]: " << __arg << std::endl; fflush(stderr);
                }
        };

        CachedField()
            : _updated(false)
        {}

        void update( StoredT const& content )
        {
            _content = content;
            _updated = true;
        }

        bool isUpdated() const { return _updated; }
        void outdate  ()       { _updated = false; }

        StoredT const& get() const
        {
            if ( !_updated )
                throw new CachedFieldException("Cached field outdated");

            return _content;
        }

        protected:
            bool    _updated;
            StoredT _content;
    };

     /*! \brief                 Primitive wrapper base class. Requires implementation of #pos() and #dir() functions.
      *                         It's also good to implement constructor from default input (i.e. from point and normal for planes)
      *  \tparam _EmbedSpaceDim Hack to know, if the embedding space is 2D or 3D. TODO: don't use smartgeometry::fitlinearprimitve in \ref Segment::fitLocal()
      *  \tparam _Dim           Length of internal vector storage. Should be able to fit any primitive.
      *  \tparam _Scalar        Internal data type. Concept: float.
      */
    template <int _EmbedSpaceDim, int _Dim, typename _Scalar = float>
    class Primitive : public ::globopt::Taggable<_Scalar>
    {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            typedef ::globopt::Taggable<_Scalar> TaggableT;
            typedef TaggableT                ParentT;
            enum { EmbedSpaceDim = _EmbedSpaceDim };
#if 0
            enum STATUS_VALUES { ACTIVE = 1 //!< Used to show, that this primitive is active. Also used as starting point in \ref GF2::Solver.
                               , SMALL  = 2 //!< Used to show, that this primitive is *not* active, but needs to be kept for later use, so copied to output everywhere.
                               //, UNSET  = ParentT::TAG_UNSET
                               , INVALID = -2 //!< Used to show, that the status value has not been read yet.
                               };
#endif
            //! A weak attempt to get strongly typed enums. Won't work because of UNSET us used across types.
            struct STATUS_VALUES
            {
                static const char ACTIVE  = 1; //!< Used to show, that this primitive is active. Also used as starting point in \ref GF2::Solver.
                static const char SMALL   = 2; //!< Used to show, that this primitive is *not* active, but needs to be kept for later use, so copied to output everywhere.
                static const char UNSET   = ParentT::TAG_UNSET;
                static const char INVALID = -2; //!< Used to show, that the status value has not been read yet.
                static const char FIXED   = 3;
            };

            //! \brief Defines the tags (ids) that this primitive can manage using setTag and getTag functions.
            struct TAGS {
                static const GidT    GID       = 0;  //!< group id             - which group this primitive is supposed to explain
                static const DidT    DIR_GID   = 1;  //!< direction group id   - which group this primitive got it's direction from
                static const char    STATUS    = 2;  //!< an additional flag to store, if this is part of a solution.
                static const _Scalar GEN_ANGLE/* = 3*/;  //!< flag to show which angle the direction is commited to
            };//...TAGS

            struct GEN_ANGLE_VALUES {
                static const _Scalar UNSET /*= ParentT::TAG_UNSET*/;
            };

            struct LONG_VALUES {
                static const GidT UNSET /*= ParentT::TAG_UNSET*/;
            };

            // ____________________TYPEDEFS____________________
            typedef _Scalar                      Scalar;     //!< Scalar typedef to reach from outside.
            enum {  Dim                        = _Dim };     //!< Standard typedef for _Dim template parameter.
            typedef Eigen::Matrix<Scalar,_Dim,1> VectorType; //!< Standard typedef of internal storage Eigen::Matrix.
            typedef Eigen::Matrix<Scalar,3,1>    Position;   //!< \brief Standard 3D point typedef.
            typedef std::vector<Eigen::Matrix<Scalar,3,1> > ExtremaT; //!< Stores extrema of finite primitive.

            // ____________________CONSTRUCTORS____________________
            //! \brief Default constructor, initializing coeffs to zeros.
            Primitive()
                : _coeffs( Eigen::Matrix<Scalar,_Dim,1>::Zero()  )
            {}

            //! \brief Constructor that takes raw data in Eigen format as input.
            Primitive( Eigen::Matrix<Scalar,_Dim,1> coeffs )
                : _coeffs( coeffs )
            {}

            //! \brief Constructor that takes raw data in std::vector format as input.
            Primitive( std::vector<Scalar> const& coeffs )
                : _coeffs( coeffs.size() )
            {
                std::copy( coeffs.data(), coeffs.data() + coeffs.size(), _coeffs.data() );
            }

            // ____________________GETTERS____________________
            //! \brief Const reference getter for internal storage. Used from pcl::SampleConsensusModelLine and pcl::SampleConsensusModelPlane.
            Eigen::Matrix<Scalar,_Dim,1> const& coeffs    () const { return _coeffs; }
            //! \brief Reference getter for internal storage. Used from pcl::SampleConsensusModelLine and pcl::SampleConsensusModelPlane.
            Eigen::Matrix<Scalar,_Dim,1>      & coeffs    ()       { return _coeffs; }
            //! \brief Const reference getter alias for internal storage.
            Eigen::Matrix<Scalar,_Dim,1> const& operator()() const { return _coeffs; }
            //! \brief Reference getter alias for internal storage.
            Eigen::Matrix<Scalar,_Dim,1>      & operator()()       { return _coeffs; }
            /*! \brief Converts to Eigen::Matrix<...> to avoid operator() as getter.
             * \return Copy of internally stored data.
             */
            /*explicit */ operator VectorType()          { return _coeffs; }
            /*! \brief Converts to Eigen::Matrix<...> to avoid operator() as getter. Const version.
             *  \return Copy of internally stored data.
             */
            /*explicit */ operator VectorType() const    { return _coeffs; }

            // ____________________VIRTUALS____________________
            //! \brief  Pure virtual function to return the location of the primitive.
            //! \return Location of the primitive as a 3D Eigen::Vector.
            virtual Eigen::Matrix<Scalar,3,1>   pos() const = 0;
            //! \brief  Pure virtual function to return the direction of the primitive.
            //! \return Direction of the primitive as a 3D Eigen::Vector.
            virtual Eigen::Matrix<Scalar,3,1>   dir() const = 0;

            // ____________________CONSOLE____________________
            //! \brief Convenience function for plotting.
            inline std::string toString() const
            {
                std::stringstream ss;
                for ( size_t d = 0; d != _Dim; ++d )
                    ss << _coeffs(d) << ", ";

                return ss.str();
            } //...toString()
            inline ExtremaT const& getExtentCached() const { return _extents.get(); }
            inline void setExtentOutdated() const { _extents.outdate(); }

        protected:
            // ____________________FIELDS____________________
            Eigen::Matrix<Scalar,_Dim,1>    _coeffs;  //!< \brief Internal storage of arbitrary length vector representing the primitive.
            mutable CachedField<ExtremaT>   _extents; //!< \brief Stores extents of finite primitive.

            // ____________________DEPRECATED____________________
#if 0
        public:
            //! \deprecated Convenience getter for older implementations, should be deprecated.
            std::vector<Scalar>
            coeffsVector() const
            {
                std::vector<Scalar> coeffs( _coeffs.rows() );
                std::copy( _coeffs.data(), _coeffs.data()+_coeffs.rows(), coeffs.begin() );

                return coeffs;
            }

            // FILE operations
            //! \deprecated Saves primitive to file. Deprecated, since it does not save tags.
            //! \todo The functionality from io::writePrimitives could be moved here, if taggable would have a dump function as well.
            template <class PrimitiveT> static int
            dump( std::string                   const& path
                  , std::vector<PrimitiveT>     const& lines
                  , std::vector<int>            const* lines_clusters
                  , std::vector<int>            const* mask
                  )
            {
                std::ofstream f( path );
                if ( !f.is_open() )
                {
                    std::cerr << "[" << __func__ << "] couldn't open file" << std::endl;
                    return EXIT_FAILURE;
                }

                for ( size_t i = 0; i != lines.size(); ++i )
                {
                    if ( mask && !(*mask)[i] ) continue;

                    std::vector<Scalar> coeffsVector = lines[i].coeffsVector();
                    for ( size_t d = 0; d != coeffsVector.size(); ++d )
                        f << coeffsVector[d] << ((d == coeffsVector.size()-1) ? "" : ", "); // << ((d!=coeffsVector.size()-1)?", ":"");
                    if ( lines_clusters ) f << ", " << (*lines_clusters)[i];
                    f << "\n";
                }

                f.close();

                return EXIT_SUCCESS;
            }

            //! \deprecated Reads primitive from file. Deprecated, since it does not read tags.
            //! \todo The functionality from io::readPrimitives could be moved here, if taggable would have a dump function as well.
            template <class PrimitiveT> static std::vector<PrimitiveT>
            read( std::string const& path
                  , std::vector<int> *lines_clusters = NULL )
            {
                std::vector<PrimitiveT> lines;

                std::ifstream f( path );
                if ( !f.is_open() )
                {
                    std::cerr << "[" << __func__ << "] couldn't open file" << std::endl;
                    return lines;
                }

                std::string line;
                while ( getline(f, line) )
                {
                    std::vector<Scalar> floats;
                    std::istringstream iss( line );
                    std::string        tmp_str;
                    while ( std::getline(iss, tmp_str, ',') )
                    {
                        floats.push_back( atof(tmp_str.c_str()) );
                    }

                    if ( lines_clusters )
                    {
                        if ( (lines.size() == 0) && (floats.size() != _Dim+1) )
                        {
                            std::cerr << "[" << __func__ << "]: " << "floats.size() != " << _Dim+1 << "...setting cluster to 0\\n";
                            //return lines;
                            lines_clusters->push_back( 0 );
                        }
                        else
                            lines_clusters->push_back( floats.back() );
                    }

                    if ( floats.size() > _Dim ) floats.resize( _Dim );

                    lines.emplace_back( PrimitiveT(floats) );
                }

                return lines;
            } //...read()
#endif
    }; //...Primitive

    template <class _Primitive>
    class PrimitiveVectorT : public std::vector< std::vector<_Primitive> >
    {
        public:
            typedef _Primitive                                   PrimitiveT;
            typedef std::vector<PrimitiveT>                      InnerPrimitiveContainerT;
            typedef std::vector<InnerPrimitiveContainerT>        ParentT;

            using ParentT::ParentT;
    };

//    template <class _Primitive>
//    class PrimitiveMapT : public std::map< PidT, std::vector<_Primitive> >
//    {
//        public:
//            typedef _Primitive                                   PrimitiveT;
//            typedef std::vector<PrimitiveT>                      InnerPrimitiveContainerT;
//            typedef std::map   < PidT, InnerPrimitiveContainerT> ParentT;

//            using ParentT::ParentT;
//    };


    template <int _EmbedSpaceDim, int _Dim, typename _Scalar>
    const _Scalar Primitive<_EmbedSpaceDim, _Dim, _Scalar>::TAGS::GEN_ANGLE = _Scalar( 3. );
    template <int _EmbedSpaceDim, int _Dim, typename _Scalar>
    const _Scalar Primitive<_EmbedSpaceDim, _Dim, _Scalar>::GEN_ANGLE_VALUES::UNSET = Primitive<_EmbedSpaceDim, _Dim, _Scalar>::ParentT::TAG_UNSET;
    template <int _EmbedSpaceDim, int _Dim, typename _Scalar>
    const GidT Primitive<_EmbedSpaceDim, _Dim, _Scalar>::LONG_VALUES::UNSET = Primitive<_EmbedSpaceDim, _Dim, _Scalar>::ParentT::TAG_UNSET;
} //...ns GF2

#endif // __AM__PRIMITIVE_HPP__
