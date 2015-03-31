#ifndef __GF2_POINTPRIMITIVE_H__
#define __GF2_POINTPRIMITIVE_H__

#include "globfit2/primitives/primitive.h"
#include "globfit2/simple_types.h"
//#include "globfit2/primitives/taggable.h"

namespace GF2
{
    //! \brief Wrapper class for oriented 3D points.
    class PointPrimitive : public ::GF2::Primitive<3,6>
    {
        typedef ::GF2::Primitive<3,6> ParentT;
        public:
            //! \brief custom tags for taggable
            struct TAGS {
                static const PidT PID  = 0;   //!< \brief Point ID \warning Currently unused.
                // WARNING: don't ever change GID=1, historical data won't work!
                static const GidT GID  = 1;   //!< \brief Group ID. Should match LinePrimitive::GID and PlanePrimitive::GID (not the actual enum value, but PointPrimitive.getTag(PointPrimitive::GID) == LinePrimitive.getTag(LinePrimitive::GID)).
                static const LidT LID  = 2;   //!< \brief Primitive ID \warning Currently UNUSED, GID is used instead.
                static const LidT LID0 = 3;  //!< \brief Primitive ID \warning Currently UNUSED, GID is used instead.
            }; //...enum TAGS

            typedef ParentT::Scalar Scalar;

            // ____________________CONSTRUCTORS____________________
#if 0 //__cplusplus > 199711L
            //! \brief Inherited constructors from parent class.
            using ::GF2::Primitive<Dim>::Primitive;
#else
            PointPrimitive() : ParentT() {}

            //! \brief Constructor that takes raw data in Eigen format as input.
            PointPrimitive( Eigen::Matrix<Scalar,Dim,1> coeffs ) : ParentT( coeffs ) {}

            //! \brief Constructor that takes raw data in std::vector format as input.
            PointPrimitive( std::vector<Scalar> const& coeffs ) : ParentT( coeffs ) {}
#endif
            //! \brief          Creates PointPrimitive from point position and direction.
            //! \param[in] p0   Point position.
            //! \param[in] dir  Point orientation.
            PointPrimitive( Eigen::Matrix<Scalar,3,1> const& p0 )
            {
                _coeffs.head   <3>( ) = p0;
                _coeffs.segment<3>(3) = Scalar(-1.) * Eigen::Matrix<Scalar,3,1>::Ones();
            }

            //! \brief          Creates PointPrimitive from point position and direction.
            //! \param[in] p0   Point position.
            //! \param[in] dir  Point orientation.
            PointPrimitive( Eigen::Matrix<Scalar,3,1> const& p0, Eigen::Matrix<Scalar,3,1> const& dir )
            {
                _coeffs.head   <3>( ) = p0;
                _coeffs.segment<3>(3) = dir.normalized();
            }

            //! \brief ::GF2::Primitive<Dim>::operator() are inherited convenience getters from parent class. \todo Use explicit operator VectorType() instead.
            using ParentT::operator();

            //! \brief Constructor functor. Used in Solver::show(), should be deprecated.
            struct Allocator
            {
                //! \param[in] pnt  3D point to store.
                //! \return         A PointPrimitive object containing the \p pnt and a default {1,1,1} orientation.
                static inline PointPrimitive eval( Eigen::Matrix<Scalar,3,1> pnt ) { return PointPrimitive( (VectorType() << pnt, 1, 1, 1).finished() ); }
            };
            //! \brief Constructor functor. Used in Solver::show(), should be deprecated.
            struct RawAllocator
            {
                //! \param[in] pnt  3D point to store.
                //! \return         An Eigen::Vector object containing the \p pnt and a default {1,1,1} orientation.
                static inline VectorType     eval( Eigen::Matrix<Scalar,3,1> pnt ) { return (VectorType() << pnt, 1, 1, 1).finished(); }
            };

            // ____________________VIRTUALS____________________
            //! \brief  Compulsory virtual overload of position getter. The position of the point is the location stored at the first three coordinates of #_coeffs.
            //! \return The position of the point as a 3D Eigen::Vector.
            virtual Eigen::Matrix<Scalar,3,1> pos() const
            {
                return _coeffs.template head<3>();
            }

            //! \brief  Compulsory virtual overload of orientation getter. The orientation of the point is the direction stored at the second three coordinates of #_coeffs.
            //! \return The orientation of the point as a 3D Eigen::Vector.
            virtual Eigen::Matrix<Scalar,3,1> dir() const
            {
                return _coeffs.template segment<3>(3) ;
            }

            // ____________________GETTERS____________________
            //! \brief  Additional convenience getter to convert oriented point to 3D point.
            //! \return The position of the point as a 3D Eigen::Vector.
            /*explicit*/ inline operator Eigen::Matrix<Scalar,3,1>() const { return this->pos(); }
            //! \brief  Additional convenience const getter to convert oriented point to 3D point.
            //! \return The position of the point as a 3D Eigen::Vector.
            /*explicit*/ inline operator Eigen::Matrix<Scalar,3,1>()       { return this->pos(); }

            // ____________________STATICS____________________
            //! \brief                      Convenience conversion of an STL container of this class' objects to another container using the provided allocator functor and the push_back() function of \p cloud.
            //! \tparam CloudPtrT           A pointer to a container that implements reserve and push_back. Concept: pcl::PointCloud<pcl::PointXYZ>::Ptr.
            //! \tparam PointsContainerT    An STL container holding instances of this (PointPrimitive) class. Concept: std::vector<PointPrimitive>.
            //! \tparam AllocatorFunctorT   A functor implementing a create function that takes an Eigen::Vector and returns an instance that can be pushed back to the \p cloud. Concept: PCLPointAllocator<PointPrimitive::Dim>.
            //! \param[out] cloud           Output cloud filled from the input.
            //! \param[in]  points          Input stl container of this class' instances.
            //! \return                     EXIT_SUCCESS
            template <class CloudPtrT, class PointsContainerT, class AllocatorFunctorT>
            static inline int toCloud( CloudPtrT               & cloud
                                       , PointsContainerT const& points )
            {
                typedef typename PointsContainerT::const_iterator const_iterator;
                typedef typename CloudPtrT::element_type::PointType CloudPointT;
                typedef typename PointsContainerT::value_type::VectorType VectorType;

                if ( !cloud )
                {
                    std::cerr << "[" << __func__ << "]: " << "cloud needs to be initialized...exiting\n";
                    return EXIT_FAILURE;
                }

                cloud->reserve( points.size() );
                const const_iterator end_it = points.end();
                for ( const_iterator it = points.begin(); it != end_it; ++it )
                {
                    cloud->push_back(
                                AllocatorFunctorT::template create<CloudPointT>( (*it)() )
                                );
                }

                return EXIT_SUCCESS;
            } //...toCloud()

    }; //...class PointPrimitive

    class PointPrimitiveVector : public std::vector<PointPrimitive>
    {
        public:
            typedef std::vector<PointPrimitive> ParentT;
            typedef PointPrimitive PrimitiveT;
    };

} //...ns GF2

#endif // __GF2_POINTPRIMITIVE_H__
