#ifndef GF2_PCL_UTIL_HPP
#define GF2_PCL_UTIL_HPP

#ifdef GF2_USE_PCL

namespace GF2
{
namespace pclutil
{

template <class PointAllocatorFunctorT, class PointContainerT, class CloudPtrT> inline int
cloudToVector( PointContainerT &container, CloudPtrT const& cloud )
{
    const size_t count = cloud->size();
    container.reserve( count );
    for ( size_t pid = 0; pid != count; ++pid )
    {
        container.push_back( PointAllocatorFunctorT::eval(cloud->at(pid).getVector3fMap()) );
    }

    return EXIT_SUCCESS;
}

} //...ns pclutil
} //...ns GF2

#endif //GF2_USE_PCL

#endif // GF2_PCL_UTIL_HPP
