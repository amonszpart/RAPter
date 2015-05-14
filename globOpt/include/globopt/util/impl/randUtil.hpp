#ifndef GO_RANDUTIL_HPP
#define GO_RANDUTIL_HPP

namespace globopt
{
    template <typename _Scalar>
    inline _Scalar randf() { return rand() / _Scalar(RAND_MAX); }

    template <typename _Scalar>
    inline _Scalar randf(_Scalar scale) { return scale * rand() / _Scalar(RAND_MAX); }
}

#endif // GO_RANDUTIL_HPP
