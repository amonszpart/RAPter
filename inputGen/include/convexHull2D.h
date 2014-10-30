#ifndef CONVEXHULL2D_H
#define CONVEXHULL2D_H

/*!
 * \brief 2D convex hull and utility functions stored in a STL-compatible
 * container
 *
 */
template <class _Point,      //!< \brief 2D Vector type
          class Container = std::vector<_Point> > //!< \brief STL-compatible container
class ConvexHull2D{
public:
    typedef typename _Point::Scalar Scalar;
    typedef _Point Point;

private:
    /*!
     * \brief Internal function to compute the area and the centroid
     *
     * Use a template option to let the compiler optimize when we don't want to
     * compute the centroid, so we can have 1 single code for the two operations
     * and still better performances when we only need to get the area.
     */
    template <bool _computeCentroid>
    Scalar _areaAndCentroid(Point& centroid) const;

public:
    inline ConvexHull2D(const Container& inputSet)
    { compute(inputSet); }


    /*!
     * Returns a list of points on the convex hull in counter-clockwise order.
     * \note The last point in the returned list is the same as the first one.
     *
     * \note The container is duplicated for processing, that why we don't use
     *       a reference.
     *
     * Source:
     * http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
     */
    inline void compute(Container inputSet);

    /*!
     * \brief Inclusion test for an arbitrary point p
     * \param p Evaluation point
     * \param conservative Indicates if the points on the hull must be
     * explicitly included or not
     *
     * \see distanceTo
     */
    inline bool isInside( const Point& p, bool conservative = true) const;

    /*!
     * \brief Compute the unsigned distance to the hull
     * \param p Evaluation point
     * \return Euclidean distance to the closest point on the hull
     */
    inline Scalar distanceTo( const Point& p) const;

    /*!
     * \brief Compute the centroid of the polygon
     *
     * The centroid is computed by:
     *   - Pick a point p inside the hull
     *   - For each triangle formed by p and two consecutive points on the hull
     *      - Compute the centroid and the area and store the weighted sum
     *   - Normalize and return
     *
     */
    inline Point computeCentroid() const;

    /*!
     * \brief Compute the area of the hull
     * \return
     */
    inline Scalar area() const
    {
        Point p;
        return _areaAndCentroid<false>(p);
    }

    inline const Point& at(int i) const { return m_points.at(i); }

    /*!
     * \internal The first point is duplicated at the end of the list
     */
    inline int size() const { return m_points.size()-1; }

private:
    Container m_points; //!< \brief Internal hull storage
};

#include "impl/convexHull2D.hpp"

#endif // CONVEXHULL2D_H
