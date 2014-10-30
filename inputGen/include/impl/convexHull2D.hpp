#ifndef CONVEXHULL_HPP
#define CONVEXHULL_HPP

// Uncomment only to get code completion in your favourite IDE
#include "convexHull2D.h"
// End Uncomment

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Internals


namespace ConvexHullInternal{

/*!
 * \brief  2D cross product of OA and OB vectors, i.e. z-component of their 3D
 * cross product.
 * \return a positive value, if OAB makes a counter-clockwise turn, negative for
 * clockwise turn, and zero if the points are collinear.
 */
template <class Point>
inline
typename Point::Scalar
cross(const Point &O, const Point &A, const Point &B)
{
    return (A(0) - O(0)) * (B(1) - O(1)) - (A(1) - O(1)) * (B(0) - O(0));
}

/*!
 * \brief Lexicographical point comparison operator
 */
template <class Point>
inline
bool PointCompOp(const Point &p1, const Point &p2) {
    return p1(0) < p2(0) || (p1(0) == p2(0) && p1(1) < p2(1));
}

//! Check if a point p is on the left hand side of the segment [prev, curr]
template <class Point>
inline
bool leftHand(const Point&p, Point prev, Point curr){
    typedef typename Point::Scalar Scalar;

    // first check if p is on the left hand side of prev of curr
    if ( p(0) <= prev(0) || p(0) <= curr(0) ){
        // project the current point horizontally to [prev,curr]
        Scalar alpha = (p(1)-prev(1)) / (curr(1)-prev(1));
        if( alpha > 0 && alpha < 1 ){
            Scalar xproj = prev(0) + (curr(0)-prev(0)) * alpha;
            return (p(0) <= xproj);
        }
        return false;
    }
    else
        return false;
}

} // namespace ConvexHullInternal


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ConvexHull2D Implementation

template <class Point,
          class Container>
void
ConvexHull2D<Point, Container>::compute(Container inputSet)
{
    int n = inputSet.size(), k = 0;
    m_points.clear();
    m_points.resize(2*n);

    // Sort points lexicographically
    sort(inputSet.begin(),
         inputSet.end(),
         ConvexHullInternal::PointCompOp<Point>);

    // Build lower hull
    for (int i = 0; i < n; ++i) {
        while (k >= 2 &&
               ConvexHullInternal::cross(m_points[k-2],
                                         m_points[k-1],
                                         inputSet[i]) <= 0)
            k--;
        m_points[k++] = inputSet[i];
    }

    // Build upper hull
    for (int i = n-2, t = k+1; i >= 0; i--) {
            while (k >= t &&
                   ConvexHullInternal::cross(m_points[k-2],
                                             m_points[k-1],
                                             inputSet[i]) <= 0)
                k--;
        m_points[k++] = inputSet[i];
    }

    m_points.resize(k);
}


template <class Point,
          class Container>
bool
ConvexHull2D<Point, Container>::isInside(const Point &p,
                                                 bool conservative) const
{
    if (size() < 3)
        return true;

    // Count intersection with an horizontal ray going through p and hull
    int count = 0;

    // Check between first and last element
    if( m_points.back() != m_points.front() &&
        ConvexHullInternal::leftHand(p, m_points.back(), m_points.front()))
        count++;

    for (int i = 1; i < m_points.size(); i++){
        if (p == m_points[i]) return conservative;
        if( ConvexHullInternal::leftHand(p, m_points[i-1], m_points[i]))
            count++;
    }

    // we know that the polygon is convex, inside points have only 1 hit
    return count == 1;
}



template <class Point,
          class Container>
typename ConvexHull2D<Point, Container>::Scalar
ConvexHull2D<Point, Container>::distanceTo( const Point& p) const
{
    const int nbUniquePoint = size();
    if (nbUniquePoint == 0)
        return std::numeric_limits<Scalar>::infinity();
    if (nbUniquePoint == 1)
        return - (p-this->at(0)).norm(); // negative because out of the hull

    Scalar dmin = std::numeric_limits<Scalar>::max();

    // Iterate over all segments and measure point to segment distance
    for(unsigned int i = 0; i != nbUniquePoint; i++){
        // We normalize the current segment, to get a projection expressed in
        // [0:1] if we are on the segment. In that case, we compute the distance
        // to the projected point, otherwise the distance to the closest bound,
        // eg, prev if proj<=0, and next if proj>=1.
        const Point &prev = this->at(i);
        const Point &next = this->at(i+1);

        Point seg  = next - prev;
        Point cand = p - prev;
        Scalar proj = (cand / seg.norm()).dot(seg / seg.norm());

        // Point we want the distance with
        Point ref;
        if      (proj <= 0) ref = prev;
        else if (proj >= 1) ref = next;
        else                ref = prev + proj * cand;

        Scalar d = (p-ref).norm();

        if( d < dmin) dmin = d;
    }

    return dmin;
}


template <class Point,
          class Container>
template <bool _computeCentroid>
typename ConvexHull2D<Point, Container>::Scalar
ConvexHull2D<Point, Container>::_areaAndCentroid(Point& centroid) const
{
    const int nbUniquePoint = size();

    // this is duplicated, but not the following
    if (nbUniquePoint < 3){  // not enough points

        if (_computeCentroid) centroid = Point::Zero();

        return Scalar(0);
    }
    if (nbUniquePoint == 3){ // simple case: triangle
        const Point& p =  this->at(0);
        const Point& p1 = this->at(1);
        const Point& p2 = this->at(2);

        if (_computeCentroid) centroid = ( p + p1 + p2 ) / Scalar(3.);

        return (p1-p).norm() * (p1-p).dot(p2-p);
    }

    // Pick a point p inside the hull: we use the unweighted gravity center
    Point p (Point::Zero());
    for(unsigned int i = 0; i != nbUniquePoint; i++)
        p += this->at(i);
    p /= Scalar(nbUniquePoint);

    // For each triangle formed by p and two consecutive points on the hull
    Point  sum  (Point::Zero());
    Scalar sumW = Scalar(0);
    for (unsigned int i = 0; i!= nbUniquePoint; i++){
        const Point& p1 = this->at(i);
        const Point& p2 = this->at(i+1);

        // Compute the centroid and the area and store the weighted sum
        const Scalar w = (p1-p).norm() * (p1-p).dot(p2-p);
        sumW += w;

        if (_computeCentroid) sum  += w * (p + p1 + p2) / Scalar(3.);
    }

    // Normalize and return
    if (_computeCentroid) centroid = sum / sumW;

    return sumW;
}


template <class Point,
          class Container>
Point
ConvexHull2D<Point, Container>::computeCentroid() const
{
    Point centroid;
    _areaAndCentroid<true>(centroid);
    return centroid;
}

#endif // CONVEXHULL_HPP
