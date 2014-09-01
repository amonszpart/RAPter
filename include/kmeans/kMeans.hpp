#ifndef KMEANS_HPP
#define KMEANS_HPP

#include "kMeans.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "primitive.h"
#include "pointPrimitive.h"

#include "amDefines.h"
#include "AMUtil2.h"

namespace am
{
    namespace kmeans
    {

        template <class TPoints, class _3DPoint>
        typename KMPoint<_3DPoint>::Vec
        gen_xy( int count, double radius )
        {
            double ang, r;
            std::vector<KMPoint<_3DPoint> > pts(count);
            //p_point p, pt = (point_t*)malloc(sizeof(point_t) * count);

            /* note: this is not a uniform 2-d distribution */
            std::cout << "count: " << count << "radius: " << radius << std::endl;
            for ( size_t pid = 0; pid != pts.size(); ++pid )
            {
                ang = randf( 2. * M_PI );
                r   = randf( radius );
                MyPoint pnt; pnt.x = r * cos( ang ); pnt.y = r * sin( ang ); pnt.z = 0.f;
                pts[pid].pnt = pnt;
                //pts[pid].pnt.X() = r * cos( ang );
                //pts[pid].pnt.Y() = r * sin( ang );
                std::cout << "r: " << r
                          << " ang: " << ang
                          << " rcos: " << r * cos( ang )
                          << " rsin: " << r * sin( ang )
                          << " pntxy: " << pts[pid].pnt.X() << ", " << pts[pid].pnt.Y()
                          << " pnt: " << pts[pid].pnt.pos().transpose() << std::endl;
            }

            return pts;
        }

        template <class TPoints, class TPoint>
        int
        KMeans<TPoints,TPoint>::nearest_primitive( TPoint                                    const& pt
                                                   , int                                            curr_label
                                                   , typename Primitive<TPoints,TPoint>::Vec const& primitives
                                                   , typename TPoint::Scalar                      * depth_out
                                                   , VectorType                              const* dim_coeffs
                                                   , typename Primitive<TPoints,TPoint>::Vec::const_iterator *end_arg )
        {
            typename Primitive<TPoints,TPoint>::Vec::const_iterator end = end_arg ? *end_arg : primitives.end();

            int    min_i   = curr_label;
            double min_d   = HUGE_VAL, d;
            int    prim_id = 0;
            for ( typename Primitive<TPoints,TPoint>::Vec::const_iterator it = primitives.begin();
                  it != end;
                  ++it, ++prim_id )
            {
                d = (*it)->getDistance( pt, dim_coeffs );
                if ( d < min_d )
                {
                    min_d = d;
                    min_i = prim_id;
                }
            }

            if ( depth_out )   *depth_out = min_d;
            return min_i;
        }

        template <class TPoints, class TPoint>
        void KMeans<TPoints,TPoint>::kpp( typename Primitive<TPoints,TPoint>::Vec      & primitives
                                          , std::vector<int>                           & labels
                                          , TPoints                               const& points
                                          , VectorType                            const* dim_coeffs )
        {
            const ConstIterator points_begin = points.begin();
            const ConstIterator points_end   = points.end();
            const int           points_size  = std::distance( points_begin, points_end );

            Scalar *p_depths = reinterpret_cast<Scalar*>( malloc(sizeof(Scalar) * points_size) );

            const int N = primitives.size();
            //primitives.resize ( 1 );
            //primitives.reserve( N );

            // set first cluster centroid randomly
            {
                ConstIterator it = points_begin;
                std::advance( it, rand() % primitives.size() );
                primitives[0]->asPointPrimitive()->Pnt() = it->pos();
            }

            // all other clusters
            typename Primitive<TPoints,TPoint>::Vec::const_iterator last = primitives.begin(); std::advance( last, 1 );
            for ( size_t primitive_id = 1; primitive_id != N; ++primitive_id, std::advance( last, 1) )
            {
                Scalar sum = 0;
                ConstIterator it = points_begin;
                for ( int pid = 0; pid != points_size; ++pid, ++it )
                {
                    nearest_primitive( *it, labels[pid], primitives, p_depths+pid, dim_coeffs, &last );
                    sum += p_depths[ pid ];
                }

                sum = randf( sum );

                it = points_begin;
                for ( int pid = 0; pid != points_size; ++pid, ++it )
                {
                    if ( (sum -= p_depths[pid]) > 0 ) continue;

                    (*last)->asPointPrimitive()->Pnt() = it->pos();
                    break;
                }
            }

            ConstIterator it = points_begin;
            for ( int pid = 0; pid != points_size; ++pid, ++it )
            {
                labels[pid] = nearest_primitive( *it, labels[pid], primitives, 0, dim_coeffs );
            }

            free( p_depths );
        }

        template <class TPoints, class TPoint>
        int
        KMeans<TPoints,TPoint>::lloyd( typename Primitive<TPoints,TPoint>::Vec      & primitives
                                       , std::vector<int>                           & labels
                                       , TPoints                               const& points
                                       , VectorType                            const* dim_coeffs )
        {
            const ConstIterator points_begin = points.begin();
            const ConstIterator points_end   = points.end();
            const int           points_size  = std::distance( points_begin, points_end );

            int    min_i;
            size_t changed;
            int    max_iterations = 1000;

            do
            {
                for ( size_t prim_id = 0; prim_id != primitives.size(); ++prim_id )
                {
                    primitives[prim_id]->recalculate( points, labels, prim_id );
                }

                /* find closest centroid of each point */
                changed = 0;
                int pid = 0;
                for ( ConstIterator it = points_begin; it != points_end; ++pid, ++it )
                {
                    min_i = nearest_primitive( *it, labels[pid], primitives, NULL, dim_coeffs );
                    if ( min_i != labels[pid] )
                    {
                        labels[pid] = min_i;
                        ++changed;
                    }

                }
            }
            while ( /**/ (changed > (points_size >> 10))    // stop when 99.9% of points are good
                    &&   (max_iterations--             ) );

            return EXIT_SUCCESS;
        }

        template <class TPoints, class TPoint>
        void print_eps( typename KMPoint<TPoint>::Vec const& pts
                        , typename Primitive<TPoints,TPoint>::Vec const& primitives )
        {
            const int W = 400;
            const int H = 400;

            //p_point p, c;
            double min_x, max_x, min_y, max_y, scale, cx, cy;
            double *colors = (double*)malloc(sizeof(double) * primitives.size() * 3);

            //for (c = cent, i = 0; i < n_cluster; i++, c++) {
            for ( size_t i = 0; i < primitives.size(); ++i )
            {
                colors[3*i + 0] = (3 * (i + 1) % 11)/11.;
                colors[3*i + 1] = (7 * i % 11)/11.;
                colors[3*i + 2] = (9 * i % 11)/11.;
            }

            max_x = max_y = -(min_x = min_y = HUGE_VAL);
            for ( size_t pid = 0; pid != pts.size(); ++pid )
            {
                KMPoint<TPoint> const& p = pts[pid];
                if (max_x < p.pnt.X()) max_x = p.pnt.X();
                if (min_x > p.pnt.X()) min_x = p.pnt.X();
                if (max_y < p.pnt.Y()) max_y = p.pnt.Y();
                if (min_y > p.pnt.Y()) min_y = p.pnt.Y();
            }
            scale = W / (max_x - min_x);
            if (scale > H / (max_y - min_y)) scale = H / (max_y - min_y);
            cx = (max_x + min_x) / 2;
            cy = (max_y + min_y) / 2;

            printf("%%!PS-Adobe-3.0\n%%%%BoundingBox: -5 -5 %d %d\n", W + 10, H + 10);
            printf( "/l {rlineto} def /m {rmoveto} def\n"
                    "/c { .25 sub exch .25 sub exch .5 0 360 arc fill } def\n"
                    "/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath "
                    "	gsave 1 setgray fill grestore gsave 3 setlinewidth"
                    " 1 setgray stroke grestore 0 setgray stroke }def\n"
                    );
            //for (c = cent, i = 0; i < n_cluster; i++, c++)
            for ( size_t primitive_id = 0; primitive_id != primitives.size(); ++primitive_id )
            {
                printf( "%g %g %g setrgbcolor\n",
                        colors[3*primitive_id    ],
                        colors[3*primitive_id + 1],
                        colors[3*primitive_id + 2] );
                //for ( j = 0, p = pts; j < len; j++, p++ )
                for ( size_t pid = 0; pid != pts.size(); ++pid )
                {
                    KMPoint<TPoint> const& p = pts[pid];

                    if ( (size_t)p.primitive_id != primitive_id ) continue;

                    printf("%.3f %.3f c\n",
                           (p.pnt.pos()(0) - cx) * scale + W / 2,
                           (p.pnt.pos()(1) - cy) * scale + H / 2);
                }
                printf("\n0 setgray %g %g s\n",
                       (reinterpret_cast<PointPrimitive<TPoints,TPoint>*>(primitives[primitive_id])->Pnt()(0) - cx) * scale + W / 2,
                       (reinterpret_cast<PointPrimitive<TPoints,TPoint>*>(primitives[primitive_id])->Pnt()(1) - cy) * scale + H / 2);
            }
            printf("\n%%%%EOF");
            free(colors);
        }

        inline void
        indent ( std::ofstream &file, int depth, std::string character = "  " )
        {
            MY_ASSERT( depth >= 0 );

            for ( int i = 0; i != depth; ++i )
                file << character;
        }

        template <class VectorType>
        inline void
        printNode( std::ofstream & f
                   , int         & depth
                   , int           pid
                   , VectorType    pnt
                   , std::string   color   = ""
                                             , bool          is_last = false
                                                                       , std::string   label   = ""
                                                                                                 , int           size    = 1 )
        {
            indent(f,depth++); f << "{\n";
            indent(f,depth); f << "\"id\": \"n" << pid << "\",\n";
            indent(f,depth); f << "\"label\": \""; if ( label.empty() ) f << pid; else f << label; f << "\",\n";
            indent(f,depth); f << "\"x\": \"" << pnt(1) << "\",\n";
            indent(f,depth); f << "\"y\": \"" << -1.f*pnt(0) << "\",\n";
            if ( !color.empty() ) { indent(f,depth); f << "\"color\": \"" << color << "\",\n"; }
            indent(f,depth); f << "\"size\": \"" << size << "\"\n";
            indent(f,--depth); f << "}" << (is_last?"":",") << "\n";
        }

        inline void
        printEdge( std::ofstream & f, int & depth
                   , int pid0, int pid1
                   , std::string color = ""
                                         , bool is_last = false
                                                          )
        {
            indent(f,depth++); f << "{\n";
            indent(f,depth); f << "\"id\": \"e" << pid0 << "_" << pid1 << "\",\n";
            if ( !color.empty() ) { indent(f,depth); f << "\"color\": \"" << color << "\",\n"; }
            indent(f,depth); f << "\"source\": \"n" << pid0 << "\",\n";
            indent(f,depth); f << "\"target\": \"n" << pid1 << "\"\n";
            indent(f,--depth); f << "}" << (is_last?"":",") << "\n";
        }

        template <class TPoints, class TPoint>
        int
        KMeans<TPoints, TPoint>::toJson( std::string           path
                                         , TPoint     const* p_points
                                         , int                 points_size
                                         , std::vector<int>                              & labels
                                         , typename Primitive<TPoints,TPoint>::Vec   const& primitives
                                         , std::vector<Eigen::Vector3f>             const& colours
                                         )
        {
            typedef typename TPoint::VectorType VectorType;
            std::vector<VectorType> clusters( primitives.size() );
            for ( int k = 0; k != primitives.size(); ++k )
                clusters[k] = primitives[k]->asPointPrimitive()->Pnt();

            VectorType scale;
            scale << 2.f, 1.f;

            std::ofstream f;
            f.open( path );
            int depth = 0;

            int pid = 0;
            indent(f,depth++); f << "{\n";
            // print nodes
            {
                indent(f,depth++); f << "\"nodes\": [\n";
                // print points
                for ( ; pid != points_size; ++pid )
                {
                    char col[255];
                    sprintf( col, "#%s%s%s", am::util::toHexString(colours[labels[pid]](0)).c_str()
                            , am::util::toHexString(colours[labels[pid]](1)).c_str()
                            , am::util::toHexString(colours[labels[pid]](2)).c_str() );
                    printNode<VectorType>( f,depth,pid,p_points[pid].pos().array() * scale.array(), col, false );
                }

                // print clusters
                for ( ; pid != points_size+clusters.size(); ++pid )
                {
                    char col[255];
                    sprintf( col, "#%s%s%s"
                             , am::util::toHexString(colours[pid-points_size](0)*0.2).c_str()
                            , am::util::toHexString(colours[pid-points_size](1)*0.2).c_str()
                            , am::util::toHexString(colours[pid-points_size](2)*0.2).c_str() );
                    printNode<VectorType>( f,depth,pid,clusters[pid-points_size].array() * scale.array(), col, false, am::util::sprintf("C%d",pid-points_size), 2 );
                }

                // print axis nodes
                {
                    VectorType origin( VectorType::Zero() );
                    VectorType y_end( VectorType::Zero() ), x_end( VectorType::Zero() );
                    x_end(1) = M_PI * 1.1f; y_end(0) = 2.f;

                    printNode<VectorType>( f, depth, pid++, origin, "#000", false, "(0,0)" );
                    printNode<VectorType>( f, depth, pid++, x_end.array() * scale.array(), "#000", false, am::util::sprintf("%2.2f",x_end(1)) );
                    printNode<VectorType>( f, depth, pid++, y_end.array() * scale.array(), "#000",true , am::util::sprintf("%2.2f",y_end(0)) );
                }

                indent(f,--depth); f << "],\n";
            }

            // print edges
            {
                indent(f,depth++); f << "\"edges\": [\n";

                for ( int l_pid = 0; l_pid < points_size; ++l_pid )
                {
                    printEdge( f, depth, l_pid, points_size+labels[l_pid], "#000", false );
                }

                // axees
                printEdge( f, depth, pid-3, pid-2, "#000", false );
                printEdge( f, depth, pid-3, pid-1, "#000", true );

                indent(f,--depth); f << "]\n";
            }

            indent(f,--depth); f << "}\n";

            f.close();

            return EXIT_SUCCESS;
        }

        template <class TPoints, class TPoint>
        int
        KMeans<TPoints,TPoint>::cluster(
                typename Primitive<TPoints,TPoint>::Vec        & primitives
                , std::vector<int>                             & labels
                , TPoints                                 const& points
                , VectorType                              const* dim_coeffs
                , bool                                           kmeans_pp_primitives
                )
        {
            const int points_size = std::distance( points.begin(), points.end() );

            const int K = primitives.size();
            MY_ASSERT( K > 0 );

            if ( labels.size() != points_size )
                labels.resize( points_size, 0 );

            if ( kmeans_pp_primitives )
            {
                // init primitives and labels
                kpp( primitives, labels, points, dim_coeffs );
            }

            // debug
            std::cout << "BEFORE: " << std::endl;
            for ( int k = 0; k != K; ++k )
            {
                std::cout << "\tprimitive " << k << ": "
                          << (primitives)[k]->asPointPrimitive()->Pnt().transpose() << std::endl;
            }

            // work
            lloyd( primitives, labels, points, dim_coeffs );

            std::cout << "AFTER: " << std::endl;
            for ( int k = 0; k != K; ++k )
            {
                std::cout << "\tprimitive " << k << ": "
                          << (primitives)[k]->asPointPrimitive()->Pnt().transpose() << std::endl;
            }

            return EXIT_SUCCESS;
        }

    } // ns kmeans

} // ns am

#endif // KMEANS_HPP
