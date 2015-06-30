#ifndef __RAPTER_PS_HPP__
#define __RAPTER_PS_HPP__

#include <fstream>
#include <string>
#include "rapter/util/containers.hpp"
#include "rapter/processing/util.hpp"
#include "rapter/util/util.hpp"

namespace rapter
{
namespace io
{
    inline std::string didName( DidT did )
    {
        char name[128]; sprintf(name, "\"$\\\\chi_{%ld}$\"", did );
        return name;
    }

    inline std::string getPrimName( GidT gid, DidT did )
    {
        char name[128];
        sprintf(name, "\"P%ld,%ld\"", gid, did );
        return name;
    }

    template <class _PrimitiveMapT>
    inline void getColours( std::map< DidT, Eigen::Vector3f> &colourmap, _PrimitiveMapT & prims )
    {
        typedef typename _PrimitiveMapT::mapped_type::value_type PrimitiveT;
        std::vector<Eigen::Vector3f> colours = util::paletteDarkColoursEigen();

        size_t cIndex = 0;
        for ( typename containers::PrimitiveContainer<PrimitiveT>::ConstIterator it(prims); it.hasNext(); it.step() )
        {
            if ( it->getTag(PrimitiveT::TAGS::STATUS) == PrimitiveT::STATUS_VALUES::SMALL ) continue;

            if ( colourmap.find(it.getDid()) == colourmap.end() )
                colourmap[ it.getDid() ] = colours[cIndex++];

            if ( cIndex >= colours.size() )
                cIndex = 0;
        }
    }

    //template <class DerivedT>
    inline void drawLine( FILE *fp, Eigen::Vector2f const& a, Eigen::Vector2f const& b
                        , Eigen::Vector3f colour = Eigen::Vector3f::Zero()
                        , Eigen::Vector3f center = Eigen::Vector3f::Zero() )
    {
        fprintf( fp, "newpath\n" );
        fprintf( fp, "%f %f %f setrgbcolor\n", colour(0)
                                             , colour(1)
                                             , colour(2) );
        std::vector<Eigen::Vector2f> locs(2);
        locs[0] << a(0)*200.+center(0), a(1)*200.+center(1);
        locs[1] << b(0)*200.+center(0), b(1)*200.+center(1);

        fprintf( fp, ".5 setlinewidth\n" );
        fprintf( fp, "%.3f %.3f moveto %.3f %.3f lineto \nstroke \n"
               , locs[0](0), locs[0](1)
               , locs[1](0), locs[1](1) );
    }

    inline void drawCircle( FILE *fp, Eigen::Vector2f const& a, float radius
                            , Eigen::Vector3f colour = Eigen::Vector3f::Zero()
                            , Eigen::Vector3f center = Eigen::Vector3f::Zero() )
    {
        fprintf( fp, "newpath\n" );
        fprintf( fp, "%f %f %f setrgbcolor\n", colour(0)
                                             , colour(1)
                                             , colour(2) );
        std::vector<Eigen::Vector2f> locs(1);
        locs[0] << a(0)*200.+center(0), a(1)*200.+center(1);

        //fprintf( fp, ".5 setlinewidth\n" );
        fprintf( fp, "%.3f %.3f %.3f 0 360 arc closepath\n"
               , locs[0](0), locs[0](1), radius
               );
        fprintf( fp, "%f %f %f setrgbcolor fill\n", colour(0)
                                             , colour(1)
                                             , colour(2) );
        fprintf( fp, "stroke\n" );
    }

    inline void drawFrame( FILE *fp, Eigen::Vector3f center = Eigen::Vector3f::Zero() )
    {
        drawLine( fp
                , (Eigen::Vector2f() << 0., 0.).finished()
                , (Eigen::Vector2f() << 0., 1.).finished()
                , Eigen::Vector3f::Zero()
                , center
                );
        drawLine( fp
                , (Eigen::Vector2f() << 0., 1.).finished()
                , (Eigen::Vector2f() << 1., 1.).finished()
                , Eigen::Vector3f::Zero()
                , center
                );
        drawLine( fp
                , (Eigen::Vector2f() << 1., 1.).finished()
                , (Eigen::Vector2f() << 1., 0.).finished()
                , Eigen::Vector3f::Zero()
                , center
                );
        drawLine( fp
                , (Eigen::Vector2f() << 1., 0.).finished()
                , (Eigen::Vector2f() << 0., 0.).finished()
                , Eigen::Vector3f::Zero()
                , center
                );
    }

    template <class _PrimitiveMapT, class _PointContainerT>
    inline void drawPs( _PrimitiveMapT & prims, _PointContainerT const& points, std::string path, float scale
                      , bool show
                      , bool writeNames = false
                      , bool colourCloud = false
                      , std::vector<DidT> *dids = NULL
                      , float subSample = 1.f
                      , float radius = .1f
                      )
    {
        typedef typename _PrimitiveMapT::mapped_type::value_type PrimitiveT;
        typedef typename PrimitiveT::Scalar Scalar;
        typedef typename _PointContainerT::value_type PointPrimitiveT;

        GidPidVectorMap populations; // populations[patch_id] = all points with GID==patch_id
        processing::getPopulations( populations, points );

        std::map< DidT, Eigen::Vector3f> colourMap;
        getColours( colourMap, prims );

        Eigen::Vector3f center( 240, 400, 0 );

        {
            std::string cloudPath = ("cloud_" + path);
            FILE* fpCloud = fopen( cloudPath.c_str(), "w" );
            Eigen::Matrix<Scalar,2,1> pSize; pSize << 0.001, 0.001;
            Eigen::Vector3f colour = Eigen::Vector3f::Zero();
            for ( UPidT pid = 0; pid != points.size(); ++pid )
            {
                if ( rand() / float(RAND_MAX) > subSample ) continue;
                if ( colourCloud )
                {
                    auto itt = prims.find( points[pid].getTag(PointPrimitiveT::TAGS::GID) );
                    if ( itt != prims.end() )
                        colour = colourMap[ itt->second[0].getTag(PrimitiveT::TAGS::DIR_GID) ];
                    else
                        colour = Eigen::Vector3f::Zero();
                }
                //x, y, r, 0, r arc closepath
                drawCircle( fpCloud, points[pid].template pos().template head<2>().eval()
                            , radius
                            , colour / 255.f
                            , center
                            );

//                drawLine( fpCloud, points[pid].template pos().template head<2>().eval() - pSize
//                                 , points[pid].template pos().template head<2>().eval() + pSize
//                                 , Eigen::Vector3f::Zero()
//                                 , center );
            }
            std::stringstream command;
            command << "(evince " << cloudPath << " &)";
            if ( show )
            {
                system( command.str().c_str() );
                system( ("ps2pdf " + cloudPath).c_str() );
            }
            std::cout << command.str() << std::endl;

            drawFrame( fpCloud, center );

            fprintf( fpCloud, "showpage\n" );
            fclose( fpCloud );
        }


        FILE* fp = fopen( path.c_str(), "w" );

        std::vector<Eigen::Matrix<Scalar,2,1> > drawingMinMax(2);
        drawingMinMax[0] << 1000., 1000.;
        drawingMinMax[1] << 0., 0.;
        std::vector<Eigen::Matrix<Scalar,3,1> > extents;
        for ( typename containers::PrimitiveContainer<PrimitiveT>::ConstIterator it(prims); it.hasNext(); it.step() )
        {
            if ( it->getTag(PrimitiveT::TAGS::STATUS) == PrimitiveT::STATUS_VALUES::SMALL ) continue;

            if ( dids
                 && dids->size()
                 && std::find( dids->begin(), dids->end(), it.getDid()) == (*dids).end() ) continue;

            it->template getExtent<PointPrimitiveT>( extents
                                                   , points
                                                   , scale
                                                   , &(populations[it.getGid()]) );

            fprintf( fp, "newpath\n" );
            fprintf( fp, "%f %f %f setrgbcolor\n", colourMap[it.getDid()](0) / 255.
                                                 , colourMap[it.getDid()](1) / 255.
                                                 , colourMap[it.getDid()](2) / 255. );
            std::vector<Eigen::Vector2f> locs(2);
            locs[0] << extents[0](0)*200.+center(0), extents[0](1)*200.+center(1);
            locs[1] << extents[1](0)*200.+center(0), extents[1](1)*200.+center(1);

            fprintf( fp, ".5 setlinewidth\n" );
            fprintf( fp, "%.3f %.3f moveto %.3f %.3f lineto \nstroke \n"
                   , locs[0](0), locs[0](1)
                   , locs[1](0), locs[1](1) );

            if ( writeNames )
            {
                Eigen::Vector2f mid = locs[0] + (locs[1] - locs[0])/2.f;
                mid(0) += 1.f;
                fprintf( fp, "/Times-Roman findfont\n6 scalefont\nsetfont\n%.3f %.3f moveto\n(P%ld->%ld) show\n"
                         , mid(0), mid(1), it.getDid(), it.getGid() );

                for ( int i = 0; i != 2; ++i )
                {
                    for ( int c = 0; c != 2; ++c )
                    {
                        if ( locs[i](c) < drawingMinMax[0](c) ) drawingMinMax[0](c) = locs[i](c);
                        if ( locs[i](c) > drawingMinMax[1](c) ) drawingMinMax[1](c) = locs[i](c);
                    }
                }
            }
        }

        drawFrame( fp, center );

        fprintf( fp, "showpage\n" );
        fclose( fp );

        std::stringstream command;
        command << "(evince " << path << " &)";
        if ( show )
        {
            system( command.str().c_str() );
            system( ("ps2pdf " + path).c_str() );
        }
        std::cout << command.str() << std::endl;
    }

    template <class _PrimitiveMapT, class _PointContainerT, typename _Scalar>
    inline void drawNew( std::ofstream &f, _PrimitiveMapT & prims, _PointContainerT const& points, _Scalar pw = 0.f, std::vector<DidT>* allowedDids = NULL )
    {
        typedef typename _PrimitiveMapT::mapped_type::value_type PrimitiveT;

        std::map< DidT, Eigen::Vector3f> colourMap;
        getColours( colourMap, prims );

        std::set<DidT> dids;
        for ( typename containers::PrimitiveContainer<PrimitiveT>::ConstIterator it(prims); it.hasNext(); it.step() )
        {
            if ( it->getTag( PrimitiveT::TAGS::STATUS ) == PrimitiveT::STATUS_VALUES::SMALL ) continue;
            if ( allowedDids
                 && allowedDids->size()
                 && std::find( allowedDids->begin(), allowedDids->end(), it.getDid()) == (*allowedDids).end() ) continue;

            std::string primName( getPrimName(it.getGid(), it.getDid()) );
            f << "\t" << primName << " [ pos=\"" << it->template pos()(0) * 50. << ","
                                                 << it->template pos()(1) * 50. << "!\"";
            char cStr[128];
            sprintf( cStr, "#%02x%02x%02x"
                   , static_cast<int>( colourMap[it.getDid()](0) )
                   , static_cast<int>( colourMap[it.getDid()](1) )
                   , static_cast<int>( colourMap[it.getDid()](2) ) );
            f << ", color=\"" << cStr << "\"";
            f << ", label=\"$P_{" << it.getDid() << "\\\\rightarrow " << it.getGid() << "}$\"";
            f << "];\n";

            f << "\t" << primName
              << " -- "
              << didName( it.getDid() )
              << "[ style=\"dashed\", color=\"gray\"]"
              << "\n";

            dids.insert( it.getDid() );
        }

        for ( auto it = dids.begin(); it != dids.end(); ++it )
        {
            f << "\t" << didName( *it ) << "[ "
              << "shape=\"polygon\", sides=\"6\", regular=\"true\", fontsize=\"6\""
              << "]\n";

            for ( auto it2 = dids.begin(); it2 != dids.end(); ++it2 )
            {
                if ( (it2 == it) || (*it >= *it2) )
                    continue;

                f << "\t" << didName( *it )
                  << " -- "
                  << didName( *it2 )
                  << "[ fontcolor=\"red\", penwidth=\"5\", fontsize=\"28\"";
                  if ( pw > 0. )
                      f << ", label=\"" << pw << "\"";
                  f << "]\n";
            }
        }
    }

    template <class _PrimitiveMapT, class _PointContainerT, typename _Scalar>
    inline void drawOld( std::ofstream &f, _PrimitiveMapT & prims, _PointContainerT const& points, _Scalar pw = 0.f, bool showClusters = true
            , std::vector<DidT>* dids = NULL
            , std::vector<std::pair<GidT,DidT> >* edgeSources = NULL
            )
    {
        typedef typename _PrimitiveMapT::mapped_type::value_type PrimitiveT;

        std::map< DidT, Eigen::Vector3f> colourMap;
        getColours( colourMap, prims );

        // check, if more cands / gid
        bool needPatches = false;
        GidT prevGid = -2;
        for ( typename containers::PrimitiveContainer<PrimitiveT>::ConstIterator it(prims); it.hasNext() && !needPatches; it.step() )
        {
            if ( it->getTag( PrimitiveT::TAGS::STATUS ) == PrimitiveT::STATUS_VALUES::SMALL ) continue;
            if ( prevGid == -2 ) prevGid = it.getGid();
            else if ( prevGid == it.getGid() ) needPatches = true;
        }

        prevGid = -2;
        for ( typename containers::PrimitiveContainer<PrimitiveT>::ConstIterator it(prims); it.hasNext(); it.step() )
        {
            if ( it->getTag( PrimitiveT::TAGS::STATUS ) == PrimitiveT::STATUS_VALUES::SMALL ) continue;

            if ( dids
                 && dids->size()
                 && std::find( dids->begin(), dids->end(), it.getDid()) == (*dids).end() ) continue;

            if ( needPatches && (it.getGid() != prevGid) )
            {
                if ( prevGid != -2 )
                    f << "} // cluster" << prevGid << "\n";
                prevGid = it.getGid();

                f << "subgraph cluster" << it.getGid() << " {\n";
            }

            std::string primName( getPrimName(it.getGid(), it.getDid()) );

            f << "\t" << primName << " [ pos=\""
              << it->template pos()(0) * 50. << ","
              << it->template pos()(1) * 50.;
            if ( !needPatches )
                f << "!";
            f << "\"";

            char cStr[128];
            sprintf( cStr, "#%02x%02x%02x"
                   , static_cast<int>( colourMap[it.getDid()](0) )
                   , static_cast<int>( colourMap[it.getDid()](1) )
                   , static_cast<int>( colourMap[it.getDid()](2) ) );
            f << ", color=\"" << cStr << "\"";
            f << ", label=\"$P_{" << it.getDid() << "\\\\rightarrow " << it.getGid() << "}$\"";
            f << "];\n";
        }
        if ( needPatches )
            f << "} // cluster" << prevGid << "\n";

        // end adding nodes
        // start adding edges

        for ( typename containers::PrimitiveContainer<PrimitiveT>::ConstIterator it(prims); it.hasNext(); it.step() )
        {
            if ( it->getTag( PrimitiveT::TAGS::STATUS ) == PrimitiveT::STATUS_VALUES::SMALL ) continue;

            if ( dids
                 && dids->size()
                 && std::find( dids->begin(), dids->end(), it.getDid()) == (*dids).end() ) continue;

            if ( edgeSources && (std::find(edgeSources->begin(),edgeSources->end(), std::pair<GidT,DidT>(it.getGid(),it.getDid())) == edgeSources->end()) )
                continue;

            std::string primName( getPrimName(it.getGid(), it.getDid()) );
            for ( typename containers::PrimitiveContainer<PrimitiveT>::ConstIterator it2(prims);
                  it2.hasNext(); it2.step() )
            {
                if ( it2->getTag( PrimitiveT::TAGS::STATUS ) == PrimitiveT::STATUS_VALUES::SMALL ) continue;
                if ( it.getGid() >= it2.getGid () ) continue;

                if ( dids
                     && dids->size()
                     && std::find( dids->begin(), dids->end(), it2.getDid()) == (*dids).end() ) continue;

                if ( showClusters != (it.getDid() == it2.getDid()) ) continue;

//                if ( edgeSources && (std::find(edgeSources->begin(),edgeSources->end(), std::pair<GidT,DidT>(it2.getGid(),it2.getDid())) == edgeSources->end()) )
//                    continue;

                f << "\t" << primName
                  << " -- "
                  << getPrimName( it2.getGid(), it2.getDid() );
                if ( showClusters ) f << " [ style=\"dashed\", color=\"gray\"]";
                f  << "\n";
            }
        } //...for outer


#if 0
        for ( auto it = dids.begin(); it != dids.end(); ++it )
        {
            f << "\t" << didName( *it ) << "[ "
              << "shape=\"polygon\", sides=\"6\", regular=\"true\", fontsize=\"32\""
              << "]\n";

            for ( auto it2 = dids.begin(); it2 != dids.end(); ++it2 )
            {
                if ( (it2 == it) || (*it >= *it2) )
                    continue;

                f << "\t" << didName( *it )
                  << " -- "
                  << didName( *it2 )
                  << "[ label=\"" << 99 << "\", fontcolor=\"red\", penwidth=\"5\", fontsize=\"28\" ]\n";
            }
        }
#endif
    }

    template <class _PrimitiveMapT, class _PointContainerT, typename _Scalar>
    inline void drawGraph( _PrimitiveMapT  &prims, _PointContainerT const& points, std::string path
                          , bool old = false
                         , bool     show    = false
                         , std::vector<DidT> *dids = NULL
                         , _Scalar  pw      = 1.
                         , bool showClusters = false
                         , std::vector< std::pair<GidT,DidT> > *edgeSources = NULL
                         )
    {
        std::string pPath( "prims_" + path );
        std::string gPath( pPath + ".gv" );
        std::ofstream f( gPath );
        if ( !f.is_open() ) { std::cerr << "[" << __func__ << "]: " << "could not open " << path << std::endl; return; }

        f << "graph {\n";
        f << "splines = \"spline\"\n";
        f << "node [fontsize=\"4\", sep=\"+0.1,0.1\", shape=\"ellipse\"]\n";
        f << "edge [fontsize=\"20\", sep=\"+0.1,0.1\"]\n";
        f << "subgraph cluster {\n";

        if ( old )
            drawOld( f, prims, points, pw, showClusters, dids, edgeSources );
        else
            drawNew( f, prims, points, pw, dids );

        f << "}\n}\n";

        f.close();

        std::stringstream command, cm2;
        cm2 << "inkscape -A " << pPath << ".pdf --export-latex " << pPath << ".svg";
        command << "(fdp -Tsvg -o " << pPath << ".svg" << " " << gPath << " -s1.44 && eog " << pPath << ".svg";
        command << " && " << cm2.str() << "&)";
        if ( show )
        {
            //system( command.str().c_str() );
            //system( cm2.str().c_str() );
        }
        std::cout << command.str() << std::endl;
        std::cout << cm2.str() << std::endl;

    } //...drawGraph
} //...ns io
} //...ns rapter


#endif // __RAPTER_PS_HPP__
