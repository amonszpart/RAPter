#ifndef __GF2_LINECLUSTERING_HPP__
#define __GF2_LINECLUSTERING_HPP__

#include <vector>

namespace am
{

    template <class PrimitivesT>
    struct PrimitiveClustering
    {
            PrimitiveClustering()
                : directions_in_cluster( -1 )
                , first_c( -1 )
            {}

            PrimitivesT                     lines;                   // contains all lines
            std::vector< std::vector<int> > clusters_lines;          // contains a vector of line_ids for each cluster_id
            std::vector< int              > lines_clusters;          // contains one cluster_id for each line in lines
            std::vector< int              > cluster_representatives; // contains one line_id for each cluster that has a representative direction for that cluster
            int                             directions_in_cluster;
            int                             first_c;                 // index of first line in lines, that is a consensus line direction, and not a local fit

            PrimitivesT      & planes()       { return lines; }
            PrimitivesT const& planes() const { return lines; }

            typename PrimitivesT::value_type const&
            getClusterRepresentative( int cid ) const { return lines[ cluster_representatives[cid] ]; }

            int
            getClusterID( int lid ) const { return lines_clusters[ lid ]; }

            typename PrimitivesT::value_type const&
            getFromCluster( int cid, int i ) const { return lines[ clusters_lines[cid][i] ]; }

            int
            add( typename PrimitivesT::value_type const& primitive, int cid, int origin = -1 )
            {
                lines         .push_back( primitive );
                lines_clusters.push_back( cid );

                if ( clusters_lines.size() <= cid )
                    clusters_lines.resize( cid+1 );
                clusters_lines[cid].push_back( lines.size() - 1 );

                return EXIT_SUCCESS;
            }

    };

    /* let lid : 0 < line_id < lines.size()
     * if lid < first_c, then
     *  cid0 = lines_clusters[lid]
     *  and lid is incompatible with all other lid1 >= first_c, where lines_clusters[lid1] == cid0
     */

}

#endif // __GF2_LINECLUSTERING_HPP__
