#ifndef __GF2_GROUND_TRUTH_H__
#define __GF2_GROUND_TRUTH_H__

#include <vector>
#include "opencv2/core/core.hpp"

namespace GF2
{
    class LinePrimitive;
}

namespace GF2
{
    extern int gtLine( cv::Mat                           & img
                     , std::vector<LinePrimitive>  & lines
                     , cv::Point2i                const& tl
                     , cv::Point2i                const& br
                     , float                      const  scale );
    extern int gtRect( cv::Mat                           & img
                , std::vector<LinePrimitive>      & lines
                , cv::Point2i                const& tl
                , cv::Point2i                const& br
                , float                      const  angle       = 0.f
                , cv::Point2f                const* center_arg  = NULL
                , float                      const  scale       = 1.f );

    int gtImg_stratified( cv::Mat                      &img
                          , std::vector<LinePrimitive> &lines
                          , std::vector<float> &angles
                          , int n_rects, int n_angles, float scale );

    int gtImg( cv::Mat &img
               , std::vector<LinePrimitive> &lines
               , std::vector<float> const* angles = NULL
               , std::vector<cv::Point2i> const* points_arg = NULL);

    extern int gtPearl( cv::Mat &img, std::vector<LinePrimitive> &lines, float spacing, int n_lines
                        , float const scale = 1.f );
    int gtImg2( cv::Mat &img
                , std::vector<LinePrimitive> &lines
                , int rectCount  );

    extern int gtRandomLines( cv::Mat &img
                              , std::vector<LinePrimitive> &lines
                              , int count
                              , float scale
                              , float min_angle = 0.f, float max_angle = M_PI
                              );

    extern int gtStairs( cv::Mat                          &img
                         , std::vector<LinePrimitive> &lines
                         , int                             steps = 4
                         , float line_width_percentage_of_column = 1.4f
                         , float decay                           = .1f
                         , float scale                           = 1.f );
} // ... ns GF2

#endif // __GF2_GROUND_TRUTH_H__
