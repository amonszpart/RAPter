#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include "globfit2/primitives/linePrimitive.h"

namespace GF2
{

    using ::GF2::LinePrimitive;

    int rotatePoints( cv::Point2f &p1
                      , cv::Point2f &p2
                      , cv::Point2f const& center
                      , float const angle )
    {
        cv::Point2f p1_rotated, p2_rotated;
        p1 -= center; p2-= center;
        p1_rotated.x = p1.x * cos(angle) - p1.y * sin( angle );
        p1_rotated.y = p1.x * sin(angle) + p1.y * cos( angle );
        p2_rotated.x = p2.x * cos(angle) - p2.y * sin( angle );
        p2_rotated.y = p2.x * sin(angle) + p2.y * cos( angle );
        p1 = p1_rotated + center; p2 = p2_rotated + center;

        return EXIT_SUCCESS;
    }

    // ==================== PRIMITIVES ==================== //
    int gtLine( cv::Mat                           & img
                , std::vector<LinePrimitive>      & lines
                , cv::Point2i                const& tl
                , cv::Point2i                const& br
                , float                      const  scale )
    {
        // output line primitive
        Eigen::Vector3f p0 = (Eigen::Vector3f()<<        tl.x/(float)img.cols  * scale,
                                                  (1.f - tl.y/(float)img.rows) * scale,
                                                   0.f ).finished();
        Eigen::Vector3f p1 = (Eigen::Vector3f()<<        br.x/(float)img.cols * scale,
                                                  (1.f-  br.y/(float)img.rows) * scale,
                                                   0.f ).finished();
        lines.emplace_back( LinePrimitive::fromEndPoints(p0, p1) );

        // draw line on image
        cv::line( img, tl, br, cv::Scalar::all(0) );

        return EXIT_SUCCESS;
    }

    int gtRect( cv::Mat                           & img
                , std::vector<LinePrimitive>      & lines
                , cv::Point2i                const& tl
                , cv::Point2i                const& br
                , float                      const  angle       = 0.f
                , cv::Point2f                const* center_arg  = NULL
                , float                      const  scale       = 1.f )
    {
        const std::vector<int> indices = { 0, 0, 1, 0,
                                           0, 0, 0, 1,
                                           0, 1, 1, 1,
                                           1, 0, 1, 1 };

        const cv::Point2f center = center_arg ? (*center_arg) : cv::Point2f(img.cols/2.f, img.rows/2.f);

        // check if all lines fit
        if ( angle != 0.f )
            for ( size_t cid = 0; cid != indices.size(); cid+=4 )
            {
                // endpoints
                cv::Point2f p1( (indices[cid+0]?br:tl).x, (indices[cid+1]?br:tl).y );
                cv::Point2f p2( (indices[cid+2]?br:tl).x, (indices[cid+3]?br:tl).y );

                rotatePoints( p1, p2, center, angle );

                if ( (std::min(p1.x,p2.x) < 0.f) || (std::max(p1.x,p2.x) > img.cols) ||
                     (std::min(p1.y,p2.y) < 0.f) || (std::max(p1.y,p2.y) > img.rows)) return EXIT_FAILURE;

            }
        std::cout << "rect at " << tl.x << ", " << tl.y << " and " << br.x << ", " << br.y << std::endl;

        // four lines of a  rectangle
        for ( size_t cid = 0; cid != indices.size(); cid+=4 )
        {
            // endpoints
            cv::Point2f p1( (indices[cid+0]?br:tl).x, (indices[cid+1]?br:tl).y );
            cv::Point2f p2( (indices[cid+2]?br:tl).x, (indices[cid+3]?br:tl).y );

            // rotated endpoints
            cv::Point2f p1_rotated = p1, p2_rotated = p2;

            // rotate if necessary
            if ( angle != 0.f )
            {
                rotatePoints( p1, p2, center, angle );

                //            p1_rotated.x = std::min( static_cast<float>(img.cols), std::max( 0.f, p1.x));
                //            p1_rotated.y = std::min( static_cast<float>(img.rows), std::max( 0.f, p1.y));
                //            p2_rotated.x = std::min( static_cast<float>(img.cols), std::max( 0.f, p2.x));
                //            p2_rotated.y = std::min( static_cast<float>(img.rows), std::max( 0.f, p2.y));
                p1_rotated = p1;
                p2_rotated = p2;
                if ( (std::min(p1.x,p2.x) < 0.f) || (std::max(p1.x,p2.x) > img.cols) ||
                     (std::min(p1.y,p2.y) < 0.f) || (std::max(p1.y,p2.y) > img.rows)) return EXIT_FAILURE;

            }

            // output line primitive
            lines.emplace_back( LinePrimitive::fromEndPoints((Eigen::Vector3f()<<      p1_rotated.x/(float)img.cols * scale,
                                                              (1.f-  p1_rotated.y/(float)img.rows) * scale,
                                                              0.f ).finished(),
                                                             (Eigen::Vector3f()<<      p2_rotated.x/(float)img.cols * scale,
                                                              (1.f-  p2_rotated.y/(float)img.rows) * scale,
                                                              0.f ).finished()) );
            // draw line on image
            cv::line( img, p1_rotated, p2_rotated, cv::Scalar::all(0) );
        }

        return EXIT_SUCCESS;
    }

    // ==================== IMAGES ===================== //


    /**
 * @brief gtImg Two perpendicular rectangles rotated.
 * @param img
 * @param lines
 * @param angles
 * @param points_arg
 * @return
 */
    int gtImg( cv::Mat &img
               , std::vector<LinePrimitive> &lines
               , std::vector<float> const* angles
               , std::vector<cv::Point2i> const* points_arg )
    {
        img.create( 340, 340, CV_8UC1 );
        img.setTo( 255 );
        std::vector<cv::Point2i> points;
        if ( points_arg ) points = *points_arg;
        else              points =
        {
            cv::Point2i(100,40), cv::Point2i(225,290),
            cv::Point2i(5  ,85), cv::Point2i(295,230)
        };

        if ( angles )
            if ( angles->size() != points.size()/2 )
            {
                std::cerr << "[" << __func__ << "] " << " need exactly " << points.size()/2 << " angles for this GTimg" << std::endl;
                return EXIT_FAILURE;
            }

        for ( size_t pid = 0; pid != points.size(); pid += 2 )
            gtRect( img, lines, points[pid], points[pid+1], angles ? angles->at(pid/2)*M_PI/180.f : -45.f*M_PI/180.f );

        cv::imshow( "img", img );
        cv::waitKey(100);
        return 0;
    }

    int gtRandomLines( cv::Mat &img
                       , std::vector<LinePrimitive> &lines
                       , int count
                       , float scene_size
                       , float min_angle, float max_angle
                       )
    {
        cv::Point2i center( img.cols/2, img.rows/2 );
        float diag_length = sqrtf( img.cols*img.cols+img.rows*img.rows )/2.f * .5f;
        Eigen::Vector3f dir( 0.f, 1.f, 0.f );
        int err = EXIT_SUCCESS;
        for ( int lid = 0; lid != count; ++lid )
        {
            float angle = min_angle + (rand()/static_cast<float>(RAND_MAX)) * (max_angle - min_angle);
            dir = Eigen::AngleAxisf( angle, Eigen::Vector3f::UnitZ() ) * dir;
            cv::Point2i a( center.x + dir(0) * diag_length, center.y + dir(1) * diag_length );
            cv::Point2i b( center.x - dir(0) * diag_length, center.y - dir(1) * diag_length );
            err += gtLine( img, lines, a, b, scene_size );
        }

        return err;
    }

    int gtPearl( cv::Mat &img, std::vector<LinePrimitive> &lines
                 , float spacing          = 50
                 , int   n_lines          = 5
                 , float const scene_size = 1.f )
    {
        img.create( 640, 640, CV_8UC1 );
        img.setTo( 255 );

        cv::Point2f p0( img.cols*.4, img.rows * .1), p1( img.cols * .8, img.rows * .9 );
        cv::Point2f cvdir = p1-p0;
        float norm = sqrt(cvdir.x*cvdir.x+cvdir.y*cvdir.y);
        cvdir.x /= norm;
        cvdir.y /= norm;

        const float perp_length = 150;

        // long single line
        std::vector<cv::Point2i> points =
        {
            cv::Point2i( round(p0.x + cvdir.y*perp_length),
                         round(p0.y - cvdir.x*perp_length) ),
            cv::Point2i( round(p1.x + cvdir.y*perp_length),
                         round(p1.y - cvdir.x*perp_length) ),
        };

        // perpendicular lines
        for ( int i = 1; i != n_lines+1; ++i )
        {
            cv::Point2f start = p0 + cvdir * spacing * i;
            if ( cvdir.y == 0 )
                start.y += spacing / 2.f;

            points.push_back( cv::Point2i(round(start.x),
                                          round(start.y)) );
            points.push_back( cv::Point2i(round(start.x - cvdir.y  * perp_length),
                                          round(start.y + cvdir.x  * perp_length) ) );
        }

        for ( size_t pid = 0; pid != points.size(); pid += 2 )
            gtLine( img, lines, points[pid], points[pid+1], scene_size );

        cv::imshow( "img", img );
        cv::waitKey(100);

        return 0;
    }

    int gtStairs( cv::Mat &img
                , std::vector<::GF2::LinePrimitive> &lines
                , int steps
                , float line_width_percentage_of_column
                , float decay
                , float scale )
    {
        img.create( 640, 640, CV_8UC1 );
        img.setTo( 255 );

        float col_width = img.cols / (steps * 2 + 1);
        float y_spacing = img.cols / (steps * 2 + 1);
        float line_width_offset = col_width * (line_width_percentage_of_column - 1.f) / 2.f; // will be added twice to the line width

        for ( int step = 0; step != steps; ++step )
        {
            cv::Point2i p0( (1+(step*2)  ) * col_width -              line_width_offset                  , img.cols - (1+(step*2)) * y_spacing );
            cv::Point2i p1( (1+(step*2)  ) * col_width + (col_width + line_width_offset)*(1.f-step*decay), img.cols - (1+(step*2)) * y_spacing );
            gtLine( img, lines, p0, p1, scale );
        }
        cv::imshow( "img" , img );
        cv::waitKey(100);

        return EXIT_SUCCESS;
    }


    /**
 * @brief gtImg2 Many rectangles rotated arbitrarily.
 * @param img Output image
 * @param lines Output lines
 * @param rectCount Input number of rectangles
 * @return EXIT_SUCCESS
 */
    int gtImg2( cv::Mat &img
                , std::vector<LinePrimitive> &lines
                , int rectCount  )
    {
        img.create( 340, 340, CV_8UC1 );
        img.setTo( 255 );
        const cv::Point2i xlim( 100, img.cols-100); // x0, width
        const cv::Point2i ylim( 100, img.rows-100); // y0, height

        std::vector<cv::Point2i> points;

        for ( int rect_id = 0; rect_id != rectCount*2; ++rect_id )
        {
            points.emplace_back( cv::Point2i( xlim.x + rand() % xlim.y,
                                              ylim.x + rand() % ylim.y ) );

            while ( abs(points.back().x - points[points.size()-2].x) < 7 ) // 7
            {
                points.back().x           = xlim.x + rand() % xlim.y;
                points[points.size()-2].x = xlim.x + rand() % xlim.y;
            }
            while ( abs(points.back().y - points[points.size()-2].y) < 7 )
            {
                points.back().y           = ylim.x + rand() % ylim.y;
                points[points.size()-2].y = ylim.x + rand() % ylim.y;
            }
        };

        for ( size_t pid = 0; pid != points.size(); pid += 2 )
        {
            // draw rectangle
            gtRect( img, lines, points[pid], points[pid+1], rand()/(float)RAND_MAX*90.f*M_PI/180.f );
        }

        cv::imshow( "img", img );
        cv::waitKey(100);
        return 0;
    }


    int fill_cell( cv::Mat &img
                   , std::vector<LinePrimitive> &lines
                   , cv::Rect2i roi
                   , float angle
                   , float scale )
    {
        std::cout << "[" << __func__ << "]: " << roi.x << "-" << roi.x+roi.width << "," << roi.y << "-" << roi.y+roi.height << std::endl;

        const int min_side_length = //round(sqrt( roi.width * roi.width + roi.height * roi.height ) / 2.f);
                                    round( std::min( roi.width, roi.height ) / 3.f );
        std::cout << "min_side_length: " << min_side_length << std::endl;
        const cv::Point2f center( roi.x + roi.width/2.f, roi.y + roi.height/2.f);

        //cv::Point2f p1( roi.x + rand() % roi.width, roi.y + rand() % roi.height ),
        //p2( roi.x + rand() % roi.width, roi.y + rand() % roi.height );
        //rotatePoints( p1, p2, center, angle );

        bool success = false;
        do
        {
            success = false;

            cv::Point2f p1,p2;
            bool retry_x = true, retry_y = true;
            do
            {
                if ( retry_x ) { p1.x = roi.x + rand() % roi.width ; p2.x = roi.x + rand() % roi.width ; }
                if ( retry_y ) { p1.y = roi.y + rand() % roi.height; p2.y = roi.y + rand() % roi.height; }

                rotatePoints( p1, p2, center, angle );
            }
            while (     ( retry_x  = (abs(p1.x - p2.x) < min_side_length)          )
                        || ( retry_y  = (abs(p1.y - p2.y) < min_side_length)          )
                        || ( retry_x |= (    (std::min(p1.x,p2.x) < roi.x)
                                             || (std::max(p1.x,p2.x) >= roi.x+roi.width)     )           )
                        || ( retry_y |= (    (std::min(p1.y,p2.y) < roi.y)
                                             || (std::max(p1.y,p2.y) >= roi.y+roi.height)    )           ));

            success = (EXIT_SUCCESS == gtRect(img, lines, p1, p2, angle, &center, scale ) );
        } while ( !success );

        return EXIT_SUCCESS;
    }

    int gtImg_stratified( cv::Mat &img
                          , std::vector<LinePrimitive> &lines
                          , std::vector<float> &angles
                          , int n_rects
                          , int n_angles
                          , const float scale )
    {
        if ( img.empty() )
        {
            img.create( 640, 640, CV_8UC1 );
            img.setTo( 255 );
        }

        int cells_vertical = 1, cells_horizontal = 1;
        while ( cells_vertical * cells_horizontal < n_rects )
        {
            if ( cells_vertical < cells_horizontal ) ++cells_vertical;
            else                                     ++cells_horizontal;
        }
        std::cout << cells_vertical << "x" << cells_horizontal << "(M_PI_2 / (float)n_angles): " << (M_PI_2 / (float)n_angles) << std::endl;

        for ( int i = 0; i != n_angles; ++i )
        {
            float ang = 0.f;
            do
            {
                ang = rand() / static_cast<float>(RAND_MAX);
                ang = (static_cast<int>(ang * 90.f) / 5 * 5) * M_PI / 180.f;
                std::cout << "ang: " << ang << ", " << ang * 180.f / M_PI;
                if ( angles.size() ) std::cout << ", diff: " << fabs(ang - angles.back()) << " <? " << (M_PI_2 / (float)n_angles / 2.f);
                std::cout << std::endl;
            } while (    ( angles.size()                                          )
                         && ( fabs(ang - angles.back()) < (M_PI_2 / (float)n_angles / 2.f) ) );

            angles.push_back( ang );
        }
        if ( !angles.size() ) angles.push_back(0.f);
        std::cout<<"angles:";for(size_t vi=0;vi!=angles.size();++vi)std::cout<<angles[vi] << "(" << angles[vi] * 180.f / M_PI << ")" <<" ";std::cout << "\n";

        const int cell_width  = img.cols / cells_horizontal;
        const int cell_height = img.rows / cells_vertical;
        std::cout << "cell_width: " << cell_width << std::endl;
        std::cout << "cell_height: " << cell_height << std::endl;
        const int frame = std::min( 10, std::min(cell_width,cell_height) / 2);

        int count = 0;
        for ( int y = 0; y != cells_vertical; ++y )
            for ( int x = 0; x != cells_horizontal; ++x )
            {
                int angle_id = rand() % angles.size();
                if ( count < n_rects )
                {
                    fill_cell( img, lines, cv::Rect2i(x*cell_width+frame,y*cell_height+frame,cell_width-2*frame,cell_height-2*frame), angles[angle_id], scale );
                    ++count;
                }
                std::cout << y << "," << x << std::endl;
            }

        return EXIT_SUCCESS;
    }

} // ... GF2
