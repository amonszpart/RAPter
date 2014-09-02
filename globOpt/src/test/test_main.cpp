#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/console/parse.h>
#include <boost/thread/thread.hpp>
#include <pcl/filters/extract_indices.h>
#include <pcl/sample_consensus/ransac.h>
//#include <pcl/visualization/pcl_visualizer.h>  // deprecated, treba nahradit
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/sample_consensus/sac_model_sphere.h>

int main(int argc, char *argv[])
{
  pcl::console::setVerbosityLevel(pcl::console::L_DEBUG);

  std::string file_name = "/media/Data2/Data/GF2/cloud.pcd";
  pcl::PointCloud<pcl::PointXYZRGBNormal> my_cloud;
  pcl::PCDReader input;
  if (input.read(file_name, my_cloud, 0) == -1)
  {
    std::cout << "Failed." << std::endl;
    return -1;
  }

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud;

  return 0;
}

