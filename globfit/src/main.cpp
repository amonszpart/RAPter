#include <string>
#include <iostream>
//#include <boost/thread.hpp>
#include "boost/thread.hpp"
#include "boost/program_options.hpp"
#include <osg/Group>

#include "Viewer.h"
#include "GlobFit.h"

int main(int argc, char *argv[])
{
    double paraOrthThreshold, equalAngleThreshold, coaxialThreshold, coplanarThreshold, equalLengthThreshold, equalRadiusThreshold;

    namespace po = boost::program_options;
    po::options_description options("Allowed options");
    options.add_options()
            ("help,h", "produce help message")
            ("verbose,v", "turn on verbose mode")
            ("input,i", po::value<std::string>(), "input file")
            ("paraOrthThreshold,o", po::value<double>(&paraOrthThreshold)->default_value(10.00, "10.00"), "parallel/orthogonal threshold")
            ("equalAngleThreshold,g", po::value<double>(&equalAngleThreshold)->default_value(10.00, "10.00"), "equal angle threshold")
            ("coaxialThreshold,a", po::value<double>(&coaxialThreshold)->default_value(0.02, "0.02"), "coaxial threshold")
            ("coplanarThreshold,p", po::value<double>(&coplanarThreshold)->default_value(0.02, "0.02"), "coplanar threshold")
            ("equalLengthThreshold,l", po::value<double>(&equalLengthThreshold)->default_value(0.02, "0.02"), "equal length threshold")
            ("equalRadiusThreshold,r", po::value<double>(&equalRadiusThreshold)->default_value(0.02, "0.02"), "equal radius threshold")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    if (!vm.count("input")) {
        std::cout << options << "\n";
        system("read -p \'Press any key...\'");
        return 1;
    }

    bool verbose = vm.count("verbose");

    std::vector<osg::Node*> vecViewData;
    for (size_t i = 0; i < 6; ++ i) {
        vecViewData.push_back(new osg::Group);
        vecViewData.back()->setDataVariance(osg::Object::DYNAMIC);
    }
    Viewer viewer;
    boost::thread viewerThread(viewer, vecViewData);

    std::string inputFilename = vm["input"].as<std::string>();
    size_t dotPos = inputFilename.find_last_of('.');
    std::string base = inputFilename.substr(0, dotPos);
    std::string ext = inputFilename.substr(dotPos);

    // read input file
    GlobFit globFit;
    if (!globFit.load(inputFilename)) {
        system("read -p \'Press any key...\'");
        return 1;
    }
    std::pair<osg::Node*, osg::Node*> pointsGeometry = globFit.convertPointsToGeometry();
    dynamic_cast<osg::Group*>(vecViewData[0])->addChild(pointsGeometry.first);
    dynamic_cast<osg::Group*>(vecViewData[1])->addChild(pointsGeometry.second);
    dynamic_cast<osg::Group*>(vecViewData[2])->addChild(globFit.convertPrimitivesToGeometry("Initial Primitives"));
    viewerThread.detach();

    if (!globFit.createMatlabArraies()) {
        std::cout << "[" << __func__ << "]: " << "createMatlabArrays failed..." << std::endl;
        system("read -p \'Press any key...\'");
        return 1;
    }

    std::cout << "["<< __FUNCTION__ << "]: " << "Orientation Alignment" << std::endl;
    // Orientation Alignment
    if (!globFit.orientationAlignment(paraOrthThreshold, equalAngleThreshold)) {
        std::cout << "[" << __func__ << "]: " << "orienationAlignment failed..." << std::endl;
        globFit.destoryMatlabArraies();
        system("read -p \'Press any key...\'");
        return 1;
    }
    dynamic_cast<osg::Group*>(vecViewData[3])->addChild(globFit.convertPrimitivesToGeometry("Orientation Alignment"));
    if (verbose) {
        std::string oaFilename = base+"_oa"+ext;
        globFit.save(oaFilename);
    }

    std::cout << "["<< __FUNCTION__ << "]: " << "Placement Alignment" << std::endl;
    // Placement Alignment
    if (!globFit.placementAlignment(coaxialThreshold, coplanarThreshold)) {
        std::cout << "[" << __func__ << "]: " << "placementAlignment failed..." << std::endl;
        globFit.destoryMatlabArraies();
        system("read -p \"Pressanykey...\"");
        return 1;
    }
    dynamic_cast<osg::Group*>(vecViewData[4])->addChild(globFit.convertPrimitivesToGeometry("Placement Alignment"));
    if (verbose) {
        std::string paFilename = base+"_pa"+ext;
        globFit.save(paFilename);
    }

    std::cout << "["<< __FUNCTION__ << "]: " << "Equality Alignment" << std::endl;
    // Equality Alignment
    if (!globFit.equalityAlignment(equalLengthThreshold, equalRadiusThreshold)) {
        std::cout << "[" << __func__ << "]: " << "equalityAlignment failed..." << std::endl;
        globFit.destoryMatlabArraies();
        system("read -p \'Press any key...\'");
        return 1;
    }
    dynamic_cast<osg::Group*>(vecViewData[5])->addChild(globFit.convertPrimitivesToGeometry("Equality Alignment"));
    std::string eaFilename = base+"_ea"+ext;
    globFit.save(eaFilename);

    globFit.destoryMatlabArraies();
    //std::cout << "press any key" << std::endl;
    //char key;
    //std::cin >> key;
    //system("echo \"finished\"; read -p \"Press any key to continue...\"");
    return 0;
}
