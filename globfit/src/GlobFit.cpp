#include <fstream>
#include <sstream>

#include "Types.h"
#include "Cone.h"
#include "Plane.h"
#include "Sphere.h"
#include "Cylinder.h"

#include "GlobFit.h"

bool GlobFit::stoppingAtError = false;

GlobFit::GlobFit(void)
{
}


GlobFit::~GlobFit(void)
{
  for (size_t i = 0, iEnd = _vecPointSet.size(); i < iEnd; ++ i) {
    delete _vecPointSet[i];
  }

  for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
    delete _vecPrimitive[i];
  }
}


bool GlobFit::load(const std::string& filename)
{
    std::string line;
    std::ifstream fin(filename.c_str());

    if ( !fin.good() )
    {
        std::cout << "[" << __func__ << "]: " << "!fin.good() -> return false;" << std::endl;
        return false;
    }

    size_t numPoint = 0;
    while ( fin )
    {
        getline(fin, line);
        if (line.empty() || line[0]=='#' || line[0] == 13)
        {
            continue;
        }

        std::stringstream sin(line);
        sin >> numPoint;
        std::cout << "[" << __func__ << "]: " << "numPoints: " << numPoint << std::endl;
        break;
    }

    while (fin && numPoint != 0)
    {
        getline(fin, line);

        if (line.empty() || line[0]=='#' || line[0] == 13)
        {
            continue;
        }

        --numPoint;

        RichPoint* p = new RichPoint();
        _vecPointSet.push_back(p);

        double x, y, z, nx, ny, nz, confidence;

        std::stringstream sin(line);
        sin >> x >> y >> z >> nx >> ny >> nz >> confidence;

        p->point        = Point(x, y, z);
        p->normal       = Vector(nx, ny, nz)/sqrt(nx*nx+ny*ny+nz*nz);
        p->confidence   = confidence;
    }

    size_t numPrimitive = 0;
    while (fin)
    {
        getline(fin, line);
        if ( line.empty() || line[0]=='#' || line[0] == 13)
        {
            std::cout << "[" << __func__ << "]: " << "skipping line " << line << std::endl;
            continue;
        }

        std::cout << "[" << __func__ << "]: " << "parsing line \"" << line << "\" == " << (int)line[0] << std::endl;
        std::stringstream sin(line);
        sin >> numPrimitive;
        std::cout << "[" << __func__ << "]: " << "numPrimitive: " << numPrimitive << std::endl;
        break;
    }

    while (fin && numPrimitive != 0){
        getline(fin, line);

        if (line.empty() || line[0]=='#' || line[0] == 13)
        {
            continue;
        }

        std::string indication = line.substr(0, line.find_first_of(" \t"));
        Primitive* pPrimitive = NULL;
        if (indication == "plane") {
            pPrimitive = new Plane(_vecPointSet);
        } else if (indication == "cylinder") {
            pPrimitive = new Cylinder(_vecPointSet);
        } else if (indication == "cone") {
            pPrimitive = new Cone(_vecPointSet);
        } else if (indication == "sphere") {
            pPrimitive = new Sphere(_vecPointSet);
        }

        if (pPrimitive == NULL) {
            std::cerr << "Error: bad file format!" << std::endl;
            continue;
        }

        fin.seekg(-((int)(line.size()+1)), std::ios::cur);
        if (pPrimitive->load(fin))
        {
            pPrimitive->setIdx(_vecPrimitive.size());
            _vecPrimitive.push_back(pPrimitive);
            --numPrimitive;
            //std::cout << "[" << __func__ << "]: " << "read primitive " << std::endl;
        } else {
            delete pPrimitive;
            std::cerr << "Error: bad file format!" << std::endl;
        }
    }

    return true;
}

bool GlobFit::save(const std::string& filename) const
{
    std::ofstream fout(filename.c_str());
    if (!fout.good()) {
        return false;
    }

    fout.precision(16);

    fout << "# Number of Points" << std::endl;
    fout << _vecPointSet.size() << std::endl;
    fout << "# Here comes the " << _vecPointSet.size() << " Points" << std::endl;
    fout << "# point_x point_y point_z normal_x normal_y normal_z confidence" << std::endl;
    for (std::vector<RichPoint*>::const_iterator it = _vecPointSet.begin();
         it != _vecPointSet.end();
         ++ it) {
        fout << (*it)->point << " " << (*it)->normal << " " << (*it)->confidence << std::endl;
    }
    fout << "# End of Points" << std::endl;

    fout << std::endl;
    fout << "# Number of Primitives" << std::endl;
    fout << _vecPrimitive.size() << std::endl;
    fout << "# Here comes the " << _vecPrimitive.size() << " Primitives" << std::endl;
    for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i ) {
        fout << "# Primitive " << i << std::endl;
        _vecPrimitive[i]->save(fout);
        fout << std::endl;
    }
    fout << "# End of Primitives" << std::endl;

    std::cout << "saved " << filename << std::endl;

    return true;
}

