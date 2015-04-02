#include <sstream>
#include "boost/algorithm/string.hpp"

#include "Types.h"
#include "Primitive.h"

Primitive::Primitive(const std::vector<RichPoint*>& vecPointSet, PrimitiveType primitiveType) :
    _vecPointSet(vecPointSet),
    _primitiveType(primitiveType)
{
    _precision = 0.0;
}


Primitive::~Primitive(void)
{
}

bool Primitive::load(std::ifstream& fin)
{
    bool bPointsLoaded = false;
    bool bPrimitiveLoaded = false;

    while (fin && !(bPointsLoaded && bPrimitiveLoaded)){
        std::string line;
        getline(fin, line);
        boost::algorithm::trim(line);

        if (line.empty() || line[0]=='#') {
            continue;
        }

        std::string indication = line.substr(0, line.find_first_of(" \t"));
        if (indication == "points") {
            if (bPointsLoaded) {
                return false;
            } else {
                bPointsLoaded = true;
                std::stringstream sin(line.substr(line.find_first_of(" \t")+1));
                size_t index = 0;
                sin >> index;
                while (sin) {
                    _vecPointIdx.push_back(index);
                    sin >> index;
                }
            }
        } else {
            if(bPrimitiveLoaded) {
                return false;
            } else {
                if (loadParameters(line)) {
                    bPrimitiveLoaded = true;
                    //std::cout << "[" << __func__ << "]: " << "loaded params" << std::endl;
                } else {
                    //std::cout << "[" << __func__ << "]: " << "!! did not load params" << std::endl;
                    return false;
                }
            }
        }
    }

    computePrecision();
    return true;
}

bool Primitive::save(std::ofstream& fout) const
{
    if (!saveParameters(fout)) {
        return false;
    }

    fout << "# points idx_1 idx_2 idx_3 ... " << std::endl;
    fout << "points";
    for (std::vector<size_t>::const_iterator it = _vecPointIdx.begin();
        it != _vecPointIdx.end();
        ++ it) {
            fout << " " << *it;
    }
    fout << std::endl;

    return true;
}
