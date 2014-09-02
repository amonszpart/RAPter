#ifndef SAMPLEGENERATORFACTORY_H
#define SAMPLEGENERATORFACTORY_H

#include <QDockWidget>
#include "samplegenerator.h"
#include "types.h"
#include "typesGL.h"

namespace Ui {
class SampleGeneratorFactory;
}

class SampleGeneratorFactory : public QDockWidget
{
    Q_OBJECT

public:
    explicit SampleGeneratorFactory(QWidget *parent = 0);
    ~SampleGeneratorFactory();

    inline void setPrimitives(
            std::vector< InputGen::Application::Primitive >*s)
    {
        _pSet = s;
        updateGenerator();
    }

    inline void setPoints( InputGen::Application::PointSet*s)
    {
        _pointSet = s;
    }

public slots:
    void updateGenerator();

signals:
    //! \warning Pointers are invalid outside of the connected slots
    void samplesChanged(InputGen::Application::PointSet *,
                        InputGen::Application::SampleGenerator*);

private:
    typedef InputGen::Application::SampleGenerator SampleGenerator;

    enum GENERATOR_TYPE{
        GEN_FROM_PRIMITIVE = 0,
        GEN_FROM_PONCTUAL  = 1
    };

    Ui::SampleGeneratorFactory *ui;
    std::vector< InputGen::Application::Primitive > *_pSet;
    InputGen::Application::PointSet * _pointSet;
};

#endif // SAMPLEGENERATORFACTORY_H
