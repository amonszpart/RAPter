#ifndef SAMPLERFACTORY_H
#define SAMPLERFACTORY_H

#include <QDockWidget>
#include "sampler.h"
#include "types.h"
#include "typesGL.h"

namespace Ui {
class SamplerFactory;
}

class SamplerFactory : public QDockWidget
{
    Q_OBJECT

public:
    explicit SamplerFactory(QWidget *parent = 0);
    ~SamplerFactory();

    inline void setPrimitives(
            std::vector< InputGen::Application::Primitive >*s)
    {
        _pSet = s;
        updateSampler();
    }

    inline void setPoints( InputGen::Application::PointSet*s)
    {
        _pointSet = s;
    }

public slots:
    void updateSampler();

signals:
    //! \warning Pointers are invalid outside of the connected slots
    void samplesChanged(InputGen::Application::PointSet *,
                        InputGen::Application::Sampler*);

private:
    enum SAMPLER_TYPE{
        GEN_FROM_PRIMITIVE = 0,
        GEN_FROM_PONCTUAL  = 1
    };

    Ui::SamplerFactory *ui;
    std::vector< InputGen::Application::Primitive > *_pSet;
    InputGen::Application::PointSet * _pointSet;
};

#endif // SAMPLERFACTORY_H
