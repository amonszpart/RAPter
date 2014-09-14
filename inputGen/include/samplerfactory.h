#ifndef SAMPLERFACTORY_H
#define SAMPLERFACTORY_H

#include <QDockWidget>
#include <QtXml>

#include "sampler.h"
#include "types.h"
#include "typesGL.h"
#include "project.h"

namespace Ui {
class SamplerFactory;
}

class SamplerFactory : public QDockWidget
{
    Q_OBJECT

public:
    explicit SamplerFactory(QWidget *parent = 0);
    ~SamplerFactory();

    inline void setProject(InputGen::Application::Project*p, bool updateSamples = false)
    {
        _project = p;
        if (updateSamples) updateSampler();
    }

    void saveSamples(QDomDocument& doc, QDomElement& root) const;
    void loadSamples(QDomElement& root);

public slots:
    void updateSampler();

signals:
    void samplerUpdated();

private:
    enum SAMPLER_TYPE{
        GEN_FROM_PRIMITIVE = 0,
        GEN_FROM_PUNCTUAL  = 1
    };

    Ui::SamplerFactory *ui;
    InputGen::Application::Project * _project;
};

#endif // SAMPLERFACTORY_H
