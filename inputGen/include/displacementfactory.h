#ifndef DISPLACEMENTFACTORY_H
#define DISPLACEMENTFACTORY_H

#include <QDockWidget>

#include "project.h"

namespace Ui {
class DisplacementFactory;
}

class DisplacementFactory : public QDockWidget
{
    Q_OBJECT

public:
    explicit DisplacementFactory(QWidget *parent = 0);
    ~DisplacementFactory();

    void setProject(InputGen::Application::Project* p);

signals:
    void projectUpdated();

public slots:
    void recomputeDisplacement();
    void recomputeDisplacementLayer(int layerId, bool triggerSignal=true);

private slots:
    //! \brief Add a new layer according to the
    void addLayerTriggerred();
    //! \brief Refresh displacement values wrt to UI
    void refreshFromView();

    //! \brief Triggered by UI
    void currentLayerChanged();

private:
    int getSelectedLayerFromUI();
    void configureFromUI(InputGen::Application::Project::DisplacementKernel* kernel);

private:
    Ui::DisplacementFactory *ui;

    InputGen::Application::Project* _project;
};

#endif // DISPLACEMENTFACTORY_H
