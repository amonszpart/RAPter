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

private slots:
    //! \brief Add a new layer according to the
    void addLayerTriggerred();
    //! \brief Refresh displacement values wrt to UI
    void refreshFromView();

private:
    Ui::DisplacementFactory *ui;

    InputGen::Application::Project* _project;
};

#endif // DISPLACEMENTFACTORY_H
