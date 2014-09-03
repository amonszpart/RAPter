#ifndef MYVIEW_H
#define MYVIEW_H

#include <QGraphicsView>

#include "primitive.h"
#include "types.h"
#include "typesGL.h"
#include "sampler.h"
#include "project.h"

class MyScene;

class MyView : public QGraphicsView
{
    Q_OBJECT
public:

    explicit MyView(QWidget *parent = 0);
    void setProject(InputGen::Application::Project* p);

signals:
    void projectUpdated();

public slots:

};

#endif // MYVIEW_H
