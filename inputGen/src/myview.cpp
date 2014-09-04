#include "myview.h"
#include "myscene.h"

#include <iostream>

#include <QGLWidget>

using InputGen::Application::Primitive;

MyView::MyView(QWidget *parent) :
    QGraphicsView(parent)
{
    setViewport(new QGLWidget);
    setViewportUpdateMode(FullViewportUpdate);

    MyScene* scene = new MyScene(this);
    setScene(scene);

    connect(this, SIGNAL(projectUpdated()),
            scene, SLOT(projectUpdated()));
}


void
MyView::setProject(InputGen::Application::Project* p)
{
    MyScene* sc = dynamic_cast<MyScene*> (scene());
    if (sc)
        sc->setProject(p);
    else
        std::cerr << "Unsupported scene type" << std::endl;
}
