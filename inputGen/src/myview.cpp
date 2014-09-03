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

    connect(this, SIGNAL(samplesChanged(InputGen::Application::PointSet*,
                                        InputGen::Application::Sampler*)),
              scene, SLOT(updateSamples(InputGen::Application::PointSet*,
                                        InputGen::Application::Sampler*)));
}


void
MyView::setPrimitives(std::vector<InputGen::Application::Primitive> *s)
{
    MyScene* sc = dynamic_cast<MyScene*> (scene());
    if (sc)
        sc->setPrimitives(s);
    else
        std::cerr << "Unsupported scene type" << std::endl;
}
