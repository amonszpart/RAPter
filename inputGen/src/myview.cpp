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
    setScene(new MyScene(this));
}


void
MyView::setPrimitives(const std::vector<InputGen::Application::Primitive> &s)
{
    MyScene* sc = dynamic_cast<MyScene*> (scene());
    if (sc)
        sc->setPrimitives(s);
    else
        std::cerr << "Unsupported scene type" << std::endl;
}
