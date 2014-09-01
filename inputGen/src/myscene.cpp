#include "myscene.h"


//#include "Eigen/OpenGLSupport"
#include <QtOpenGL>

#include <iostream>



void qgluPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar)
{
    const GLdouble ymax = zNear * tan(fovy * M_PI / 360.0);
    const GLdouble ymin = -ymax;
    const GLdouble xmin = ymin * aspect;
    const GLdouble xmax = ymax * aspect;
    glFrustum(xmin, xmax, ymin, ymax, zNear, zFar);
}



MyScene::MyScene(QObject *parent) :
    QGraphicsScene(parent)
{
    //setStates();
}

void
MyScene::setStates(){
    glClearColor(1.f, 1.f, 1.f, 1.0f);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_LIGHTING);
    //glEnable(GL_COLOR_MATERIAL);
    //glEnable(GL_TEXTURE_2D);
    glEnable(GL_NORMALIZE);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
}

void
MyScene::drawBackground(QPainter *painter, const QRectF &rect){
    float width = float(painter->device()->width());
    float height = float(painter->device()->height());

    painter->beginNativePainting();
    setStates();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    //qgluPerspective(60.0, width / height, 0.01, 15.0);

    glMatrixMode(GL_MODELVIEW);

    glBegin(GL_LINES);
    for(unsigned int i = 0; i != _pSet.size(); i++){
        InputGen::Application::Primitive& p = _pSet[i];
        glVertex3dv((p.coord() - p.dir() * 1000.f).eval().data());
        glVertex3dv((p.coord() + p.dir() * 1000.f).eval().data());

    }
    glEnd();


    painter->endNativePainting();
}
