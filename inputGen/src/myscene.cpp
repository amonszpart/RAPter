#include "myscene.h"


//#include "Eigen/OpenGLSupport"
#include <QtOpenGL>

#include <iostream>


// Utility functions
void qgluPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar)
{
    const GLdouble ymax = zNear * tan(fovy * M_PI / 360.0);
    const GLdouble ymin = -ymax;
    const GLdouble xmin = ymin * aspect;
    const GLdouble xmax = ymax * aspect;
    glFrustum(xmin, xmax, ymin, ymax, zNear, zFar);
}


//! Convenience wrapper to call OpenGL command with compile-time vertex definitions
template <typename _Scalar>
struct GLDisplayFunctor{
    static inline void displayVertex(const _Scalar *) {}
};

template <>
void
GLDisplayFunctor<double>::displayVertex(const double* data){
    glVertex3dv(data);
}

template <>
void
GLDisplayFunctor<float>::displayVertex(const float* data){
    glVertex3fv(data);
}


///////////////////////////////////////////////////////////////////////////////



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
    //float width = float(painter->device()->width());
    //float height = float(painter->device()->height());

    painter->beginNativePainting();
    setStates();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    //qgluPerspective(60.0, width / height, 0.01, 15.0);

    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    static GLdouble invertY [16]  = {
        1.0, 0.0, 0.0, 0.0,
        0.0,-1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    };

    glMultMatrixd(invertY);

    glBegin(GL_LINES);
    for(unsigned int i = 0; i != _pSet.size(); i++)
        _pSet[i].displayAsLine<GLDisplayFunctor>();
    glEnd();


    painter->endNativePainting();
}
