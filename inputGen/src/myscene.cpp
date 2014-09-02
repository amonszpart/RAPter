#include "myscene.h"


//#include "Eigen/OpenGLSupport"
#include <QtOpenGL>

#include <iostream>

using namespace std;

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
    QGraphicsScene(parent),
    _pSet(NULL),
    _pointSet(NULL),
    _zoom(1.)
{
    //setStates();
}

void
MyScene::wheelEvent(QGraphicsSceneWheelEvent *event){

    int numDegrees = event->delta() / 8;
    int numSteps = numDegrees / 15;

    if (event->orientation() == Qt::Vertical) {
        _zoom += _zoom*0.1*float(numSteps);

        event->accept();
        update();
    }
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

    using InputGen::Application::Scalar;

    painter->beginNativePainting();
    setStates();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if(_pSet==NULL) return

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
    glScalef(_zoom, _zoom, _zoom);

    // display lines
    {
        glBegin(GL_LINES);
        std::vector< InputGen::Application::Primitive >::const_iterator it;
        for(it = _pSet->begin(); it != _pSet->end(); it++)
            (*it).displayAsLine<GLDisplayFunctor>();
        glEnd();
    }

    // display samples
    if (_pointSet != NULL){
        glPointSize(2.f);
        glBegin(GL_POINTS);
        InputGen::Application::PointSet::const_iterator it;
        for(it = _pointSet->begin(); it != _pointSet->end(); it++){
            GLDisplayFunctor<Scalar>::displayVertex((*it).data());
        }
        glEnd();
    }


    painter->endNativePainting();
}
