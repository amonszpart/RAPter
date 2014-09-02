#ifndef MYSCENE_H
#define MYSCENE_H

#include <QGraphicsScene>
#include <QWheelEvent>

#include "primitive.h"
#include "samplegenerator.h"
#include "typesGL.h"
#include "types.h"

class MyScene : public QGraphicsScene
{
    Q_OBJECT
public:
    explicit MyScene(QObject *parent = 0);
    inline void setPrimitives(
            std::vector< InputGen::Application::Primitive >*s)
    {
        _pSet = s;
        update();
    }

protected:
    virtual void drawBackground(QPainter*, const QRectF &);
    virtual void wheelEvent    (QGraphicsSceneWheelEvent * );

private:
    void setStates();

signals:

public slots:
    inline void updateSamples(InputGen::Application::PointSet *set,
                              InputGen::Application::SampleGenerator* generator){
        delete (_generator);

        // both input pointer are invalid out of this function
        _pointSet  = set;

        _generator = generator!=NULL ? generator->copy() : NULL;
        update();
    }

private:
    std::vector< InputGen::Application::Primitive > *_pSet;
    InputGen::Application::PointSet *_pointSet;
    InputGen::Application::SampleGenerator *_generator;
    float _zoom;

};

#endif // MYSCENE_H
