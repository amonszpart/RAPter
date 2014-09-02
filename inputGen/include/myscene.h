#ifndef MYSCENE_H
#define MYSCENE_H

#include <QGraphicsScene>
#include <QWheelEvent>

#include "primitive.h"
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
    inline void updateSamples(InputGen::Application::PointSet *set){
        _pointSet = set;
        update();
    }

private:
    std::vector< InputGen::Application::Primitive > *_pSet;
    InputGen::Application::PointSet *_pointSet;
    float _zoom;

};

#endif // MYSCENE_H
