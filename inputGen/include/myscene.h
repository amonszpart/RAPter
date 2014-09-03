#ifndef MYSCENE_H
#define MYSCENE_H

#include <QGraphicsScene>
#include <QWheelEvent>

#include "primitive.h"
#include "sampler.h"
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
                              InputGen::Application::Sampler* sampler){
        delete (_sampler);

        // both input pointer are invalid out of this function
        _pointSet  = set;

        _sampler = sampler!=NULL ? sampler->copy() : NULL;
        update();
    }

private:
    std::vector< InputGen::Application::Primitive > *_pSet;
    InputGen::Application::PointSet *_pointSet;
    InputGen::Application::Sampler *_sampler;
    float _zoom;

};

#endif // MYSCENE_H
