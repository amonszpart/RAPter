#ifndef MYSCENE_H
#define MYSCENE_H

#include <QGraphicsScene>
#include <QWheelEvent>

#include "project.h"
#include "typesGL.h"
#include "types.h"

class MyScene : public QGraphicsScene
{
    Q_OBJECT
public:
    explicit MyScene(QObject *parent = 0);

    inline void setProject(InputGen::Application::Project* p)
    {
        _project = p;
        update();
    }

protected:
    virtual void drawBackground(QPainter*, const QRectF &);
    virtual void wheelEvent    (QGraphicsSceneWheelEvent * );

private:
    void setStates();

signals:

public slots:
    inline void projectUpdated(){
        update();
    }

private:
    InputGen::Application::Project *_project;
    float _zoom;

};

#endif // MYSCENE_H
