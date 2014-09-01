#ifndef MYSCENE_H
#define MYSCENE_H

#include <QGraphicsScene>

#include "primitive.h"
#include "types.h"

class MyScene : public QGraphicsScene
{
    Q_OBJECT
public:
    explicit MyScene(QObject *parent = 0);
    inline void setPrimitives(
            const std::vector< InputGen::Application::Primitive >&s)
    {
        _pSet = s;
        update();
    }

protected:
    virtual void drawBackground(QPainter*, const QRectF &);

private:
    void setStates();

signals:

public slots:

private:
    std::vector< InputGen::Application::Primitive > _pSet;

};

#endif // MYSCENE_H
