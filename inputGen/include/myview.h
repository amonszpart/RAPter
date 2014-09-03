#ifndef MYVIEW_H
#define MYVIEW_H

#include <QGraphicsView>

#include "primitive.h"
#include "types.h"
#include "typesGL.h"
#include "sampler.h"

class MyScene;

class MyView : public QGraphicsView
{
    Q_OBJECT
public:

    explicit MyView(QWidget *parent = 0);
    void setPrimitives(std::vector<InputGen::Application::Primitive> *s);

signals:
    void samplesChanged(InputGen::Application::PointSet *,
                        InputGen::Application::Sampler*);

public slots:

};

#endif // MYVIEW_H
