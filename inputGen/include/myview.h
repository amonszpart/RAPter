#ifndef MYVIEW_H
#define MYVIEW_H

#include <QGraphicsView>

#include "primitive.h"
#include "types.h"

class MyView : public QGraphicsView
{
    Q_OBJECT
public:

    explicit MyView(QWidget *parent = 0);
    inline void setPrimitives(
            const std::vector< InputGen::Application::Primitive >&s)
    {
        _pSet = s;
    }

signals:

public slots:


private:
    std::vector< InputGen::Application::Primitive > _pSet;

};

#endif // MYVIEW_H
