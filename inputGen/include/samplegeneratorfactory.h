#ifndef SAMPLEGENERATORFACTORY_H
#define SAMPLEGENERATORFACTORY_H

#include <QDockWidget>

namespace Ui {
class SampleGeneratorFactory;
}

class SampleGeneratorFactory : public QDockWidget
{
    Q_OBJECT

public:
    explicit SampleGeneratorFactory(QWidget *parent = 0);
    ~SampleGeneratorFactory();

private:
    Ui::SampleGeneratorFactory *ui;
};

#endif // SAMPLEGENERATORFACTORY_H
