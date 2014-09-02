#include "samplegeneratorfactory.h"
#include "ui_samplegeneratorfactory.h"

SampleGeneratorFactory::SampleGeneratorFactory(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::SampleGeneratorFactory)
{
    ui->setupUi(this);
}

SampleGeneratorFactory::~SampleGeneratorFactory()
{
    delete ui;
}
