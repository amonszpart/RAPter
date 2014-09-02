#include "samplegeneratorfactory.h"
#include "ui_samplegeneratorfactory.h"

#include <iostream>

using namespace std;

SampleGeneratorFactory::SampleGeneratorFactory(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::SampleGeneratorFactory),
    _pSet(NULL),
    _pointSet(new InputGen::Application::PointSet)
{
    ui->setupUi(this);
}

SampleGeneratorFactory::~SampleGeneratorFactory()
{
    delete ui;
}

void
SampleGeneratorFactory::updateGenerator()
{
    cout << "UpdateGenerator" << endl;
    typedef InputGen::Application::Primitive::vec vec;

    _pointSet->clear();

    if (_pSet != NULL){

        switch (ui->toolBox->currentIndex () ){
        case GEN_FROM_PRIMTIVE:
            _pointSet->push_back(vec(0.,0.,0.));
            break;
        default:
            break;
        };

    }

    // send the generated set
    emit samplesChanged(_pointSet);
}
