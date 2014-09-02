#include "samplegeneratorfactory.h"
#include "ui_samplegeneratorfactory.h"

#include "samplegenerator.h"
#include "types.h"

#include <iostream>

using namespace std;

SampleGeneratorFactory::SampleGeneratorFactory(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::SampleGeneratorFactory),
    _pSet(NULL),
    _pointSet(NULL)
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
    typedef InputGen::Application::Primitive::vec vec;

    if (_pSet != NULL && _pointSet != NULL){

        _pointSet->clear();

        switch (ui->toolBox->currentIndex () ){
        case GEN_FROM_PRIMITIVE:
        {
            // set generator parameters
            InputGen::PrimitiveSampleGenerator<InputGen::Application::Scalar> generator;
            generator.spacing = ui->_generatorPrimitivesParamSpacing->value();

            // generate points
            generator.generateSamples(*_pointSet, *_pSet);

            break;
        }
        default:
            break;
        };

    }

    // send the generated set
    emit samplesChanged(_pointSet);
}
