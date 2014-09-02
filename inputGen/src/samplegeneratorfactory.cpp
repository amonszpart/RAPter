#include "samplegeneratorfactory.h"
#include "ui_samplegeneratorfactory.h"

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
            InputGen::PrimitiveSampleGenerator<InputGen::Application::Scalar,
                    InputGen::Application::GLDisplayFunctor > lgen;

            // set generator parameters
            lgen.spacing = ui->_generatorPrimitivesParamSpacing->value();

            // generate points
            lgen.generateSamples(*_pointSet, *_pSet);

            emit samplesChanged(_pointSet, &lgen);

            break;
        }
        default:
            emit samplesChanged(_pointSet, NULL);
        };

    }
    //emit samplesChanged(_pointSet, _generator);
}
