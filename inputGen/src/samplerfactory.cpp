#include "samplerfactory.h"
#include "ui_samplerfactory.h"

#include "types.h"

#include <iostream>

using namespace std;

SamplerFactory::SamplerFactory(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::SamplerFactory),
    _pSet(NULL),
    _pointSet(NULL)
{
    ui->setupUi(this);
}

SamplerFactory::~SamplerFactory()
{
    delete ui;
}


void
SamplerFactory::updateSampler()
{
    typedef InputGen::Application::Primitive::vec vec;

    if (_pSet != NULL && _pointSet != NULL){

        _pointSet->clear();

        switch (ui->toolBox->currentIndex () ){
        case GEN_FROM_PRIMITIVE:
        {
            InputGen::PrimitiveSampler<InputGen::Application::Scalar,
                    InputGen::Application::GLDisplayFunctor > lgen;

            // set sampler parameters
            lgen.spacing = ui->_samplerPrimitivesParamSpacing->value();

            // generate points
            lgen.generateSamples(*_pointSet, *_pSet);

            emit samplesChanged(_pointSet, &lgen);

            break;
        }
        default:
            emit samplesChanged(_pointSet, NULL);
        };

    }
    //emit samplesChanged(_pointSet, _sampler);
}
