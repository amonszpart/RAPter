#include "samplerfactory.h"
#include "ui_samplerfactory.h"

#include "types.h"

#include <iostream>

using namespace std;

SamplerFactory::SamplerFactory(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::SamplerFactory),
    _project(NULL)
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

    if (_project != NULL){

        _project->samples.clear();

        switch (ui->toolBox->currentIndex () ){
        case GEN_FROM_PRIMITIVE:
        {
            InputGen::PrimitiveSampler<InputGen::Application::Scalar,
                    InputGen::Application::GLDisplayFunctor > lgen;

            // set sampler parameters
            lgen.spacing = ui->_samplerPrimitivesParamSpacing->value();

            // generate points
            lgen.generateSamples(_project->samples, _project->primitives);

            _project->copySampler(&lgen);

            break;
        }
        default:
            ;
        };

    }
    emit projectUpdated();
}
