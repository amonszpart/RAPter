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
    emit samplerUpdated();
}

void
SamplerFactory::saveSamples(QDomDocument& doc, QDomElement& root) const
{
    typedef InputGen::Application::Primitive::vec vec;

    if (_project != NULL){

        QDomElement samplerElement = doc.createElement( "sampler" );
        samplerElement.setAttribute("typeId", QString::number(int(ui->toolBox->currentIndex ())) );

        switch (ui->toolBox->currentIndex () ){
        case GEN_FROM_PRIMITIVE:
        {
            samplerElement.setAttribute("spacing", QString::number(ui->_samplerPrimitivesParamSpacing->value()));
            break;
        }
        default:
            return;
        };

        root.appendChild(samplerElement);
    }
}

void
SamplerFactory::loadSamples(QDomElement &root){
    std::cout << "Loading samples (load only first sampler)" << std::endl;

    QDomNode samplerNode = root.firstChild();
    while(! samplerNode.isNull()){
        QDomElement samplerElement = samplerNode.toElement(); // try to convert the node to an element.
        if(!samplerElement.isNull() && samplerElement.tagName().compare(QString("sampler")) == 0) {
            int samplerId = samplerElement.attribute("typeId").toInt();
            ui->toolBox->setCurrentIndex(samplerId);

            switch (samplerId){
            case GEN_FROM_PRIMITIVE:
            {
                ui->_samplerPrimitivesParamSpacing->setValue(samplerElement.attribute("spacing").toDouble());
                break;
            }
            default:
                return;
            };
        }
        samplerNode = samplerNode.nextSibling();
    }
}
