#include "displacementfactory.h"
#include "ui_displacementfactory.h"

#include "types.h"

#include <iostream>

DisplacementFactory::DisplacementFactory(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::DisplacementFactory),
    _project(NULL)
{
    ui->setupUi(this);
}

DisplacementFactory::~DisplacementFactory()
{
    delete ui;
}

void
DisplacementFactory::setProject(InputGen::Application::Project *p){
    _project = p;

    //todo: clear UI

    if ( p == NULL ){
        ;
    }else{
        // todo: load layers from the project

    }
}

void
DisplacementFactory::addLayerTriggerred(){

    typedef InputGen::Application::Scalar Scalar;
    typedef InputGen::Application::Primitive::vec vec;

    int ktype = ui->_displacementLayerAddCombo->currentIndex();

    if (ktype >= InputGen::DISPLACEMENT_KERNEL_TYPE::INVALID_KERNEL){
        std::cerr << "This should nerver happen... "
                  << __FILE__ << " "
                  << __LINE__ << std::endl;
        return;
    }

    int layerId = _project->nbDisplacementLayers();

    switch(ktype){
    case InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM:
    {
        InputGen::RandomDisplacementKernel<Scalar>* lkernel =
                new InputGen::RandomDisplacementKernel<Scalar>;

        _project->addDisplacementLayer(lkernel);

        lkernel->generateDisplacement( _project->displacementLayerPtr(layerId),
                                      _project->samples,
                                      _project->primitives);
        break;
    }
    case InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_BIAS:
    {
        InputGen::BiasDisplacementKernel<Scalar>* lkernel =
                new InputGen::BiasDisplacementKernel<Scalar>;
        lkernel->bias = 0.2;

        _project->addDisplacementLayer(lkernel);

        lkernel->generateDisplacement( _project->displacementLayerPtr(layerId),
                                      _project->samples,
                                      _project->primitives);
        break;
    }
    }



}

void
DisplacementFactory::refreshFromView(){
    std::cout << "refresh from view" << std::endl;
}
