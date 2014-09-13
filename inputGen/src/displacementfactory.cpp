#include "displacementfactory.h"
#include "ui_displacementfactory.h"

#include "types.h"

#include <iostream>

#include <QTableWidgetItem>

DisplacementFactory::DisplacementFactory(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::DisplacementFactory),
    _project(NULL)
{
    ui->setupUi(this);
    currentLayerChanged();
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
DisplacementFactory::savekernels(QDomDocument &doc, QDomElement &root) const
{
    typedef InputGen::Application::Scalar S;
    typedef InputGen::Application::Project::SampleContainer    SC;
    typedef InputGen::Application::Project::PrimitiveContainer PC;

    if (_project != NULL){
        for (unsigned int i = 0; i != _project->nbDisplacementLayers(); i++)
        {
            InputGen::Application::Project::DisplacementKernel* kernel = _project->displacementKernel(i);

            QDomElement kernelElement = doc.createElement( "kernel" );
            kernelElement.setAttribute("typeId", int( kernel->type ));
            kernelElement.setAttribute("enabled", int( _project->isDisplacementLayerEnabled(i) ) );

            switch(kernel->type){
            case::InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_BIAS:
            {

                InputGen::BiasDisplacementKernel<S,SC,PC>* lkernel =
                        dynamic_cast<InputGen::BiasDisplacementKernel<S,SC,PC>*> (kernel);
                if(lkernel == NULL) {
                    std::cerr << "This should nerver happen... "
                              << __FILE__ << " "
                              << __LINE__ << std::endl;
                    break;
                }

                kernelElement.setAttribute("bias", QString::number( lkernel->bias ));
                break;
            }
            case::InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM_UNIFORM:
            {
                InputGen::UniformRandomDisplacementKernel<S,SC,PC>* lkernel =
                        dynamic_cast<InputGen::UniformRandomDisplacementKernel<S,SC,PC>*> (kernel);
                if(lkernel == NULL) {
                    std::cerr << "This should nerver happen... "
                              << __FILE__ << " "
                              << __LINE__ << std::endl;
                    break;
                }

                kernelElement.setAttribute("distributionMin", QString::number( lkernel->distributionMin() ));
                kernelElement.setAttribute("distributionMax", QString::number( lkernel->distributionMax() ));
                break;
            }
            case::InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM_NORMAL:
            {
                InputGen::NormalRandomDisplacementKernel<S,SC,PC>* lkernel =
                        dynamic_cast<InputGen::NormalRandomDisplacementKernel<S,SC,PC>*> (kernel);
                if(lkernel == NULL) {
                    std::cerr << "This should nerver happen... "
                              << __FILE__ << " "
                              << __LINE__ << std::endl;
                    break;
                }

                kernelElement.setAttribute("distributionMean", QString::number( lkernel->distributionMean() ));
                kernelElement.setAttribute("distributionStdDev", QString::number( lkernel->distributionStdDev() ));
                break;
            }
            }

            root.appendChild(kernelElement);
        }
    }
}


void
DisplacementFactory::addLayerTriggerred(){
    std::cout << "addLayerTriggerred" << std::endl;

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

    typedef InputGen::Application::Project::SampleContainer    SampleContainer;
    typedef InputGen::Application::Project::PrimitiveContainer PrimitiveContainer;

    InputGen::Application::Project::DisplacementKernel *kernel = NULL;

    switch(ktype){
    case InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM_UNIFORM:
        kernel = new InputGen::UniformRandomDisplacementKernel<Scalar,SampleContainer,PrimitiveContainer>;
        break;
    case InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM_NORMAL:
        kernel = new InputGen::NormalRandomDisplacementKernel<Scalar,SampleContainer,PrimitiveContainer>;
        break;
    case InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_BIAS:
        kernel = new InputGen::BiasDisplacementKernel<Scalar,SampleContainer,PrimitiveContainer>;
        break;
    }


    // Update UI
    if (kernel != NULL)
    {
        // add layer and compute it
        _project->addDisplacementLayer(kernel);

        configureFromUI(kernel);

        kernel->generateDisplacement( _project->displacementLayerPtr(layerId),
                                      _project->samples,
                                      _project->primitives);
        // add new row in the table
        int rowId = ui->_displacementLayerTable->rowCount();
        ui->_displacementLayerTable->setRowCount(rowId+1);

        // Row title
        QTableWidgetItem *_qtablewidgetitem1 = new QTableWidgetItem();
        _qtablewidgetitem1->setText(QString::number(rowId));
        ui->_displacementLayerTable->setVerticalHeaderItem(rowId, _qtablewidgetitem1);

        QTableWidgetItem *_qtablewidgetitem2 = new QTableWidgetItem();
        _qtablewidgetitem2->setText(QString::fromStdString(kernel->name));
        ui->_displacementLayerTable->setItem(rowId, 0, _qtablewidgetitem2);

        QTableWidgetItem *_qtablewidgetitem3 = new QTableWidgetItem();
        _qtablewidgetitem3->setCheckState(Qt::Checked);
        ui->_displacementLayerTable->setItem(rowId, 1, _qtablewidgetitem3);

        ui->_displacementLayerTable->selectRow(rowId);

        // implicitely call during the row selection (above)
        //emit projectUpdated();
    }


}

void
DisplacementFactory::itemChanged(QTableWidgetItem *item){
    std::cout << "itemChanged()" << std::endl;

    // the current cell is in the "enable" column
    if (ui->_displacementLayerTable->column(item) != 1 || _project == NULL)
        return;

    int layerId = ui->_displacementLayerTable->row(item);
    if (layerId < 0) return;

    _project->enableDisplacementLayer(layerId, item->checkState()==Qt::Checked);

    emit projectUpdated();
}

int
DisplacementFactory::getSelectedLayerFromUI(){

    QList<QTableWidgetItem*> selectedItems = ui->_displacementLayerTable->selectedItems();
    if (selectedItems.size() == 0)
        return -1;

    return ui->_displacementLayerTable->row(selectedItems.front());
}

bool
DisplacementFactory::configureFromUI(InputGen::Application::Project::DisplacementKernel *kernel){
    typedef InputGen::Application::Scalar S;
    typedef InputGen::Application::Project::SampleContainer    SC;
    typedef InputGen::Application::Project::PrimitiveContainer PC;

    bool needUpdate = false;

    switch(kernel->type){
    case::InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_BIAS:
    {

        InputGen::BiasDisplacementKernel<S,SC,PC>* lkernel =
                dynamic_cast<InputGen::BiasDisplacementKernel<S,SC,PC>*> (kernel);
        if(lkernel == NULL) {
            std::cerr << "This should nerver happen... "
                      << __FILE__ << " "
                      << __LINE__ << std::endl;
            break;
        }

        // set parameters
        if (lkernel->bias != ui->_displacementParamBiasValue->value()){
            needUpdate = true;
            lkernel->bias = ui->_displacementParamBiasValue->value();
        }

        break;
    }
    case::InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM_UNIFORM:
    {
        InputGen::UniformRandomDisplacementKernel<S,SC,PC>* lkernel =
                dynamic_cast<InputGen::UniformRandomDisplacementKernel<S,SC,PC>*> (kernel);
        if(lkernel == NULL) {
            std::cerr << "This should nerver happen... "
                      << __FILE__ << " "
                      << __LINE__ << std::endl;
            break;
        }

        // set parameters
        if (lkernel->distributionMin() != ui->_displacementParamRandomUniformMinValue->value() ||
            lkernel->distributionMax() != ui->_displacementParamRandomUniformMaxValue->value()){
            needUpdate = true;
            lkernel->setDistributionRange(ui->_displacementParamRandomUniformMinValue->value(),
                                          ui->_displacementParamRandomUniformMaxValue->value());
        }

        break;
    }
    case::InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM_NORMAL:
    {
        InputGen::NormalRandomDisplacementKernel<S,SC,PC>* lkernel =
                dynamic_cast<InputGen::NormalRandomDisplacementKernel<S,SC,PC>*> (kernel);
        if(lkernel == NULL) {
            std::cerr << "This should nerver happen... "
                      << __FILE__ << " "
                      << __LINE__ << std::endl;
            break;
        }

        // set parameters
        if (lkernel->distributionMean()   != ui->_displacementParamRandomNormalMeanValue->value() ||
            lkernel->distributionStdDev() != ui->_displacementParamRandomNormalStddevValue->value()){
            needUpdate = true;
            lkernel->setDistributionProperties(ui->_displacementParamRandomNormalMeanValue->value(),
                                               ui->_displacementParamRandomNormalStddevValue->value());
        }

        break;
    }
    }

    return needUpdate;
}

// Hide show the correct parameter groups
void
DisplacementFactory::currentLayerChanged(){
    typedef InputGen::Application::Scalar S;
    typedef InputGen::Application::Project::SampleContainer    SC;
    typedef InputGen::Application::Project::PrimitiveContainer PC;

    std::cout << "currentLayerChanged()" << std::endl;

    int layerId = getSelectedLayerFromUI();

    // hide all param groups
    if (layerId == -1){
        ui->_displacementParamBiasGroup->hide();
        ui->_displacementParamRandomUniformGroup->hide();
        ui->_displacementParamRandomNormalGroup->hide();
    }

    if (_project == NULL) return;

    InputGen::Application::Project::DisplacementKernel *kernel = _project->displacementKernel(layerId);
    if(kernel == NULL) return;


    switch(kernel->type){
    case::InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_BIAS:
    {
        InputGen::BiasDisplacementKernel<S,SC,PC>* lkernel =
                dynamic_cast<InputGen::BiasDisplacementKernel<S,SC,PC>*> (kernel);
        if(lkernel == NULL) {
            std::cerr << "This should nerver happen... "
                      << __FILE__ << " "
                      << __LINE__ << std::endl;
            return;
        }

        ui->_displacementParamBiasValue->setValue(lkernel->bias);

        ui->_displacementParamBiasGroup->show();
        ui->_displacementParamRandomUniformGroup->hide();
        ui->_displacementParamRandomNormalGroup->hide();
        break;
    }
    case::InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM_UNIFORM:
    {
        InputGen::UniformRandomDisplacementKernel<S,SC,PC>* lkernel =
                dynamic_cast<InputGen::UniformRandomDisplacementKernel<S,SC,PC>*> (kernel);
        if(lkernel == NULL) {
            std::cerr << "This should nerver happen... "
                      << __FILE__ << " "
                      << __LINE__ << std::endl;
            return;
        }

        ui->_displacementParamRandomUniformMinValue->setValue(lkernel->distributionMin());
        ui->_displacementParamRandomUniformMaxValue->setValue(lkernel->distributionMax());

        ui->_displacementParamBiasGroup->hide();
        ui->_displacementParamRandomUniformGroup->show();
        ui->_displacementParamRandomNormalGroup->hide();
        break;
    }
    case::InputGen::DISPLACEMENT_KERNEL_TYPE::DISPLACEMENT_RANDOM_NORMAL:
    {
        InputGen::NormalRandomDisplacementKernel<S,SC,PC>* lkernel =
                dynamic_cast<InputGen::NormalRandomDisplacementKernel<S,SC,PC>*> (kernel);
        if(lkernel == NULL) {
            std::cerr << "This should nerver happen... "
                      << __FILE__ << " "
                      << __LINE__ << std::endl;
            return;
        }

        ui->_displacementParamRandomNormalMeanValue->setValue(lkernel->distributionMean());
        ui->_displacementParamRandomNormalStddevValue->setValue(lkernel->distributionStdDev());

        ui->_displacementParamBiasGroup->hide();
        ui->_displacementParamRandomUniformGroup->hide();
        ui->_displacementParamRandomNormalGroup->show();
        break;
    }
    }
}

void
DisplacementFactory::refreshFromView(){
    std::cout << "refreshFromView()" << std::endl;

    int layerId = getSelectedLayerFromUI();

    if (_project == NULL || layerId < 0) return;

    InputGen::Application::Project::DisplacementKernel *kernel = _project->displacementKernel(layerId);
    if(kernel == NULL) return;

    // recompute only when required
    if (configureFromUI(kernel))
        recomputeDisplacementLayer(layerId, false);

    emit projectUpdated();
}

void
DisplacementFactory::recomputeDisplacement(){
    std::cout << "recomputeDisplacement()" << std::endl;

    if(_project == NULL) return;

    for(unsigned int i = 0; i<_project->nbDisplacementLayers(); i++)
        recomputeDisplacementLayer(i, false);

    // this signal must be emitted even without displacement layers,
    // to transmit the update information to other application elements
    emit projectUpdated();
}

void
DisplacementFactory::recomputeDisplacementLayer(int layerId, bool triggerSignal){
    std::cout << "recomputeDisplacementLayer()" << std::endl;

    if(_project == NULL) return;

    InputGen::Application::Project::DisplacementKernel *kernel = _project->displacementKernel(layerId);
    if(kernel == NULL) return;


    kernel->generateDisplacement( _project->displacementLayerPtr(layerId),
                                  _project->samples,
                                  _project->primitives);

    if(triggerSignal)
        emit projectUpdated();

}
