#include "mergedialog.h"
#include "ui_mergedialog.h"

MergeDialog::MergeDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MergeDialog)
{
    ui->setupUi(this);
}

MergeDialog::~MergeDialog()
{
    delete ui;
}


void
MergeDialog::getParams(InputGen::MergeParam<InputGen::Application::Scalar> &params){
    params.useAngular = ui->_angularGroup->isChecked();
    params.angleRef   = ui->_angularGroupAngleValue->value()/InputGen::Application::Scalar(180.)*M_PI;
    params.usePeriodicAngles = ui->_angularGroupPeriodic->isChecked();
}
