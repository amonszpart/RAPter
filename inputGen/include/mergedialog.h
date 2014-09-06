#ifndef MERGEDIALOG_H
#define MERGEDIALOG_H

#include <QDialog>

#include "types.h"

namespace Ui {
class MergeDialog;
}

class MergeDialog : public QDialog
{
    Q_OBJECT

public:
    explicit MergeDialog(QWidget *parent = 0);
    ~MergeDialog();

    void getParams(InputGen::MergeParam<InputGen::Application::Scalar>& params);

private:
    Ui::MergeDialog *ui;
};

#endif // MERGEDIALOG_H
