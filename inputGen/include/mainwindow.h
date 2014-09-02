#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "ui_mainwindow.h"
#include <QMainWindow>


class MainWindow : public QMainWindow, private Ui::MainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_actionLoad_SVG_triggered();

protected:
    virtual void closeEvent(QCloseEvent*);

private:
    void readSettings();
};

#endif // MAINWINDOW_H
