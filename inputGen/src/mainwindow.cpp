#include "mainwindow.h"

#include <iostream>

#include <QFile>
#include <QFileDialog>
#include <QSettings>
#include <QXmlStreamReader>

using std::cout;
using std::endl;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{
    setupUi(this);
    readSettings();
}


void
MainWindow::readSettings()
{
    QSettings settings("SmartGeometryProcessingGroup", "GlobFitInputGen");
    restoreGeometry(settings.value("MainWindow/geometry").toByteArray());
    restoreState(settings.value("MainWindow/windowState").toByteArray());
}

void
MainWindow::closeEvent(QCloseEvent *event)
{
    QSettings settings("SmartGeometryProcessingGroup", "GlobFitInputGen");
    settings.setValue("MainWindow/geometry", saveGeometry());
    settings.setValue("MainWindow/windowState", saveState());
    QMainWindow::closeEvent(event);
}

MainWindow::~MainWindow()
{
}

void MainWindow::on_actionLoad_SVG_triggered()
{
    QSettings settings;
    QString defaultPath = settings.value("Path/svgOpen").toString();

    QString path =
    QFileDialog::getOpenFileName(this,
                                 tr("Load svg file"),
                                 defaultPath,
                                 "Scalable Vector Graphics (*.svg)");
    if (! path.isNull()){
        QFile input(path);
        if (input.open(QIODevice::ReadOnly)) {
            settings.setValue("Path/svgOpen", path);

            QByteArray ba = input.readAll();
            input.close();

            // Read the svg by extracting the paths
            QXmlStreamReader reader(ba);
            while (reader.readNextStartElement()) {

                while (! ( reader.hasError() || reader.atEnd() ) ) {
                    if (reader.name().compare("path") == 0){

                        // Extract path coordinates (attribute d), and split them using ' '
                        QStringList attrList (reader.attributes().value("d").toString().split(" "));

                        // The path is not empty and contains coordinates
                        if( attrList.size() > 1 &&
                                attrList.front().compare("M") == 0){

                            cout << "Reading new path" << endl;

                            attrList.pop_front();

                            QStringList::const_iterator constIterator;
                            for (constIterator = attrList.constBegin(); constIterator != attrList.constEnd();
                                 ++constIterator)
                                if ((*constIterator).contains(','))
                                    cout << (*constIterator).toLocal8Bit().constData() << endl;
                        }
                    }
                    reader.readNext();
                }
            }
        }
    }
}
