#include "mainwindow.h"

#include <iostream>
#include <math.h>       /* fmod */

#include <QFile>
#include <QFileDialog>
#include <QSettings>
#include <QXmlStreamReader>
#include <QTextStream>
#include <QDir>
#include <QtXml>

#include "mergedialog.h"

#include "primitive.h"
#include "types.h"

using std::cout;
using std::endl;
using InputGen::Application::Primitive;
using InputGen::Application::Project;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    _project(new Project)
{
    setupUi(this);
    readSettings();

    connect( this,             SIGNAL(currentProjectUpdated()),
             _samplerDoc,      SLOT(updateSampler()));

    connect( _samplerDoc,      SIGNAL(samplerUpdated()),
             _displacementDoc, SLOT  (recomputeDisplacement()));

    connect( _displacementDoc, SIGNAL(projectUpdated()),
             _graphicsView,    SIGNAL(projectUpdated()));

    _samplerDoc->setProject(_project);
    _graphicsView->setProject(_project);
    _displacementDoc->setProject(_project);

    //emit currentProjectUpdated();
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
    //using InputGen::LinearPrimitive<Scalar>/*::LINE_2D*/;
    //using InputGen::LinearPrimitive<Scalar>::vec;
    using InputGen::Application::Scalar;


    QSettings settings;
    QString defaultPath = settings.value("Path/svgOpen").toString();

    QString path =
    QFileDialog::getOpenFileName(this,
                                 tr("Load svg file"),
                                 defaultPath,
                                 "Scalable Vector Graphics (*.svg)");

    // can be replaced by a new project to handle multi-view
    _project->clear();

    if (! path.isNull()){
        QFile input(path);
        if (input.open(QIODevice::ReadOnly)) {
            settings.setValue("Path/svgOpen", path);

            QByteArray ba = input.readAll();
            input.close();

            // Read the svg by extracting the paths
            QXmlStreamReader reader(ba);
            while (reader.readNextStartElement()) {
                Primitive::vec globalTranslation(Primitive::vec::Zero());

                while (! ( reader.hasError() || reader.atEnd() ) ) {
                    if (reader.name().compare("g") == 0){
                        QStringList attrTrans (reader.attributes().value("transform").toString().split(QRegExp("[()]")));


                        globalTranslation = Primitive::vec::Zero();

                        while( attrTrans.size() > 0){
                            if (attrTrans.front().compare("translate") == 0 && attrTrans.size() > 1){
                                QStringList coordTr = attrTrans.at(1).split(',');
                                globalTranslation << coordTr.at(0).toDouble(),
                                                     coordTr.at(1).toDouble(),
                                                     Scalar(0.);
                                attrTrans.pop_front();
                                attrTrans.pop_front();
                                cout << "Detect new global translation" << globalTranslation.transpose() << endl;
                            }else
                                attrTrans.pop_front();
                        }

                    }
                    if (reader.name().compare("path") == 0){

                        // Extract path coordinates (attribute d), and split them using ' '
                        QStringList attrList  (reader.attributes().value("d").toString().split(" "));
                        QStringList attrTrans (reader.attributes().value("transform").toString().split(QRegExp("[()]")));

                        // The path is not empty and contains coordinates
                        if( attrList.size() > 1){
                            std::vector< Primitive > lines;
                            Primitive::vec localTranslation(Primitive::vec::Zero());

                            while( attrTrans.size() > 0){
                                if (attrTrans.front().compare("translate") == 0 && attrTrans.size() > 1){
                                    QStringList coordTr = attrTrans.at(1).split(',');
                                    localTranslation << coordTr.at(0).toDouble(),
                                                        coordTr.at(1).toDouble(),
                                                        Scalar(0.);
                                    attrTrans.pop_front();
                                    attrTrans.pop_front();
                                }else
                                    attrTrans.pop_front();
                            }

                            //cout << "Path translation = " << localTranslation.transpose() << endl;

                            if( attrList.front().compare("M") == 0  ||
                                    attrList.front().compare("m") == 0 ){

                                cout << "Reading new set of lines" << endl;

                                bool relativeCoord = attrList.front().compare("M");

                                attrList.pop_front();


                                // read a list of lines from the file.
                                // At this stage, only the position are extracted
                                QStringList::const_iterator constIterator;
                                for (constIterator = attrList.constBegin(); constIterator != attrList.constEnd();
                                     ++constIterator)
                                    if ((*constIterator).contains(',')){
                                        QStringList coordLists = (*constIterator).split(',');
                                        if (coordLists.size() == 2){
                                            Primitive line (Primitive::LINE_2D);
                                            line.setCoord(Primitive::vec(coordLists.at(0).toDouble(),
                                                                         coordLists.at(1).toDouble(),
                                                                         0));
                                            
                                            if ( relativeCoord && lines.size() != 0)
                                                line.setCoord(line.coord() + lines.back().coord());
                                            lines.push_back(line);
                                        }
                                    }

                                // compute direction of the n-1 lines
                                for (unsigned int i = 0; i < lines.size()-1; i++){
                                    Primitive& l1 = lines[i];
                                    const Primitive& l2 = lines[i+1];
                                    const Primitive::vec dortho= (l2.coord() - l1.coord()).normalized();
                                    const Primitive::vec normal ( dortho(1), - dortho(0), Scalar(0.) );
                                    const Primitive::vec2 dim (-(l2.coord() - l1.coord()).norm(), Scalar(0.));

                                    l1.setDim(dim);
                                    l1.setNormal(normal);

                                    cout << "l1=" << l1.coord().transpose() << endl;
                                    cout << "l2=" << l2.coord().transpose() << endl;

//                                    cout << l1.coord().transpose() << " "
//                                         << l1.dir().transpose() << " "
//                                         << l1.dim()(0)
//                                         << endl;
                                }

                                // second step of the extraction, we now compute the directions
                                // the path is considered as closed if the last character of the path='z'
                                if (attrList.last().compare("z") != 0) // open path, last line must be removed
                                    lines.pop_back();
                                else{
                                    Primitive& l1 = lines.back();
                                    const Primitive& l2 = lines.front();
                                    const Primitive::vec dortho= (l2.coord() - l1.coord()).normalized();
                                    const Primitive::vec normal ( dortho(1), - dortho(0), Scalar(0.) );
                                    const Primitive::vec2 dim (-(l2.coord() - l1.coord()).norm(), Scalar(0.));

                                    l1.setDim(dim);
                                    l1.setNormal(normal);

                                    cout << "l1=" << l1.coord().transpose() << endl;
                                    cout << "l2=" << l2.coord().transpose() << endl;

                                    //cout << l1.coord().transpose() << " "
                                    //     << l1.dir().transpose() << " "
                                    //     << l1.dim()(0)
                                    //     << endl;
                                }
                            }
                            // update with transformations
                            for (unsigned int i = 0; i < lines.size(); i++){
                                lines[i].setCoord(lines[i].coord() + localTranslation + globalTranslation);
                            }


                            _project->primitives.insert(_project->primitives.end(),
                                                        lines.begin(),
                                                        lines.end());
                        }
                    }
                    reader.readNext();
                }
            }
        }
    }

    emit currentProjectUpdated();

}


void MainWindow::on_actionLoad_Project_triggered()
{
    using InputGen::Application::Scalar;

    QSettings settings;
    QString defaultPath = settings.value("Path/prjOpen").toString();

    QString path =
    QFileDialog::getOpenFileName(this,
                                 tr("Load project file"),
                                 defaultPath,
                                 "InputGen project (*.prj)");

    // can be replaced by a new project to handle multi-view
    _project->clear();

    if (! path.isNull()){
        QFile input(path);
        if (input.open(QIODevice::ReadOnly)) {
            settings.setValue("Path/prjOpen", path);

            QDomDocument doc("mydocument");
            if (!doc.setContent(&input)) {
                input.close();
                return;
            }
            input.close();

            // print out the element names of all elements that are direct children
            // of the outermost element.
            QDomElement docElem = doc.documentElement();

            QDomNode n = docElem.firstChild();
            while(!n.isNull()) {
                QDomElement e = n.toElement(); // try to convert the node to an element.
                if(!e.isNull()) {
                    if (e.tagName().compare(QString("primitives")) == 0){
                        cout << "Loading primitives ..." << endl;

                        QDomNode primNode = e.firstChild();
                        while(! primNode.isNull()){
                            QDomElement primElement = primNode.toElement(); // try to convert the node to an element.
                            if(!primElement.isNull() && primElement.tagName().compare(QString("primitive")) == 0) {
                                // create a new primitive
                                Primitive line (Primitive::LINE_2D,
                                                primElement.attribute("uid").toInt(),
                                                primElement.attribute("did").toInt());

                                QStringList coordLists = primElement.attribute("pos").split(' ');
                                if (coordLists.size() == 3){
                                    line.setCoord(Primitive::vec(coordLists.at(0).toDouble(),
                                                                 coordLists.at(1).toDouble(),
                                                                 coordLists.at(2).toDouble()));
                                } else{
                                    std::cerr << "Unexpected error while reading primitive position" << endl;
                                    continue;
                                }

                                coordLists = primElement.attribute("dir").split(' ');
                                if (coordLists.size() == 3){
                                    line.setNormal(Primitive::vec(coordLists.at(0).toDouble(),
                                                                  coordLists.at(1).toDouble(),
                                                                  coordLists.at(2).toDouble()));
                                } else{
                                    std::cerr << "Unexpected error while reading primitive direction" << endl;
                                    continue;
                                }

                                coordLists = primElement.attribute("dim").split(' ');
                                if (coordLists.size() == 2){
                                    line.setDim(Primitive::vec2(coordLists.at(0).toDouble(),
                                                                coordLists.at(1).toDouble()));
                                } else{
                                    std::cerr << "Unexpected error while reading primitive dimension" << endl;
                                    continue;
                                }

                                _project->primitives.push_back(line);
                            }else
                                std::cerr << "Unsupported primitive type" << endl;

                            primNode = primNode.nextSibling();
                        }
                    }
                    else if (e.tagName().compare(QString("samplers")) == 0){
                        _samplerDoc->loadSamples(e);
                    }
                    else if (e.tagName().compare(QString("displacements")) == 0){
                        _displacementDoc->loadkernels(e);
                    }
                }
                n = n.nextSibling();
            }
        }
    }
    emit currentProjectUpdated();
}



void MainWindow::on_actionMerge_primitives_triggered()
{
    if (_project == NULL)
        return;

    using InputGen::Application::Scalar;

    MergeDialog* dialog = new MergeDialog;
    if (dialog->exec() == QDialog::Accepted){
        InputGen::MergeParam<Scalar> p;
        dialog->getParams(p);

        const Scalar angleTh = 10.e-3; //in rad

        // merge primitive according to dialog properties
        // do it here now, should be moved to dedicated class
        for(Project::PrimitiveContainer::iterator it1 = _project->primitives.begin();
            it1 != _project->primitives.end(); it1++){

            // we process only for the primitives that have not been merged previously
            if ((*it1).did() == (*it1).uid()){

                const Primitive::vec& normal = (*it1).normal();
                for(Project::PrimitiveContainer::iterator it2 = it1+1;
                    it2 != _project->primitives.end(); it2++){

                    if (p.useAngular){
                        Scalar angle = std::acos(normal.dot((*it2).normal()));
                        Scalar fm = std::fmod(angle, p.angleRef);
                        if (p.usePeriodicAngles){
                            if (std::abs(std::fmod(angle, p.angleRef)) <= angleTh)
                                (*it2).did() = (*it1).did();
                        }else if (std::abs(angle-p.angleRef) <= angleTh) // merge
                            (*it2).did() = (*it1).did();
                    }
                }
            }
        }
    }

    // clean
    delete dialog;
}

void MainWindow::on_actionSave_points_triggered()
{
    if (_project == NULL)
        return;

    QSettings settings;
    QString defaultPath = settings.value("Path/plySave").toString();

    QString path =
    QFileDialog::getSaveFileName(this,
                                 tr("Save cloud as PLY file"),
                                 defaultPath,
                                 "Pointcloud (*.ply)");


    if (! path.isNull()){
        writeSamples(path);
    }
}

void MainWindow::on_actionSave_primitives_triggered()
{
    if (_project == NULL)
        return;

    QSettings settings;
    QString defaultPath = settings.value("Path/priSave").toString();

    QString path =
    QFileDialog::getSaveFileName(this,
                                 tr("Save primitives as CSV file"),
                                 defaultPath,
                                 "CSV (*.csv)");


    if (! path.isNull()){
        writePrimitives(path);
    }
}

void MainWindow::on_actionSave_assignement_triggered()
{
    if (_project == NULL)
        return;

    QSettings settings;
    QString defaultPath = settings.value("Path/assPath").toString();

    QString path =
    QFileDialog::getSaveFileName(this,
                                 tr("Save assignements as CSV file"),
                                 defaultPath,
                                 "CSV (*.csv)");


    if (! path.isNull()){
        writeAssignement(path);
    }
}

void MainWindow::on_actionSave_all_triggered()
{
    if (_project == NULL)
        return;

    QSettings settings;
    QString defaultPath = settings.value("Path/allPath").toString();

    QString path =
            QFileDialog::getExistingDirectory(
                this,
                tr("Select output folder"),
                defaultPath);

    QStringList st;
    st << "primitives.csv"
       << "points_primitives.csv"
       << "cloud.ply";


    if (! path.isNull()){
        // check if files exists
        QDir dirPath(path);
        bool override = dirPath.entryList(st).size() != 0;
        if (override && QMessageBox::question(this,
                                              "Existing directory",
                                              "Existing files will be overwritten. Proceed anyway? ",
                                              QMessageBox::Save,
                                              QMessageBox::Cancel)
                == QMessageBox::Cancel)
            return;

        dirPath.mkdir("gt");

        QString projectName = dirPath.dirName();

        writeProject    (path+QString("/gt/")+projectName+QString(".prj"), projectName);
        writePrimitives (path+QString("/gt/primitives.csv"));
        writeAssignement(path+QString("/gt/points_primitives.csv"));
        writeSamples    (path+QString("/cloud.ply"));

        settings.setValue("Path/allPath", path);
    }
}

void MainWindow::writePrimitives(QString path){
    if (_project == NULL)
        return;

    using InputGen::Application::Scalar;

    QFile outfile(path);
    if (outfile.open(QIODevice::WriteOnly |
                   QIODevice::Truncate  |
                   QIODevice::Text)) {

        QSettings settings;
        settings.setValue("Path/priSave", path);

        QTextStream out(&outfile);

        out << "#Describes primitives of the scene" << endl;
        out << "#x,y,z,nx,ny,nz,primitiveId,orientationId,used" << endl;

        for(InputGen::Application::Project::PrimitiveContainer::const_iterator it = _project->primitives.begin();
            it != _project->primitives.end(); it++){
            const InputGen::Application::Primitive::vec& coord  = (*it).coord();
            const InputGen::Application::Primitive::vec& normal = (*it).normal();
            out /*<< (*it).coord()(0)     << ","
                << (*it).coord()(1)     << ","
                << (*it).coord()(2)     << ","*/
                << coord(0)             << ","
                << Scalar(1.)-coord(1)  << ","
                << coord(2)             << ","
                <<  normal(0)           << ","
                << -normal(1)           << ","
                <<  normal(2)           << ","
                << (*it).uid()          << ","
                << (*it).did()          << ","
                << "1"               << endl; //1 means used
        }
        outfile.close();
    }
}

void
MainWindow::writeProject(QString path, QString projectName){
    if (_project == NULL)
        return;

    QFile file( path );
    if( !file.open( QIODevice::WriteOnly ) )
        return;

    QDomDocument doc (projectName);

    QDomElement root = doc.createElement( "scene" );
    root.setAttribute( "name", projectName);
    doc.appendChild( root );

    // Write primitives
    {
        QDomElement rootPrimitives = doc.createElement( "primitives" );
        root.appendChild( rootPrimitives );


        std::for_each(_project->primitives.begin(),
                      _project->primitives.end(),
                      [&doc, &rootPrimitives] (typename Project::PrimitiveContainer::const_reference p)
        {
            QDomElement primElement = doc.createElement( "primitive" );
            primElement.setAttribute( "pos",
                                      QString::number(p.coord()(0)) + QString(" ") +
                                      QString::number(p.coord()(1)) + QString(" ") +
                                      QString::number(p.coord()(2)));
            primElement.setAttribute( "dir",
                                      QString::number(p.normal()(0)) + QString(" ") +
                                      QString::number(p.normal()(1)) + QString(" ") +
                                      QString::number(p.normal()(2)));
            primElement.setAttribute( "dim",
                                      QString::number(p.dim()(0)) + QString(" ") +
                                      QString::number(p.dim()(1)));
            primElement.setAttribute( "uid", QString::number(p.uid()) );
            primElement.setAttribute( "did", QString::number(p.did()) );

            rootPrimitives.appendChild( primElement );
        });
    }

    // Write samplers
    QDomElement rootSamplers = doc.createElement( "samplers" );
    root.appendChild( rootSamplers );
    _samplerDoc->saveSamples(doc, rootSamplers);

    // Write displacement
    QDomElement rootDisplacement = doc.createElement( "displacements" );
    root.appendChild( rootDisplacement );
    _displacementDoc->savekernels(doc, rootDisplacement);

    QTextStream ts( &file );
    ts << doc.toString();

    file.close();
}

void MainWindow::writeAssignement(QString path){
    if (_project == NULL)
        return;

    QFile outfile(path);
    if (outfile.open(QIODevice::WriteOnly |
                   QIODevice::Truncate  |
                   QIODevice::Text)) {

        QSettings settings;
        settings.setValue("Path/assPath", path);

        QTextStream out(&outfile);

        out << "#Describes point to primitive assignation" << endl;
        out << "#pointId,primitiveId,orientationId" << endl;

        unsigned int sampleId = 0;
        for(InputGen::Application::SampleSet::const_iterator it = _project->samples.begin();
            it != _project->samples.end(); it++, sampleId++){
            out << sampleId << "," << (*it).primitiveId << ",-1" << endl;
        }
        outfile.close();
    }
}

void MainWindow::writeSamples(QString path){

    if (_project == NULL)
        return;

    QFile outfile(path);
    if (outfile.open(QIODevice::WriteOnly |
                   QIODevice::Truncate  |
                   QIODevice::Text)) {

        QSettings settings;
        settings.setValue("Path/plySave", path);

        QTextStream out(&outfile);

        out << "ply\n"
            << "format ascii 1.0\n"
            << "comment Generated by InputGen\n"
            << "element vertex " << _project->samples.size() << "\n"
            << "property float x\n"
            << "property float y\n"
            << "property float z\n"
            << "property float nx\n"
            << "property float ny\n"
            << "property float nz\n"
            << "end_header\n";

        unsigned int sampleId = 0;
        for(InputGen::Application::SampleSet::const_iterator it = _project->samples.begin();
            it != _project->samples.end(); it++, sampleId++){
            InputGen::Application::Primitive::vec pos =
                    ((*it) + _project->computeTotalDisplacement(sampleId));
            out << pos(0) << " "
                << InputGen::Application::Scalar(1.) - pos(1) << " "     // flip Y
                << pos(2) << " 0 0 0" << endl;
        }

        outfile.close();
    }
}
