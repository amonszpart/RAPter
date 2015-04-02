#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osgViewer/CompositeViewer>
#include <osgViewer/ViewerEventHandlers>

#include "Viewer.h"

Viewer::Viewer()
{
}

Viewer::~Viewer(void)
{
}

void Viewer::operator()(std::vector<osg::Node*> vecViewData)
{
  if (vecViewData.size() == 0) {
    return;
  }

  osg::GraphicsContext::WindowingSystemInterface* wsi = osg::GraphicsContext::getWindowingSystemInterface();
  if (!wsi) {
    osg::notify(osg::NOTICE)<<"Error, no WindowSystemInterface available, cannot create windows."<<std::endl;
    return;
  }

  unsigned int width, height;
  wsi->getScreenResolution(osg::GraphicsContext::ScreenIdentifier(0), width, height);
  width -= 200;
  height -= 200;

  size_t dataNum = vecViewData.size();
  size_t row = (size_t)(std::ceil(std::sqrt((double)(dataNum))));
  size_t col = (size_t)(std::ceil((double)(dataNum)/row));

  width/= row;
  height/=col;

  osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
  traits->x = 100;
  traits->y = 100;
  traits->width = row*width;
  traits->height = col*height;
  traits->windowDecoration = true;
  traits->doubleBuffer = true;
  traits->sharedContext = 0;

  osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits.get());
  if (gc.valid()) {
    gc->setClearColor(osg::Vec4f(1.0f, 1.0f, 1.0f, 1.0f));
    gc->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  } else {
    osg::notify(osg::NOTICE)<<"  GraphicsWindow has not been created successfully."<<std::endl;
  }

  osg::ref_ptr<osgViewer::CompositeViewer> viewer = new osgViewer::CompositeViewer;
  osgGA::TrackballManipulator* tbManipulator = new osgGA::TrackballManipulator;

  for (size_t i = 0; i < row; ++ i) {
    for (size_t j = 0; j < col; ++ j) {
      size_t idx = j + i*col;
      if (idx >= dataNum) {
        break;
      }

      osgViewer::View* view = new osgViewer::View;
      viewer->addView(view);

      view->setSceneData(vecViewData[idx]);
      view->getCamera()->setGraphicsContext(gc.get());
      view->getCamera()->setClearColor(osg::Vec4f(1.0f, 1.0f, 1.0f, 1.0f));
      view->getCamera()->setViewport(new osg::Viewport(i*width, (col-1-j)*height, width, height));
      view->getCamera()->setProjectionMatrixAsPerspective(30.0, double(width)/height, 1.0, 1000.0);

      view->setCameraManipulator(tbManipulator);
      view->addEventHandler(new osgGA::StateSetManipulator(vecViewData[idx]->getOrCreateStateSet()));
    }
  }

  viewer->run();
  return;
}
