#include <osg/Depth>
#include <osg/Group>
#include <osg/Geode>
#include <osg/Geometry>
#include <osgText/Text>

#include "Types.h"
#include "ColorMap.h"
#include "Primitive.h"

#include "GlobFit.h"

static osg::Group* generateLabel(const std::string& title)
{
    osg::Camera* hudCamera = new osg::Camera;
    osgText::Text* labelText = new osgText::Text();
    labelText->setText(title);
    labelText->setPosition(osg::Vec3(-1.0f, 0.9f, 0.0f));
    labelText->setColor(osg::Vec4(0.2f,0.2f,0.6f,1.0f));
    labelText->setDataVariance(osg::Object::DYNAMIC);
    labelText->setCharacterSize(0.1f);

    std::string font("fonts/times.ttf");
    labelText->setFont(font);

    hudCamera->setReferenceFrame(osg::Transform::ABSOLUTE_RF_INHERIT_VIEWPOINT);
    hudCamera->setViewMatrix(osg::Matrix::identity());
    hudCamera->setRenderOrder(osg::Camera::POST_RENDER);
    hudCamera->setClearMask(GL_DEPTH_BUFFER_BIT);

    osg::Geode* geode = new osg::Geode();
    osg::StateSet* stateset = geode->getOrCreateStateSet();
    stateset->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
    stateset->setMode(GL_DEPTH_TEST,osg::StateAttribute::OFF);
    geode->addDrawable( labelText );
    hudCamera->addChild(geode);

    osg::Group* group = new osg::Group;
    group->addChild(hudCamera);

    return group;
}

static void generatePartition(const std::vector<size_t>& vecIndex,
    const std::vector<RichPoint*>& vecPointSet,
    osg::Geode* confidenceGeode,
    osg::Geode* partitionGeode,
    const osg::Vec4* partitionColor)
{
    if (!confidenceGeode && !partitionGeode) {
        return;
    }

    size_t partitionSize = 10000;
    size_t indexNum = vecIndex.size();

    std::vector<size_t> indexRange;
    indexRange.push_back(0);
    for (size_t i = partitionSize; i < indexNum; i += partitionSize) {
        indexRange.push_back(i);
    }
    if (indexNum%partitionSize != 0) {
        indexRange.push_back(indexNum);
    }

    for (size_t i = 0, iEnd = indexRange.size()-1; i < iEnd; ++ i) {
        size_t length = indexRange[i+1]-indexRange[i];

        osg::Vec3Array* vertices = new osg::Vec3Array;
        osg::Vec4Array* confidenceColors = new osg::Vec4Array;
        vertices->reserve(length);
        if (confidenceGeode) {
            confidenceColors->reserve(length);
        }

        for (size_t j = indexRange[i]; j < indexRange[i+1]; ++ j) {
            const RichPoint* richPoint = vecPointSet[vecIndex[j]];
            vertices->push_back(Vec3Caster<Point>(richPoint->point));
            if (confidenceGeode) {
                confidenceColors->push_back(ColorMap::Instance().getColor(ColorMap::JET, (1-richPoint->confidence), 0, 1.0));
            }
        }

        if (confidenceGeode) {
            osg::Geometry* geometry = new osg::Geometry;
            geometry->setUseDisplayList(true);
            geometry->setUseVertexBufferObjects(true);
            geometry->setVertexArray(vertices);
            geometry->setColorArray(confidenceColors);
            geometry->setColorBinding(osg::Geometry::BIND_PER_VERTEX);
            geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, length));
            confidenceGeode->addDrawable(geometry);
        }

        if (partitionGeode) {
            osg::Vec4Array* partitionColors = new osg::Vec4Array;
            partitionColors->push_back(*partitionColor);

            osg::Geometry* geometry = new osg::Geometry;
            geometry->setUseDisplayList(true);
            geometry->setUseVertexBufferObjects(true);
            geometry->setVertexArray(vertices);
            geometry->setColorArray(partitionColors);
            geometry->setColorBinding(osg::Geometry::BIND_OVERALL);
            geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, length));
            partitionGeode->addDrawable(geometry);
        }
    }

    return;
}

std::pair<osg::Node*, osg::Node*> GlobFit::convertPointsToGeometry() const
{
    osg::Geode* confidenceGeode = new osg::Geode;
    osg::Geode* partitionGeode = new osg::Geode;

    std::vector<bool> vecFlag(_vecPointSet.size(), true);
    for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
        const std::vector<size_t>& vecPointIdx = _vecPrimitive[i]->getPointIdx();
        osg::Vec4 partitionColor = ColorMap::Instance().getColor(ColorMap::JET, i, 0, iEnd);
        generatePartition(vecPointIdx, _vecPointSet, confidenceGeode, partitionGeode, &partitionColor);
        for (size_t j = 0, jEnd = vecPointIdx.size(); j < jEnd; ++ j) {
            vecFlag[vecPointIdx[j]] = false;
        }
    }

    std::vector<size_t> vecUnassociatedPointIdx;
    for (size_t i = 0, iEnd = _vecPointSet.size(); i < iEnd; ++ i) {
        if (vecFlag[i]) {
            vecUnassociatedPointIdx.push_back(i);
        }
    }
    generatePartition(vecUnassociatedPointIdx, _vecPointSet, confidenceGeode, NULL, NULL);

    osg::Group* confidenceGroup = new osg::Group;
    osg::Group* partitionGroup = new osg::Group;

    confidenceGroup->addChild(confidenceGeode);
    partitionGroup->addChild(partitionGeode);

    confidenceGroup->addChild(generateLabel("All Points"));
    partitionGroup->addChild(generateLabel("RANSAC Partitions"));

    confidenceGroup->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OVERRIDE|osg::StateAttribute::OFF);
    partitionGroup->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OVERRIDE|osg::StateAttribute::OFF);

    return std::make_pair(confidenceGroup, partitionGroup);
}

osg::Node* GlobFit::convertPrimitivesToGeometry(const std::string& title) const
{
    osg::Group* root = new osg::Group;

    for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
        const osg::Vec4& primitiveColor = ColorMap::Instance().getColor(ColorMap::JET, i, 0, iEnd);
        root->addChild(_vecPrimitive[i]->toGeometry(primitiveColor));
    }

    root->addChild(generateLabel(title));

    osg::StateSet* stateSet = root->getOrCreateStateSet();
    stateSet->setMode(GL_BLEND, osg::StateAttribute::ON);
    stateSet->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);
    stateSet->setMode(GL_DEPTH_TEST, osg::StateAttribute::ON);
    osg::Depth* depth = new osg::Depth;
    depth->setWriteMask(false);
    stateSet->setAttributeAndModes(depth, osg::StateAttribute::ON);
    stateSet->setMode(GL_LIGHTING, osg::StateAttribute::OVERRIDE|osg::StateAttribute::OFF);

    return root;
}