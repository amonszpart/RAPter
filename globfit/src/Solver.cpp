#include <boost/filesystem.hpp>

#include <engine.h>
#pragma comment( lib, "libmx.lib" )
#pragma comment( lib, "libeng.lib" )

#include "Types.h"
#include "Primitive.h"
#include "RelationEdge.h"

#include "GlobFit.h"

#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds

namespace {
  Engine* matlabEngine = NULL;
  mxArray* numVertices = NULL;
  mxArray* maxIterNum = NULL;
  mxArray* primitiveType = NULL;
  mxArray* coordX = NULL;
  mxArray* coordY = NULL;
  mxArray* coordZ = NULL;
  mxArray* normalX = NULL;
  mxArray* normalY = NULL;
  mxArray* normalZ = NULL;
  mxArray* confVertices = NULL;
}

    // by Nicolas
std::string getExecPath()
{
    std::string fullFileName = "";

    // Code taken from: http://www.gamedev.net/community/forums/topic.asp?topic_id=459511
    std::string path = "";
    pid_t pid = getpid();
    char buf[20] = {0};
    sprintf(buf,"%d",pid);
    std::string _link = "/proc/";
    _link.append( buf );
    _link.append( "/exe");
    char proc[512];
    int ch = readlink(_link.c_str(),proc,512);
    if (ch != -1) {
        proc[ch] = 0;
        path = proc;
        std::string::size_type t = path.find_last_of("/");
        path = path.substr(0,t);
    }

    fullFileName = path + std::string("/");
    return fullFileName;
}

bool GlobFit::createMatlabArraies()
{
    //matlabEngine = engOpen("\0");
    matlabEngine = engOpen("matlab -nodesktop");
    if (NULL == matlabEngine) {
        fprintf(stderr, "Could not initialize the engine.\n");
        return false;
    }

    // by Nicolas
    {
        std::string pathCommand( "addpath('");
        pathCommand.append(getExecPath());
        pathCommand.append("../../matlab')");
        std::cout << "[" << __func__ << "]: " << "calling " << pathCommand << std::endl;
        engEvalString(matlabEngine, pathCommand.c_str());

    }

    size_t numPrimitives = _vecPrimitive.size();
    for ( size_t i = 0; i < numPrimitives; ++i )
    {
        const Primitive* pPrimitive = _vecPrimitive[i];
        if (pPrimitive->getType() != Primitive::PT_CONE)
        {
            //std::cout << "[" << __func__ << "]: " << "skipping non-cone" << std::endl;
            continue;
        }

        Vector normal;
        pPrimitive->getNormal( normal );
        const std::vector<size_t>& vecVertexIdx = pPrimitive->getPointIdx();
        for (size_t j = 0, jEnd = vecVertexIdx.size(); j < jEnd; ++j) {
            size_t idx = j*numPrimitives+i;
            RichPoint* pRichPoint = _vecPointSet[vecVertexIdx[j]];
            if (pRichPoint->normal*normal < 0) {
                pRichPoint->normal = -pRichPoint->normal;
            }
        }
    }
    if ( !_vecPrimitive.size() )
        std::cout << "[" << __func__ << "]: " << "no primitives..." << std::endl;

    numVertices = mxCreateNumericMatrix(numPrimitives, 1, mxINT32_CLASS, mxREAL);
    int *pNumVertices = (int*)mxGetData(numVertices);

    maxIterNum = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    int* pMaxIterNum = (int*)mxGetData(maxIterNum);
    *pMaxIterNum = 500;

    primitiveType = mxCreateNumericMatrix(numPrimitives, 1, mxINT32_CLASS, mxREAL);
    int *pPrimitiveType = (int*)mxGetData(primitiveType);

    int maxNumVertices = 0;
    for (size_t i = 0; i < numPrimitives; ++i) {
        const Primitive* pPrimitive = _vecPrimitive[i];
        pNumVertices[i] = std::max(1, (int)(pPrimitive->getPointIdx().size()));
        pPrimitiveType[i] = pPrimitive->getType();
        if (maxNumVertices < pNumVertices[i])
            maxNumVertices = pNumVertices[i];
    }

    coordX = mxCreateDoubleMatrix(numPrimitives, maxNumVertices, mxREAL);
    double* pCoordX = mxGetPr(coordX);

    coordY = mxCreateDoubleMatrix(numPrimitives, maxNumVertices, mxREAL);
    double* pCoordY = mxGetPr(coordY);

    coordZ = mxCreateDoubleMatrix(numPrimitives, maxNumVertices, mxREAL);
    double* pCoordZ = mxGetPr(coordZ);

    normalX = mxCreateDoubleMatrix(numPrimitives, maxNumVertices, mxREAL);
    double* pNormalX = mxGetPr(normalX);

    normalY = mxCreateDoubleMatrix(numPrimitives, maxNumVertices, mxREAL);
    double* pNormalY = mxGetPr(normalY);

    normalZ = mxCreateDoubleMatrix(numPrimitives, maxNumVertices, mxREAL);
    double* pNormalZ = mxGetPr(normalZ);

    confVertices = mxCreateDoubleMatrix(numPrimitives, maxNumVertices, mxREAL);
    double* pConfVertices = mxGetPr(confVertices);

    memset(pCoordX, 0, maxNumVertices*sizeof(double));
    memset(pCoordY, 0, maxNumVertices*sizeof(double));
    memset(pCoordZ, 0, maxNumVertices*sizeof(double));
    memset(pNormalX, 0, maxNumVertices*sizeof(double));
    memset(pNormalY, 0, maxNumVertices*sizeof(double));
    memset(pNormalZ, 0, maxNumVertices*sizeof(double));
    memset(pConfVertices, 0, maxNumVertices*sizeof(double));
    for (size_t i = 0; i < numPrimitives; ++i) {
        const Primitive* pPrimitive = _vecPrimitive[i];
        const std::vector<size_t>& vecVertexIdx = pPrimitive->getPointIdx();
        for (size_t j = 0, jEnd = vecVertexIdx.size(); j < jEnd; ++j) {
            size_t idx = j*numPrimitives+i;
            const RichPoint* pRichPoint = _vecPointSet[vecVertexIdx[j]];
            pCoordX[idx] = pRichPoint->point.x();
            pCoordY[idx] = pRichPoint->point.y();
            pCoordZ[idx] = pRichPoint->point.z();
            pNormalX[idx] = pRichPoint->normal.x();
            pNormalY[idx] = pRichPoint->normal.y();
            pNormalZ[idx] = pRichPoint->normal.z();
            pConfVertices[idx] = pRichPoint->confidence;
        }
    }

    engPutVariable(matlabEngine, "numVertices", numVertices);
    engPutVariable(matlabEngine, "primitiveType", primitiveType);
    engPutVariable(matlabEngine, "coordX", coordX);
    engPutVariable(matlabEngine, "coordY", coordY);
    engPutVariable(matlabEngine, "coordZ", coordZ);
    engPutVariable(matlabEngine, "normalX", normalX);
    engPutVariable(matlabEngine, "normalY", normalY);
    engPutVariable(matlabEngine, "normalZ", normalZ);
    engPutVariable(matlabEngine, "confVertices", confVertices);
    engPutVariable(matlabEngine, "maxIterNum", maxIterNum);

    return true;
}

void GlobFit::destoryMatlabArraies()
{
    mxDestroyArray(numVertices);
    mxDestroyArray(primitiveType);
    mxDestroyArray(coordX);
    mxDestroyArray(coordY);
    mxDestroyArray(coordZ);
    mxDestroyArray(normalX);
    mxDestroyArray(normalY);
    mxDestroyArray(normalZ);
    mxDestroyArray(confVertices);
    mxDestroyArray(maxIterNum);
}

void GlobFit::dumpData(const std::vector<RelationEdge>& vecRelationEdge, const std::string& stageName)
{
    size_t maxVerticesNum = 0;
    for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
        Primitive* pPrimitive = _vecPrimitive[i];
        pPrimitive->prepareParameters();
        maxVerticesNum = std::max(pPrimitive->getPointIdx().size(), maxVerticesNum);
    }

    std::locale loc("C");

    std::string path = boost::filesystem::current_path().string();
    path = path+"/matlab/data/"+stageName+"/";
    boost::filesystem::create_directories(path);
    std::ofstream constraints((path+"constraints.dat").c_str());
    std::ofstream primitiveType((path+"primitiveType.dat").c_str());
    std::ofstream inputParameters((path+"inputParameters.dat").c_str());
    std::ofstream numVertices((path+"numVertices.dat").c_str());
    std::ofstream coordX((path+"coordX.dat").c_str());
    std::ofstream coordY((path+"coordY.dat").c_str());
    std::ofstream coordZ((path+"coordZ.dat").c_str());
    std::ofstream normalX((path+"normalX.dat").c_str());
    std::ofstream normalY((path+"normalY.dat").c_str());
    std::ofstream normalZ((path+"normalZ.dat").c_str());
    std::ofstream confVertices((path+"confVertices.dat").c_str());

    constraints.imbue(loc);
    primitiveType.imbue(loc);
    inputParameters.imbue(loc);
    inputParameters.precision(16);
    numVertices.imbue(loc);
    coordX.imbue(loc);
    coordX.precision(16);
    coordY.imbue(loc);
    coordY.precision(16);
    coordZ.imbue(loc);
    coordZ.precision(16);
    normalX.imbue(loc);
    normalX.precision(16);
    normalY.imbue(loc);
    normalY.precision(16);
    normalZ.imbue(loc);
    normalZ.precision(16);
    confVertices.imbue(loc);
    confVertices.precision(16);

    for (size_t i = 0, iEnd = vecRelationEdge.size(); i < iEnd; ++ i) {
        const RelationEdge& relationEdge = vecRelationEdge[i];
        constraints << relationEdge << std::endl;
    }

    for (size_t i = 0, iEnd = _vecPrimitive.size(); i < iEnd; ++ i) {
        Primitive* pPrimitive = _vecPrimitive[i];

        primitiveType << pPrimitive->getType() << std::endl;

        pPrimitive->prepareParameters();
        for (size_t j = 0, jEnd = Primitive::getNumParameter(); j < jEnd; ++ j) {
            inputParameters << pPrimitive->getParameter(j) << " ";
        }
        inputParameters << std::endl;

        numVertices << pPrimitive->getPointIdx().size() << std::endl;
        for (size_t j = 0, jEnd = pPrimitive->getPointIdx().size(); j < jEnd; ++ j) {
            const RichPoint* richPoint = _vecPointSet[pPrimitive->getPointIdx()[j]];
            const Point& point = richPoint->point;
            const Vector& normal = richPoint->normal;
            coordX << point.x() << " ";
            coordY << point.y() << " ";
            coordZ << point.z() << " ";
            normalX << normal.x() << " ";
            normalY << normal.y() << " ";
            normalZ << normal.z() << " ";
            confVertices << richPoint->confidence << " ";
        }

        for (size_t j = pPrimitive->getPointIdx().size(); j < maxVerticesNum; ++ j) {
            coordX << 0 << " ";
            coordY << 0 << " ";
            coordZ << 0 << " ";
            normalX << 0 << " ";
            normalY << 0 << " ";
            normalZ << 0 << " ";
            confVertices << 0 << " ";
        }
        coordX << std::endl;
        coordY << std::endl;
        coordZ << std::endl;
        normalX << std::endl;
        normalY << std::endl;
        normalZ << std::endl;
        confVertices << std::endl;
    }

    return;
}

bool GlobFit::solve(std::vector<RelationEdge>& vecRelationEdge, RelationEdge::RelationEdgeType currentStage, const std::string& stageName)
{


    // dump data to file for debugging in matlab
    std::cout << "[" << __func__ << "]: " << "wrote to " << stageName << std::endl;
    dumpData(vecRelationEdge, stageName);

    size_t nConstraintNum = vecRelationEdge.size();
    std::string optimization;
    if (currentStage < RelationEdge::RET_COAXIAL) {
        optimization = "OptimizeNormal";
        std::cout << "Optimize Normal..." << std::endl;
    } else if (currentStage < RelationEdge::RET_COPLANAR) {
        optimization = "OptimizePoint";
        std::cout << "Optimize Point..." << std::endl;
    } else if (currentStage < RelationEdge::RET_EQUAL_RADIUS) {
        optimization = "OptimizeDistance";
        std::cout << "Optimize Distance..." << std::endl;
    } else {
        optimization = "OptimizeRadius";
        std::cout << "Optimize Radius..." << std::endl;
    }

    if (nConstraintNum == 0)
    {
        std::cout << "Empty constraint set." << std::endl;
        return true;
    }

    size_t numPrimitives = _vecPrimitive.size();
    mxArray* inputParameters = mxCreateDoubleMatrix(numPrimitives, Primitive::getNumParameter(), mxREAL);
    double* pInputParameters = mxGetPr(inputParameters);
    for (size_t i = 0; i < numPrimitives; ++i) {
        Primitive* pPrimitive = _vecPrimitive[i];
        pPrimitive->prepareParameters();
        for (size_t j = 0; j < Primitive::getNumParameter(); ++ j) {
            pInputParameters[j*numPrimitives+i] = pPrimitive->getParameter(j);
        }
    }
    engPutVariable(matlabEngine, "inputParameters", inputParameters);

    mxArray* constraints = mxCreateNumericMatrix(nConstraintNum, RelationEdge::getNumParameter(), mxINT32_CLASS, mxREAL);
    int* pConstraints = (int*)mxGetData(constraints);
    for (size_t i = 0; i < nConstraintNum; ++i) {
        vecRelationEdge[i].dumpData(pConstraints, nConstraintNum, i);
    }
    engPutVariable(matlabEngine, "constraints", constraints);

    std::string path = boost::filesystem::current_path().string();
    path = "cd "+path+"/matlab;";
    engEvalString(matlabEngine, path.c_str());
    size_t szOutputBuffer = 65536;
    
    char* matlabOutputBuffer = new char[szOutputBuffer];
    engOutputBuffer(matlabEngine, matlabOutputBuffer, szOutputBuffer);

    std::string output = "[outputParameters, initialFittingError, exitFittingError, exitFlag]";
    std::string input = "(inputParameters, maxIterNum, numVertices, primitiveType, coordX, coordY, coordZ, normalX, normalY, normalZ, confVertices, constraints);";

    std::string command = output+"="+optimization+input;
    std::cout << "[" << __func__ << "]: " << " calling matlab: " << command << std::endl;
    
       
    // by Nicolas
    bool running = true;
    // from http://www.bogotobogo.com/cplusplus/C11/3_C11_Threading_Lambda_Functions.php
    std::thread printOutBuffer([&running, &matlabOutputBuffer](){
        
        std::ofstream logfile;
        while(running){
            logfile.open ("run_globOpt.log", std::ofstream::trunc);
            logfile << matlabOutputBuffer << std::endl; 
            logfile.close ();
            std::this_thread::sleep_for (std::chrono::seconds(1));            
        }
        }); 

    
    engEvalString(matlabEngine, command.c_str());
    
    running = false;
    printOutBuffer.join();
    
    std::cout << "Optimisation done" << std::endl;
    

    matlabOutputBuffer[szOutputBuffer - 1] = '\0';
    printf("%s\n", matlabOutputBuffer);
    engOutputBuffer(matlabEngine, NULL, 0);
    delete[] matlabOutputBuffer;

    mxArray* outputParameters = engGetVariable(matlabEngine, "outputParameters");
    double *pOutputParameters = mxGetPr(outputParameters);
    mxArray* initialFittingError = engGetVariable(matlabEngine, "initialFittingError");
    double *pInitialFittingError = mxGetPr(initialFittingError);
    mxArray* exitFittingError = engGetVariable(matlabEngine, "exitFittingError");
    double *pExitFittingError = mxGetPr(exitFittingError);
    mxArray* exitFlag = engGetVariable(matlabEngine, "exitFlag");
    double *pExitFlag = mxGetPr(exitFlag);

    bool bValidOptimization = (*pExitFlag >= 0);
    // posterior check: consider invalid if fitting error increased too much
    // however, if the threshold is very big, the fitting error may increase a lot
    // so, be careful with this
    bValidOptimization &= (*pExitFittingError < 10*(*pInitialFittingError));
    if (!bValidOptimization) {
        mxDestroyArray(constraints);
        mxDestroyArray(inputParameters);
        mxDestroyArray(outputParameters);
        mxDestroyArray(initialFittingError);
        mxDestroyArray(exitFittingError);
        mxDestroyArray(exitFlag);

        std::cout << "No feasible solution found." << std::endl;
        return false;
    }

    // update primitives
    for (size_t i = 0; i < numPrimitives; ++ i) {
        Primitive* pPrimitive = _vecPrimitive[i];
        for (size_t j = 0; j < Primitive::getNumParameter(); ++ j) {
            pPrimitive->setParameter(j, pOutputParameters[j*numPrimitives+i]);
        }
        pPrimitive->applyParameters();
    }

    // destroy matrix
    mxDestroyArray(constraints);
    mxDestroyArray(inputParameters);
    mxDestroyArray(outputParameters);
    mxDestroyArray(initialFittingError);
    mxDestroyArray(exitFittingError);
    mxDestroyArray(exitFlag);

    return true;
}
