#include <iostream>
#include "CLogger/CLogger.h"
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include "CSensorManager.h"
#include "CObjectCatalogManager.h"
#include "CLabelManager.h"
#include <filesystem>
#include <boost/histogram.hpp>
#include <pcl/visualization/pcl_plotter.h>
#include <pcl/visualization/cloud_viewer.h>
#include "pcd_viewer.h"
#include "CAMImage.h"
#include "pcl/io/ply_io.h"
#include "pcl/io/auto_io.h"
#include "PCCommonFun.h"
#include "CActivePanel.h"
#include "RplyLib\CPlyIO.h"
#include "pcl/sample_consensus/ransac.h"
#include "pcl/sample_consensus/impl/ransac.hpp"
#include "pcl/sample_consensus/sac_model_plane.h"
#include "pcl/sample_consensus/impl/sac_model_plane.hpp"
#include "pcl/filters/extract_indices.h"
#include "pcl/sample_consensus/sac_model_perpendicular_plane.h"
#include "pcl/sample_consensus/impl/sac_model_perpendicular_plane.hpp"
//#include "pcl/segmentation/conditional_euclidean_clustering.h"
//#include "pcl/segmentation/impl/conditional_euclidean_clustering.hpp"
#include "pcl/segmentation/extract_clusters.h"
#include "pcl/segmentation/impl/extract_clusters.hpp"
#include "CRSTManager.h"
#include "CCmdOptions.h"
#include "E57CloudReader.h"
#include <complex>
#include <cmath>
#include <vector>
#include <iostream>
#include "pocketfft_hdronly.h"


using namespace pcl;
using namespace std;
using namespace pocketfft;

//#define PointType PointXYZRGBNormal
//#define PointType PointXYZRGB

typedef struct
{
    ON_Plane oPlane;
    PointUV vUVMin;
    PointUV vUVMax;
} SSurface;

typedef struct
{
    ON_Plane oFragmentPlane;
    vector<SSurface> aCandidate;
} SFragmentFamily;


#define pp
#if defined ppp

SCandidatePlaneInfo RANSACPlane(PointCloud<PointType>::Ptr& pOrigCloud, ON_Plane & oViewPlane, bool bForceOrientation, int iPlaneOrientation, SObjCatalogInfo& sObjInfo)
{
    SCandidatePlaneInfo sPlaneInfo;

    // Make a copy of original cloud
    pcl::PointCloud<PointType>::Ptr pCloud(new pcl::PointCloud<PointType>);
    pcl::copyPointCloud<PointType>(*pOrigCloud, *pCloud);

    // Init local variables
    vector< ON_Plane> aCandidatePlane;
    vector< PointCloud<PointType>::Ptr> aCandidateClouds;
    ON_3dVector vVertical(0, 0, 1);
    float fMinDist = 1e+30;
    ON_3dVector vCameraAt = oViewPlane.Normal();  // Store initial camera at direction
    PointUV vUVMin, vUVMax;

    Eigen::Vector3f vPlaneConstrain;
    float fAngularTol;
    if (iPlaneOrientation == 0)
    {
        vPlaneConstrain.x() = 0;
        vPlaneConstrain.y() = 0;
        vPlaneConstrain.z() = 1;
        fAngularTol = 10;
    }
    else
    {
        vPlaneConstrain.x() = oViewPlane.zaxis.x;
        vPlaneConstrain.y() = oViewPlane.zaxis.y;
        vPlaneConstrain.z() = oViewPlane.zaxis.z;
        fAngularTol = 60;
    }

    int nPoint = pCloud->size();
    bool bFound = 0;

    float fEccentricityH = sObjInfo.fHSizeMax / sObjInfo.fHSizeMin;

    for (int i = 0; i < 5; i++)
    {
        // Create "PerpendicularPlane" model
        SampleConsensusModelPerpendicularPlane<PointType>::Ptr model_p(new SampleConsensusModelPerpendicularPlane<PointType>(pCloud));
        model_p->setAxis(vPlaneConstrain);
        model_p->setEpsAngle(fAngularTol *3.1415926/180);  // angle tol=30°

        // Created RandomSampleConsensus object for that model
        RandomSampleConsensus<PointType> ransac(model_p);
        ransac.setDistanceThreshold(.01);
        ransac.setMaxIterations(5000);

        // Compute 
        bool bResult =ransac.computeModel();
        if (!bResult && !bFound)  // Computation failed and no previous result available
        {
            sPlaneInfo.iError = 1;
            return sPlaneInfo;
        }

        // Get inliers
        Indices inliers;
        ransac.getInliers(inliers);

        if (inliers.size() > nPoint / 100)  // Check if num inliers > 1/100 of all cloud point
        {
            // Get plane coeff
            Eigen::VectorXf model_coefficients;
            ransac.getModelCoefficients(model_coefficients);

            // Build candidate plane
            ON_PlaneEquation oEquation(model_coefficients[0], model_coefficients[1], model_coefficients[2], model_coefficients[3]);
            ON_Plane oPlane(oEquation);


            float fOrientCheck = fabs(oPlane.Normal() * vVertical);
            float fTolH = cos(5 * 3.1415926 / 180.);  // 5 deegree tollerance
            float fTolV = sin(5 * 3.1415926 / 180.);  // 5 deegree tollerance
            if ((iPlaneOrientation == 0 && fOrientCheck > fTolH) || (iPlaneOrientation == 1 && fOrientCheck < fTolV))  // Good candidate
            {
                PointCloud<PointType>::Ptr pTmp(new PointCloud<PointType>);
                copyPointCloud(*pCloud, inliers, *pTmp);
                
                if (1)  // Esco al primo buono
                {
                    sPlaneInfo.sPlane = oPlane;
                    sPlaneInfo.pInlierCloud = pTmp;
                    bFound = 1;

                    break; // ESCO AL PRIMO BUONO
                }
                else  // Cerco il piano più vicino  NOTA: questo ramo potrebbe essere abilitato in automatico se il bbox (orientato come la camera) è molto lungo rispetto alla larghezza
                {
                    float fTmp = CloudPlaneMinDistance(oViewPlane, pCloud);
                    if (fTmp < fMinDist)
                    {
                        fMinDist = fTmp;
                        sPlaneInfo.sPlane = oPlane;
                        sPlaneInfo.pInlierCloud = pTmp;
                        bFound = 1;

                        //string sName = "D:\\Andrea\\Dev\\C++\\PointCloud\\Package\\Project2\\SegmentedPC\\Ransac";
                        //sName += to_string(i) + string(".pcd");
                        //io::save(sName.c_str(), *pTmp);
                    }
                }
            }
        }

        // Calcolo elenco indici rimanenti 
        if (inliers.size() > 0)
        {
            Indices iAll;
            Indices iDiff;
            iAll.resize(pCloud->size());
            for (int k = 0; k < pCloud->size(); k++)
                iAll[k] = k;   // iAll == elenco completo di indici

            set_difference(iAll.begin(), iAll.end(), inliers.begin(), inliers.end(), inserter(iDiff, iDiff.begin()));  // Rimuovo gli inliers

            if (iDiff.size() < nPoint / 20)  // Remaining point fewer than 5% of initial cloud
                break;

            // Creo nuvola temporanea con i soli ounti restanti
            PointCloud<PointType>::Ptr pCloudTmp(new pcl::PointCloud<PointType>);
            copyPointCloud(*pCloud, iDiff, *pCloudTmp);

            // Eseguo swap dei punti
            pCloud->points.swap(pCloudTmp->points);
        }
        else  // No more valid point can be found
            break;

    }

    static int iCount = 0;
    if (bFound)
    {
    //+++++++++++++++++++++++
     //   string sName = "D:\\Andrea\\Dev\\C++\\PointCloud\\Package\\Project2\\SegmentedPC\\Ransac";
     //   sName += to_string(iCount) + string(".pcd");
     //   io::save(sName.c_str(), *(sPlaneInfo.pInlierCloud));
     //   iCount++;
    //+++++++++++++++++++++++

        
        if (sPlaneInfo.sPlane.DistanceTo(oViewPlane.Origin()) < 0) // Normale di oPlane è sbagliata --> la inverto
        {
            ON_3dPoint vTmpOrigin = sPlaneInfo.sPlane.Origin();
            ON_3dVector vTmpNormal = -sPlaneInfo.sPlane.Normal();
            sPlaneInfo.sPlane.CreateFromNormal(vTmpOrigin, vTmpNormal);
        }
        // Rotate plane
        if (iPlaneOrientation == 1)  // Piano verticale
        {
            LocalY2GlobalZ(sPlaneInfo.sPlane); // Ruoto il TargetPlane in modo da portare l'asse Y il più possibile allineato allo Z globale
            if (!bForceOrientation)
                MinimizeBBox(sPlaneInfo.pInlierCloud, sPlaneInfo.sPlane, vUVMin, vUVMax, 10, 1); // Provo a  ruotare sino a 10°, con step iniziale di 1
        }
        else  // Piano orizzontale
            MinimizeBBox(sPlaneInfo.pInlierCloud, sPlaneInfo.sPlane, vUVMin, vUVMax, 45, 3);  // Provo a  ruotare sino a 45°, con step iniziale di 3


        // Build Result struct to be returned
        sPlaneInfo.fConvergenceMeanErr = 0;
        sPlaneInfo.nConvergenceIter = 1;
        sPlaneInfo.vUVMin = vUVMin;
        sPlaneInfo.vUVMax = vUVMax;
        sPlaneInfo.fOrtoFactor = -(sPlaneInfo.sPlane.Normal() * vCameraAt);
        sPlaneInfo.iError = 0;

        ON_3dPoint P0 = ON_3dPoint(0, 0, 0);
        ON_3dVector X0 = ON_3dVector(1, 0, 0);
        ON_3dVector Y0 = ON_3dVector(0, 1, 0);
        ON_3dVector Z0 = ON_3dVector(0, 0, 1);
        ON_3dPoint P1 = sPlaneInfo.sPlane.Origin() + vUVMin.u * sPlaneInfo.sPlane.xaxis + vUVMin.v * sPlaneInfo.sPlane.yaxis;

        // Calcolo rotazione che porta da sistema riferimento globale a sistema locale del piano (com origine in Umin,Vmin)
        sPlaneInfo.xform.Rotation(P0, X0, Y0, Z0, P1, sPlaneInfo.sPlane.zaxis, sPlaneInfo.sPlane.xaxis, sPlaneInfo.sPlane.yaxis);
    }
    else
        sPlaneInfo.iError = 1;

    return sPlaneInfo;


}

void SearchAndRemoveFarClusters(PointCloud<PointType>::Ptr& pCloud, SCameraInfo& sCameraInfo)
{
    // Creating the KdTree object for the search method of the extraction
    pcl::search::KdTree<pcl::PointType>::Ptr tree(new pcl::search::KdTree<pcl::PointType>);
    tree->setInputCloud(pCloud);

    std::vector<pcl::PointIndices> cluster_indices;

    // Identify clusters
    int nClusterMin = pCloud->size() / 20; // 5% of original cloud
    int nClusterMax = pCloud->size();       

    pcl::EuclideanClusterExtraction<pcl::PointType> ec;
    ec.setClusterTolerance(0.01); // 50cm
    ec.setMinClusterSize(nClusterMin);
    ec.setMaxClusterSize(nClusterMax);
    ec.setSearchMethod(tree);
    ec.setInputCloud(pCloud);
    ec.extract(cluster_indices);

    if (cluster_indices.size() <= 1)
        return;


    // Compute min dist
    ON_3dPoint vCameraPos(sCameraInfo.vPos.x(), sCameraInfo.vPos.y(), sCameraInfo.vPos.z());

    double fMinDist = 1e+30;
    int iClust = -1;
    for (int i = 0; i < cluster_indices.size(); i++)  //Loop over culsters
    {
        ON_3dPoint vClusterPos( 0,0,0);
        for (const auto& idx : cluster_indices[i].indices)  // Loop over indices of i-th cluster
            vClusterPos += ON_3dPoint(pCloud->points[idx].x, pCloud->points[idx].y, pCloud->points[idx].z);

        vClusterPos /= cluster_indices[i].indices.size();
        double fTmpDist=vClusterPos.DistanceTo(vCameraPos);
        if (fTmpDist < fMinDist)  // If current cluster is the closest --> store its index
        {
            fMinDist = fTmpDist;
            iClust = i;
        }

    }

    // Creo nuvola temporanea con i soli ounti restanti
    PointCloud<PointType>::Ptr pCloudTmp(new pcl::PointCloud<PointType>);
    copyPointCloud(*pCloud, cluster_indices[iClust].indices, *pCloudTmp);

    // Eseguo swap dei punti
    pCloud->points.swap(pCloudTmp->points);

}

SCandidatePlaneInfo DetectPlane(PointCloud<PointType>::Ptr& pCloud, SCameraInfo& sCameraInfo, bool bForceOrientation, int iViewFace, int iMaxIteration, SObjCatalogInfo& sObjInfo)
{
    // Recupero orientamento camera
    Eigen::Matrix3d R = sCameraInfo.mRot.inverse();

    CActivePanel oActivePanel;
    int iPlaneOrientation;
    if (iViewFace == ESideType::eTop || iViewFace == ESideType::eBottom)
        iPlaneOrientation = 0;  // Horizzontal
    else
        iPlaneOrientation = 1;  // Vertical

    // Metodo standard

    // Costruisco piano camera
    ON_3dPoint vOrigin(sCameraInfo.vPos[0], sCameraInfo.vPos[1], sCameraInfo.vPos[2]);
    ON_3dVector vNormal(R(2, 0), R(2, 1), R(2, 2));
    ON_Plane oPlane(vOrigin, vNormal);

    // Detect plane using RANSAC method
    SCandidatePlaneInfo sInfo = RANSACPlane(pCloud, oPlane, bForceOrientation, iPlaneOrientation, sObjInfo);

    if(sInfo.iError==0) 
        return sInfo;

    // If code is here RANSAC has failed --> try with ActivePanel
    cout << endl << "RANSAC plane estimation failed!!!!!!   I'll try with AdactivePanel method" << endl;

    // Detect plane using ActivePanel method
    sInfo = oActivePanel.SimpleTargetPlaneEstimation(oPlane, pCloud, 3, 3, bForceOrientation, iPlaneOrientation, iMaxIteration);

    if (sInfo.nConvergenceIter < iMaxIteration && sInfo.fOrtoFactor > 0.707)  // Ottenuta convergenza 
        return sInfo;

    // **************** Provo metodo avvicinamento progressivo ***************
    cout << endl<<"Standard plane estimation failed!!!!!!   I'll try with progressive method" << endl;

    // Ricostruisco piano camera iniziale 
    vOrigin.Set(sCameraInfo.vPos[0], sCameraInfo.vPos[1], sCameraInfo.vPos[2]);
    vNormal.Set(R(2, 0), R(2, 1), R(2, 2));
    oPlane.CreateFromNormal(vOrigin, vNormal);

    sInfo = oActivePanel.ProgressiveTargetPlaneEstimation(oPlane, pCloud, 3, 3, bForceOrientation, iPlaneOrientation, iMaxIteration);

    if (sInfo.nConvergenceIter < iMaxIteration && sInfo.fOrtoFactor > 0.707)  // Ottenuta convergenza 
        return sInfo;


    return sInfo;
}

void ComputeTransform(SObjCatalogInfo& sObjInfo, SLabelInfo& sLabelInfo, SCandidatePlaneInfo & sPlaneInfo, ON_Xform & mXform, ON_Xform & mScale, ON_3dPoint vBBoxSide)
{
    // Calcolo trasformazione che porta:
    
    // Calcolo traslazione
    ON_3dPoint vFaceCenter = sPlaneInfo.sPlane.origin+ 0.5 * (sPlaneInfo.vUVMax.u + sPlaneInfo.vUVMin.u) * sPlaneInfo.sPlane.xaxis + 0.5 * (sPlaneInfo.vUVMax.v + sPlaneInfo.vUVMin.v) * sPlaneInfo.sPlane.yaxis;
    if (sObjInfo.iFrontAxis == 0 && sObjInfo.iVerticalAxis == 2) 
    {
        if (sLabelInfo.iViewFace == ESideType::eFront) // X ---> oPlane.Z, Y ---> oPlane.X, Z ---> oPlane.Y
        {
            // Calcolo rotazione 
            mXform[0][0] = sPlaneInfo.sPlane.zaxis.x;
            mXform[1][0] = sPlaneInfo.sPlane.zaxis.y;
            mXform[2][0] = sPlaneInfo.sPlane.zaxis.z;
            mXform[0][1] = sPlaneInfo.sPlane.xaxis.x;
            mXform[1][1] = sPlaneInfo.sPlane.xaxis.y;
            mXform[2][1] = sPlaneInfo.sPlane.xaxis.z;
            mXform[0][2] = sPlaneInfo.sPlane.yaxis.x;
            mXform[1][2] = sPlaneInfo.sPlane.yaxis.y;
            mXform[2][2] = sPlaneInfo.sPlane.yaxis.z;

            // Aggiungo correzione distanza faccia da centro oggetto
            vFaceCenter = vFaceCenter - sPlaneInfo.sPlane.zaxis * mScale[0][0] * vBBoxSide.x/2;

            // Aggiungo offset
            vFaceCenter = vFaceCenter - sPlaneInfo.sPlane.xaxis * sObjInfo.fOffsetY * mScale[1][1] - sPlaneInfo.sPlane.yaxis * sObjInfo.fOffsetZ * mScale[2][2] - sPlaneInfo.sPlane.zaxis * sObjInfo.fOffsetX * mScale[0][0];
        }
        else if (sLabelInfo.iViewFace == ESideType::eTop) // Z ---> oPlane.Z, X ---> oPlane.X, Y ---> oPlane.Y
        {
            // Calcolo rotazione 
            mXform[0][0] = sPlaneInfo.sPlane.xaxis.x;
            mXform[1][0] = sPlaneInfo.sPlane.xaxis.y;
            mXform[2][0] = sPlaneInfo.sPlane.xaxis.z;
            mXform[0][1] = sPlaneInfo.sPlane.yaxis.x;
            mXform[1][1] = sPlaneInfo.sPlane.yaxis.y;
            mXform[2][1] = sPlaneInfo.sPlane.yaxis.z;
            mXform[0][2] = sPlaneInfo.sPlane.zaxis.x;
            mXform[1][2] = sPlaneInfo.sPlane.zaxis.y;
            mXform[2][2] = sPlaneInfo.sPlane.zaxis.z;

            // Aggiungo correzione distanza faccia da centro oggetto
            vFaceCenter = vFaceCenter - sPlaneInfo.sPlane.zaxis  * mScale[2][2] * vBBoxSide.z/2;

            // Aggiungo offset
            vFaceCenter = vFaceCenter - sPlaneInfo.sPlane.xaxis * sObjInfo.fOffsetX * mScale[0][0] - sPlaneInfo.sPlane.yaxis * sObjInfo.fOffsetY * mScale[1][1] - sPlaneInfo.sPlane.zaxis * sObjInfo.fOffsetZ * mScale[2][2];
        }
        else  if (sLabelInfo.iViewFace == ESideType::eBottom)
        {
            // Calcolo rotazione 
            mXform[0][0] = -sPlaneInfo.sPlane.xaxis.x;
            mXform[1][0] = -sPlaneInfo.sPlane.xaxis.y;
            mXform[2][0] = -sPlaneInfo.sPlane.xaxis.z;
            mXform[0][1] = -sPlaneInfo.sPlane.yaxis.x;
            mXform[1][1] = -sPlaneInfo.sPlane.yaxis.y;
            mXform[2][1] = -sPlaneInfo.sPlane.yaxis.z;
            mXform[0][2] = -sPlaneInfo.sPlane.zaxis.x;
            mXform[1][2] = -sPlaneInfo.sPlane.zaxis.y;
            mXform[2][2] = -sPlaneInfo.sPlane.zaxis.z;

            // Aggiungo correzione distanza faccia da centro oggetto
            vFaceCenter = vFaceCenter - sPlaneInfo.sPlane.zaxis * mScale[2][2] * vBBoxSide.z / 2;

            // Aggiungo offset
            vFaceCenter = vFaceCenter - sPlaneInfo.sPlane.xaxis * sObjInfo.fOffsetX * mScale[0][0] - sPlaneInfo.sPlane.yaxis * sObjInfo.fOffsetY * mScale[1][1] - sPlaneInfo.sPlane.zaxis * sObjInfo.fOffsetZ * mScale[2][2];
        }
    }

    mXform[0][3] = vFaceCenter.x;
    mXform[1][3] = vFaceCenter.y;
    mXform[2][3] = vFaceCenter.z;

    mXform[3][0] = mXform[3][1] = mXform[3][2] = 0;
    mXform[3][3] = 1;
}
int DetectPlanesParallelToSegment(PointCloud<PointType>::Ptr pCloud, SFragmentFamily & SFagmentFamily, CCmdOptions & oCmdOptions)
{
    Eigen::Vector3f vPlaneConstrain(SFagmentFamily.oFragmentPlane.zaxis.x, SFagmentFamily.oFragmentPlane.zaxis.y, SFagmentFamily.oFragmentPlane.zaxis.z);

    bool bContinue = true;
    float fAngularTol = 4;
    while (bContinue)
    {
        // Create "PerpendicularPlane" model
        SampleConsensusModelPerpendicularPlane<PointType>::Ptr model_p(new SampleConsensusModelPerpendicularPlane<PointType>(pCloud));
        model_p->setAxis(vPlaneConstrain);
        model_p->setEpsAngle(fAngularTol * 3.1415926 / 180);  // angle tol=1°

        // Created RandomSampleConsensus object for that model
        RandomSampleConsensus<PointType> ransac(model_p);
        ransac.setDistanceThreshold(.02);
        ransac.setMaxIterations(5000);

        // Compute 
        bool bResult = ransac.computeModel();
        if (!bResult)  // Computation failed 
            break;

        // Refine plane parameters
        //ransac.refineModel();

        // Get inliers
        Indices inliers;
        ransac.getInliers(inliers);

        if (inliers.size() > 10000)
        {
            // Get plane coeff
            Eigen::VectorXf model_coefficients;
            ransac.getModelCoefficients(model_coefficients);

            // Build candidate plane
            ON_PlaneEquation oEquation(model_coefficients[0], model_coefficients[1], model_coefficients[2], model_coefficients[3]);
            ON_Plane oPlane(oEquation);

//+++++++++++++++++++++++++++++++++++++
//            float fTmp1 = abs(SFagmentFamily.oFragmentPlane.zaxis * oPlane.zaxis);
//            float fTresh = cos(2 * 3.1415926 / 180);
//            if (fTmp1 < fTresh)
//                continue;
//+++++++++++++++++++++++++++++++++++++


            LocalY2GlobalZ(oPlane); // Ruoto il TargetPlane in modo da portare l'asse Y il più possibile allineato allo Z globale

            // Creo nuvola temporanea con gli inliers
            PointCloud<PointType>::Ptr pTmp(new PointCloud<PointType>);
            copyPointCloud(*pCloud, inliers, *pTmp);

            //+++++++++++++++++++++++++++++++
            static int iCount = 0;
            string WallFilename = oCmdOptions.sSegmentePath + string("WallFinder") + to_string(iCount) + string(".pcd");
            io::save(WallFilename, *pTmp);
            iCount++;
            //+++++++++++++++++++++++++++++++


            PointUV vUVMin, vUVMax;
            ProjectPointOnPlane(oPlane, pTmp, vUVMin, vUVMax);

            //+++++++++++++++++++++++++++++++
            float fDensity = pTmp->size() / ((vUVMax.u - vUVMin.u) * (vUVMax.v - vUVMin.v));
            cout << " segment :"<<WallFilename << " fDensity=" << fDensity << endl;
                //+++++++++++++++++++++++++++++++

            if (vUVMax.u - vUVMin.u > oCmdOptions.fMinWallLength && vUVMax.v - vUVMin.v > oCmdOptions.fMinWallHeight)
            {
                // Add current surface to candidate list
                SSurface sSurface;
                sSurface.oPlane = oPlane;
                sSurface.vUVMin = vUVMin;
                sSurface.vUVMax = vUVMax;
                SFagmentFamily.aCandidate.push_back(sSurface);

                // Remove current inliers from global cloud
                Indices iAll;
                Indices iDiff;
                iAll.resize(pCloud->size());
                for (int k = 0; k < pCloud->size(); k++)
                    iAll[k] = k;   // iAll == elenco completo di indici

                set_difference(iAll.begin(), iAll.end(), inliers.begin(), inliers.end(), inserter(iDiff, iDiff.begin()));  // Rimuovo gli inliers

                if (iDiff.size() < pCloud->size() / 20)  // Remaining point fewer than 5% of initial cloud
                    break;

                // Creo nuvola temporanea con i soli ounti restanti
                PointCloud<PointType>::Ptr pCloudTmp(new pcl::PointCloud<PointType>);
                copyPointCloud(*pCloud, iDiff, *pCloudTmp);

                // Eseguo swap dei punti
                pCloud->points.swap(pCloudTmp->points);
            }

        }
        else
            bResult = false;

    }


    return 0;
}

#endif


int LoadAndSliceCloud(PointCloud<PointType>::Ptr pCloud, ON_Plane& oFloorPlane, ON_Plane& oCeilingPlane, string& sPCFileName)
{
    int iErr;
    // Read fragment cloud
    PointCloud<PointType>::Ptr pFullCloud(new PointCloud<PointType>);
    if (sPCFileName.rfind(string(".e57")) != string::npos || sPCFileName.rfind(string(".E57")) != string::npos)
        iErr = ReadA57_d((const char*)(sPCFileName.c_str()), *pFullCloud);
    else
        iErr = pcl::io::load((const char*)(sPCFileName.c_str()), *pFullCloud);
    if (iErr == -1)     if (iErr == -1)
    {
        cout << "Cannot read file " << sPCFileName << endl;
        return -1;
    }

    float fXmin = 1e+30;
    float fYmin = 1e+30;
    float fXmax = -1e+30;
    float fYmax = -1e+30;
    float fZZmin = 1e+30;
    float fZZmax = -1e+30;

    float zMin = oFloorPlane.origin.z;
    float zMax = oCeilingPlane.origin.z;
    for (auto& point : *pFullCloud)
    {
///        if (point.z > zMin + 0.5 && point.z < zMax - 0.5 && point.x> -27 && point.x< 23 && point.y>-38.5 && point.y < 21.5)   // NOTA: LIMITI CABLATI SU ESEMPIO!!!!!
          if (point.z > zMin + 0.5 && point.z < zMax - 0.5)
        {
            fXmin = min(fXmin, point.x);
            fYmin = min(fYmin, point.y);
            fXmax = max(fXmax, point.x);
            fYmax = max(fYmax, point.y);

            pCloud->points.push_back(point);
        }

        fZZmin = min(fZZmin, point.z);
        fZZmax = max(fZZmax, point.z);
    }
    pCloud->width = pCloud->size();
    pCloud->height = 1;
    pCloud->is_dense = false;

    //++++++++++++++++++++++++++++++++++
    iErr = io::save("D:\\Andrea\\Dev\\C++\\PointCloud\\Package - 2022\\Project3\\SegmentedPC\\SlicedCloud.pcd", *pCloud);
    //++++++++++++++++++++++++++++++++++


    //------------------------------------------------------

    float fDeltaX = fXmax - fXmin;
    float fDeltaY = fYmax - fYmin;
    float fDs = 0.01;   // 1 pixel==1cm

    int nX = fDeltaX / fDs + 1;
    int nY = fDeltaY / fDs + 1;

    CAMImage imgPng(nX, nY, 1);
    unsigned char* pImgBuffer = imgPng.GetRawData();

    float* pfBuffer = new float[nX * nY];
    //    float* pfGradBuffer = new float[nX * nY];

    memset(pfBuffer, 0, nX * nY*sizeof(float));
    //    memset(pfGradBuffer, 0, nX * nY);

    float fMax = 0;
    float fMean = 0;
    int nTot = 0;

    // Riempio pfBuffer
    int k, kX, kY;
    for (auto& point : *pCloud)
    {
        kX = (point.x - fXmin) / fDs;
        kY = nY - 1 - (point.y - fYmin) / fDs;

        if (point.y < 80)
            point.y = point.y;

        if (kX > nX - 1)
            kX = nX - 1;
        else if (kX < 0)
            kX = 0;

        if (kY > nY - 1)
            kY = nY - 1;
        else if (kY < 0)
            kY = 0;

        k = kX + kY * nX;
        pfBuffer[k]++;
        fMax = max(fMax, pfBuffer[k]);
    }

    vector<float> aDen;
    // Calcolo valor medio (su valori non nulli)
    float fPositiveMin = 1e+30;
    for (int k = 0; k < nX * nY; k++)
    {
        if (pfBuffer[k] > 0)
        {
            fPositiveMin = min(fPositiveMin, pfBuffer[k]);
            aDen.push_back(pfBuffer[k]);
            fMean += pfBuffer[k];
            nTot++;
        }
    }
    fMean /= nTot;


//+++++++++++++++++++++++++++++++++++
    memset(pImgBuffer, 255, nX * nY);
    for (int k = 0; k < nX * nY; k++)
    {
        if(pfBuffer[k]>0)
            pImgBuffer[k] = 255*(1.f- pfBuffer[k]/ fMax);
    }
    string sTmp = string("D:\\Andrea\\Dev\\C++\\PointCloud\\Package - 2022\\Project3\\Prova.Png");
    imgPng.WritePng(sTmp.c_str());

//+++++++++++++++++++++++++++++++++++





    using namespace boost::histogram;
    int nStep = 256;
    double fDZ = (fMax - fPositiveMin) / nStep;
    if (fDZ < 1)
    {
        fDZ = 1;
        nStep = int(fMax - fPositiveMin) + 1;

    }

    auto Hist = make_histogram(axis::regular<>(nStep, fPositiveMin, fMax));
    Hist.fill(aDen);

    static bool s_bDebugPanel = false;
    if (s_bDebugPanel)
    {
        vector<double> aX;
        vector<double> aY;

        aX.resize(Hist.size() - 1);
        aY.resize(Hist.size() - 1);
        aX[0] = fPositiveMin;
        aY[0] = Hist[0];
        for (int i = 1; i < Hist.size() - 1; i++)
        {
            aX[i] = aX[i - 1] + fDZ;
            aY[i] = Hist[i];
        }

        visualization::PCLPlotter oPlotter1;
        oPlotter1.addPlotData(aX, aY);
        oPlotter1.plot();
    }

    float afThrash[6] = { .75,.8,.85,.9,.95,.98 };
    float fT=0;
    for (int kk = 0; kk < 6; kk++)
    {
        float fTmp = 0;
        for (int i = 0; i < Hist.size(); i++)
        {
            fTmp += Hist[i];
            if (fTmp / nTot > afThrash[kk])
            {
                fT = fPositiveMin + i * fDZ;
                break;
            }
        }

        memset(pImgBuffer, 255, nX* nY);
        for (int k = 0; k < nX * nY; k++)
        {
            if (pfBuffer[k] > fT)
                pImgBuffer[k] = 0;
        }
        string sTmp = string("D:\\Andrea\\Dev\\C++\\PointCloud\\Package - 2022\\Project3\\ProvaBis") + to_string(kk) + string(".png");
        imgPng.WritePng(sTmp.c_str());


    }



    // Calcolo densità corretta
    /*
    float fNorm;
    float fDeltafNorm = (fMax - fMean) / 20;
    for (int iTest0 = 1; iTest0 < 2; iTest0++)
    {
        fNorm = fMean+ iTest0* fDeltafNorm;

        for (int iTest1=1; iTest1 < 2; iTest1++)
        {

            cout << "fMax=" << fMax << "  fMean=" << fMean << " iTest0=" << iTest0 << " iTest1=" << iTest1 << " fNorm=" << fNorm << endl;

            for (int k = 0; k < nX * nY; k++)
            {
                float fBuf = pfBuffer[k] / fNorm;
                if (fBuf > 1)
                {
                    pImgBuffer[k] = 255;
                    pfBuffer[k] = 1;
                }
                else if (fBuf > iTest1)
                {
                    pImgBuffer[k] = 255 * fBuf;
                    pfBuffer[k] = fBuf;
                }
                else if (pfBuffer[k] > 0)
                {
                    pImgBuffer[k] = 0;
                    pfBuffer[k] = 0;
                }

            }
            string sTmp = string("D:\\Andrea\\Dev\\C++\\PointCloud\\Package\\Project3\\Prova") + to_string(iTest0) +string("_")+to_string(iTest1)+ string(".png");
            imgPng.WritePng(sTmp.c_str());
            memset(pImgBuffer, 0, nX * nY);
        }
    }

*/

    return 0;
}


int DetectFragmentPlane(ON_Plane& oPlane, string sFragmentFile)
{
    // Read fragment cloud
    PointCloud<PointType>::Ptr pCloud(new PointCloud<PointType>);
    int iErr = pcl::io::load((const char*)(sFragmentFile.c_str()), *pCloud);
    if (iErr == -1)
    {
        cout << "Cannot read file " << sFragmentFile << endl;
        return -1;
    }


    // Create "PerpendicularPlane" model
    SampleConsensusModelPlane<PointType>::Ptr model_p(new SampleConsensusModelPlane<PointType>(pCloud));

    // Created RandomSampleConsensus object for that model
    RandomSampleConsensus<PointType> ransac(model_p);
    ransac.setDistanceThreshold(.01);
    ransac.setMaxIterations(5000);

    // Compute 
    bool bResult = ransac.computeModel();
    if (!bResult)  // Computation failed
    {
        cout << "Cannot extract reference plane from fragment " << sFragmentFile << endl;
        return -2;
    }

    // Optimize plane detection
    ransac.refineModel();

    // Get inliers
    Indices inliers;
    ransac.getInliers(inliers);


    if (inliers.size() < pCloud->size() * .8)  // Check if num inliers > 80% of all cloud point
        cout << "WARNING: plane detected from fragment " << sFragmentFile << " might suffer of low precision." << endl;

    // Get plane coeff
    Eigen::VectorXf model_coefficients;
    ransac.getModelCoefficients(model_coefficients);

    // Build candidate plane
    ON_PlaneEquation oEquation(model_coefficients[0], model_coefficients[1], model_coefficients[2], model_coefficients[3]);
    oPlane = ON_Plane(oEquation);

    return 0;
}

using json = nlohmann::json;
int GetWallDirFromJsonFile(const char* sFileName, vector<Eigen::Vector3f> & aWallNormal)
{

    FILE* pFile = fopen(sFileName, "r");
    if (pFile == NULL)
        return -1;

    auto jDoc = json::parse(pFile);
    fclose(pFile);

    // Get shapes
    auto aShapes = jDoc["shapes"];
    int ind = 0;
    for (auto shape = aShapes.begin(); shape != aShapes.end(); shape++)
    {
        // Skip shape if not "polygon"
        string sShapeType = (*shape)["shape_type"];
        if (sShapeType != string("polygon"))
            continue;

        auto aPoints = (*shape)["points"];
        vector<float> aPolyPointX, aPolyPointY;

        for (auto point = aPoints.begin(); point != aPoints.end(); point++)
        {
            aPolyPointX.push_back((*point)[0]);
            aPolyPointY.push_back((*point)[1]);
        }
        if (aPolyPointY.size() != 2)  // Invalid polygon
            continue;

        Eigen::Vector3f vNormal(aPolyPointY[1]- aPolyPointY[0], aPolyPointX[0] - aPolyPointX[1],0);  // 
        vNormal.normalize();
        aWallNormal.push_back(vNormal);

    }


    return 0;
}

void GetBorderIntersect(Eigen::Vector3f& vCentre, Eigen::Vector3f& vDir, int nX, int nY, float & xx, float & yy, int& iBord)
{
    float tX=1e+30, tY = 1e+30;

    float vX = vDir.x();
    int iBorderX = -1;
    if (vX > 0)
    {
        tX = (nX-1 - vCentre.x()) * 1.0 / vX;
        iBorderX = 1; // Right
    }
    else if (vX < 0)
    {
        tX = -vCentre.x() * 1.0 / vX;
        iBorderX = 3; // Left
    }

    float vY = vDir.y();
    int iBorderY = -1;
    if (vY > 0)
    {
        tY = (nY-1 - vCentre.y()) * 1.0 / vY;
        iBorderY = 2; // Bottom
    }
    else if (vY < 0)
    {
        tY = -vCentre.y() * 1.0 / vY;
        iBorderY = 0; // Top
    }

    if (tX < tY)
    {
        iBord = iBorderX;
        if (iBord == 1)
            xx = nX-1;
        else
            xx = 0;

        yy = vCentre.y() + vY * tX;
        if (yy < 0)
            yy = 0;
        else if (yy >= nY)
            yy = nY - 1;
    }
    else
    {
        iBord = iBorderY;
        if (iBord == 2)
            yy = nY-1;
        else
            yy = 0;

        xx = vCentre.x() + vX * tY;
        if (xx < 0)
            xx = 0;
        else if (xx >= nX)
            xx = nX - 1;
    }
}

void CreateFftMask(CAMImage& imgFFTMask, Eigen::Vector3f& vDir, int nBandWidth, int nCenterRadius)
{
    // Reset filter
    imgFFTMask.Reset();
    Eigen::Vector3f vNegDir(-vDir.x(), -vDir.y(), 0);
    Eigen::Vector3f vNorm(vDir.y(), -vDir.x(), 0);
    int nX = imgFFTMask.GetSizeX();
    int nY = imgFFTMask.GetSizeY();

    Eigen::Vector3f vCentre(nX / 2, nY / 2, 0);

    if (nBandWidth > 0)
    {
        Eigen::Vector3f vCentreDown = vCentre - vNorm * nBandWidth;  // Mi sposto dal centro di bandwidth (nel verso opposto alla normale)
        Eigen::Vector3f vCentreUp = vCentre + vNorm * nBandWidth; // Mi sposto dal centro di bandwidth (nel verso della normale)

        int iBord[4]; // 0==>Top, 1==>Right, 2==>Bottom, 3==>Left
        vector<float> xx, yy;
        xx.resize(4);
        yy.resize(4);

        GetBorderIntersect(vCentreDown, vDir, nX, nY, xx[0], yy[0], iBord[0]);
        GetBorderIntersect(vCentreDown, vNegDir, nX, nY, xx[1], yy[1], iBord[1]);
        GetBorderIntersect(vCentreUp, vNegDir, nX, nY, xx[2], yy[2], iBord[2]);
        GetBorderIntersect(vCentreUp, vDir, nX, nY, xx[3], yy[3], iBord[3]);

        int iTopY = nY + 1;
        int iBottomY = -1;
        int iLeftX = nX + 1;
        int iRightX = -1;

        for (int i = 0; i < 4; i++)
        {
            iTopY = min(iTopY, (int)yy[i]);
            iBottomY = max(iBottomY, (int)yy[i]);
            iLeftX = min(iLeftX, (int)xx[i]);
            iRightX = max(iRightX, (int)xx[i]);
        }

        // Create white band
        imgFFTMask.PolyRasterize(xx, yy, iTopY, iBottomY, iLeftX, iRightX, 255);
    }
    else
        imgFFTMask.Reset(255);

    if (nCenterRadius > 0)
    {
        // Add black circle at centre
        int nRad2 = nCenterRadius * nCenterRadius;
        for (int i = -nCenterRadius; i <= nCenterRadius; i++)
        {
            for (int j = -nCenterRadius; j <= nCenterRadius; j++)
            {
                if (i * i + j * j <= nRad2)
                    imgFFTMask.SetR(i + vCentre.x(), j + vCentre.y(), 0);
            }

        }
    }
    else
    {
        // Add black circle out of centre
        int nRad2 = nCenterRadius * nCenterRadius;
        for (int i = nCenterRadius; i <= -nCenterRadius; i++)
        {
            for (int j = nCenterRadius; j <= -nCenterRadius; j++)
            {
                if (i * i + j * j >= nRad2)
                    imgFFTMask.SetR(i + vCentre.x(), j + vCentre.y(), 0);
            }

        }
    }


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    imgFFTMask.WritePng("D:\\Andrea\\Dev\\C++\\PointCloud\\Package\\Project3\\MASK_ProvaBis1.png");
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++


}


int main(int argc, char** argv)
//***********************************************************************************************************************************************
// 0. Load project settings                                 [ReadCommandFile()]
// 1. Load xml file containing camera infos                 [CObjectCatalogManager::LoadFromXmlFile()]
// 2. Loop over all lbi objects
//    2a. if .lbi is older than .ply --> just skip
//    2b. Lettura file .lbi                                 [ReadLabelInfoFile()]
//    3c. Loop sulle label associate al segmento corrente (tipicamente 1!)
//        3ca.  Lettura del primo LabelMAnager valido       [CLabelManager::Load()]
//        3cb.  Get ObjInfo                                 [GetObjectInfo()]
//        3cc.  Get 3D model file name (.ply)      
//    3d. Load point cloud associated to current segment
//    3e. Loop over all cameras used for current object segmentation (tipicamente 1!)
//        3ea.  Get camera infos
//        3eb.  Get Label manager
//        3ec.  
//
//
//
//
//***********************************************************************************************************************************************
{
    if (argc < 2)
    {
        cout << "Missing command file name on command line!";
        cout << " Usage: WallFinder.exe CmdOptions.xml" << endl;
        return -1;
    }

    CCmdOptions oCmdOptions;
    int iErr = oCmdOptions.ReadCommandFile(argv[1]);

    if (iErr < 0)
        return -1;

    // Get Cloud with point included between floor and ceiling
    PointCloud<PointType>::Ptr pCloud(new PointCloud<PointType>);
    ON_Plane oFloorPlane;
    ON_Plane oCeilingPlane;
//    double zCentre= 106.5;
    double zCentre = 0;
        oFloorPlane.origin.z = -1.3 - zCentre;
    oCeilingPlane.origin.z = 7.42- zCentre;
    LoadAndSliceCloud(pCloud, oFloorPlane, oCeilingPlane, oCmdOptions.sPCFileName);


    ON_3dVector vGlobalZ( double(0), double(0), double(1));

    vector< SFragmentFamily> aWallSurfaces;
    if (0)   // Temporaneamente escluso!!!!!!!!!!!!!!!!!!!!!!!
    {
        // Read floor cloud segment and detect floor plane
        ON_Plane oFloorPlane;
        iErr = DetectFragmentPlane(oFloorPlane, oCmdOptions.sFloorSegmentFile);
        if (oFloorPlane.Normal() * vGlobalZ < 0)
            oFloorPlane.Flip();

        // Read ceiling cloud segment and detect ceiling plane
        ON_Plane oCeilingPlane;
        iErr = DetectFragmentPlane(oCeilingPlane, oCmdOptions.sCeilingSegmentFile);
        if (oCeilingPlane.Normal() * vGlobalZ < 0)
            oCeilingPlane.Flip();

        float fDelta = 2;  // Soglia di tolleranza (in gradi) perchè 2 segment siano considerati paralleli
        float fThresh = cos(fDelta * 3.1415926 / 180);
        SFragmentFamily sFragmentFamily;
        // Read wall fragments and remove duplicates
        for (int i = 0; i < oCmdOptions.asWallSegment.size(); i++)
        {
            iErr = DetectFragmentPlane(sFragmentFamily.oFragmentPlane, oCmdOptions.asWallSegment[i]);
            if (iErr)
                continue;

            // Check if duplicates
            bool bOK = true;
            for (int k = 0; k < aWallSurfaces.size(); k++)
            {
                float fTmp = fabs(aWallSurfaces[k].oFragmentPlane.Normal() * sFragmentFamily.oFragmentPlane.Normal());
                if (fTmp > fThresh)
                {
                    bOK = false;
                    break;
                }
            }

            if (bOK)
                aWallSurfaces.push_back(sFragmentFamily);
        }

        // Get Cloud with point included between floor and ceiling
        PointCloud<PointType>::Ptr pCloud(new PointCloud<PointType>);
        LoadAndSliceCloud(pCloud, oFloorPlane, oCeilingPlane, oCmdOptions.sPCFileName);
    }

    // Read walss normals from Json file
    vector<Eigen::Vector3f> aWallNormal;
    GetWallDirFromJsonFile(oCmdOptions.sWallDirectionsFile.c_str(), aWallNormal);

    //++++++++++++++++++++++++++++++++++
    CAMImage imgPng;
    imgPng.ReadFile("D:\\Andrea\\Dev\\C++\\PointCloud\\Package\\Project3\\ProvaBis1.png");

    int nX = imgPng.GetSizeX();
    int nY = imgPng.GetSizeY();
    int n = nX * nY;

//+++++++++++++++++++++++++++++++++++++++++++++++++++
    /*
    CAMImage imgFFTMask1;
    imgFFTMask1.Init(nX, nY, 1);
    int nBandWidth = 30;
    int nCenterRad = 10;
    CreateFftMask(imgFFTMask1, aWallNormal[0], nBandWidth, nCenterRad);
    */
//++++++++++++++++++++++++++++++++++++++++++++++++++++++


    // ---------------- Compute FFT of density image ---------------- 

    // Prepare fft data
    shape_t iShape; iShape.push_back(nX); iShape.push_back(nY);
    stride_t iStrideIn{ sizeof(complex<float>), sizeof(complex<float>) * nX };
    stride_t iStrideOut{ sizeof(complex<float>), sizeof(complex<float>) * nX };
    shape_t iAxes; iAxes.push_back(0); iAxes.push_back(1);
    vector<complex<float>> cImgData;
    vector<complex<float>> cDataOut; cDataOut.resize(n);

    // Copy image uint buffer in complex buffer
    cImgData.resize(n);
    int k = 0;
    for (int j = 0; j < nY; j++)
    {
        for (int i = 0; i < nX; i++)
        {
            cImgData[k] = imgPng.GetR(i, j);
            k++;
        }
    }

    // Perform fft computation
    c2c(iShape, iStrideIn, iStrideOut, iAxes, FORWARD, cImgData.data(), cDataOut.data(), 1.f);

    // Shift high freq to contre
    vector<complex<float>> cShiftedFft;
    cShiftedFft.resize(n);
    fftshift(cShiftedFft.data(), cDataOut.data(), nY, nX);


    if (0)  // Scrittura png della fft Temporaneamente esclusa
    {
        CAMImage imgFFT;
        float fMaxIntensity = 0;
        for (int i = 0; i < n; i++)
            fMaxIntensity = max(fMaxIntensity, powf(abs(cShiftedFft[i]), 0.35));

        cout << "fMaxIntensity=" << fMaxIntensity << endl;

        imgFFT.Init(nX, nY, 1);
        unsigned char* aucTmp = imgFFT.GetRawData();

        //    vector<unsigned char> aucTmp;
        //    aucTmp.resize(n);

        for (int i = 0; i < n; i++)
            aucTmp[i] = unsigned char(255.0 * (powf(abs(cShiftedFft[i]), 0.35) / fMaxIntensity));

        //    fftshift(imgFFT.GetRawData(), aucTmp.data(), nY, nX);

        imgFFT.WritePng("D:\\Andrea\\Dev\\C++\\PointCloud\\Package\\Project3\\FFT_ProvaBis1.png");
    }



    // Loop over all wall fragments
    CAMImage imgFFTMask;
    imgFFTMask.Init(nX, nY, 1);
    int nBandWidth = -1;
    int nCenterRad = -500;// 50;
    for(int i=0; i< aWallNormal.size(); i++)
    {

        // Create mask for current direction
        CreateFftMask(imgFFTMask, aWallNormal[i], nBandWidth, nCenterRad);

        // Apply mask to FFT
        unsigned char* aucMask = imgFFTMask.GetRawData();
        for (k = 0; k < n; k++)
            cShiftedFft[k] *= float(aucMask[k]) / 255;

        //++++++++++++++++++++++++
        CAMImage imgFFT;
        imgFFT.Init(nX, nY, 1);
        unsigned char* aucTmp = imgFFT.GetRawData();

        float fMaxIntensity = 0;
        for (int i = 0; i < n; i++)
            fMaxIntensity = max(fMaxIntensity, powf(abs(cShiftedFft[i]), 0.35));

        for (int i = 0; i < n; i++)
            aucTmp[i] = unsigned char(255.0 * (powf(abs(cShiftedFft[i]), 0.35) / fMaxIntensity));

        imgFFT.WritePng("D:\\Andrea\\Dev\\C++\\PointCloud\\Package\\Project3\\MASKED_FFT_ProvaBis1.png");
        //++++++++++++++++++++++++

        // Shift back
        ifftshift(cDataOut.data(), cShiftedFft.data(), nY, nX);

        // Perform inverse fft computation
        c2c(iShape, iStrideIn, iStrideOut, iAxes, BACKWARD, cDataOut.data(), cImgData.data(), 1.f/n);

        //+++++++++++++++++++++++++++++++++
        int k = 0;
        unsigned char ucTmp;
        for (int j = 0; j < nY; j++)
        {
            for (int i = 0; i < nX; i++)
            {
                ucTmp = unsigned char(cImgData[k].real());
                imgPng.SetR(i, j, ucTmp);
                k++;
            }
        }
        imgPng.WritePng("D:\\Andrea\\Dev\\C++\\PointCloud\\Package\\Project3\\ProvaBis1_FILTERED.png");



        // Detect all valid planes parallel to current wall fragment
//        iErr = DetectPlanesParallelToSegment(pCloud, aWallSurfaces[i], oCmdOptions);

        // Identify facing planes --> 3D wall
        //.....................

        // Build .ply model
        //................

        // Update oCRSTM
    }

    // Identify perimeter

    // Save RevitSummaryTable

}

#if defined(pippo)
{
    // Load xml file containing camera infos
    CSensorManager oSM;
    int iEr=oSM.LoadSensorAndCamera(oCmdOptions);
    if (iEr != 0)
    {
        cout << "Cannot read file " << oCmdOptions.sInputCameraFile << endl;
        return -1;
    }

    //  Load xml file containing object catalog
    CObjectCatalogManager oCM;
    XMLError xmlIer = oCM.LoadFromXmlFile(oCmdOptions.sInputCatalogFile.c_str());
    if (xmlIer != XML_SUCCESS)
    {
        cout << "Cannot read file " << oCmdOptions.sInputCatalogFile << endl;
        return -1;
    }

    CRSTManager oCRSTM;

    string sTmp;
    int i = 0;
    string sExtFilter=".lbi";

    // ----------------- Loop over all lbi objects ---------------------------
    for (const auto& file : filesystem::directory_iterator(oCmdOptions.sSegmentePath))
    {
        sTmp = filesystem::path(file).extension().string();  // Get current file
        if(_stricmp((const char*)(sTmp.c_str()), (const char*)(sExtFilter.c_str())))  // If not .lbi just skip it
            continue;
        // Now file is current file with label infos
        filesystem::path oPath(file);

        // Check dependency: if .lbi is older than .ply --> just skip
        struct stat resultLbi,resultPly;
        stat(oPath.string().c_str(), &resultLbi);  // Get .lbi infos

        if (stat(oPath.replace_extension(".ply").string().c_str(), &resultPly) == 0)  // File exist
        {
            if (resultPly.st_mtime > resultLbi.st_mtime)  // .lbi is older than .ply --> no need to recalculate
                continue;
        }

        cout << "Processing segment " << filesystem::path(file).filename().string() << endl;

        // Lettura file lbi associato
        vector< SSegmentLabelInfo> aSegmentLabelInfo;
        sTmp = filesystem::path(file).string();
        ReadLabelInfoFile(sTmp, aSegmentLabelInfo);

        CLabelManager oLM;
        SCandidatePlaneInfo sPlaneInfo, sPlaneInfoTmp;
        SObjCatalogInfo sObjInfo;
        SLabelInfo sLabelInfo, sLabelInfoTmp;
        string sModelFileName;
        bool bSkip = true;
        for (int j = 0; j < aSegmentLabelInfo.size(); j++)
        {
            // Init label manager for current camera
            filesystem::path oPath(aSegmentLabelInfo[j].sLabelFileName);
            int iErr = oLM.Load(oPath, aSegmentLabelInfo[j].iLabelFormat, oCM);
            if (iErr < 0)
                continue;

            // Get LabelInfo and ObjInfo
            sLabelInfoTmp = oLM.GetLabelInfoFromInd(aSegmentLabelInfo[j].iLabelInd);
            sObjInfo = oCM.GetObjectInfo(sLabelInfoTmp.iObjCatalogID);

            // Get 3D model file name (.ply)
            sModelFileName = string(oCM.GetModelLibraryPath() + "\\" + sObjInfo.sModelName);
            if (sObjInfo.sModelName != "" && filesystem::exists(sModelFileName))  // valid name AND file exist
            {
                bSkip = false;
                break;
            }
        }
       
        if (bSkip)  // No 3d model available --> skip current object
        {
            cout << "Invalid model name, or 3d model file is missing" << endl;
            cout << "This segment will be skipped" << endl << endl;
            continue;
        }

        cout << "  ... loading point cloud...";
        // Load point cloud
        PointCloud<PointType>::Ptr pCloud(new PointCloud<PointType>);
        sTmp = filesystem::path(file).replace_extension(".pcd").string();
        iErr = pcl::io::load((const char*)(sTmp.c_str()), *pCloud);
        if (iErr == -1)
        {
            cout << "Cannot read file " << filesystem::path(file).string().c_str() << endl;
            cout << "This segment will be skipped" << endl << endl;
            continue;
        }

        cout << "... loaded!"<< endl;

        string sSegmentFileName = sTmp;

        cout << "  ... computing position ...";

        // Loop over all cameras used for current object segmentation
        bool bPositionFound = true;
        for (int j = 0; j < aSegmentLabelInfo.size(); j++)
        {

            // ------------------ In questo blocco viene costruita la matrice di proiezione associata alla camera ---------------------------------
            // Get index of camera with same name of current mask image
            int iRotIndex = 0;
            int iIndex= oSM.GetCameraIndex((char*)(aSegmentLabelInfo[j].sCameraName.c_str()),&iRotIndex);
            if (iIndex < 0)
            {
                cout << endl << "Cannot find camera named "<< aSegmentLabelInfo[j].sCameraName << endl;
                continue;
            }

            // Get camera infos
            SCameraInfo sCameraInfo = oSM.GetCameraInfo(iIndex, iRotIndex);  // [R, cameraPos] = SM.GetCameraRT(i);

            // Init label manager for current camera
            filesystem::path oPath(aSegmentLabelInfo[j].sLabelFileName);
            int iErr = oLM.Load(oPath, aSegmentLabelInfo[j].iLabelFormat, oCM);
            if (iErr < 0)
                continue;

//+++++++++++++++++++++++++++++++++++++
///            // Pulisco l'oggetto usando le normali (calcolate al volo)
///            CleanCloudUsingNormal(pCloud, sCameraInfo);
//+++++++++++++++++++++++++++++++++++++


//+++++++++++++++++++++++++++++++++++++
//            // Verifico la presenza di cluster e nel caso rimuovo quelli più lontani dalla camera
//            SearchAndRemoveFarClusters(pCloud, sCameraInfo);
//+++++++++++++++++++++++++++++++++++++


            int iMaxIter = 15;
            if(sObjInfo.sShapeType == "Cube")
            {
                sLabelInfo = oLM.GetLabelInfoFromInd(aSegmentLabelInfo[j].iLabelInd);

                sPlaneInfoTmp = DetectPlane(pCloud, sCameraInfo, oCmdOptions.bForceVertical, sLabelInfo.iViewFace, iMaxIter,sObjInfo);
                if (sPlaneInfoTmp.iError > 0)  // Errore non gestibile nell'algoritmo di identificazione
                {
                    cout << "Cannot succesfully process segment " << sSegmentFileName <<". It will be skipped"<< endl << endl;
                    bPositionFound = false;
                    continue;
                }
                else if (sPlaneInfoTmp.nConvergenceIter >= iMaxIter)  // Fallita convergenza --> proviamo a pulire nuvola
                {
                    CleanCloudUsingNormal(pCloud, sCameraInfo);

                    //++++++++++++++++++++++++++
                    //iErr = io::save("D:\\Andrea\\Dev\\C++\\PointCloud\\Package\\Project2\\SegmentedPC\\PillarCleaned.pcd", *pNewCloud);
                    //++++++++++++++++++++++++++

                    // Ritento calcolo con nuvola ripulita 
                    sPlaneInfoTmp = DetectPlane(pCloud, sCameraInfo, oCmdOptions.bForceVertical, sLabelInfo.iViewFace, iMaxIter, sObjInfo);

                }

                if (j == 0)  // Prima camera
                    sPlaneInfo = sPlaneInfoTmp;
                else if (sPlaneInfoTmp.fConvergenceMeanErr < sPlaneInfo.fConvergenceMeanErr)
                    sPlaneInfo = sPlaneInfoTmp;
            }
            else if (sObjInfo.sShapeType == "Sphere")
            {
                continue; // TODO TODO
            }

        }

        // If location is failed skip to next segment
        if (bPositionFound)
            cout << "... done!" << endl;
        else
        {
            cout << "...computation failed!!!! This segment will be skipped." << endl << endl;
                continue;
        }


        cout << "  ... read 3d model file ...";
        // Leggo modello di oggetto
        CPlyIO< ON_3dPoint, ON_3dVector, ON_Xform> oPlyIO;
        oPlyIO.ReadFile(sModelFileName.c_str());
        oPlyIO.SetOrigInBBoxCentre();  // Necessary to correctly manage rotations
        SRSTInfo sRSTInfo;

        cout << "... done!" << endl;

        if (sObjInfo.sShapeType == "Cube")
        {
            cout << "  ... computing transformation ...";

            // Compute model transformation (scale,rotation and translation)
            ON_Xform mXScale = ON_Xform::IdentityTransformation;
            ON_3dPoint vBBoxSide = oPlyIO.GetBBoxSide();
            sRSTInfo.vOriginalBBoxSize = { (float)vBBoxSide.x,(float)vBBoxSide.y,(float)vBBoxSide.z};
            if (sObjInfo.iVerticalAxis == 2 && sObjInfo.iFrontAxis == 0)  // Axis: Z oggetto --> Y piano, Y oggetto --> X piano,  X oggetto --> Normale piano
            {
                if (sLabelInfo.iViewFace == ESideType::eFront)  // Scaling: Y oggetto --> U piano, Z oggetto --> V piano, X oggetto --> (Sy+Sz)/2 [or fixed]
                {
                    if (sObjInfo.iScaleTypeY == 1 && vBBoxSide.y > 0)
                        mXScale[1][1] = (sPlaneInfo.vUVMax.u - sPlaneInfo.vUVMin.u) / vBBoxSide.y;
                    else if (sObjInfo.iScaleTypeY == 2)
                        mXScale[1][1] = sObjInfo.fFixedScaleY;

                    if (sObjInfo.iScaleTypeZ == 1 && vBBoxSide.z > 0)
                        mXScale[2][2] = (sPlaneInfo.vUVMax.v - sPlaneInfo.vUVMin.v) / vBBoxSide.z;
                    else if (sObjInfo.iScaleTypeZ == 2)
                        mXScale[2][2] = sObjInfo.fFixedScaleZ;

                    if (sObjInfo.iScaleTypeX == 1 && vBBoxSide.x > 0)
                        mXScale[0][0] = (mXScale[1][1] + mXScale[2][2]) / 2;  // Lungo X uso valor medi di scale Y e Z
                    else if (sObjInfo.iScaleTypeX == 2)
                        mXScale[0][0] = sObjInfo.fFixedScaleX;
                }
                else if (sLabelInfo.iViewFace == ESideType::eTop || sLabelInfo.iViewFace == ESideType::eBottom)  // Scaling: X oggetto --> U piano, Y oggetto --> V piano, Z oggetto --> (Sy+Sz)/2 [or fixed]
                {
                    if (sObjInfo.iScaleTypeX == 1 && vBBoxSide.x > 0)
                        mXScale[0][0] = (sPlaneInfo.vUVMax.u - sPlaneInfo.vUVMin.u) / vBBoxSide.x;
                    else if (sObjInfo.iScaleTypeX == 2)
                        mXScale[0][0] = sObjInfo.fFixedScaleX;

                    if (sObjInfo.iScaleTypeY == 1 && vBBoxSide.y > 0)
                        mXScale[1][1] = (sPlaneInfo.vUVMax.v - sPlaneInfo.vUVMin.v) / vBBoxSide.y;
                    else if (sObjInfo.iScaleTypeY == 2)
                        mXScale[1][1] = sObjInfo.fFixedScaleY;

                    if (sObjInfo.iScaleTypeZ == 1 && vBBoxSide.z > 0)
                        mXScale[2][2] = (mXScale[0][0] + mXScale[1][1]) / 2;  // Lungo X uso valor medi di scale Y e Z
                    else if (sObjInfo.iScaleTypeZ == 2)
                        mXScale[2][2] = sObjInfo.fFixedScaleZ;

                }

            }
            else
                ; // TODO

            sRSTInfo.vScale = { (float)mXScale[0][0],(float)mXScale[1][1],(float)mXScale[2][2] };

            ON_Xform mXform;
            ComputeTransform(sObjInfo, sLabelInfo, sPlaneInfo, mXform, mXScale, oPlyIO.GetBBoxSide());

            sRSTInfo.vCenterPos = { (float)mXform[0][3],(float)mXform[1][3],(float)mXform[2][3] };
            sRSTInfo.mRotation = mXform;

            mXform = mXform * mXScale;
            ON_3dPoint v0(0, 0, 0);

            // Apply transformation
            oPlyIO.ApplyXform(mXform, v0);

            cout << "... done!" << endl;

            // If requested compute texture mapping
            if(sObjInfo.bCreateTxt)
            {
                cout << "  ... computing texture mapping ...";

                bool bFace[6];
                memset(bFace, 0, 6 * sizeof(bool));
                for (int j = 0; j < aSegmentLabelInfo.size(); j++)
                {

                    // ------------------ In questo blocco viene costruita la matrice di proiezione associata alla camera ---------------------------------
                    // Get index of camera with same name of current mask image
                    int iRotIndex = 0;
                    int iIndex = oSM.GetCameraIndex((char*)(aSegmentLabelInfo[j].sCameraName.c_str()),&iRotIndex);
                    if (iIndex < 0)
                    {
                        cout << endl << "Cannot find camera named " << aSegmentLabelInfo[j].sCameraName << endl;
                        continue;
                    }

                    // Init label manager for current camera
                    filesystem::path oPath(aSegmentLabelInfo[j].sLabelFileName);
                    int iErr = oLM.Load(oPath, aSegmentLabelInfo[j].iLabelFormat, oCM);
                    if (iErr < 0)
                        continue;

                    // Get camera infos
                    SCameraInfo sCameraInfo = oSM.GetCameraInfo(iIndex, iRotIndex);  // [R, cameraPos] = SM.GetCameraRT(i);
                    // Get label info
                    sLabelInfoTmp = oLM.GetLabelInfoFromInd(aSegmentLabelInfo[j].iLabelInd);

                    if (!bFace[sLabelInfoTmp.iViewFace])  // Current side not yet mapped
                    {
                        // Get sensor infos
                        SSensorInfo sSensorInfo = oSM.GetSensorInfo(sCameraInfo.iSensorID);

                        CreateFullTxtMapping(oPlyIO, sLabelInfoTmp, sCameraInfo, sSensorInfo, oCmdOptions);

                        bFace[sLabelInfoTmp.iViewFace] = true;

                        // ++++++++ TODO TODO TODO ++++++++++++++++++
                        // Dato che il formato .ply supporta 1 sola texture per modello al momento non si può mappare più di una faccia, dato che
                        // lati diversi utilizzerebbero texture (==foto) diverse. L'unicomodo per ovviare alproblema è di costruire una texture ad hoc 
                        // componendo contributi provenienti dalle varie foto.
                        break;
                        // ++++++++++++++++++++++++++++++++++++++++++
                    }


                }
                cout << "... done!" << endl;

            }


            // Salvo modello
            cout << "  ... saving model ...";
            sTmp = filesystem::path(file).replace_extension(".ply").string();
            oPlyIO.WriteFileTransformed(sTmp.c_str());


            sRSTInfo.iCatalogObjID = sObjInfo.iObjID;
            sRSTInfo.sCatalogObjName = sObjInfo.sObjName;
            sRSTInfo.s3DObjFileName = filesystem::path(sTmp).filename().string();
            sRSTInfo.sLabelName = sLabelInfo.sLabelName;
            sRSTInfo.sPicName = filesystem::path(sLabelInfo.sImageFileName).filename().string();
            sRSTInfo.sTxtName = filesystem::path(oPlyIO.GetTextureName()).filename().string();

            oCRSTM.AddEntry(sRSTInfo);


            cout << "... done!" << endl << endl;
        }

    } // Fine ciclo su oggetti segmentati

    if (oCmdOptions.sOutputRSTFileName != "")
        oCRSTM.WriteTable(oCmdOptions.sOutputRSTFileName);

}
#endif




