#ifndef _2DHoughTransform2D_
#define _2DHoughTransform2D_

// #include <TF1.h>
#include <Rtypes.h>
#include <TVector3.h>
#include <TH2.h>
#include "TMath.h"

#include <iostream>
#include <vector>
using namespace std;

class HoughTransform2D
{
public:
    static const Int_t maxlength;
    
    HoughTransform2D();
    ~HoughTransform2D();

    /*
    vector<TVector3> Point 笛卡尔坐标系下3D坐标点
    std::vector<std::pair<double theta,double rho>> HoughPara  霍夫空间下参数半径r以及极角theta
    Bool_t XYPlane  flag 2维平面选择XZ  电场方向值作为Z坐标，束流方向值作为X坐标（原点选在方形左下角）  电场方向值作为Z坐标，束流方向值作为Y坐标（原点选在方形左上角）
    Bool_t YZPlane  
    Bool_t XZPlane
    */
    // void HoughSpace(std::vector<TVector3> Point,std::vector<std::pair<Double_t ,Double_t>>& HoughPara,TH2D *HoughHist,
    //                 Bool_t XYPlane,Bool_t YZPlane,Bool_t XZPlane,Int_t NumOfLines,Double_t bin_theta=1,Double_t bin_r=1);

    // Int_t FindPeaks(TH2* hist,std::vector<std::pair<Double_t ,Double_t>>& HoughPara,Int_t RebinNum,Int_t npeaks=3);

    void HoughSpace(std::vector<TVector3> Point,std::vector<std::pair<Double_t ,Double_t>>& HoughPara,TH2D *HoughHist,
                    Bool_t XYPlane,Bool_t YZPlane,Bool_t XZPlane,Double_t bin_theta=1.0,Double_t bin_r=1.,Int_t threshold=3);
    void HoughSpace(vector<TVector3> Point,vector<pair<double ,double>>& HoughPara,
                Bool_t XYPlane,Bool_t YZPlane,Bool_t XZPlane,Double_t bin_theta=1.,Double_t bin_r=1.,Int_t threshold=3);

ClassDef(HoughTransform2D,0);
};


#endif