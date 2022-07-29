#include <TH2.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TMath.h>
#include <TF1.h>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TStopwatch.h>

#include "include/HoughTransform2D.hpp"

#include <iostream>
#include <vector>
using namespace std;

#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"

//轨迹空间下方程 y=k*x+b
double LineFunc(Double_t *x,Double_t *par){
    double bb = par[1]/sin((par[0])*TMath::DegToRad());
    double kk = - 1/TMath::Tan(par[0]*TMath::DegToRad());

    return kk*x[0]+bb;
}

int main(int argc,char **argv){
    // gStyle->SetOptStat(0);

    Double_t bin_theta;  //大于10度 参数 0.22 1.0，小于10° 参数0.1，2.0
    Double_t bin_r;
    cin>>bin_theta;
    cin>>bin_r;
    cout<<"bin_theta: "<<bin_theta<<"   bin_r :"<<bin_r<<endl;

    TApplication app("app", &argc, argv);
    //霍夫变换的初始化声明
    HoughTransform2D *Rec2D = new HoughTransform2D();
    TH2D *HoughHist =new TH2D("Hough Space in 2D","Hough Space in 2D",round(180./bin_theta),0,180,2*round(HoughTransform2D::maxlength/bin_r),-HoughTransform2D::maxlength,HoughTransform2D::maxlength);
    HoughHist->GetXaxis()->SetNameTitle("#theta(deg)","#theta(deg)");
    HoughHist->GetYaxis()->SetNameTitle("r(mm)","r(mm)");
    vector<pair<double ,double>> HoughXYPara_dig;//霍夫空间变换参数，第一个为角度theta，第二个为距离r
    vector<pair<double ,double>> HoughXZPara_dig;
    vector<pair<double ,double>> HoughYZPara_dig;

    //假设存在一条3D直线(x-4)/2 =(y+1)/1 =(z-3)/5
    //产生拟合数据
    TVector3 pos;
    vector<TVector3> data;
    for (int i = 0; i < 100; i++)
    {
        //x在-10至10之间随机产生数据
        pos[0] = -10 +(10-(-10))*gRandom->Uniform();
        pos[1] = ((pos[0]-4)/2)*1-1;
        pos[2] = ((pos[0]-4)/2)*5+3;

        //smeared data
        pos[0] = gRandom->Gaus(pos[0],0.5);
        pos[1] = gRandom->Gaus(pos[1],0.5);
        pos[2] = gRandom->Gaus(pos[2],0.5);

        data.push_back(pos);
    }

    //Do the fit  划分大bin做粗略拟合
    Rec2D->HoughSpace(data,HoughXYPara_dig,HoughHist,true,false,false,bin_theta,bin_r,10);  //XY  bin_theta=.2 bin_r=1.划分大bin做粗略拟合
    Rec2D->HoughSpace(data,HoughYZPara_dig,HoughHist,false,true,false,bin_theta,bin_r,10);
    Rec2D->HoughSpace(data,HoughXZPara_dig,HoughHist,false,false,true,bin_theta,bin_r,10);
    cout<<"霍夫变换后XY: theta "<<HoughXYPara_dig[0].first<<"   r "<<HoughXYPara_dig[0].second<<endl;
    cout<<"霍夫变换后YZ: theta "<<HoughYZPara_dig[0].first<<"   r "<<HoughYZPara_dig[0].second<<endl;
    cout<<"霍夫变换后XZ: theta "<<HoughXZPara_dig[0].first<<"   r "<<HoughXZPara_dig[0].second<<endl;

    TGraph *XY = new TGraph(100);
    TGraph *YZ = new TGraph(100);
    TGraph *XZ = new TGraph(100);
    for (int i = 0; i < 100; i++)
    {
        XY->SetPoint(i,data[i][0],data[i][1]);
        YZ->SetPoint(i,data[i][1],data[i][2]);
        XZ->SetPoint(i,data[i][0],data[i][2]);
    }

    TF1 *f= new TF1("line in track space",LineFunc,-10,10,2);
    f->SetParameters(HoughXYPara_dig[0].first,HoughXYPara_dig[0].second);
    TCanvas *XYCan = new TCanvas("XYCan","XYCan");
    XY->SetNameTitle("XY plane","XY plane");
    XY->Draw("AP*");
    XY->Fit(f,"RQ+");
    HoughXYPara_dig[0].first = f->GetParameter(0);
    HoughXYPara_dig[0].second = f->GetParameter(1);
    cout<<"拟合后XY:  theta "<<HoughXYPara_dig[0].first<<"   r "<<HoughXYPara_dig[0].second<<endl;
    cout<<"拟合后XY:  k "<<- 1/TMath::Tan(HoughXYPara_dig[0].first*TMath::DegToRad())<<"   b "<<HoughXYPara_dig[0].second/sin((HoughXYPara_dig[0].first)*TMath::DegToRad())<<endl;

    TF1 *f1= new TF1("line in track space",LineFunc,-10,10,2);
    f1->SetParameters(HoughYZPara_dig[0].first,HoughYZPara_dig[0].second);
    TCanvas *YZCan = new TCanvas("YZCan","YZCan");
    YZ->SetNameTitle("YZ plane","YZ plane");
    YZ->Draw("AP*");
    YZ->Fit(f1,"RQ+");
    HoughYZPara_dig[0].first = f1->GetParameter(0);
    HoughYZPara_dig[0].second = f1->GetParameter(1);
    cout<<"拟合后YZ:  theta "<<HoughYZPara_dig[0].first<<"   r "<<HoughYZPara_dig[0].second<<endl;
    cout<<"拟合后YZ:  k "<<- 1/TMath::Tan(HoughYZPara_dig[0].first*TMath::DegToRad())<<"   b "<<HoughYZPara_dig[0].second/sin((HoughYZPara_dig[0].first)*TMath::DegToRad())<<endl;

    TF1 *f2= new TF1("line in track space",LineFunc,-10,10,2);
    f2->SetParameters(HoughXZPara_dig[0].first,HoughXZPara_dig[0].second);
    TCanvas *XZCan = new TCanvas("XZCan","XZCan");
    XZ->SetNameTitle("XZ plane","XZ plane");
    XZ->Draw("AP*");
    XZ->Fit(f2,"RQ+");
    HoughXZPara_dig[0].first = f2->GetParameter(0);
    HoughXZPara_dig[0].second = f2->GetParameter(1);
    cout<<"拟合后XZ:  theta "<<HoughXZPara_dig[0].first<<"   r "<<HoughXZPara_dig[0].second<<endl;
    cout<<"拟合后XZ:  k "<<- 1/TMath::Tan(HoughXZPara_dig[0].first*TMath::DegToRad())<<"   b "<<HoughXZPara_dig[0].second/sin((HoughXZPara_dig[0].first)*TMath::DegToRad())<<endl;
    
    delete Rec2D;

    app.Run();
    return 0;
}
// 从刻度杆入射 alpha源刻度  长面