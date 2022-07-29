#include "HoughTransform2D.hpp"

#include <iostream>
#include <vector>
using namespace std;


//霍夫空间下参数方程 r=x*cos(theta)+y*sin(theta)
Double_t HoughSpaceFunc(Double_t *x,Double_t *par){
    return par[0]*TMath::Cos(x[0]*TMath::DegToRad())+par[1]*TMath::Sin(x[0]*TMath::DegToRad());
}

bool judge(const pair<pair<double,double>,int> a, pair<pair<double,double>,int> b) {
    return a.second>b.second;
}

const Int_t HoughTransform2D::maxlength = floor(sqrt(300*300+150*150));

HoughTransform2D::HoughTransform2D()
{
}

HoughTransform2D::~HoughTransform2D()
{
}

/*
vector<TVector3> Point 笛卡尔坐标系下3D坐标点
std::vector<std::pair<double theta,double rho>> HoughPara  霍夫空间下参数半径r以及极角theta
Bool_t XYPlane  flag 2维平面选择XZ  电场方向值作为Z坐标，束流方向值作为X坐标（原点选在方形左下角）  电场方向值作为Z坐标，束流方向值作为Y坐标（原点选在方形左上角）
Bool_t YZPlane  
Bool_t XZPlane
Double_t bin_theta    theta [0，180）划分地bin大小,默认值为 1°
Double_t bin_r        r划分地bin大小,默认值为1
*/
void HoughTransform2D::HoughSpace(vector<TVector3> Point,vector<pair<double ,double>>& HoughPara,TH2D *HoughHist,
                Bool_t XYPlane,Bool_t YZPlane,Bool_t XZPlane,Double_t bin_theta,Double_t bin_r,Int_t threshold){
    HoughPara.clear();
    HoughHist->Reset();
    HoughHist->SetStats(0);
    
    int numangle, numrho;

    //由角度和距离的分辨率得到角度和距离的数量，即霍夫变换后角度和距离的个数  
    numangle = round(180. / bin_theta);  
    numrho = round(2*(maxlength/bin_r));
    
    //为累加器数组分配内存空间  
    //该累加器数组其实就是霍夫空间，它是用一维数组表示二维空间 
    vector<int> accumulator(((numangle+2) * (numrho+2)),0);//多分配一行一列，主要是方便 stage2 中4邻域的比较，否则比较时会溢出
    //为排序数组分配内存空间 
    vector<int> sort_theta(numangle * numrho);
    vector<int> sort_rho(numangle * numrho);

    // stage 1. fill accumulator  
    //执行步骤1，逐点进行霍夫空间变换，并把结果放入累加器数组内  
    for(int i=0;i<Point.size();i++){  
        Double_t par[2];
        if(XYPlane){par[0]=Point[i][0];par[1]=Point[i][1];}
        else if(YZPlane){par[0]=Point[i][1];par[1]=Point[i][2];}
        else if(XZPlane){par[0]=Point[i][0];par[1]=Point[i][2];}
        else {cout<<"Wrong flag......EXit!"<<endl;exit(-1);}
        
        for(double deg=0;deg<180.0;deg+=bin_theta){
            double rho = HoughSpaceFunc(&deg,par);

            if(abs(rho)>maxlength)continue;
            HoughHist->Fill(deg,rho);
            // cout<<deg<<"\t"<<rho<<endl;
            rho += maxlength;
            // cout<<(numangle+2) * (numrho+2)<<"\t"<<(Int_t(rho/bin_r)+1)*(numangle+2)+(Int_t(deg/bin_theta)+1)<<endl;
            accumulator[(Int_t(rho/bin_r)+1)*(numangle+2)+(Int_t(deg/bin_theta)+1)]++;//Int_t(rho/bin_r)+1是为了第一行空出来
                                                                                     //numangle+2 是总共的列数
                                                                                     //Int_t(deg/bin_theta)+1把第一列空出来，stage 2需要比较4邻域累加器中值的大小
        }
        // cout<<i<<endl;
    }
    
    // stage 2. find local maximums  
    //执行步骤2，找到局部极大值，即非极大值抑制  
    pair<double,double> LinePolar;//first为极角参数，second为rho参数
    pair<pair<double,double>,int> sort_buf;
    vector<pair<pair<double,double>,int>> Sort;
    
    for(int r = 0; r < numrho; r++ )  
        for(int n = 0; n < numangle; n++ )  
        {  
            //得到当前值在累加器数组的位置  
            int base = (r+1) * (numangle+2) + n+1;   
            if( accumulator[base] > threshold &&    //必须大于所设置的阈值  
                //在10邻域内进行非极大值抑制  
                accumulator[base] > accumulator[base - 1] && accumulator[base] > accumulator[base + 1] &&  
                accumulator[base] > accumulator[base - 2] && accumulator[base] > accumulator[base + 2] &&
                accumulator[base] > accumulator[base - numangle] && accumulator[base] > accumulator[base + numangle] &&
                accumulator[base] > accumulator[base - numangle - 1] && accumulator[base] > accumulator[base + numangle + 1] &&
                accumulator[base] > accumulator[base - numangle - 2] && accumulator[base] > accumulator[base + numangle + 2]
                ) {
                    LinePolar.first = n*bin_theta;
                    LinePolar.second = (r-numrho/2.-1)*bin_r;
                    // cout<<threshold<<"  "<<"bin_r:"<<bin_r<<"   (r-r/2):"<<(r-r/2)<<"   "<<(r-r/2)*bin_r<<"  bin_theta:"<<bin_theta<<"  (n*bin_theta):"<<(n*bin_theta)<<" "<<accumulator[base]<<endl;
                    sort_buf.first=LinePolar;
                    sort_buf.second=accumulator[base];
                    Sort.push_back(sort_buf);
                    } 
        }  
    
    //对Houghpeak按计数从大到小排序
    sort(Sort.begin(),Sort.end(),judge);
    
    for(auto i=0;i<Sort.size();i++){
        HoughPara.push_back(Sort[i].first);
    }

    for (vector<std::pair<double,double>>::iterator it = HoughPara.begin(); it != HoughPara.end();) {
        if(it == HoughPara.end()-1)break;
        if (((*it).first <= (*(it+1)).first+5 && ((*it).first >= (*(it+1)).first-5)) && ((*it).second <= (*(it+1)).second+5 && ((*it).second >= (*(it+1)).second-5))) {
            it = HoughPara.erase(it+1);
        } else {
            ++it;
        }
    }

}

void HoughTransform2D::HoughSpace(vector<TVector3> Point,vector<pair<double ,double>>& HoughPara,
                Bool_t XYPlane,Bool_t YZPlane,Bool_t XZPlane,Double_t bin_theta,Double_t bin_r,Int_t threshold){
    HoughPara.clear();
    
    int numangle, numrho;

    //由角度和距离的分辨率得到角度和距离的数量，即霍夫变换后角度和距离的个数  
    numangle = round(180. / bin_theta);  
    numrho = round(2*(maxlength/bin_r));
    
    //为累加器数组分配内存空间  
    //该累加器数组其实就是霍夫空间，它是用一维数组表示二维空间 
    vector<int> accumulator(((numangle+2) * (numrho+2)),0);//多分配一行一列，主要是方便 stage2 中4邻域的比较，否则比较时会溢出
    //为排序数组分配内存空间 
    vector<int> sort_theta(numangle * numrho);
    vector<int> sort_rho(numangle * numrho);

    // stage 1. fill accumulator  
    //执行步骤1，逐点进行霍夫空间变换，并把结果放入累加器数组内  
    for(int i=0;i<Point.size();i++){  
        Double_t par[2];
        if(XYPlane){par[0]=Point[i][0];par[1]=Point[i][1];}
        else if(YZPlane){par[0]=Point[i][1];par[1]=Point[i][2];}
        else if(XZPlane){par[0]=Point[i][0];par[1]=Point[i][2];}
        else {cout<<"Wrong flag......EXit!"<<endl;exit(-1);}
        
        for(double deg=0;deg<180.0;deg+=bin_theta){
            double rho = HoughSpaceFunc(&deg,par);

            // cout<<deg<<"\t"<<rho<<endl;
            rho += maxlength;
            // cout<<(numangle+2) * (numrho+2)<<"\t"<<(Int_t(rho/bin_r)+1)*(numangle+2)+(Int_t(deg/bin_theta)+1)<<endl;
            accumulator[(Int_t(rho/bin_r)+1)*(numangle+2)+(Int_t(deg/bin_theta)+1)]++;//Int_t(rho/bin_r)+1是为了第一行空出来
                                                                                     //numangle+2 是总共的列数
                                                                                     //Int_t(deg/bin_theta)+1把第一列空出来，stage 2需要比较4邻域累加器中值的大小
        }
    }
    
    // stage 2. find local maximums  
    //执行步骤2，找到局部极大值，即非极大值抑制 
    pair<double,double> LinePolar;//first为极角参数，second为rho参数
    pair<pair<double,double>,int> sort_buf;
    vector<pair<pair<double,double>,int>> Sort;
    
    for(int r = 0; r < numrho; r++ )  
        for(int n = 0; n < numangle; n++ )  
        {  
            //得到当前值在累加器数组的位置  
            int base = (r+1) * (numangle+2) + n+1;   
            if( accumulator[base] > threshold &&    //必须大于所设置的阈值  
                //在10邻域内进行非极大值抑制  
                accumulator[base] > accumulator[base - 1] && accumulator[base] > accumulator[base + 1] &&  
                accumulator[base] > accumulator[base - 2] && accumulator[base] > accumulator[base + 2] &&
                accumulator[base] > accumulator[base - numangle] && accumulator[base] > accumulator[base + numangle] &&
                accumulator[base] > accumulator[base - numangle - 1] && accumulator[base] > accumulator[base + numangle + 1] &&
                accumulator[base] > accumulator[base - numangle - 2] && accumulator[base] > accumulator[base + numangle + 2]
                ) {
                    LinePolar.first = n*bin_theta;
                    LinePolar.second = (r-numrho/2.-1)*bin_r;
                    // cout<<threshold<<"  "<<"bin_r:"<<bin_r<<"   (r-r/2):"<<(r-r/2)<<"   "<<(r-r/2)*bin_r<<"  bin_theta:"<<bin_theta<<"  (n*bin_theta):"<<(n*bin_theta)<<" "<<accumulator[base]<<endl;
                    sort_buf.first=LinePolar;
                    sort_buf.second=accumulator[base];
                    Sort.push_back(sort_buf);
                    } 
        } 
    
    //对Houghpeak按计数从大到小排序
    sort(Sort.begin(),Sort.end(),judge);
    
    for(auto i=0;i<Sort.size();i++){
        HoughPara.push_back(Sort[i].first);
    } 

    for (vector<std::pair<double,double>>::iterator it = HoughPara.begin(); it != HoughPara.end();) {
        if(it == HoughPara.end()-1)break;
        if (((*it).first <= (*(it+1)).first+5 && ((*it).first >= (*(it+1)).first-5)) && ((*it).second <= (*(it+1)).second+5 && ((*it).second >= (*(it+1)).second-5))) {
            it = HoughPara.erase(it+1);
        } else {
            ++it;
        }
    }
}