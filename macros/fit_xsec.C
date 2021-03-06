// fit_xsec: Fits the gluino and squark cross-sections. It does it with a polynomial m^-n to see
//           roughly the power of the mass with which it decreases, but the good fit includes
//           the exponential (expfit)

#include <fstream>
#include <iostream>
#include <vector>
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TF1.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TGraphErrors.h"

using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double num, int decimals, double denom=1.);
void setCanvas(TCanvas &can, float lMargin, float tMargin, float rMargin, float bMargin);

TF1 fit_xsec(TString particle="gluino", float minx=600, float maxx=2300){
  //// Creating canvas
  gStyle->SetOptStat(0);  
  float lMargin(0.14), tMargin(0.08), rMargin(0.04), bMargin(0.14);
  TCanvas can("canvas","", 700, 600);
  setCanvas(can, lMargin, tMargin, rMargin, bMargin);
  float Xmin(700), Xmax(1750), Ymin(0), Ymax(1800);

  ifstream filexsec("txt/"+particle+"_xsec.txt");
  vector<float> mass, emass, xsec, exsec;
  TString dummy, smass, sxsec, sexsec;
  for(int ind(0); ind < 9; ind++) filexsec >> dummy;
  while(filexsec){
    filexsec >> dummy >> smass >> dummy >> dummy >> sxsec >> dummy >> sexsec >> dummy;
    if(smass.Atof()<=0) break;
    mass.push_back(smass.Atof());
    emass.push_back(0);
    xsec.push_back(sxsec.Atof()*1000);
    exsec.push_back(sexsec.Atof()*xsec[xsec.size()-1]/100.);
  }

  TString fit_func = "[0]*pow(x, -[1]*exp(x*[2]))";
  TF1 expfit("expfit",fit_func, minx,maxx);
  expfit.SetParameter(0, 5e17);
  expfit.SetParameter(1, 5);
  expfit.SetParameter(2, 6e-5);
  expfit.SetLineWidth(4);


  TGraphErrors graph(mass.size(), &mass[0], &xsec[0], &emass[0], &exsec[0]);
  TGraphErrors graph2(mass.size(), &mass[0], &xsec[0], &emass[0], &emass[0]);
  graph.SetTitle("Fit from "+RoundNumber(minx,0)+" to "+RoundNumber(maxx,0));
  graph.SetLineColor(1); graph.SetFillColor(kOrange); graph.SetLineWidth(1); graph.SetLineStyle(2);
  graph.Draw("a3");
  graph2.SetLineColor(1); graph2.SetLineWidth(1); graph2.SetLineStyle(2);
  graph2.Draw("L same");
  graph.Fit(&expfit,"M Q N","",minx,maxx);

  //// Important to find good initial values due to the correlation between par0 and par1
  TF1 polyfit("polyfit","[0]*pow(x, -[1])",minx,maxx);
  float expon = (log(expfit.Eval(maxx)) - log(expfit.Eval(minx)))/(log(minx)-log(maxx));
  float sigma1 = expfit.Eval(maxx)/pow(maxx, -expon);

  polyfit.SetParameter(0, sigma1);
  polyfit.SetParameter(1, expon);
  polyfit.SetLineWidth(4);
  polyfit.SetLineColor(4);

  graph.Fit(&polyfit,"M Q N","",minx,maxx);
  expfit.Draw("same");
  polyfit.Draw("same");

  graph.GetXaxis()->SetLabelFont(42);
  graph.GetXaxis()->SetLabelSize(0.035);
  graph.GetXaxis()->SetTitleFont(42);
  graph.GetXaxis()->SetTitleSize(0.05);
  graph.GetXaxis()->SetTitleOffset(1.2);
  graph.GetXaxis()->SetLabelOffset(0.009);
  graph.GetXaxis()->SetTitle(particle+" mass [GeV]");
  graph.GetYaxis()->SetLabelFont(42);
  graph.GetYaxis()->SetLabelSize(0.035);
  graph.GetYaxis()->SetTitleFont(42);
  graph.GetYaxis()->SetTitleSize(0.05);
  graph.GetYaxis()->SetTitleOffset(1.35);
  graph.GetYaxis()->SetTitle(particle+" pair-production cross section [fb]");

  double legX = 0.37, legY = 0.87;
  double legW = 0.12, legH = 0.12*2;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.056); leg.SetFillColor(0); leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(132);

  TString leglabel = "1-parameter fit: m_{"+particle+"}^{-"+RoundNumber(polyfit.GetParameter(1),1)+"}";
  leg.AddEntry(&polyfit, leglabel,"l");
  leglabel = "2-parameter fit: m_{"+particle+"}^{-"+RoundNumber(expfit.GetParameter(1),1)+"#timese^{" 
    +RoundNumber(expfit.GetParameter(2),6)+"#timesm_{"+particle+"}}}";
  leg.AddEntry(&expfit, leglabel,"l");
  leg.Draw();

  //// Printing fit to xsec
  for(int par=0; par<expfit.GetNpar(); par++){
    TString val=""; val+=expfit.GetParameter(par);
    TString par_s = "["; par_s += par; par_s += "]";
    fit_func.ReplaceAll(par_s, val);
    fit_func.ReplaceAll("x", "mass");
    fit_func.ReplaceAll("emassp", "exp");
  }
  cout<<endl<<"Fit to "<<particle<<" cross-section is "<<fit_func<<endl<<endl;

  if(false){
  } // if debug

  can.SetLogy(1);
  can.SaveAs(particle+"_xsec_"+RoundNumber(minx,0)+"_"+RoundNumber(maxx,0)+".pdf");

  return polyfit;
}


void gluino_fit(float minx=600, float maxx=2300){
  TF1 expfit = fit_xsec("gluino", minx, maxx);
  vector<float> masses = {600, 700, 1000, 1400, 1800, 2000, 2100, 2300};
  vector<float> xsecs = {9.20353, 3.5251, 0.325388, 0.0252977, 0.00276133, 0.000981077, 0.000591918, 0.000219049};
  cout<<endl;
  for(size_t im=0; im<masses.size(); im++){
    float xs = expfit.Eval(masses[im])/1000;
    cout<<"Fitted xsec for mass "<<setw(4)<<masses[im]<<" is "<<setw(8)<<RoundNumber(xs*1000,3)<<" fb, off by "
	<<setw(5)<<RoundNumber((xs-xsecs[im])/xsecs[im]*100,2)<<" %"<<endl;
  }
  cout<<endl;
}

void stop_fit(float minx=100, float maxx=300){
  TF1 expfit = fit_xsec("stop", minx, maxx);
  vector<float> masses = {100, 200, 300, 500, 700, 1000, 1200};
  vector<float> xsecs = {1521.11, 64.5085, 8.51615, 0.51848, 0.0670476, 0.00615134, 0.00159844 };
  cout<<endl;
  for(size_t im=0; im<masses.size(); im++){
    float xs = expfit.Eval(masses[im])/1000;
    cout<<"Fitted xsec for mass "<<setw(4)<<masses[im]<<" is "<<setw(11)<<RoundNumber(xs*1000,3)<<" fb, off by "
	<<setw(5)<<RoundNumber((xs-xsecs[im])/xsecs[im]*100,2)<<" %"<<endl;
  }
  cout<<endl;

  expfit = fit_xsec("stop", 300, 1200);
  cout<<endl;
  for(size_t im=0; im<masses.size(); im++){
    float xs = expfit.Eval(masses[im])/1000;
    cout<<"Fitted xsec for mass "<<setw(4)<<masses[im]<<" is "<<setw(11)<<RoundNumber(xs*1000,3)<<" fb, off by "
	<<setw(5)<<RoundNumber((xs-xsecs[im])/xsecs[im]*100,2)<<" %"<<endl;
  }
  cout<<endl;
}

void setCanvas(TCanvas &can, float lMargin, float tMargin, float rMargin, float bMargin){
  can.SetTickx(1);
  can.SetTicky(1);
  can.SetLeftMargin(lMargin);
  can.SetTopMargin(tMargin);
  can.SetRightMargin(rMargin);
  can.SetBottomMargin(bMargin);
}

TString RoundNumber(double num, int decimals, double denom){
  if(denom==0) return " - ";
  double neg = 1; if(num*denom<0) neg = -1;
  num /= neg*denom; num += 0.5*pow(10.,-decimals);
  long num_int = static_cast<long>(num);
  long num_dec = static_cast<long>((1+num-num_int)*pow(10.,decimals));
  TString s_dec = ""; s_dec += num_dec; s_dec.Remove(0,1);
  TString result="";
  if(neg<0) result+="-";
  result+= num_int;
  if(decimals>0) {
    result+="."; result+=s_dec;
  }

  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<decimals-afterdot.Length(); i++)
    result += "0";
  return result;
}

