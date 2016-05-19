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
#include "TGraphErrors.h"

using namespace std;
using std::cout;
using std::endl;

TString RoundNumber(double num, int decimals, double denom=1.);


TF1 fit_xsec(TString particle="gluino", float minx=600, float maxx=2300){
  TCanvas can;

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

  double legX = 0.55, legY = 0.87;
  double legW = 0.12, legH = 0.074*2;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(0.056); leg.SetFillColor(0); leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(132);

  TString fit_func = "[0]*pow(x, -[1]*exp(x*[2]))";
  TF1 expfit("expfit",fit_func, minx,maxx);
  expfit.SetParameter(0, 5e17);
  expfit.SetParameter(1, 5);
  expfit.SetParameter(2, 6e-5);
  expfit.SetLineWidth(4);

  TF1 polyfit("polyfit","[0]*pow(x, -[1])",minx,maxx);
  polyfit.SetParameter(0, 3e21);
  polyfit.SetParameter(1, 7);
  polyfit.SetLineWidth(4);
  polyfit.SetLineColor(4);

  TGraphErrors graph(mass.size(), &mass[0], &xsec[0], &emass[0], &exsec[0]);
  graph.SetTitle("");
  graph.Draw("ALP");
  graph.Fit(&expfit,"M Q N","",minx,maxx);
  //graph.Fit(&polyfit,"M Q N","",minx,maxx);
  polyfit.Draw("same");
  expfit.Draw("same");
  graph.GetXaxis()->SetTitle(particle+" mass (GeV)");
  graph.GetYaxis()->SetTitle("Cross section (fb)"); 

  TString leglabel = "m_{"+particle+"}^{-"+RoundNumber(polyfit.GetParameter(1),1)+"}";
  leg.AddEntry(&polyfit, leglabel);
  leglabel = "m_{"+particle+"}^{-"+RoundNumber(expfit.GetParameter(1),1)+" e^{" 
    +RoundNumber(expfit.GetParameter(2),6)+"m_{"+particle+"}}}";
  leg.AddEntry(&expfit, leglabel);
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

  return expfit;
}


void gluino_fit(float minx=600, float maxx=2300){
  TF1 expfit = fit_xsec("gluino", minx, maxx);
  vector<float> masses = {600, 700, 1000, 1400, 1800, 2000, 2300};
  vector<float> xsecs = {9.20353, 3.5251, 0.325388, 0.0252977, 0.00276133, 0.000981077, 0.000219049};
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

