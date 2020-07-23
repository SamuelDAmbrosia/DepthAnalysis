#include "TFile.h"
#include "TH2F.h"
#include "TRandom3.h"
#include <math.h>
#include <iostream>

void to1d(){

  TFile *file = new TFile("img_0.root", "READ");
  TH2F *img_raw = (TH2F*)file->Get("image_raw");
  //img_raw->Draw();

  //std::cout << img_raw->GetXaxis()->GetNbins();

  TH1F *img_1d = new TH1F("img_1d", "img_1d", img_raw->GetXaxis()->GetNbins(),0,img_raw->GetXaxis()->GetNbins());

  for (int i = 0; i < img_raw->GetXaxis()->GetNbins(); i++){
    int summer = 0;
    for (int j = 0; j < img_raw->GetYaxis()->GetNbins(); j++){
      summer = summer + img_raw->GetBin(i,j);
    }

    //std::cout << summer;
    img_1d->SetBinContent(i, summer);
  }

  img_1d->Draw();
  //img_raw->RebinX(5); img_raw->RebinY(5);
  //img_raw->Draw("COLZ");

}

double mkideal(double energy, double depth, bool img){

  //Energy (keV), depth (0-675 micron)

  //Get deposits from energy
  int ndeposits = (int) ((energy * 1000)/(3.77));
  
  //Get sigma from depth
  if(depth < 0){
    cout << "Outside CCD! Too low." << endl;
    return 0;
  } if(depth > 675){
    cout << "Outside CCD! Too high." << endl;
    return 0;
  }
  
  
  double sigma = sqrt(215 * -log(abs(1-(1.3*pow(10, -3)) * depth)));
  
  TH2F* event = new TH2F("event", "Event", 100, -80, 80, 100, -80, 80);
  TRandom3 r(0);

  //Gen distribution
  for(int i = 0; i<ndeposits; i++){

    double depositx = r.Gaus(0, sigma);
    double deposity = r.Gaus(0, sigma);
    event->Fill(depositx, deposity);
    
  }


  //Calc variance... mean = 0

  TH1D* px = event->ProjectionX("px", 0, 100);
  TH1D* py = event->ProjectionY("py", 0, 100);

  //cout << (py->GetRMS()*py->GetRMS()) << endl;

  double var = (py->GetRMS()*py->GetRMS());

  if(img){
    event->Draw();
    return var;
  } else {
    event->~TH2F();
    return var;
  }
}

double getVarVar(int nvars, double energy, double depth){

  TH1F* vars = new TH1F("vars", "Variances", 100, 0, 808);
  
  for(int i = 0; i < nvars; i++){
    vars->Fill(mkideal(energy, depth, false));

  }

  double varvar = (vars->GetRMS() * vars->GetRMS());

  vars->~TH1F();

  return varvar;

}

double oneovern(double energy, double scalar){

  return (scalar/energy);

}

void energyVarVar(double emin, double emax, int ebins, double depth){
  TH1F* evarvars = new TH1F("evarvars", "Energy to VarVars", ebins, emin, emax);

  double einc = (emax - emin) / ebins;
  
  for(int i=0; i<ebins; i++){
    evarvars->AddBinContent(i, getVarVar(100, emin + i * einc, depth));

  }

  evarvars->Draw();

}
