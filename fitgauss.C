#include "TFile.h"
#include "TH2F.h"
#include "TRandom3.h"
#include <math.h>
#include <iostream>
#include <cstdlib>

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

double mkideal(double energy, double depth, int nbins, bool img){

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
  
  TH2F* event = new TH2F("event", "Event", nbins, -75, 75, nbins, -80, 80);
  TRandom3 r(0);

  //Gen random position wrt pixels

  //cout << RAND_MAX << endl;
  int xmean = 15*((double) rand())/((double) RAND_MAX);
  int ymean = 15*((double) rand())/((double) RAND_MAX);
  //cout << rand() << endl;
  //cout << 15*((double) rand())/((double) RAND_MAX) << endl;
  //cout << 15*((double) rand())/((double) RAND_MAX) << endl;
  
  //Gen distribution
  for(int i = 0; i<ndeposits; i++){

    double depositx = r.Gaus(xmean, sigma);
    double deposity = r.Gaus(ymean, sigma);
    event->Fill(depositx, deposity);
    
  }

  //event->RebinX(10); event->RebinY(10);

  //Calc variance...

  TH1D* px = event->ProjectionX("px", 0, nbins);
  TH1D* py = event->ProjectionY("py", 0, nbins);

  double var = (py->GetRMS()*py->GetRMS());

  if(img){
    event->Draw("COLZ");
    return var;
  } else {
    event->~TH2F();
    return var;
  }
}

double getVarVar(int nvars, double energy, double depth){

  TH1F* vars = new TH1F("vars", "Variances", 100, 0, 808);
  
  for(int i = 0; i < nvars; i++){
    vars->Fill(sqrt(mkideal(energy, depth, 10, false)));

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

void edvarvar(double emin, double emax, int ebins, double dmin, double dmax, int dbins, int ntrials){
  TH2F* edf = new TH2F("edf", "Depth and Energy to VarVar", ebins, emin, emax, dbins, dmin, dmax);

  double einc = (emax - emin) / ebins;
  double dinc = (dmax - dmin) / dbins;

  cout << ebins << dbins << endl;
  
  for(int i=0; i<=ebins; i++){
    for(int j=0; j<=dbins; j++){

      double varvar = getVarVar(ntrials, emin + i * einc, dmin + j * dinc);
      cout << i << j << varvar << endl;
      
      edf->SetBinContent(i, j, sqrt(varvar));
    }
  }

  edf->Draw("COLZ");

}

double getFails(int ntests, double energy, double depth){

  //see how many correct predictions at a certain energy and depth
  int failcount = 0;
  int passcount = 0;

  double realsig = sqrt(215 * -log(abs(1-(1.3*pow(10, -3)) * depth)))/15;
  //cout << "RealSig: " << realsig << endl;
  
  for(int i = 0; i < ntests; i++){
    double predsig = sqrt(mkideal(energy, depth, 10, false))/15;
    //cout << predsig << endl;

    if((realsig < 0.35 and predsig > 0.35) or (realsig > 1.25 and predsig < 1.25) or (realsig > 0.35 and predsig < 0.35) or (realsig < 1.25 and predsig > 1.25)){
      //cout << "fail" << endl;
      failcount++;
    } else {
      //cout << "pass" << endl;
    }
  }

  return ((double) failcount)/ntests;

}

void energypassfail(double emin, double emax, int ebins, double depth){
  TH1F* epassfail = new TH1F("epassfail", "Energy to Fail%", ebins, emin, emax);

  double einc = (emax - emin) / ebins;
  
  for(int i=0; i<=ebins; i++){
    epassfail->AddBinContent(i, getFails(1000, emin + i * einc, depth));

  }

  epassfail->Draw();

}

void edpassfail(double emin, double emax, int ebins, double dmin, double dmax, int dbins, int ntrials){
  TH2F* edf = new TH2F("edf", "Depth and Energy to Fail%", ebins, emin, emax, dbins, dmin, dmax);

  double einc = (emax - emin) / ebins;
  double dinc = (dmax - dmin) / dbins;

  cout << ebins << dbins << endl;
  
  for(int i=0; i<=ebins; i++){
    for(int j=0; j<=dbins; j++){

      double fails = getFails(ntrials, emin + i * einc, dmin + j * dinc);
      cout << i << j << fails << endl;
      if(fails == 0){
	edf->SetBinContent(i, j, 0.00001);
      }else{
	edf->SetBinContent(i, j, fails);
      }
    }
  }

  edf->Draw("COLZ");

}
