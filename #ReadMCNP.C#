#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TFITS.h"
#include "TTree.h"
#include "TROOT.h"
#include "TParameter.h"
#include "TList.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TVectorD.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <vector>

#include "ConfigValue.C"

TH2F* ResizeH2(TH2F* h, Int_t nbins_x, Int_t nbins_y, Int_t xlow, Int_t ylow){
    
    TString name = h->GetName();
    h->SetNameTitle("h","h");
    TH2F* t = new TH2F(name, name, nbins_x, 0, nbins_x, nbins_y, 0, nbins_y);
    Int_t nx = h->GetNbinsX();
    Int_t ny = h->GetNbinsY();
    
    for(Int_t i=1; i<=nx; i++) for(Int_t j=1; j<=ny; j++) if(i<=nbins_x && j<=nbins_y)
        t->SetBinContent(i+xlow-1,j+ylow-1,h->GetBinContent(i,j));
 
    h->Delete();
    return t;
}

void ReadMCNP( TString mcnp_fname, TString out_name, Int_t runid=0, TString cfgfile="", TString blank_name="" ){
    
    LoadConfigFile(cfgfile);
    
    //SIMULATION PARAMETERS//
    Double_t e_factor = ConfigValue("sim_e_conv_factor");
    Int_t bit_pix = ConfigValue("sim_bit_pix");
    Double_t pedestal_sigma = ConfigValue("sim_pedestal_sigma");
    Double_t pedestal_mean = ConfigValue("sim_pedestal_mean");
    Int_t nreads = ConfigValue("sim_nreads");
    Double_t diff_A = ConfigValue("sim_diff_A");
    Double_t diff_b = ConfigValue("sim_diff_b");
    Double_t diff_offset = ConfigValue("sim_diff_offset");
    Double_t fano_factor = ConfigValue("sim_fano_factor");
    Bool_t sim_cti_fast = ConfigValue("sim_cti_fast");
    Double_t sim_cti_1 = ConfigValue("sim_cti_1");
    Double_t sim_cti_2 = ConfigValue("sim_cti_2");
    Double_t sim_cti_x = ConfigValue("sim_cti_x");
    Int_t sim_cti_bound = ConfigValue("sim_cti_bound");
    Bool_t sim_flip_CCD = ConfigValue("sim_flip_CCD");
    
    //Added by KR
    Double_t sim_inject_rate = ConfigValue("sim_inject_rate"); // Dark current rate e-/pix/img
    Double_t diff_alpha = ConfigValue("sim_diff_alpha"); // Slope of sigma_max um/keV
    
    //======================================//

    //file to output data
    TFile* fin = NULL;
    
    Int_t hwid  = 0;
    
    if(mcnp_fname.Contains(".mcnp")){
    
        hwid=-1;
    }
    else if(mcnp_fname.Contains(".pds")){
        
        hwid=-2;
    }
    else{
        
        std::cout << "Error: File is neither .mcnp nor .pds!" << std::endl;
        return;
    }
    
    //some general things
    //First add the most general information
    FILE* sysout;
    char ct[80];
    int lastchar;
    sysout = popen("pwd","r");
    lastchar = fread(ct, 1, 80, sysout);
    ct[lastchar-1] = '\0';
    TString pwd(ct);
    pclose(sysout);
    sysout = popen("hostname","r");
    lastchar = fread(ct, 1, 80, sysout);
    ct[lastchar-1] = '\0';
    TString hostname(ct);
    pclose(sysout);
    sysout = popen("svnversion","r");
    lastchar = fread(ct, 1, 80, sysout);
    ct[lastchar-1] = '\0';
    TString svnversion(ct);
    pclose(sysout);
    
    //some general values
    Double_t e_ion = 3.77E-3; //Mean energy to ionize e in Si
    Double_t adc_unit = e_ion/e_factor;

    Int_t nps, ntal, nx, ny;
    Double_t bw_x, bw_y, bw_z;
    Int_t simid;
    
    Bool_t leakexists = false;
    
    TH2F* edep = NULL;
    TH2F* depth = NULL;
    TH2F* meanx = NULL;
    TH2F* meany = NULL;
    TH2F* eleak = NULL; // Added by KR
    
    //if .mcnp read from the mesh file
    
    if(hwid==-1) {
        
        //object from .fits
        ifstream mcnp_file(mcnp_fname);
        
        //This is for reading TMESH MCNPX files
        TString dummy="";
        
        //only implement readout of 1 CCD, ID 0
        simid=0;

        //At some point may want to store more of these header values
        mcnp_file>>dummy>>dummy>>dummy>>dummy>>dummy>>nps;
        
        while(dummy!="ntal") mcnp_file>>dummy;
        
        mcnp_file>>ntal;
        
        if(ntal<1 || ntal>4){
            
            std::cout << "Error: " << ntal << " is an incorrect number of tallies. Must be between 1 and 4!" << std::endl;
            return;
        }
        
        while(dummy!="f") mcnp_file>>dummy;
        
        mcnp_file>>dummy>>dummy>>nx>>ny>>dummy;
        
        Double_t x1, x2;
        mcnp_file>>x1>>x2;
        for(Int_t k=0; k<=nx-2; k++) mcnp_file>>dummy;
        bw_x = x2 - x1;
        
        Double_t y1, y2;
        mcnp_file>>y1>>y2;
        for(Int_t k=0; k<=ny-2; k++) mcnp_file>>dummy;
        bw_y = y2 - y1;
        
        Double_t z1, z2;
        mcnp_file>>z1>>z2;
        bw_z = z2 - z1;
        
        //read the values now
        while(dummy!="vals") mcnp_file>>dummy;
        
        edep = new TH2F("edep", "edep", nx, 0, nx, ny, 0, ny);
        depth = new TH2F("depth", "depth", nx, 0, nx, ny, 0, ny);
        meanx = new TH2F("meanx", "meanx", nx, 0, nx, ny, 0, ny);
        meany = new TH2F("meany", "meany", nx, 0, nx, ny, 0, ny);
        eleak = new TH2F("eleak","eleak",nx,0,nx,ny,0,ny);  // Added by KR
        
        Double_t val;
        for(Int_t j=1; j<=ny; j++)
            for(Int_t i=1; i<=nx; i++){
                
                mcnp_file>>val;
                edep->SetBinContent(i,j,val);
                mcnp_file>>val;
                
            }
        
        if(ntal>1){
            
            dummy="";
            
            while(dummy!="vals") mcnp_file>>dummy;
            
            for(Int_t j=1; j<=ny; j++)
                for(Int_t i=1; i<=nx; i++){
                    
                    mcnp_file>>val;
                    //extract z and correct for edep
                    Double_t ebc = edep->GetBinContent(i,j);
                    if(TMath::Abs(ebc)>1E-6){
                        
                        Double_t z = z2 - val/ebc;
                        depth->SetBinContent(i,j,z);
                    }
                    
                    mcnp_file>>val;
                }
        }
        
        if(ntal>2){
            
            dummy="";
            
            while(dummy!="vals") mcnp_file>>dummy;
            
            for(Int_t j=1; j<=ny; j++)
                for(Int_t i=1; i<=nx; i++){
                    
                    mcnp_file>>val;
                    //extract x and correct for edep
                    Double_t ebc = edep->GetBinContent(i,j);
                    if(TMath::Abs(ebc)>1E-6){
                        
                        Double_t xx = val/ebc-x1;
                        meanx->SetBinContent(i,j,xx);
                    }
                    
                    mcnp_file>>val;
                }
            
        }
        
        if(ntal>3){
            
            dummy="";
            
            while(dummy!="vals") mcnp_file>>dummy;
            
            for(Int_t j=1; j<=ny; j++)
                for(Int_t i=1; i<=nx; i++){
                    
                    mcnp_file>>val;
                    //extract y and correct for edep
                    Double_t ebc = edep->GetBinContent(i,j);
                    if(TMath::Abs(ebc)>1E-6){
                        
                        Double_t yy = val/ebc-y1;
                        meany->SetBinContent(i,j,yy);
                    }
                    
                    mcnp_file>>val;
                }
        }
        
        mcnp_file.close();
        
        //Done reading file
        //Change edep bin contents to energy
        
        for(Int_t j=1; j<=ny; j++)
            for(Int_t i=1; i<=nx; i++){
                
                Double_t ebc = edep->GetBinContent(i,j);
                Double_t e = ebc*bw_x*bw_y*bw_z*nps*1000.; //last factor to change to keV
                edep->SetBinContent(i,j,e);
            }
        
    }
    
    else{
        
        fin = new TFile(mcnp_fname,"READ");
        gROOT->cd();
        TTree* sinfo = (TTree*) fin->Get("sinfo");
        TTree* deposits = NULL; deposits = (TTree*) fin->Get("deposits");
        if(deposits==NULL){
            
            std::cout << ".pds file does not contain deposits tree!" << std::endl;
            return;
        }
        
        ntal=4;
        sinfo->SetBranchAddress("simid",&simid);
        sinfo->SetBranchAddress("nbins_x",&nx);
        sinfo->SetBranchAddress("nbins_y",&ny);
        sinfo->SetBranchAddress("bw_x",&bw_x);
        sinfo->SetBranchAddress("bw_y",&bw_y);
        sinfo->SetBranchAddress("bw_z",&bw_z);
        sinfo->GetEntry(0);
        
        edep = new TH2F("edep","edep",nx,0,nx,ny,0,ny);
        depth = new TH2F("depth","depth",nx,0,nx,ny,0,ny);
        meanx = new TH2F("meanx","meanx",nx,0,nx,ny,0,ny);
        meany = new TH2F("meany","meany",nx,0,nx,ny,0,ny);
        eleak = new TH2F("eleak","eleak",nx,0,nx,ny,0,ny);       // Added by KR
        
        //Now fill 2D histograms with energy deposits from tree
        Double_t dep_x, dep_y, dep_z, dep_e;
        deposits->SetBranchAddress("dep_x",&dep_x);
        deposits->SetBranchAddress("dep_y",&dep_y);
        deposits->SetBranchAddress("dep_z",&dep_z);
        deposits->SetBranchAddress("dep_e",&dep_e);
        
        Int_t nentries = deposits->GetEntries();
        
        for(int i=0; i<nentries; i++){
            
            deposits->GetEntry(i);
            int binx = dep_x / bw_x + 1;
            int biny = dep_y / bw_y + 1;
            
            Double_t pz = depth->GetBinContent(binx,biny);
            Double_t px = meanx->GetBinContent(binx,biny);
            Double_t py = meany->GetBinContent(binx,biny);
            Double_t pe = edep->GetBinContent(binx,biny);
            
            depth->SetBinContent(binx,biny, (dep_z*dep_e+pz*pe)/(dep_e+pe) );
            meanx->SetBinContent(binx,biny, (dep_x*dep_e+px*pe)/(dep_e+pe) );
            meany->SetBinContent(binx,biny, (dep_y*dep_e+py*pe)/(dep_e+pe) );
            edep->SetBinContent(binx,biny, dep_e+pe );
        }
    }
    
    //Create/load the blanks (pedestal noise templates)
    vector<TH2F*> raws;
    vector<Int_t> xid;
    vector<TTree*> vinfo;
    vector<Long_t> sat_val;
    //This only applies to extensions > 0
    Int_t xlow=1, xhigh=nx, ylow=1, yhigh=ny;
    
    TRandom3 r(0);
        
    if(blank_name==""){
    
        raws.push_back(new TH2F("image_raw", "image_raw", nx, 0, nx, ny, 0, ny));
        xid.push_back(simid);
        
        for(Int_t i=1; i<=nx; i++)
            for(Int_t j=1; j<=ny; j++){
                
                Double_t cont = r.Gaus(pedestal_mean,TMath::Sqrt(nreads*pedestal_sigma*pedestal_sigma));
                raws.at(0)->SetBinContent(i,j,cont);
            }
        
        TList* vars = new TList();
        vars->Add(new TNamed("HOSTNAME",hostname.Data()));
        vars->Add(new TNamed("PWD",pwd.Data()));
        vars->Add(new TNamed("SVNVERSION",svnversion.Data()));
        
        ifstream hinfo("header.txt");
        //Now the information from the .fits header
        TString vname, vtype;
        while(hinfo>>vname && hinfo>>vtype){
            
            if(vname=="HWID") vars->Add(new TParameter<double>("HWID", hwid));
            else if (vname=="RUNID") continue;
            else if(vname=="BITPIX"){
                
                vars->Add(new TParameter<double>("BITPIX",bit_pix));
                sat_val.push_back(TMath::Power(2,bit_pix)-1);
            }
            else if(vname=="NTAL") vars->Add(new TParameter<double>("NTAL", ntal));
            else if(vname=="NAXIS") vars->Add(new TParameter<double>("NAXIS",2));
            else if(vname=="NAXIS1") vars->Add(new TParameter<double>("NAXIS1",nx));
            else if(vname=="NAXIS2") vars->Add(new TParameter<double>("NAXIS2",ny));
            else if(vname=="TRIMSEC"){
                
                TString trimsec = "[1:"; trimsec += nx; trimsec += ",1:"; trimsec += ny; trimsec += "]";
                //in this case image was enlarged to avoid cache problems in MCNP
                if(nx==2056 && ny==4104)
                    trimsec = "[5:"; trimsec += nx-4; trimsec += ",5:"; trimsec += ny-4; trimsec += "]";
                vars->Add(new TNamed("TRIMSEC",trimsec.Data()));
            }
            
            else{
                
                if(vtype=="D" || vtype=="I") //For now store everything as double
                    vars->Add(new TParameter<double>(vname,0));
            
                else if(vtype=="S")
                    vars->Add(new TNamed(vname,TString("")));
            }
        }
        hinfo.close();
        
        TString oname = out_name; oname += "_"; oname += xid.back(); oname += ".root";
        vars->AddFirst(new TNamed(TString("ROOTFPATH"), oname));
        vars->AddFirst(new TNamed(TString("FITSFPATH"), TString("")));
        vars->AddFirst(new TNamed(TString("SIMFPATH"), mcnp_fname));
        
        vinfo.push_back(new TTree("finfo","finfo"));
        vinfo.back()->Branch("RUNID", &runid, "RUNID/I");
        vinfo.back()->Branch("EXTID", &simid, "EXTID/I");
        vinfo.back()->Branch(vars);
        vinfo.back()->Fill();
        
        vars->Delete();
    }
    
    else
        for(Int_t k=0; k<100; k++)
            for(Int_t l=0; l<100; l++){
                
                TString sex = "";
                sex += k;
                if(l>0) { sex = ""; sex += l; sex += "|"; sex += k; }
               
                if(!ConfigEnabled("extensions_enabled",sex)) continue;
       			         
                Int_t left_extension = l==0 ? k : l;
                Int_t right_extension = k;
                
                TFITSHDU *hdu_right = new TFITSHDU(blank_name, right_extension);
                TFITSHDU *hdu_left = NULL;
                if(left_extension != right_extension) hdu_left = new TFITSHDU(blank_name, left_extension);
                
                if(hdu_right->IsZombie() || (hdu_left!=NULL && hdu_left->IsZombie())){
                    
                    std::cout << "Could not open blank " << blank_name << std::endl;
                    return;
                }
                
                Bool_t compressed = ((TString) hdu_right->GetKeywordValue("XTENSION")) == "'BINTABLE'";
                
                Int_t nbins_x = compressed ? ((TString) hdu_right->GetKeywordValue("ZNAXIS1")).Atoi() : ((TString) hdu_right->GetKeywordValue("NAXIS1")).Atoi();
                Int_t nbins_y = compressed ? ((TString) hdu_right->GetKeywordValue("ZNAXIS2")).Atoi() : ((TString) hdu_right->GetKeywordValue("NAXIS2")).Atoi();
                raws.push_back(new TH2F("image_raw", "image_raw", nbins_x, 0, nbins_x, nbins_y, 0, nbins_y));
                xid.push_back(right_extension);
                
                //fill histogram
                for(Int_t i=1; i<=nbins_x; i++){
                    
                    TVectorD* currcol;
                    if(i<=nbins_x/2 && hdu_left!=NULL) currcol = hdu_left->GetArrayColumn(i-1);
                    else currcol = hdu_right->GetArrayColumn(i-1);
                    
                    for(Int_t j=1; j<=nbins_y; j++){
                        
                        Int_t cont = (Int_t) (*currcol)[j-1];
                        raws.back()->SetBinContent(i,j,cont);
                    }
                    
                    currcol->Delete();
                }
                
                TList* vars = new TList();
                vars->Add(new TNamed("HOSTNAME",hostname.Data()));
                vars->Add(new TNamed("PWD",pwd.Data()));
                vars->Add(new TNamed("SVNVERSION",svnversion.Data()));
                
                ifstream hinfo("header.txt");
                //Now the information from the .fits header
                TString vname, vtype, vval;
                while(hinfo>>vname && hinfo>>vtype){
                    
                    if(vname=="RUNID") { runid = ((TString) hdu_right->GetKeywordValue(vname)).Atoi(); continue; }
                    
                    if( (vname=="NAXIS1" || vname=="NAXIS2" || vname=="BITPIX") && compressed)
                        vval = hdu_right->GetKeywordValue((TString) "Z"+vname);
                    
                    else
                        vval = hdu_right->GetKeywordValue(vname);
                    
                    if(vval.Length()==0){
                        
                        if(vname=="HWID") vars->Add(new TParameter<double>(vname, hwid));
                        if(vname=="NTAL") vars->Add(new TParameter<double>(vname, ntal));
                    }
                    
                    else{
                        
                        if(vtype=="D" || vtype=="I") //For now store everything as double
                            vars->Add(new TParameter<double>(vname,vval.Atof()));
                        
                        else if(vtype=="S")
                            vars->Add(new TNamed(vname,vval));
                        
                        if(vname=="BITPIX")
                            sat_val.push_back(TMath::Power(2,vval.Atoi())-1);
                        
                        if(vname=="TRIMSEC"){
                            
                            xlow = ((TString) vval(vval.First('[')+1,vval.First(':')-vval.First('[')-1)).Atoi();
                            xhigh = ((TString) vval(vval.First(':')+1,vval.First(',')-vval.First(':')-1)).Atoi();
                            ylow = ((TString) vval(vval.First(',')+1,vval.Last(':')-vval.First(',')-1)).Atoi();
                            yhigh = ((TString) vval(vval.Last(':')+1,vval.First(']')-vval.Last(':')-1)).Atoi();
                            //We can't tell if half image is empty from header so far - hard code here :(
                            if( (nbins_x==4408) || (nbins_x==5408) || (nbins_x==5108) || (nbins_x==9204) || (nbins_x==8544) || (nbins_x==1062) || (nbins_x==1762) || (nbins_x==3758) || (nbins_x==9244) || (nbins_x==8534) || (nbins_x==8562) || (nbins_x==8304) || (nbins_x==8368) ) xlow = nbins_x/2 + 1;
                        }
                    }
                }
                hinfo.close();
                
                TString oname = out_name; oname += "_"; oname += xid.back(); oname += ".root";
                vars->AddFirst(new TNamed(TString("ROOTFPATH"), oname));
                vars->AddFirst(new TNamed(TString("FITSFPATH"), blank_name));
                vars->AddFirst(new TNamed(TString("SIMFPATH"), mcnp_fname));
                
                vinfo.push_back(new TTree("finfo","finfo"));
                vinfo.back()->Branch("RUNID", &runid, "RUNID/I");
                vinfo.back()->Branch("EXTID", &right_extension, "EXTID/I");
                vinfo.back()->Branch(vars);
                vinfo.back()->Fill();
                
                vars->Delete();
                hdu_right->Delete();
                if(hdu_left!=NULL) hdu_left->Delete();
            }
    
    //intermediate histogram after applying diffusion
    //entries are number of electrons in pixel
    TH2F* ediff = (TH2F*) edep->Clone("ediff");
    ediff->SetTitle("ediff");
    ediff->Reset();
    
    //add deposits to histogram
    for(Int_t i=1; i<=nx; i++)
        for(Int_t j=1; j<=ny; j++){
            
            //Added by KR
            Double_t Esim = edep->GetBinContent(i,j);
            
            Double_t ne = edep->GetBinContent(i,j)/e_ion;
			Double_t z_p=0;
            if(sim_flip_CCD==false)    //false means the depth is calculated from the top
                z_p = depth->GetBinContent(i,j)*1E4;
            else
                z_p = (bw_z-depth->GetBinContent(i,j))*1E4;
            
            //will deal with small/negative entries later
            //will deal with out of range depth later
            if(ne<1 || z_p<0 || z_p>bw_z*1E4) continue;
            
            Int_t eobs = 0;
            //number observed depends on probability of collecting one
            if(fano_factor<=0)
                eobs = ne;
            else
                eobs = r.Gaus(ne,TMath::Sqrt(fano_factor*ne)) + 0.5;

            Double_t x_p = meanx->GetBinContent(i,j)*1E4; //in um
            Double_t y_p = meany->GetBinContent(i,j)*1E4;
            
            if(x_p==0) x_p = (edep->GetXaxis()->GetBinCenter(i))*bw_x*1E4; //in um
            if(y_p==0) y_p = (edep->GetYaxis()->GetBinCenter(j))*bw_y*1E4; //in um
            
            Double_t var = 0;
            if(z_p > diff_offset) 
                var = -diff_A*TMath::Log(1. - diff_b*(z_p-diff_offset))*TMath::Power(1. + diff_alpha*Esim/(TMath::Sqrt(-diff_A*TMath::Log(1. - diff_b*(bw_z-diff_offset)))),2.);
                //var = -diff_A*TMath::Log(1. - diff_b*(z_p-diff_offset));
            
            for(Int_t k=0; k<eobs; k++){
                
                TVector3 v(0,0,0);
                if(var>0){
                    
                    v.SetX(r.Gaus(0,TMath::Sqrt(var)));
                    v.SetY(r.Gaus(0,TMath::Sqrt(var)));
                }
                
                Int_t binx = (x_p + v.X())/(bw_x*1E4) + 1;
                Int_t biny = (y_p + v.Y())/(bw_y*1E4) + 1;
               
                //consider that these clipped events may end within the final image if it is larger than ediff/edep histo
                if(binx>nx || biny>ny)
                    std::cout << "Warning: Deposit at " << binx << " " << biny << " outside ediff histo" << std::endl;
                    
                else
                    ediff->SetBinContent(binx,biny,ediff->GetBinContent(binx,biny)+1);
            }
        }
        
    //Added by KR
    //histogram of leakage current
    if (sim_inject_rate>0.){
        for(Int_t i=1; i<=nx; i++){
            for(Int_t j=1; j<=ny; j++){
                eleak->SetBinContent(i,j,r.Poisson(sim_inject_rate));
            }
        }
        ediff->Add(eleak,1);
    }
        
    //now simulate read out, considering CTI
    //the fast procedure should be very accurate (although not exact) and only works for a single CTI parameter in the CCD
    if(sim_cti_fast && sim_cti_1 > 0){
        
        if(sim_cti_2>0)
            std::cout << "Warning: fast CTI method only considers sim_cti_1 parameter" << std::endl;
        
        if(sim_cti_1>0.01)
            std::cout << "Warning: CTI is too large, fast procedure may not work" << std::endl;
        
        Double_t cti = sim_cti_1;
        
        for(Int_t i=1; i<=nx; i++){
            
            Double_t qo = 0;
            Double_t jo = 0;
            Double_t next = 0;
            Double_t sub = 0;
            Bool_t end = true;
            
            for(Int_t j=1; j<=ny; j++){
                
                Double_t bc = ediff->GetBinContent(i,j);
                Double_t nbc = bc + next;
                
                if(bc>0 || j==ny){
                
                    ediff->SetBinContent(i,jo,qo-sub);
                    jo=j;
                    qo=nbc;
                    sub=0;
                    
                    if(j==ny) ediff->SetBinContent(i,j,nbc);
                }
                
                else if(end)
                    continue;
                
                else
                    ediff->SetBinContent(i,j,nbc);
                
                //first opertion is the expected value
                next = qo*TMath::Power(1.-cti,jo)*TMath::Exp(TMath::LnGamma(j+2)-TMath::LnGamma(jo+1)-TMath::LnGamma(j-jo+2)+(j+1-jo)*TMath::Log(cti));
                //throw a random number only if expected value gets a 1 in at least once every 1000 trials
                if(next>0.001){
                    next = r.Poisson(next);
                    end = false;
                }
                else{
                    next = 0;
                    end = true;
                }
                
                sub += next;
            }
        }
    }
    
    //if not do the slow procedure. One extra loop (slower)
    else if(sim_cti_1 > 0 || sim_cti_2 > 0){
        
        TH2F* eremain = (TH2F*) ediff->Clone("eremain");
        ediff->Reset();
        
        for(Int_t i=1; i<=nx; i++){
            
            Int_t binx = i;
            
            for(Int_t j=1; j<=ny; j++){
                
                Int_t biny = j;
                
                for(Int_t jj=1; jj<=ny-j+1; jj++){
                    
                    Double_t cti = jj<sim_cti_bound ? sim_cti_1 : sim_cti_2;
                    
                    Double_t ntotal = eremain->GetBinContent(i,jj);
                    
                    if(ntotal==0) continue;
                    
                    //number of electrons that will remain
                    Double_t nremain = 0;
                    
                    if(cti>0){
                        
                        if(ntotal*cti > 10)
                            nremain = r.Gaus(ntotal*cti,TMath::Sqrt(ntotal*cti*(1.-cti)));
                        
                        else if(ntotal > 20 && cti < 0.05)
                            nremain = r.Poisson(ntotal*cti);
                        
                        else
                            nremain = r.Binomial(ntotal,cti);
                    }
                    
                    //move first row to image
                    
                    if(jj==1)
                        ediff->SetBinContent(binx,biny,ediff->GetBinContent(binx,biny) + (ntotal-nremain));
                    
                    else
                        eremain->SetBinContent(i,jj-1,eremain->GetBinContent(i,jj-1) + (ntotal-nremain));
                    
                    eremain->SetBinContent(i,jj,nremain);
                }
            }
        }
        
        eremain->Delete();
    }
    
    //do CTI in x (only fast procedure)
    if(sim_cti_x > 0){
        
        if(sim_cti_x>0.01)
            std::cout << "Warning: CTI is too large, fast procedure may not work" << std::endl;
        
        Double_t cti = sim_cti_x;
        
        for(Int_t j=1; j<=ny; j++){
            
            Double_t qo = 0;
            Double_t io = 0;
            Double_t next = 0;
            Double_t sub = 0;
            Bool_t end = true;
            
            for(Int_t i=1; i<=nx; i++){
                
                Double_t bc = ediff->GetBinContent(i,j);
                Double_t nbc = bc + next;
                
                if(bc>0 || i==nx){
                    
                    ediff->SetBinContent(io,j,qo-sub);
                    io=i;
                    qo=nbc;
                    sub=0;
                    
                    if(i==nx) ediff->SetBinContent(i,j,nbc);
                }
                
                else if(end)
                    continue;
                
                else
                    ediff->SetBinContent(i,j,nbc);
                
                //first opertion is the expected value
                next = qo*TMath::Power(1.-cti,io)*TMath::Exp(TMath::LnGamma(i+2)-TMath::LnGamma(io+1)-TMath::LnGamma(i-io+2)+(i+1-io)*TMath::Log(cti));
                //throw a random number only if expected value gets a 1 in at least once every 1000 trials
                if(next>0.001){
                    next = r.Poisson(next);
                    end = false;
                }
                else{
                    next = 0;
                    end = true;
                }
                
                sub += next;
            }
        }
    }
    
    //put the histogram on raws
    for(size_t n=0; n<raws.size(); n++){
        
        Int_t nbins_x = raws.at(n)->GetNbinsX();
        Int_t nbins_y = raws.at(n)->GetNbinsY();
        
        Int_t xofs=0, yofs=0;
        
        if(blank_name!=""){
            
            xofs = xlow-1;
            yofs = ylow-1;
        }
        
        for(Int_t i=1; i<=nx; i++){
            
            Int_t binx = i+xofs;
            if(binx>nbins_x) continue;
            
            for(Int_t j=1; j<=ny; j++){
                
                Int_t biny = j+yofs;
                if(biny>nbins_y) continue;
                
                raws.at(n)->SetBinContent(binx,biny,raws.at(n)->GetBinContent(binx,biny) + adc_unit*ediff->GetBinContent(i,j));
            }
        }
    }
    
    //round histogram content to the nearest integer
    for(size_t k=0; k<raws.size(); k++){
        
        Int_t nbins_x = raws.at(k)->GetNbinsX();
        Int_t nbins_y = raws.at(k)->GetNbinsY();
        
        for(Int_t i=1; i<=nbins_x; i++)
            for(Int_t j=1; j<=nbins_y; j++){
                
                Double_t cont = raws.at(k)->GetBinContent(i,j);
                if(cont>sat_val.at(k))
                    raws.at(k)->SetBinContent(i,j,sat_val.at(k));
                else if(cont<0)
                    raws.at(k)->SetBinContent(i,j,0);
                else
                    raws.at(k)->SetBinContent(i,j,(Int_t) cont+0.5);
            }
    }
    
    //write the image    
    for(size_t k=0; k<raws.size(); k++){
        
        gROOT->cd();
        TString oname = out_name; oname += "_"; oname += xid.at(k); oname += ".root";
        
        if(blank_name!="" && k==0){
            
            //In this case resize the simulation histograms to match the extension one
            edep = ResizeH2(edep, raws.at(k)->GetNbinsX(), raws.at(k)->GetNbinsY(), xlow, ylow);
            ediff = ResizeH2(ediff, raws.at(k)->GetNbinsX(), raws.at(k)->GetNbinsY(), xlow, ylow);
            depth = ResizeH2(depth, raws.at(k)->GetNbinsX(), raws.at(k)->GetNbinsY(), xlow, ylow);
            meanx = ResizeH2(meanx, raws.at(k)->GetNbinsX(), raws.at(k)->GetNbinsY(), xlow, ylow);
            meany = ResizeH2(meany, raws.at(k)->GetNbinsX(), raws.at(k)->GetNbinsY(), xlow, ylow);
        }
        
        TFile* f = new TFile(oname, "RECREATE");
        vinfo.at(k)->Write();
        raws.at(k)->Write();
        f->mkdir("sim");
        f->cd("sim");
        
        edep->Write();
        depth->Write();
        meanx->Write();
        meany->Write();
        ediff->SetNameTitle("nelec","nelec");
        ediff->Write();
        
        f->Close();
        f->Delete();
    }
    
    if(fin!=NULL) fin->Close();
}
