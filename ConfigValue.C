#include <fstream>
#include <iostream>
#include "TMap.h"
#include "TObjString.h"
#include "TTree.h"
#include "TList.h"
#include "TParameter.h"
#include "TNamed.h"

TMap* config = new TMap();

void LoadSingleConfigFile(TString fname){
    
    string rawline;
    ifstream in(fname);
    
    while(getline(in,rawline)){
        
        TString line = rawline;
        if(line.BeginsWith("#")) continue;
        
        TString cname = "";
        TString cval = "";
        
        if(line.Length()==0) continue;
        while(line.First(" ")==0) line.Remove(0,1);
        cname=line(0,line.First(" "));
        line.Remove(0,cname.Length()+1);
        while(line.First(" ")==0) line.Remove(0,1);
        cval=line;
        
        if(!config->Contains(cname)){
            
            TObjString* oname = new TObjString(cname);
            TObjString* oval = new TObjString(cval);
            config->Add(oname, oval);
        }
    }
    
    in.close();
}

void LoadConfigFile(TString fname){
    
    LoadSingleConfigFile(fname);
    LoadSingleConfigFile("default.cfg");
}

Bool_t ConfigEnabled(TString cname, TString cint){
    
    if(config->Contains(cname)){
        
        TString cval = " ";
        cval += cint;
        cval += " ";
        TString s = " ";
        s += ((TObjString*) config->GetValue(cname))->GetString();
        s += " ";
        return s.Contains(cval);
    }
    
    else{
        
        std::cout << cname << " is not a config parameter!" << std::endl;
        return false;
    }
}

Double_t ConfigValue(TString cname){
    
    if(config->Contains(cname)){
        
        TString s = ((TObjString*) config->GetValue(cname))->GetString();
        if(s=="true") return 1;
        else if(s=="false") return 0;
        else return s.Atof();
    }
    
    else{
        
        std::cout << cname << "is not a config parameter!" << std::endl;
        return 0;
    }
}

void AddConfigToTree(TTree* t){
    
    TList* blist = new TList();
    TString preffix = "CFG_";
    
    TIter mit(config);
    TObjString* ckey;
    TObjString* cval;
    while((ckey = (TObjString*) mit.Next())){
                   
        cval = (TObjString*) config->GetValue(ckey);
        TString skey = preffix + ckey->GetString();
        TString sval = cval->GetString();
        
        if(sval.IsFloat())
            blist->Add(new TParameter<double>(skey,sval.Atof()));
        else if(sval=="true")
            blist->Add(new TParameter<double>(skey,1));
        else if(sval=="false")
            blist->Add(new TParameter<double>(skey,0));
        else
            blist->Add(new TNamed(skey,sval));
    }
    
    //Add branches to tree with their values
    //Add the same value for all entries
    t->Branch(blist);
    for(Int_t i=0; i<blist->GetSize(); i++)
        for(Int_t j=0; j<t->GetEntries(); j++)
            (t->GetBranch(blist->At(i)->GetName()))->Fill();
    
    blist->SetOwner();
    blist->Delete();
}

