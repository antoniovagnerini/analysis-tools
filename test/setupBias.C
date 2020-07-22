/* 
   Creator: Antonio Vagnerini 
   14/01/2019 
*/
#include <iostream>
#include <memory>
#include <utility>
#include <string>
#include <tuple>

#include "Analysis/Models/src/classes.h"
#include "Analysis/Models/interface/RooDoubleGausExp.h"
#include "Analysis/Models/interface/RooQuadGausExp.h"

//RooFit includes
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooDataHist.h"
#include "RooBukinPdf.h"
#include "RooEffProd.h"
#include "RooNovosibirsk.h"
#include "RooDataSet.h"
#include "RooBukinPdf.h"

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

//string SignalModel(const vector<string>& );
//void setupSignal(  std::pair<string,string>&  );                            // freeze signal parameters
//void makeRooMultiPdfWorkspace(  std::tuple<string,string,string,string>& ); //multipdf for main and alternative models


void setupBias(){
  
   // Load the combine Library 
   gSystem->Load("libHiggsAnalysisCombinedLimit.so");

   //Signal fit [masspoint, model]
   std::vector<std::pair<string,string>> signalfit_models = { {"130","doublegausexp"}, {"140","doublegausexp"}, {"160","doublegausexp"}, {"180","doublegausexp"}, {"200","doublegausexp"}, {"250","doublegausexp"}, {"350","bukin"}, {"400","bukin"}, {"450","bukin"}, {"500","bukin"}, {"600","bukin"}, {"700","bukin"} };

   //Bkg fit [subrange ,main model, alternative model 1, alternative model2]
   std::vector<std::tuple<string,string,string,string>> bkgfit_models;
   bkgfit_models.push_back(std::make_tuple("sr1", "extnovoeffprod","superdijet0","berneffprod7"));
   //   bkgfit_models.push_back(std::make_tuple("sr1", "berneffprod7","extnovoeffprod","superdijet0"));
   bkgfit_models.push_back(std::make_tuple("sr2","supernovo0" ,"superdijetlinearprod0" ,"exppolynom3" ));


   //Setup signal workspace
   for ( auto& signalfit : signalfit_models )
     {
       // Open the signal workspaces
       string mp = signalfit.first;
       string model = signalfit.second;
 
       TFile *f_sig = TFile::Open( ("fitresults/signal_M"+mp+"_"+model+"_REGFSR/workspace/FitContainer_workspace.root").c_str());
       RooWorkspace &w_sig = *((RooWorkspace*)f_sig->Get("workspace")); //output workspace
       RooRealVar& mbb = *( (RooRealVar*) w_sig.var("mbb") );

       // Save to a new workspace
       TFile *fout = new TFile(("inputbiasfiles/Signal/signal_m"+mp+"_pdf.root").c_str(),"RECREATE");
       RooWorkspace wout("workspace");
   
       //Signal shape parameters list
       vector<string> parameters;
       vector<string> var_parameters; //parameters list that are not constant
       RooArgList vars = w_sig.allVars();
       auto iter = vars.createIterator();
       RooRealVar *v;
       while ((v = (RooRealVar*) iter->Next()) ){
	 if(string(v->GetName()) == string(mbb.GetName())) continue; // skip mbb
	 parameters.push_back( string(v->GetName()) );
       }

       //Form strings for RooFormulaVars
       map<string,string> formula;
       //   RooArgList nuisance_list; // Needed to add nuisances

       for(const auto &par : parameters){
	 RooRealVar& vCent = (RooRealVar&) *w_sig.var(par.c_str());
	 formula[par] = to_string(vCent.getValV());
       }

       //Make a RooFormulaVars:
       for(const auto &par : parameters){
	 //     if(find(var_parameters.begin(),var_parameters.end(),par) != var_parameters.end()){//If parameter is var = f(nuisances)
	 //       RooFormulaVar fPar( (par).c_str(),(par).c_str(),formula[par].c_str(),nuisance_list);
	 //       w.import(fPar);
	 //     }
	 //     else {
	 //if(par == "signal_norm") continue;
	 RooRealVar constant(("f_"+par).c_str(),("f_"+par).c_str(),stod(formula[par].c_str()));
	 constant.setConstant();
	 wout.import(constant);
	 RooFormulaVar fconstant(par.c_str(),par.c_str(),"@0",RooArgSet(constant));
	 wout.import(fconstant);
       }

       // Make signal function from parameters
       //  string model = SignalModel(parameters);
       std::unique_ptr<RooAbsPdf> func;
       if(model == "bukin"){
	 func = std::unique_ptr<RooBukinPdf>(new RooBukinPdf ("signal","signal",mbb,(RooFormulaVar&) *wout.function("Xp"),(RooFormulaVar&) *wout.function("sigp"),
							      (RooFormulaVar&)*wout.function("xi"),(RooFormulaVar&) *wout.function("rho1"),
							      (RooFormulaVar&)*wout.function("rho2")) );
       }
       else if (model == "doublegausexp"){
	 func = make_unique<analysis::models::RooDoubleGausExp>("signal","signal",mbb,(RooFormulaVar&) *wout.function("mean"),(RooFormulaVar&) *wout.function("sigmaL"),
									 (RooFormulaVar&) *wout.function("sigmaR"),(RooFormulaVar&)*wout.function("tail_shift"),
									 (RooFormulaVar&) *wout.function("tail_sigma"));
       }
       else if (model == "quadgausexp"){
	 func = make_unique<analysis::models::RooQuadGausExp>("signal","signal",mbb,(RooFormulaVar&) *wout.function("mean"),(RooFormulaVar&) *wout.function("sigmaL1"),
								       (RooFormulaVar&) *wout.function("sigmaL2"),
								       (RooFormulaVar&) *wout.function("sigmaR1"),
								       (RooFormulaVar&) *wout.function("sigmaR2"),
								       (RooFormulaVar&)*wout.function("tail_shift"),
								       (RooFormulaVar&) *wout.function("tail_sigma"),
								       (RooFormulaVar&) *wout.function("norm_g1"),
								       (RooFormulaVar&) *wout.function("norm_g2"));
       }
       else throw logic_error("ERROR");
       wout.import(*func.get());
       wout.Print();
       wout.Write();
   
     }
  

   //Setup multipdf object with background fits and data_obs
   for ( auto& bkgfit : bkgfit_models )
     {
       //Get subrange and model name
       string subrange = std::get<0>(bkgfit);
       string main     = std::get<1>(bkgfit); //models
       string a1       = std::get<2>(bkgfit);
       string a2       = std::get<3>(bkgfit);

       // !!!!!!!!!!!!!!!!!!!!!!!!!!!!! nobserved hard-coded in FITMODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       int n_observed; 
       if (subrange == "sr1") n_observed = 129643;
       else n_observed = 1932;
       int nbins = 5000;

       std::vector<std::string> models { main,a1,a2 } ; 
       int model_index = 0;

       for ( auto& model : models )
	 {
	   // Open the bkg workspaces  //USE ITERATOR ON MODEL
	   TFile *f_model = TFile::Open(("fitresults/background_"+subrange+"_"+model+"/workspace/FitContainer_workspace.root").c_str());
	   RooWorkspace *w_model = (RooWorkspace*)f_model->Get("workspace");
     
	   //Prepare rooDataHist for data_obs
	   RooRealVar *mbb =  w_model->var("mbb");
	   mbb->setBins(nbins);

	   //Turn-on Parameter Fix
	   if ( subrange == "sr1" )
	     {
	       RooRealVar & slope_model  =  *w_model->var("slope_novoeff");
	       RooRealVar & turnon_model =  *w_model->var("turnon_novoeff");
	       slope_model.setConstant();
	       turnon_model.setConstant();
	       //	   slope_model.SetName(("slope_"+ model).c_str());
	       //	   turnon_model.SetName(("turnon_"+model).c_str());
	       
	       if (  model.find("dijet") == std::string::npos) {
		 cout << "model " << model << endl;
		 RooAbsReal & effprod_model =  *w_model->function("background_eff");
		 effprod_model.SetName(("background_eff_"+ model).c_str());
	       }

	     }
	   RooAbsPdf *pdf_model = w_model->pdf("background");
	   pdf_model->SetName(model.c_str());
       
	   // Fit the functions to the data to set the "prefit" state (note this can and should be redone with combine when doing 
	   // bias studies as one typically throws toys from the "best-fit"
	   RooDataSet *data = (RooDataSet*)w_model->data("data_obs");
           pdf_model->fitTo(*data);  // index 0 uses 5000 bins fit -> improves a lot bias 
       
	   // Make a plot (data is a toy dataset)
	   RooPlot *plot = mbb->frame();   data->plotOn(plot);
	   pdf_model->plotOn(plot,RooFit::LineColor(kGreen));
	   plot->SetTitle("PDF fits to toy data");
	   plot->Draw();
       
	   // Make a RooCategory object. This will control which of the pdfs is "active"
	   RooCategory cat("pdf_index","Index of Pdf which is active");

	   // Make a RooMultiPdf object. The order of the pdfs will be the order of their index, ie for below 
	   // 0 == model	      
	   RooArgList mypdfs;
	   mypdfs.add(*pdf_model);
	   RooMultiPdf multipdf("roomultipdf","All Pdfs",cat,mypdfs);
       
	   // As usual make an extended term for the background with _norm for freely floating yield
	   //   RooRealVar bg_norm("roomultipdf_norm","Number of background events",0,n_observed);
	   RooRealVar bg_norm("roomultipdf_norm","Number of background events",0,n_observed*2);

	   TFile *fmulti = new TFile(("inputbiasfiles/Background/background_pdf_"+to_string(model_index)+"_"+subrange+".root").c_str(),"RECREATE");
	   RooWorkspace wmulti("backgrounds","backgrounds");
	   wmulti.import(cat);
	   wmulti.import(bg_norm);
	   wmulti.import(multipdf);
	   wmulti.Print();
	   wmulti.Write();
	 
	   // Save data_obs
	   bool blinded = true;
	   TFile *fdata = new TFile(("inputbiasfiles/Data/data_obs_"+subrange+".root").c_str(),"RECREATE");
	   std::cout << "data_obs from Integral = " << data->sumEntries() << " while real n_observed = " << n_observed << std::endl;                                     	   w_model->Write();
	   ++model_index;
	 }
    }
 



}



