#include "tdrstyle.C"
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "Fit/Fitter.h"
#include "TStyle.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH2D.h"
#include "TLegend.h"


#include "HbbStylesNew.cc"
using namespace std;

void biasResults()
{   
  HbbStylesNew style;
  style.SetStyle();
  gStyle->SetOptStat(1011); // display fit on plot
  setTDRStyle();

  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);

  const int mumax=0;
  const string toypdf_index = "2";
  const string fitpdf_index = "0";
  const string subrange = "sr1";
  const int ntoys = 1000; 

  std::vector<std::string> fits = {"fit_b","fit_s"};

  std::map<std::string, TH1F*> h1;
  std::map<std::string, TH2F*> h2;

  std::vector<std::string> mps;

    if ( subrange == "sr1")  mps = {"120"};
  //  if ( subrange == "sr2")  mps = {"250","350","400",};
  //  if ( subrange == "sr3")  mps = {"500","600","700"};
    
  //Parameters bounds (Novosibirsk with fixed turn-on parameters)
  std::map<std::string, std::pair<double,double>> params = {
    {"CMS_lumi_13TeV", {std::make_pair(-0.0001,0.0001)} },
    {"peak", {std::make_pair(138.,142.)}},
    {"shapeBkg_background_bbHTo4b__norm", {std::make_pair(125000.,135000.)} },
    {"tail", {std::make_pair(-0.65,-0.55)} },
    {"width", {std::make_pair(50.,60.)} },
    {"r" , {std::make_pair(-4.,4.)} },
  };

  const string pdf_indices = "toypdf_"+toypdf_index+"_fitpdf_"+ fitpdf_index;
  const string output = "PLOTS/" +subrange+ "_"+ pdf_indices +".pdf";

  //Latex table for output bias results
  ofstream myfile;
  myfile.open ( ("bias_results_"+ pdf_indices+"_"+subrange+".tex").c_str()) ; 
  myfile << "\\begin{table}"  << std::endl ; 
  myfile <<  "\\centering"  << std::endl ;
  myfile <<  "\\begin{tabular}{|c|c|c|c|c|c|}"  << std::endl ;
  myfile <<  "\\hline" << std::endl ;
  myfile <<  " $m_{A/H}$  & Toy failed $\\%$ & \\textbf{$\\mu= 0$}  & \\textbf{$\\mu=1$}  & \\textbf{$\\mu=2$}  & \\textbf{$\\mu=3$} \\\\ "  << std::endl ;

  TCanvas* c = style.MakeCanvas("c","",700,700);

  for (auto & mp : mps){

    int pull_MINOS     [mumax+1] ;   
    int pull_LIKELIHOOD[mumax+1] ; 

    int pull_err_MINOS     [mumax+1] ;   
    int pull_err_LIKELIHOOD[mumax+1] ; 

    double maxtoy_failpc_MINOS =0 ;
    double maxtoy_failpc_LIKELIHOOD =0;

    double mu_maxtoyfailpc_MINOS = -1;  //signal strength corresponding to max # toy failures
    double mu_maxtoyfailpc_LIKELIHOOD = -1;
    
    for ( int mu =0 ; mu<=mumax ; mu++)
      {

	//Title header
	c->cd();
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.9, 1, 1.0);
	pad1->SetBottomMargin(0); // Upper and lower plot are joined
	gStyle->SetOptTitle(0);
	pad1->Draw();             
	pad1->cd();               
	gStyle->SetOptTitle(0);
	TPaveLabel *title = new TPaveLabel(0.,0.,1.,1., ("#mu = "+to_string(mu)+"\t M_{A/H} ="+ mp+" GeV \t"+pdf_indices ).c_str() ,"brndc");
	title->Draw();
	
	//Histograms
	c->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.9);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->Draw();
	pad2->Divide(3,2);
      	
	//Likelihood and MINOS trees
	TFile *f1 = TFile::Open(("results/fitDiagnostics_mu"                 +to_string(mu)+"_m" + mp+  "_"+ pdf_indices+".root").c_str());
	TFile *f2 = TFile::Open(("results/higgsCombineTest.FitDiagnostics_mu"+to_string(mu)+"_m" + mp+  "_"+ pdf_indices+".root").c_str());

	//Useful params
	int injected_mu = mu;
	int low = -4;
	int high = 4;
	low+=mu;
	high+=mu;

	//Fitted mu -MINOS
	pad2->cd(1);
	TTree * tree_fit_sb = (TTree*)f1->Get("tree_fit_sb"); 
	//tree_fit_sb->Draw(("r>>h(20,"+to_string(low)+","+to_string(high)+")").c_str(), "fit_status>=0 && fabs(nll_min) < 50.");
	tree_fit_sb->Draw(("r>>h(20,"+to_string(low)+","+to_string(high)+")").c_str(), "fit_status>=0");
	TH1F *h = (TH1F*)gDirectory->Get("h");
	style.InitHist(h, "#mu_{fit}","N toys",kBlue,0);
	h->SetMaximum(ntoys/5);
	h->Draw();
	h->Fit("gaus");
	//Check fraction of toys failing the fits
	double toyfail_thismu_M = 1-(h->GetEntries()/ntoys);
	if ( toyfail_thismu_M > maxtoy_failpc_MINOS ) 
	  {
	    maxtoy_failpc_MINOS = toyfail_thismu_M;
	    mu_maxtoyfailpc_MINOS = mu;
	  }

	//Fitted mu -LIKELIHOOD
	pad2->cd(4);
	TTree * limit = (TTree*)f2->Get("limit"); 
	limit->Draw(("limit>>h2(20,"+to_string(low)+","+to_string(high)+")").c_str(),"limitErr>0.1");
	TH1F *h2 = (TH1F*)gDirectory->Get("h2");
	style.InitHist(h2, "#mu_{fit}","N toys",kBlue,0);
	h2->SetMaximum(ntoys/5);
	h2->Draw();
	h2->Fit("gaus");  
	double toyfail_thismu_L = 1-(h2->GetEntries()/ntoys);
	if ( toyfail_thismu_L > maxtoy_failpc_LIKELIHOOD )
	  { 
	    maxtoy_failpc_LIKELIHOOD = toyfail_thismu_L;
	    mu_maxtoyfailpc_LIKELIHOOD = mu;	    
	  }

	//rErr - MINOS      
	pad2->cd(2);
	//tree_fit_sb->Draw("rErr>>herr(100,0,8)", "fit_status>=0 && fabs(nll_min) < 50. ");
	tree_fit_sb->Draw("rErr>>herr(100,0,8)", "fit_status>=0");
	TH1F *herr = (TH1F*)gDirectory->Get("herr");
	style.InitHist(herr, "#sigma_{#mu_{fit}}","N toys",kBlue,0);
	gStyle->SetOptStat(1); // display
	//	herr->SetTitleSize(0.15);
	herr->SetMaximum(ntoys);
	herr->Draw();

	//rErr -LIKELIHOOD
	pad2->cd(5); 
	limit->Draw("limitErr>>h2err(100,0,8)","limitErr>0.1");
	TH1F *h2err = (TH1F*)gDirectory->Get("h2err");
	style.InitHist(h2err, "#sigma_{#mu_{fit}}","N toys",kBlue,0);
	gStyle->SetOptStat(1); // display
	h2err->SetMaximum(ntoys);
	h2err->Draw();

	low-=mu;
	high-=mu;

	//pull mu -MINOS
	gStyle->SetOptStat(1011); // display
	pad2->cd(3);
	tree_fit_sb->Draw(("(r-"+to_string(mu)+")/rErr>>hpull(20,"+to_string(low)+","+to_string(high)+")").c_str(), "fit_status>=0");
	//	tree_fit_sb->Draw(("(r-"+to_string(mu)+")/rErr>>hpull(20,"+to_string(low)+","+to_string(high)+")").c_str(), "fit_status>=0 && fabs(nll_min) < 50.");
	TH1F *hpull = (TH1F*)gDirectory->Get("hpull");
	style.InitHist(hpull, "(#mu_{fit}-#mu_{in})/#sigma_{#mu_{fit}}","N toys",kBlue,0);
	hpull->SetMaximum(ntoys/5);
	hpull->Draw();

	//Gaussian fit for bias value extraction
	TF1 *g1 = new TF1("g1","gaus",low,high);
	hpull->Fit(g1,"R");  
	pull_MINOS[mu]     = g1->GetParameter(1)*100;
        pull_err_MINOS[mu] = g1->GetParError(1)*100;

     
	//pull mu -LIKELIHOOD
	pad2->cd(6);
	limit->Draw(("(limit-"+to_string(mu)+")/limitErr>>h2pull(20,"+to_string(low)+","+to_string(high)+")").c_str(),"limitErr>0.1");
	TH1F *h2pull = (TH1F*)gDirectory->Get("h2pull");
	style.InitHist(h2pull, "(#mu_{fit}-#mu_{in})/#sigma_{#mu_{fit}}","N toys",kBlue,0);
	h2pull->SetMaximum(ntoys/5);
	h2pull->Draw();
	h2pull->Fit("gaus");  

	TF1 *g2 = new TF1("g2","gaus",low,high);
	h2pull->Fit(g2,"R");  
	pull_LIKELIHOOD[mu]     = g2->GetParameter(1)*100;
        pull_err_LIKELIHOOD[mu] = g2->GetParError(1)*100;

	c->cd();
	c-> Update();

	string page;
	if      ( mu == 0     && mp==mps[0]            ) page = "(";
	else page = "";
	c-> Print( ( output+page ).c_str(), ("Mu = "+to_string(mu)+" M_{A/H} ="+ mp+" GeV "+pdf_indices ).c_str());

      } // end signal strength loop

    myfile.precision(3);
    myfile <<  "\\hline" << std::endl;
    myfile <<  mp << " & " << maxtoy_failpc_MINOS*100  << " for $\\mu = $" << mu_maxtoyfailpc_MINOS  << " & " << pull_MINOS[0] << " $\\pm$ " << pull_err_MINOS[0] << " & " << pull_MINOS[1] << " $\\pm$ " << pull_err_MINOS[1] << " & " << pull_MINOS[2] << " $\\pm$ " << pull_err_MINOS[2] << " & " << pull_MINOS[3] << " $\\pm$ " << pull_err_MINOS[3]  << " \\\\ "  << std::endl ;

    myfile <<  "\\hline"<< std::endl ;
  

    myfile <<  mp << " & " << maxtoy_failpc_LIKELIHOOD*100  << " for $\\mu = $" << mu_maxtoyfailpc_LIKELIHOOD   << " & " << pull_LIKELIHOOD[0] << " $\\pm$ " << pull_err_LIKELIHOOD[0] << " & " << pull_LIKELIHOOD[1] << " $\\pm$ " << pull_err_LIKELIHOOD[1] << " & " << pull_LIKELIHOOD[2] << " $\\pm$ " << pull_err_LIKELIHOOD[2] << " & " << pull_LIKELIHOOD[3] << " $\\pm$ " << pull_err_LIKELIHOOD[3]  << " \\\\ "  << std::endl ;

     
  } // end masspoint loop

  myfile <<  "\\hline" << std::endl;
  myfile <<  "\\end{tabular}"  << std::endl ; 
  myfile <<  "\\end{table}"  << std::endl ; 

  myfile.close();

  //  c-> Print( (output+")").c_str(), "");                                   


  // Sanity checks for fit parameters and par errors distributions
  //Parameter loop
  TFile *fin = TFile::Open(("fit_params/fit_params_"+subrange+"_"+pdf_indices+".root").c_str());

   for (auto & mp : mps) 
     {
     for ( int mu =0 ; mu<=mumax ; mu++)
       {      
	 for (auto & fit : fits )
	   {

	     int npar= 0;
	     int nerr= 0;

	     //-----First page --- PARMATERS -----------------------
	     c->cd();
  	     TPad *pad1 = new TPad("pad1", "pad1", 0, 0.9, 1, 1.0);
	     pad1->SetBottomMargin(0); 
	     gStyle->SetOptTitle(0);
	     pad1->Draw();             
	     pad1->cd();               
	     gStyle->SetOptTitle(0);
	     TPaveLabel *title = new TPaveLabel(0.,0.,1.,1., ("#mu = "+to_string(mu)+"\t M_{A/H} ="+ mp+" GeV \t"+pdf_indices +" par_"+ fit ).c_str() ,"brndc");
	     title->Draw();

	     //Histograms
	     c->cd();
	     TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.9);
	     pad2->SetTopMargin(0);
	     pad2->SetBottomMargin(0.2);
	     pad2->Draw();
	     pad2->Divide(3,2);

	     //Parameters
	     for ( auto & param : params )
	       {		 
		 ++npar;
		 pad2-> cd(npar);
		 TH1F * h = new TH1F;
		 cout << param.first << mu <<mp<< fit << endl; 
		 h = (TH1F*)fin->Get(Form("%s_mu%i_m%s_%s",param.first.c_str(), mu, mp.c_str(), fit.c_str()));
		 style.InitHist(h, param.first.c_str() ,"N toys",kBlue,0);
		 //		 h->SetMaximum(ntoys/20);
		 h->Draw();
	       }
	   
	     //NLL 
	     pad2->cd(6);
	     TH1F * h_nll = new TH1F;
	     h_nll = (TH1F*)fin->Get(Form("nll_min_mu%i_m%s_%s", mu, mp.c_str(), fit.c_str())) ;
	     style.InitHist(h_nll, "nll_min" ,"N toys",kBlue,0);
	     //  h->SetMaximum(ntoys/20);
	     h_nll->Draw();

	     c-> Update();
	     c-> Print( output.c_str(), ("Mu = "+to_string(mu)+" M_{A/H} ="+ mp+" GeV "+pdf_indices +"par_"+ fit ).c_str());


	     //Second page -------- PARAMETERS ERRORS -----------------------------------
	     c->cd();
	     //Histograms
	     TPad *pad_2 = new TPad("pad_2", "pad_2", 0, 0.05, 1, 0.9);
	     pad_2->SetTopMargin(0);
	     pad_2->SetBottomMargin(0.2);
	     pad_2->Draw();
	     pad_2->Divide(3,2);

	     //Parameters errs
	     for ( auto & param : params )
	       {		 
		 ++nerr;
		 pad_2-> cd(nerr);
		 TH1F * herr = new TH1F;
		 herr = (TH1F*)fin->Get(Form("%s_err_mu%i_m%s_%s",param.first.c_str(), mu, mp.c_str(), fit.c_str()));
	         style.InitHist(herr, (param.first+"_error").c_str() ,"N toys",kBlue,0);
		 herr->Draw();
	       }

	     string page;
	     if      ( mu == mumax && mp==mps[mps.size()-1] && fit==fits[fits.size() -1]  ) page = ")";
	     else page = "";
	     c-> Update();
	     c-> Print( (output+page).c_str(), ("Mu = "+to_string(mu)+" M_{A/H} ="+ mp+" GeV "+pdf_indices +"parerr_"+ fit ).c_str());
       
	   } //end fit loop

       } //end mu loop

     }//end mp loop	   	 
	 
}	
       
	







