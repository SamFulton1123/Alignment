#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <string>
#include <algorithm> 
#include <stdio.h>
#include <string.h> 

vector<int>accepted_tracks;
int spill = 1;
char infile[ ] = "../data/599/SSD20.root";
double station0_zpos = -170;
double station1_zpos = 160;
double station2_zpos = 920;
double station_thickness = 120;
double ssdPitch = 0.06;
double constr = 0.07;


//.........................................................................

vector<int> find_bounds()
{
   
   ///Return Vector with starting index of each event

   //Read in data from sssd.root file
   TFile *f = new TFile(infile);
   TTree *t1 = (TTree*)f->Get("SSDtree");
   int fer, module, row;
   vector<int> result = {0};
   t1->SetBranchAddress("fer",&fer);
   t1->SetBranchAddress("module",&module);
   t1->SetBranchAddress("row",&row);
   
   //read all entries and fill the histograms
   for (int i = 0; i<t1->GetEntries(); i++) {
     t1->GetEntry(i);
     	int first_fer = fer; 
     t1->GetEntry(i+1);
        int second_fer = fer;
    if(first_fer == 3 && second_fer ==0){
            result.push_back(i);	    
            }
   }   
   return result;
}

vector<int> bounds = find_bounds(); 


//.........................................................................


vector<double> initial_xalign(vector<double> residuals){
    
   // Returns vector of averages of hits for each SSD
   // Use to center SSDs before alignment 

   vector<double> new_residuals = {0,0,0,0,0,0}; 
   vector<int>counts = {0,0,0,0,0,0};

   //Read in data from SSD.root tree 
   TFile *f = new TFile(infile);
   TTree *t1 = (TTree*)f->Get("SSDtree");
   int fer, module, row;
   t1->SetBranchAddress("fer",&fer);
   t1->SetBranchAddress("module",&module);
   t1->SetBranchAddress("row",&row);
   auto h1 = new TH2F("h1","xz-Hits",150,-200,1200,100,-35,35);//histogram with parameters
   auto event = new TH2I("event","xz-Hits",150,-200,1200,40,-45,45);//histogram with individual event
   
   

//loops over single event and fills vector
for(int j=0; j<bounds.size(); j++){

   vector<vector<double>> xpos = {{},{},{},{},{},{}};
   vector<vector<double>> zpos = {{},{},{},{},{},{}};
   
   
   for(int i=bounds[j];i<bounds[j+1];i++){
	t1->GetEntry(i);
	//Filling vectors for xz-plot
	if(fer ==0 && module==1){
		xpos[0].push_back((row-320)*ssdPitch-residuals[0]);
		zpos[0].push_back(station0_zpos);		
	}else if(fer ==0 && module==3){
		xpos[1].push_back((row-320)*ssdPitch-residuals[1]);
		zpos[1].push_back(station0_zpos+station_thickness);		
	}else if(fer ==1 && module==1){
		xpos[2].push_back((row-320)*ssdPitch-residuals[2]);
		zpos[2].push_back(station1_zpos);		
	}else if(fer ==1 && module==3){
		xpos[3].push_back((row-320)*ssdPitch-residuals[3]);
		zpos[3].push_back(station1_zpos+station_thickness);		
	}else if(fer ==2 && module==0){
		xpos[4].push_back((row)*ssdPitch-residuals[4]);
		zpos[4].push_back(station2_zpos);
	}else if(fer ==2 && module==1){
		xpos[4].push_back(-(row)*ssdPitch-residuals[4]);
		zpos[4].push_back(station2_zpos);		
	}else if(fer ==3 && module==0){
		xpos[5].push_back((row-640)*ssdPitch-residuals[5]);
		zpos[5].push_back(station2_zpos+station_thickness);		
	}else if(fer ==3 && module==1){
		xpos[5].push_back(-(row-640)*ssdPitch-residuals[5]);
		zpos[5].push_back(station2_zpos+station_thickness);
		
	}	
	
   }
  
   //Only keep events with single track
   if( (zpos[0].size()==1 || (zpos[0].size()==2 && abs(xpos[0][0]-xpos[0][1])<constr ))
    && (zpos[1].size()==1 || (zpos[1].size()==2 && abs(xpos[1][0]-xpos[1][1])<constr ))
    && (zpos[2].size()==1 || (zpos[2].size()==2 && abs(xpos[2][0]-xpos[2][1])<constr ))
    && (zpos[3].size()==1 || (zpos[3].size()==2 && abs(xpos[3][0]-xpos[3][1])<constr ))
    && (zpos[4].size()==1 || (zpos[4].size()==2 && abs(xpos[4][0]-xpos[4][1])<constr ))
    && (zpos[5].size()==1 || (zpos[5].size()==2 && abs(xpos[5][0]-xpos[5][1])<constr )) ){
	accepted_tracks.push_back(j);

	
	for(int i = 0; i<zpos.size();i++){
		if(xpos[i].size() !=0 && zpos[i].size() !=0){
			new_residuals[i] += xpos[i][0];
			h1->Fill(zpos[i][0],xpos[i][0]);
			}
		}

	event->Reset("ICESM");		
	}
  
  }
   
   auto c = new TCanvas("c","c");
   h1->SetBarWidth(10);
   h1->SetFillStyle(0);
   h1->SetFillColor(kGray);
   h1->SetLineColor(kBlue);
   h1->GetYaxis()->SetTitle("x-position mm");
   h1->GetXaxis()->SetTitle("z-position mm");
   h1->SetStats(0);
   h1->Draw("violiny(112000000)");
   
   for(int i=0;i<new_residuals.size();i++){
	new_residuals[i]/=accepted_tracks.size();
	}
   new_residuals.push_back(new_residuals[5]);
   new_residuals.push_back(new_residuals[5]);
   new_residuals[5] = 0;
   new_residuals[5]+=new_residuals[4];

  return new_residuals; 
}



//.........................................................................

vector<double> xalign(vector<double> residuals){
   
   //Input: vector with SSD offsets
   //Returns a vector with new offsets
  
   //Read in data from SSD.root tree 
   TFile *f = new TFile(infile);
   TTree *t1 = (TTree*)f->Get("SSDtree");
   int fer, module, row;
   t1->SetBranchAddress("fer",&fer);
   t1->SetBranchAddress("module",&module);
   t1->SetBranchAddress("row",&row);
   vector<double> new_residuals;
   //Histstograms filled with residuals
   TH1F** residual_array = new TH1F*[8]; 
   for (int i=0;i<8;i++) { 
	residual_array[i] = new TH1F(Form("h%d",i),Form("Residuals %d",i),100,-300,300);
	}

   auto event = new TH2I("event","xz-Hits",150,-200,1200,40,-45,45);//histogram with single event
   
   

//loops over single event and fills vector
for(int j=0; j<bounds.size(); j++){

   vector<vector<double>> xpos = {{},{},{},{},{},{},{},{}};
   vector<vector<double>> zpos = {{},{},{},{},{},{},{},{}};
   
   
   for(int i=bounds[j];i<bounds[j+1];i++){
	t1->GetEntry(i);
	//Filling vectors for xz-plot
	if(fer ==0 && module==1){
		xpos[0].push_back((row-320)*ssdPitch-residuals[0]);
		zpos[0].push_back(station0_zpos);		
	}else if(fer ==0 && module==3){
		xpos[1].push_back((row-320)*ssdPitch-residuals[1]);
		zpos[1].push_back(station0_zpos+station_thickness);		
	}else if(fer ==1 && module==1){
		xpos[2].push_back((row-320)*ssdPitch-residuals[2]);
		zpos[2].push_back(station1_zpos);		
	}else if(fer ==1 && module==3){
		xpos[3].push_back((row-320)*ssdPitch-residuals[3]);
		zpos[3].push_back(station1_zpos+station_thickness);		
	}else if(fer ==2 && module==0){
		xpos[4].push_back((row)*ssdPitch-residuals[4]);
		zpos[4].push_back(station2_zpos);
	}else if(fer ==2 && module==1){
		xpos[5].push_back(-(row)*ssdPitch-residuals[5]);
		zpos[5].push_back(station2_zpos);		
	}else if(fer ==3 && module==0){
		xpos[6].push_back((row-640)*ssdPitch-residuals[6]);
		zpos[6].push_back(station2_zpos+station_thickness);		
	}else if(fer ==3 && module==1){
		xpos[7].push_back(-(row-640)*ssdPitch-residuals[7]);
		zpos[7].push_back(station2_zpos+station_thickness);
		
	}	
	
	
   }
  
   //Only keep events with single track
   if( (zpos[0].size()==1 || (zpos[0].size()==2 && abs(xpos[0][0]-xpos[0][1])<constr ))
    && (zpos[1].size()==1 || (zpos[1].size()==2 && abs(xpos[1][0]-xpos[1][1])<constr ))
    && (zpos[2].size()==1 || (zpos[2].size()==2 && abs(xpos[2][0]-xpos[2][1])<constr ))
    && (zpos[3].size()==1 || (zpos[3].size()==2 && abs(xpos[3][0]-xpos[3][1])<constr ))
    && (((zpos[4].size()==1|| (zpos[4].size()==2 && abs(xpos[4][0]-xpos[4][1])<constr )) && zpos[5].size()==0)
    ||  ((zpos[5].size()==1|| (zpos[5].size()==2 && abs(xpos[5][0]-xpos[5][1])<constr )) && zpos[4].size()==0))
    && (((zpos[6].size()==1|| (zpos[6].size()==2 && abs(xpos[6][0]-xpos[6][1])<constr )) && zpos[7].size()==0)
    ||  ((zpos[7].size()==1|| (zpos[7].size()==2 && abs(xpos[7][0]-xpos[7][1])<constr )) && zpos[6].size()==0)) ){
	accepted_tracks.push_back(j);

	for(int i=0;i<zpos.size();i++){
		for(int j=0; j<zpos[i].size();j++){
			event->Fill(zpos[i][0],xpos[i][0]); //fill event histogram with data
			}
			
		}




   	TF1 *fit = new TF1("fit","[1]*x+[0]",-200,1200);
	event->Fit(fit,"Q"); //fit single event
	double intercept = fit->GetParameter(0); //get parameters
	double slope = fit->GetParameter(1);

	for(int i = 0; i<zpos.size();i++){
		if(xpos[i].size() !=0 && zpos[i].size() !=0){
			double delta_x = (xpos[i][0]-intercept-slope*zpos[i][0]);//calculate difference between fitted and actual data
			residual_array[i]->Fill(delta_x*sqrt(12)/ssdPitch);
			}
		}

	event->Reset("ICESM");
        	
	}
  
  }
   delete event;

   //Draw histograms of Residuals
   TCanvas *ce4 = new TCanvas("ce4","ce4",1800,700);
   ce4->Divide(4,2);
   for(int i=0;i<8;i++){
        new_residuals.push_back(residuals[i]+(residual_array[i]->GetMean())*ssdPitch/sqrt(12)); //Fill vector with offsets	

        residual_array[i]->GetYaxis()->SetTitle("counts");
        residual_array[i]->GetXaxis()->SetTitle("residual");
        ce4->cd(i+1);
   	residual_array[i]->Draw();
	}
   
 
  return new_residuals; 
}


//.........................................................................


void xhist(vector<double> residuals){

   //Input: vector with ssd offsets
   ///Plots a histogram of the beam profile

   //Read in data from SSD.root tree 
   TFile *f = new TFile(infile);
   TTree *t1 = (TTree*)f->Get("SSDtree");
   int fer, module, row;
   t1->SetBranchAddress("fer",&fer);
   t1->SetBranchAddress("module",&module);
   t1->SetBranchAddress("row",&row);
   auto h1 = new TH2F("h1",Form("Spill %d: xz-Hits",spill),150,-200,1200,40,-45,45);
   

//loops over single event and fills vector
for(int i=0; i<t1->GetEntries(); i++){
	t1->GetEntry(i);
      	
	//Filling vectors for xz-plot
	if(fer ==0 && module==1){
		h1->Fill(station0_zpos,(row-320)*ssdPitch-residuals[0]);
	}else if(fer ==0 && module==3){
		h1->Fill(station0_zpos+station_thickness,(row-320)*ssdPitch-residuals[1]);
	}else if(fer ==1 && module==1){
		h1->Fill(station1_zpos,(row-320)*ssdPitch-residuals[2]);
	}else if(fer ==1 && module==3){
		h1->Fill(station1_zpos+station_thickness,(row-320)*ssdPitch-residuals[3]);	
	}else if(fer ==2 && module==0){
		h1->Fill(station2_zpos,(row)*ssdPitch-residuals[4]);		
	}else if(fer ==2 && module==1){
		h1->Fill(station2_zpos,-(row)*ssdPitch-residuals[5]);		
	}else if(fer ==3 && module==0){
		h1->Fill(station2_zpos+station_thickness,(row-640)*ssdPitch-residuals[6]);		
	}else if(fer ==3 && module==1){
		h1->Fill(station2_zpos+station_thickness,-(row-640)*ssdPitch-residuals[7]);
		
	}	
	
   }
  
   
   //cout<<"hello";
   h1->SetBarWidth(11);
   h1->SetFillStyle(0);
   h1->SetFillColor(kGray);
   h1->SetLineColor(kBlue);
   h1->GetYaxis()->SetTitle("x-position mm");
   h1->GetXaxis()->SetTitle("z-position mm");
   h1->SetStats(0);
   h1->Draw("violiny(112000000)");

}


//.........................................................................

void xplot(vector<double> residuals, int event){

   //Input: vector of offsets and starting event
   //Plots xz-positions for six events

   vector<vector<double>> xpos = {{},{},{},{},{},{},{},{}};
   vector<vector<double>> zpos = {{},{},{},{},{},{},{},{}};//store x,z positions
   

   //Read in data from SSD.root tree 
   TFile *f = new TFile(infile);
   TTree *t1 = (TTree*)f->Get("SSDtree");
   int fer, module, row;
   t1->SetBranchAddress("fer",&fer);
   t1->SetBranchAddress("module",&module);
   t1->SetBranchAddress("row",&row);
   int start = -1;
   //loops over a single track and fills vectors
   int numgraphs = 6;
   for(int j = 0; j<numgraphs; j++){
           start ++;
	   for(int i=1+bounds[event+j];i<=bounds[event+j+1];i++){
		t1->GetEntry(i);
      	
		//Filling vectors for xz-plot
		if(fer ==0 && module==1){
			xpos[start].push_back((row-320)*ssdPitch-residuals[0]);
			zpos[start].push_back(station0_zpos);
		}else if(fer ==0 && module==3){
			xpos[start].push_back((row-320)*ssdPitch-residuals[1]);
			zpos[start].push_back(station0_zpos+station_thickness);
		}else if(fer ==1 && module==1){
			xpos[start].push_back((row-320)*ssdPitch-residuals[2]);
			zpos[start].push_back(station1_zpos);
		}else if(fer ==1 && module==3){
			xpos[start].push_back((row-320)*ssdPitch-residuals[3]);
			zpos[start].push_back(station1_zpos+station_thickness);
		}else if(fer ==2 && module==0){
			xpos[start].push_back((row)*ssdPitch-residuals[4]);
			zpos[start].push_back(station2_zpos);
		}else if(fer ==2 && module==1){
			xpos[start].push_back(-(row)*ssdPitch-residuals[5]);
			zpos[start].push_back(station2_zpos);
		}else if(fer ==3 && module==0){
			xpos[start].push_back((row-640)*ssdPitch-residuals[6]);
			zpos[start].push_back(station2_zpos+station_thickness);
		}else if(fer ==3 && module==1){
			xpos[start].push_back(-(row-640)*ssdPitch-residuals[7]);
			zpos[start].push_back(station2_zpos+station_thickness);
		}	
	
	   }
	}
       TGraph** xzPlot_array = new TGraph*[numgraphs]; 
       for (int i=0;i<numgraphs;i++) { 
	     if(xpos[i].size()!=0){
		     Double_t xlist[xpos[i].size()], zlist[xpos[i].size()];
		     for(int j = 0; j<xpos[i].size(); j++){
			   xlist[j] = xpos[i][j];
			   zlist[j] = zpos[i][j];
			   }

		     xzPlot_array[i] = new TGraph(xpos[i].size(),zlist,xlist);
	     	     }
             }
   TCanvas *ce4 = new TCanvas("ce4","ce4",1800,700);
   ce4->Divide(3,2);
   for(int i=0;i<numgraphs;i++){
	if(xpos[i].size()!=0){
	        xzPlot_array[i]->SetTitle(Form("Event %d xz-Hits",event+i));
        	xzPlot_array[i]->GetYaxis()->SetTitle("x-position mm");
        	xzPlot_array[i]->GetXaxis()->SetTitle("z-position mm");
		xzPlot_array[i]->GetYaxis()->SetRangeUser(-45,45);
        	xzPlot_array[i]->GetXaxis()->SetRangeUser(-200,1200);
        	ce4->cd(i+1);
   		xzPlot_array[i]->Draw("AP*");
		TF1 *fit = new TF1("fit","[0]*x+[1]",-200,1200);
		xzPlot_array[i]->Fit(fit,"Q");
		}
	}
   
}


//.........................................................................

double with_fit(vector<double>residuals, int event){

    //Input: vector of offsets and starting event
    // Returns chi squared of a single track

   vector<double> xpos;
   vector<double> zpos;//store x,z positions
   

   //Read in data from SSD.root tree 
   TFile *f = new TFile(infile);
   TTree *t1 = (TTree*)f->Get("SSDtree");
   int fer, module, row;
   t1->SetBranchAddress("fer",&fer);
   t1->SetBranchAddress("module",&module);
   t1->SetBranchAddress("row",&row);
   
   //loops over a single track and fills vectors
   for(int i=1+bounds[event];i<=bounds[event+1];i++){
	t1->GetEntry(i);
      	
	//Filling vectors for xz-plot
	if(fer ==0 && module==1){
		xpos.push_back((row-320)*ssdPitch-residuals[0]);
		zpos.push_back(station0_zpos);		
	}else if(fer ==0 && module==3){
		xpos.push_back((row-320)*ssdPitch-residuals[1]);
		zpos.push_back(station0_zpos+station_thickness);		
	}else if(fer ==1 && module==1){
		xpos.push_back((row-320)*ssdPitch-residuals[2]);
		zpos.push_back(station1_zpos);		
	}else if(fer ==1 && module==3){
		xpos.push_back((row-320)*ssdPitch-residuals[3]);
		zpos.push_back(station1_zpos+station_thickness);		
	}else if(fer ==2 && module==0){
		xpos.push_back((row)*ssdPitch-residuals[4]);
		zpos.push_back(station2_zpos);
	}else if(fer ==2 && module==1){
		xpos.push_back(-(row)*ssdPitch-residuals[4]);
		zpos.push_back(station2_zpos);		
	}else if(fer ==3 && module==0){
		xpos.push_back((row-640)*ssdPitch-residuals[5]);
		zpos.push_back(station2_zpos+station_thickness);		
	}else if(fer ==3 && module==1){
		xpos.push_back(-(row-640)*ssdPitch-residuals[5]);
		zpos.push_back(station2_zpos+station_thickness);
		
	}	
	


   }
        if(xpos.size()==0) return 0;
        //Plotting xz-graph
        Double_t xlist[xpos.size()], zlist[xpos.size()];
	for(int i = 0; i<xpos.size(); i++){
		xlist[i] = xpos[i];
		zlist[i] = zpos[i];
	}
	TGraph* xzPlot = new TGraph(xpos.size(),zlist,xlist);
	TF1 *fit = new TF1("fit","[0]*x+[1]",-200,1200);
	xzPlot->Fit(fit,"Q");
	TFitResultPtr r = xzPlot->Fit(fit,"S");
	double chi2 = r->Chi2();
	f->Close();
	return chi2;


	
}


//.........................................................................


vector<vector<double>> alignments;

void xwith_constr(){
   
  //xhist();
  /*
  alignments.push_back(initial_xalign({0,0,0,0,0,0,0,0}));
  
  for(int i=0;i<6;i++){
	alignments.push_back(xalign(alignments[i]));
	}
  */	 
   /*
   xalign({0,0,0,0,0,0,0,0});
   
   TH1* h1 = new TH1I("h1", "#chi^{2} with alignment", 104, 0, 8000);
   TH1* h2 = new TH1I("h2", "#chi^{2} without alignment", 104, 0,8000);
 
   for(int i=0; i<accepted_tracks.size(); i++){
        cout<<i<<", ";
   	double with = with_fit({ -1.2695828, -0.33297059, 0.95427081, 1.6136993, -1.8158913, -2.7851771, 1.6198023, 1.7130082 },accepted_tracks[i]);
	double without =with_fit({0,0,0,0,0,0,0,0},accepted_tracks[i]);
        
	h1->Fill(with/0.0036);
	h2->Fill(without/0.0036);
        
   }
  

  THStack* hstack = new THStack("hstack", "Chi Squared Before and After Alignment");

  
  h1->SetLineColor(kRed);
  h1->SetFillStyle(3354);
  h1->SetLineWidth(2);
  h1->Draw();

  h2->SetLineColor(kBlue);
  h2->SetFillStyle(3354);
  h2->SetLineWidth(2);



  hstack->Add(h1);
  hstack->Add(h2);
  hstack->Draw();
  hstack->GetXaxis()->SetTitle("Chi Squared");
  hstack->GetYaxis()->SetTitle("counts");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  //gPad->SetLogx();
  //{ -0.45951492, 0.27765722, 1.2212789, 1.6770144, -3.0814708, -3.5427020, 0.14455361, 0.68543828} start with average
*/
}










