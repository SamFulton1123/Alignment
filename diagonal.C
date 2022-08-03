#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <string>
#include <algorithm> 
#include <stdio.h>
#include <string.h> 
#include <numeric>
#include "printing.C"


char infile[ ] = "../data/599/SSD1.root";
double station0_zpos = -170;
double station1_zpos = 194;
double station2_zpos = 768;
double station_thickness = 124;
double station2_thickness = 226;
double ssdPitch = 0.06;
double constr = 0.07;

vector<double> xresiduals =  { -0.47521848, -0.20327105, 0.51013523, 0.53933468, -2.0499423, 0.37452182, -2.4160758, -0.21869782 };//300, { -0.48217263, -0.21014961, 0.50731036, 0.53800422, -2.0419823, 0.37994331, -2.4031378, -0.21027787 };
vector<double> yresiduals =  { -1.4594831, -0.24913770, 1.6139949, 2.5456382, 0.52258020, -1.7611825, 5.4433118, 3.2463359 };//300 { -1.4549245, -0.24419160, 1.6180551, 2.5482706, 0.52150012, -1.7596827, 5.4426099, 3.2489963 };
vector<double> xcenters={xresiduals[4]/2+xresiduals[5]/2,xresiduals[6]/2+xresiduals[7]/2};
vector<double> ycenters={yresiduals[4]/2+yresiduals[5]/2,yresiduals[6]/2+yresiduals[7]/2};

//{-1.4853	,   -0.499804	,   0.87247	,   1.57393	,   -0.815158	,   -3.13893	,   3.4796	,   1.23058};


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


vector<int> find_accepted_tracks(){

   //Returns a vector with single particle events
   vector<int> result;
  
   //Read in data from SSD.root tree 
   TFile *f = new TFile(infile);
   TTree *t1 = (TTree*)f->Get("SSDtree");
   int fer, module, row;
   t1->SetBranchAddress("fer",&fer);
   t1->SetBranchAddress("module",&module);
   t1->SetBranchAddress("row",&row);
   

   //loops over single event and fills vector
   for(int j=0; j<bounds.size(); j++){

      vector<vector<double>> ypos = {{},{},{},{},{},{},{},{}};
      vector<vector<double>> yzpos = {{},{},{},{},{},{},{},{}};
      vector<vector<double>> xpos = {{},{},{},{},{},{},{},{}};
      vector<vector<double>> xzpos = {{},{},{},{},{},{},{},{}};
      vector<vector<double>> wpos = {{},{},{},{},{},{}};
      vector<vector<double>> wzpos = {{},{},{},{},{},{}};

      for(int i=bounds[j];i<bounds[j+1];i++){
   	t1->GetEntry(i);

   	//Filling vectors for xz-plot
      if(fer ==0 && module==1){
         xpos[0].push_back((-row+320)*ssdPitch-xresiduals[0]);
         xzpos[0].push_back(station0_zpos);
      }else if(fer ==0 && module==3){
         xpos[1].push_back((-row+320)*ssdPitch-xresiduals[1]);
         xzpos[1].push_back(station0_zpos+station_thickness);
      }else if(fer ==1 && module==1){
         xpos[2].push_back((-row+320)*ssdPitch-xresiduals[2]);
         xzpos[2].push_back(station1_zpos);
      }else if(fer ==1 && module==3){
         xpos[3].push_back((-row+320)*ssdPitch-xresiduals[3]);
         xzpos[3].push_back(station1_zpos+station_thickness);
      }else if(fer ==2 && module==0){
         xpos[4].push_back((row)*ssdPitch-xresiduals[4]);
         xzpos[4].push_back(station2_zpos);
      }else if(fer ==2 && module==1){
         xpos[5].push_back((-row)*ssdPitch-xresiduals[5]);
         xzpos[5].push_back(station2_zpos);
      }else if(fer ==3 && module==0){
         xpos[6].push_back((-row+640)*ssdPitch-xresiduals[6]);
         xzpos[6].push_back(station2_zpos+station2_thickness);
      }else if(fer ==3 && module==1){
         xpos[7].push_back((row-640)*ssdPitch-xresiduals[7]);
         xzpos[7].push_back(station2_zpos+station2_thickness);
      }


   	//Filling vectors for yz-plot
      else if(fer ==0 && module==0){
         ypos[0].push_back((row-320)*ssdPitch-yresiduals[0]);
         yzpos[0].push_back(station0_zpos);
      }else if(fer ==0 && module==2){
         ypos[1].push_back((row-320)*ssdPitch-yresiduals[1]);
         yzpos[1].push_back(station0_zpos+station_thickness);      
      }else if(fer ==1 && module==0){
         ypos[2].push_back((row-320)*ssdPitch-yresiduals[2]);
         yzpos[2].push_back(station1_zpos);      
      }else if(fer ==1 && module==2){
         ypos[3].push_back((row-320)*ssdPitch-yresiduals[3]);
         yzpos[3].push_back(station1_zpos+station_thickness);      
      }else if(fer ==2 && module==2){
         ypos[4].push_back((row-640)*ssdPitch-yresiduals[4]);
         yzpos[4].push_back(station2_zpos);      
      }else if(fer ==2 && module==3){
         ypos[5].push_back(-(row-640)*ssdPitch-yresiduals[5]);
         yzpos[5].push_back(station2_zpos);      
      }else if(fer ==3 && module==2){
         ypos[6].push_back((row-640)*ssdPitch-yresiduals[6]);
         yzpos[6].push_back(station2_zpos+station2_thickness);      
      }else if(fer ==3 && module==3){
         ypos[7].push_back(-(row-640)*ssdPitch-yresiduals[7]);
         yzpos[7].push_back(station2_zpos+station2_thickness);
      } 

      //w positions 
      else if(fer ==1 && module==4){
         wpos[0].push_back((row-320)*ssdPitch);
         wzpos[0].push_back(station1_zpos);
      }else if(fer ==1 && module==5){
         wpos[1].push_back((row-320)*ssdPitch);
         wzpos[1].push_back(station1_zpos+station_thickness);
      }else if(fer ==2 && module==4){
         wpos[2].push_back((row)*ssdPitch);
         wzpos[2].push_back(station2_zpos);
      }else if(fer ==2 && module==5){
         wpos[3].push_back((row)*ssdPitch);
         wzpos[3].push_back(station2_zpos);
      }else if(fer ==3 && module==4){
         wpos[4].push_back((row)*ssdPitch);
         wzpos[4].push_back(station2_zpos+station2_thickness);
      }else if(fer ==3 && module==5){
         wpos[5].push_back((row)*ssdPitch);
         wzpos[5].push_back(station2_zpos+station2_thickness);
      }
   }
  
   //Only keep events with single track
   if(  //ypart
       (  yzpos[0].size()==1 || (yzpos[0].size()==2 && abs(ypos[0][0]-ypos[0][1])<constr ))
    && (  yzpos[1].size()==1 || (yzpos[1].size()==2 && abs(ypos[1][0]-ypos[1][1])<constr ))
    && (  yzpos[2].size()==1 || (yzpos[2].size()==2 && abs(ypos[2][0]-ypos[2][1])<constr ))
    && (  yzpos[3].size()==1 || (yzpos[3].size()==2 && abs(ypos[3][0]-ypos[3][1])<constr ))
    &&((( yzpos[4].size()==1 || (yzpos[4].size()==2 && abs(ypos[4][0]-ypos[4][1])<constr )) && yzpos[5].size()==0)
    ||((  yzpos[5].size()==1 || (yzpos[5].size()==2 && abs(ypos[5][0]-ypos[5][1])<constr )) && yzpos[4].size()==0))
    &&((( yzpos[6].size()==1 || (yzpos[6].size()==2 && abs(ypos[6][0]-ypos[6][1])<constr )) && yzpos[7].size()==0)
    ||((  yzpos[7].size()==1 || (yzpos[7].size()==2 && abs(ypos[7][0]-ypos[7][1])<constr )) && yzpos[6].size()==0))
      
    //xpart
    && (  xzpos[0].size()==1 || (xzpos[0].size()==2 && abs(xpos[0][0]-xpos[0][1])<constr ))
    && (  xzpos[1].size()==1 || (xzpos[1].size()==2 && abs(xpos[1][0]-xpos[1][1])<constr ))
    && (  xzpos[2].size()==1 || (xzpos[2].size()==2 && abs(xpos[2][0]-xpos[2][1])<constr ))
    && (  xzpos[3].size()==1 || (xzpos[3].size()==2 && abs(xpos[3][0]-xpos[3][1])<constr ))
    && (((xzpos[4].size()==1 || (xzpos[4].size()==2 && abs(xpos[4][0]-xpos[4][1])<constr )) && xzpos[5].size()==0)
    || (( xzpos[5].size()==1 || (xzpos[5].size()==2 && abs(xpos[5][0]-xpos[5][1])<constr )) && xzpos[4].size()==0))
    && (((xzpos[6].size()==1 || (xzpos[6].size()==2 && abs(xpos[6][0]-xpos[6][1])<constr )) && xzpos[7].size()==0)
    || ((  xpos[7].size()==0 || xzpos[7].size()==1  || (xzpos[7].size()==2 && abs(xpos[7][0]-xpos[7][1])<constr )) && xzpos[6].size()==0))

    //wpart
    
    &&(  wzpos[0].size()==1 || (wzpos[0].size()==2 && abs(wpos[0][0]-wpos[0][1])<constr ))
    &&(  wzpos[1].size()==1 || (wzpos[1].size()==2 && abs(wpos[1][0]-wpos[1][1])<constr ))
    
    &&(((wzpos[2].size()==1 || (wzpos[2].size()==2 && abs(wpos[2][0]-wpos[2][1])<constr )) && wzpos[3].size()==0)
    ||(( wzpos[3].size()==1 || (wzpos[3].size()==2 && abs(wpos[3][0]-wpos[3][1])<constr )) && wzpos[2].size()==0))
    
    &&(((wzpos[4].size()==1 || (wzpos[4].size()==2 && abs(wpos[4][0]-wpos[4][1])<constr )) && wzpos[5].size()==0)
    ||(( wzpos[5].size()==1 || (wzpos[5].size()==2 && abs(wpos[5][0]-wpos[5][1])<constr )) && wzpos[4].size()==0))
    
    ){
      
      
      // Only take events with reasonably small chi squared
      int xcount = 0;
      int ycount = 0; 
      for(int i=0;i<yzpos.size();i++)  if( ypos[i].size() != 0)  ycount ++; 
      for(int i=0;i<xzpos.size();i++)  if( xpos[i].size() != 0)  xcount ++;  
      cout<<endl<<xcount<<", "<<ycount<<endl;
      
      Double_t xlist[xcount], ylist[ycount], zlist[ycount];

      int numy = 0;
      for(int i=0;i<ypos.size();i++){
         if( yzpos[i].size()!=0 ){
            ylist[numy] = getAverage(ypos[i]);
            zlist[numy] = yzpos[i][0];
            numy++;
         }      
      }
      int numx = 0;
      for(int i=0;i<ypos.size();i++){
         if( xzpos[i].size()!=0 ){
            xlist[numx] = getAverage(xpos[i]);
            numx++;
         }      
      }




      TGraph* yevent = new TGraph(ycount,zlist,ylist);
      TGraph* xevent = new TGraph(xcount,zlist,xlist);

     //drop events when predicted hit is in dead region but we still have a hit
      TF1 *first_x_points = new TF1("first_x_points","[5]*x+[4]",-200,station2_zpos+100);//fit a function to the first five points
      xevent->Fit(first_x_points,"Quiet","R"); //Use range option to only fit first five points
      if( ((first_x_points->Eval(station2_zpos+station2_thickness,0,0)) > (-256*ssdPitch-xresiduals[7])) 
      && ((first_x_points->Eval(station2_zpos+station2_thickness,0,0)) < (-128*ssdPitch-xresiduals[7])) 
      && (xpos[6].size() + xpos[7].size() != 0) ){
         //cout<< accepted_tracks[j]<< endl;
         continue;
         }
     //drop events when predicted hit is outside dead region and we don't have a hit
     if( (((first_x_points->Eval(station2_zpos+station2_thickness,0,0)) < (-256*ssdPitch-xresiduals[7])) 
      || ((first_x_points->Eval(station2_zpos+station2_thickness,0,0)) > (1-28*ssdPitch-xresiduals[7])) )
      && (xpos[6].size() + xpos[7].size() == 0) ){
           //cout<< accepted_tracks[j]<< endl;
           continue;
           }

      TF1 *yfit = new TF1("yfit","[1]*x+[0]",-200,1200);
      TF1 *xfit = new TF1("xfit","[3]*x+[2]",-200,1200);

      TFitResultPtr ry = yevent->Fit(yfit,"S");
      TFitResultPtr rx = xevent->Fit(xfit,"S");

      double chi2y = ry->Chi2();
      double chi2x = rx->Chi2();
      

	   if(chi2y<0.3 && chi2x<0.3)result.push_back(j);

	   }
   }

   return result; 
}

vector<int> accepted_tracks = find_accepted_tracks();



//.........................................................................

vector<double> theta_residuals(){
   
   //Input: vector with SSD offsets
   //Returns a vector with new offsets

   vector<double> Variables = {0 ,  -0.0308491  ,  0.00150159  ,  -0.0273688  ,  0.00356327  ,  -0.0192158  ,  0.00319301  ,  -0.0139853  ,0,0,0,0,0,0,0,0, 3.40833  ,  3.86717 };
   double xtot_shift = Variables[16];
   double ytot_shift = Variables[17];

   double chix_with =0;
   double chiy_with =0;
   double chix_without =0;
   double chiy_without =0;
   
   //Histstograms filled with residuals with and without rotations
   TH1F** residualx_array = new TH1F*[6]; 
   TH1F** residualy_array = new TH1F*[6];
   TH1F** baselinex_array = new TH1F*[6];
   TH1F** baseliney_array = new TH1F*[6];
   double size = 10;
   int bin = 1000;
   baselinex_array[0] = new TH1F("b0","First Station Baseline x",bin,-size,size);
   baselinex_array[1] = new TH1F("b1","Second Station Baseline x",bin,-size,size);
   baselinex_array[2] = new TH1F("b2","Third Station Baseline x",bin,-size,size);
   baselinex_array[3] = new TH1F("b3","Fourth Station Baseline x",bin,-size,size);
   baselinex_array[4] = new TH1F("b4","Fith Station Baseline x",bin,-size,size);
   baselinex_array[5] = new TH1F("b5","Fith Station Baseline x",bin,-size,size);

   baseliney_array[0] = new TH1F("b6","First Station Baseline y",bin,-size,size);
   baseliney_array[1] = new TH1F("b7","Second Station Baseline y",bin,-size,size);
   baseliney_array[2] = new TH1F("b8","Third Station Baseline y",bin,-size,size);
   baseliney_array[3] = new TH1F("b9","Fourth Station Baseline y",bin,-size,size);
   baseliney_array[4] = new TH1F("b10","Fith Station Baseline y",bin,-size,size);
   baseliney_array[5] = new TH1F("b11","Fith Station Baseline y",bin,-size,size);

   residualx_array[0] = new TH1F("r0","First Station x",bin,-size,size);
   residualx_array[1] = new TH1F("r1","Second Station x",bin,-size,size);
   residualx_array[2] = new TH1F("r2","Third Station x",bin,-size,size);
   residualx_array[3] = new TH1F("r3","Fourth Station x",bin,-size,size);
   residualx_array[4] = new TH1F("r4","Fith Station x",bin,-size,size);
   residualx_array[5] = new TH1F("r5","Fith Station x",bin,-size,size);

   residualy_array[0] = new TH1F("r6","First Station y",bin,-size,size);
   residualy_array[1] = new TH1F("r7","Second Station y",bin,-size,size);
   residualy_array[2] = new TH1F("r8","Third Station y",bin,-size,size);
   residualy_array[3] = new TH1F("r9","Fourth Station y",bin,-size,size);
   residualy_array[4] = new TH1F("r10","Fith Station y",bin,-size,size);
   residualy_array[5] = new TH1F("r11","Fith Station y",bin,-size,size);


   //Read in data from SSD.root tree 
   TFile *f = new TFile(infile);
   TTree *t1 = (TTree*)f->Get("SSDtree");
   int fer, module, row;
   t1->SetBranchAddress("fer",&fer);
   t1->SetBranchAddress("module",&module);
   t1->SetBranchAddress("row",&row);
  
   //loops over multiple events and fills vector
   for(int j=0; j<200; j++){

      vector<vector<double>> ypos = {{},{},{},{},{},{}};
      vector<vector<double>> yzpos = {{},{},{},{},{},{}};
      vector<vector<double>> xpos = {{},{},{},{},{},{}};
      vector<vector<double>> xzpos = {{},{},{},{},{},{}};
      vector<vector<double>> wpos = {{},{},{},{},{},{}};
      vector<vector<double>> wzpos = {{},{},{},{},{},{}};


      for(int i=bounds[accepted_tracks[j]]+1;i<=bounds[accepted_tracks[j]+1];i++){
			t1->GetEntry(i);

         //Filling vectors for xz-plot
         if(fer ==0 && module==1){
            xpos[0].push_back((-row+320)*ssdPitch-xresiduals[0]-xtot_shift);
            xzpos[0].push_back(station0_zpos);
         }else if(fer ==0 && module==3){
            xpos[1].push_back((-row+320)*ssdPitch-xresiduals[1]-xtot_shift);
            xzpos[1].push_back(station0_zpos+station_thickness);
         }else if(fer ==1 && module==1){
            xpos[2].push_back((-row+320)*ssdPitch-xresiduals[2]-xtot_shift);
            xzpos[2].push_back(station1_zpos);
         }else if(fer ==1 && module==3){
            xpos[3].push_back((-row+320)*ssdPitch-xresiduals[3]-xtot_shift);
            xzpos[3].push_back(station1_zpos+station_thickness);
         }else if(fer ==2 && module==0){
            xpos[4].push_back((row)*ssdPitch-xresiduals[4]-xtot_shift);
            xzpos[4].push_back(station2_zpos);
         }else if(fer ==2 && module==1){
            xpos[4].push_back((-row)*ssdPitch-xresiduals[5]-xtot_shift);
            xzpos[4].push_back(station2_zpos);
         }else if(fer ==3 && module==0){
            xpos[5].push_back((-row+640)*ssdPitch-xresiduals[6]-xtot_shift);
            xzpos[5].push_back(station2_zpos+station2_thickness);
         }else if(fer ==3 && module==1){
            xpos[5].push_back((row-640)*ssdPitch-xresiduals[7]-xtot_shift);
            xzpos[5].push_back(station2_zpos+station2_thickness);

         //Filling vectors for yz-plot
         }else if(fer ==0 && module==0){
            ypos[0].push_back((row-320)*ssdPitch-yresiduals[0]-ytot_shift);
            yzpos[0].push_back(station0_zpos);
         }else if(fer ==0 && module==2){
            ypos[1].push_back((row-320)*ssdPitch-yresiduals[1]-ytot_shift);
            yzpos[1].push_back(station0_zpos+station_thickness);     
         }else if(fer ==1 && module==0){
            ypos[2].push_back((row-320)*ssdPitch-yresiduals[2]-ytot_shift);
            yzpos[2].push_back(station1_zpos);     
         }else if(fer ==1 && module==2){
            ypos[3].push_back((row-320)*ssdPitch-yresiduals[3]-ytot_shift);
            yzpos[3].push_back(station1_zpos+station_thickness);     
         }else if(fer ==2 && module==2){
            ypos[4].push_back((row-640)*ssdPitch-yresiduals[4]-ytot_shift);
            yzpos[4].push_back(station2_zpos);     
         }else if(fer ==2 && module==3){
            ypos[4].push_back(-(row-640)*ssdPitch-yresiduals[5]-ytot_shift);
            yzpos[4].push_back(station2_zpos);     
         }else if(fer ==3 && module==2){
            ypos[5].push_back((row-640)*ssdPitch-yresiduals[6]-ytot_shift);
            yzpos[5].push_back(station2_zpos+station2_thickness);     
         }else if(fer ==3 && module==3){
            ypos[5].push_back(-(row-640)*ssdPitch-yresiduals[7]-ytot_shift);
            yzpos[5].push_back(station2_zpos+station2_thickness); 
         }

         //w positions 
         else if(fer ==1 && module==4){
            wpos[0].push_back((row-320)*ssdPitch);
            wzpos[0].push_back(station1_zpos);
         }else if(fer ==1 && module==5){
            wpos[1].push_back((row-320)*ssdPitch);
            wzpos[1].push_back(station1_zpos+station_thickness);
         }else if(fer ==2 && module==4){
            wpos[2].push_back((row)*ssdPitch);
            wzpos[2].push_back(station2_zpos);
         }else if(fer ==2 && module==5){
            wpos[3].push_back((row)*ssdPitch);
            wzpos[3].push_back(station2_zpos);
         }else if(fer ==1 && module==5){
            wpos[4].push_back((row)*ssdPitch);
            wzpos[4].push_back(station2_zpos+station2_thickness);
         }else if(fer ==1 && module==5){
            wpos[5].push_back((row)*ssdPitch);
            wzpos[5].push_back(station2_zpos+station2_thickness);
         }

      }
     
  
   
   // get x,y,z positions
   int count = 0; 
   for(int i=0;i<xpos.size();i++)  if(xpos[i].size()!=0 && xzpos[i].size()!=0 && ypos[i].size()!=0 && yzpos[i].size()!=0)  count ++;
   //cout<<endl<<count<<endl;
   Double_t xwithout_rot[count], ywithout_rot[count], zfinal[count],xfinal[count],yfinal[count];
   int num = 0;
   

   //Fill lists to be used for plotting
   for(int i=0;i<xpos.size();i++){
      if(xpos[i].size()!=0 && xzpos[i].size()!=0 && ypos[i].size()!=0 && yzpos[i].size()!=0){
         if(yzpos[i][0]!=xzpos[i][0])continue;
         double x = getAverage(xpos[i]);
         double y = getAverage(ypos[i]);
         xwithout_rot[num] = x;
         ywithout_rot[num] = y;
         xfinal[num] = x;
         yfinal[num] = y;
         zfinal[num] = getAverage(yzpos[i]);
         num++;

      }        
   }
   

   //Rotation of the first 8 SSDs
   //Rotating the first 4 stations
   for(int i=0;i<4;i++){
      double xc = 0 - xresiduals[i]-xtot_shift;
      double yc = 0 - yresiduals[i]-ytot_shift;
      double x_angle = Variables[2*i];
      double y_angle = Variables[2*i+1];
      xfinal[i] = (xc + (xwithout_rot[i]-xc)/cos(x_angle) - yc*tan(x_angle) - ((ywithout_rot[i]-yc)/cos(y_angle))*tan(x_angle))/(1+tan(x_angle)*tan(y_angle));
      yfinal[i] = (yc + (ywithout_rot[i]-yc)/cos(y_angle) + xc*tan(y_angle) + ((xwithout_rot[i]-xc)/cos(x_angle))*tan(y_angle))/(1+tan(x_angle)*tan(y_angle));
   }  


   //Rotate SSDs in the last 2 stations
   //Only rotate last station if there are a total of 6 hits 
   int last_station = 0;
   if(count==6)last_station++;
   for(int i=0;i<1+last_station;i++){
      if(xpos[i+4][0] > xcenters[i]-xtot_shift && ypos[i+4][0] > ycenters[i]-ytot_shift){
         //first quandrant
         double xc = 19.2 - xresiduals[5+2*i]-xtot_shift;
         double yc = 19.2 - yresiduals[5+2*i]-ytot_shift;
         double x_angle = Variables[8+4*i];
         double y_angle = Variables[9+4*i];
         xfinal[i+4] = (xc + (xwithout_rot[i+4]-xc)/cos(x_angle) - yc*tan(x_angle) - ((ywithout_rot[i+4]-yc)/cos(y_angle))*tan(x_angle))/(1+tan(x_angle)*tan(y_angle));
         yfinal[i+4] = (yc + (ywithout_rot[i+4]-yc)/cos(y_angle) + xc*tan(y_angle) + ((xwithout_rot[i+4]-xc)/cos(x_angle))*tan(y_angle))/(1+tan(x_angle)*tan(y_angle));
     
         
      }else if(xpos[i+4][0] < xcenters[i]-xtot_shift && ypos[i+4][0] > ycenters[i]-ytot_shift){
         //second quandrant
         double xc = -19.2 - xresiduals[4+2*i]-xtot_shift;
         double yc = 19.2 - yresiduals[5+2*i]-ytot_shift;
         double x_angle = Variables[10+4*i];
         double y_angle = Variables[9+4*i];
         xfinal[i+4] = (xc + (xwithout_rot[i+4]-xc)/cos(x_angle) - yc*tan(x_angle) - ((ywithout_rot[i+4]-yc)/cos(y_angle))*tan(x_angle))/(1+tan(x_angle)*tan(y_angle));
         yfinal[i+4] = (yc + (ywithout_rot[i+4]-yc)/cos(y_angle) + xc*tan(y_angle) + ((xwithout_rot[i+4]-xc)/cos(x_angle))*tan(y_angle))/(1+tan(x_angle)*tan(y_angle));
     
      }else if(xpos[i+4][0] < xcenters[i]-xtot_shift && ypos[i+4][0] < ycenters[i]-ytot_shift){
         //third quandrant
         double xc = -19.2 - xresiduals[4+2*i]-xtot_shift;
         double yc = -19.2 - yresiduals[4+2*i]-ytot_shift;
         double x_angle = Variables[10+4*i];
         double y_angle = Variables[11+4*i];
         xfinal[i+4] = (xc + (xwithout_rot[i+4]-xc)/cos(x_angle) - yc*tan(x_angle) - ((ywithout_rot[i+4]-yc)/cos(y_angle))*tan(x_angle))/(1+tan(x_angle)*tan(y_angle));
         yfinal[i+4] = (yc + (ywithout_rot[i+4]-yc)/cos(y_angle) + xc*tan(y_angle) + ((xwithout_rot[i+4]-xc)/cos(x_angle))*tan(y_angle))/(1+tan(x_angle)*tan(y_angle));
          
      }else if(xpos[i+4][0] > xcenters[i]-xtot_shift && ypos[i+4][0] < ycenters[i]-ytot_shift){
         //fourth quandrant
         double xc =  19.2 - xresiduals[5+2*i]-xtot_shift;
         double yc = -19.2 - yresiduals[4+2*i]-ytot_shift;
         double x_angle = Variables[8+4*i];
         double y_angle = Variables[11+4*i];
         xfinal[i+4] = (xc + (xwithout_rot[i+4]-xc)/cos(x_angle) - yc*tan(x_angle) - ((ywithout_rot[i+4]-yc)/cos(y_angle))*tan(x_angle))/(1+tan(x_angle)*tan(y_angle));
         yfinal[i+4] = (yc + (ywithout_rot[i+4]-yc)/cos(y_angle) + xc*tan(y_angle) + ((xwithout_rot[i+4]-xc)/cos(x_angle))*tan(y_angle))/(1+tan(x_angle)*tan(y_angle));
     
      }
   }


   TGraph* x_event = new TGraph(count,zfinal,xfinal);
	TGraph* y_event = new TGraph(count,zfinal,yfinal);
   TGraph* xwithout_rot_event = new TGraph(count,zfinal,xwithout_rot);
   TGraph* ywithout_rot_event = new TGraph(count,zfinal,ywithout_rot);
	

   TF1 *x_fit = new TF1("x_fit","[1]*x+[0]",-200,1200);
   TF1 *y_fit = new TF1("y_fit","[3]*x+[2]",-200,1200);

   TF1 *xwithout_rot_fit = new TF1("xwithout_rot_fit","[5]*x+[4]",-200,1200);
   TF1 *ywithout_rot_fit = new TF1("ywithout_rot_fit","[7]*x+[6]",-200,1200);

   x_event->Fit(x_fit,"Q");
   y_event->Fit(y_fit,"Q"); //fit single event

   xwithout_rot_event->Fit(xwithout_rot_fit,"Q"); 
   ywithout_rot_event->Fit(ywithout_rot_fit,"Q"); 


   int counting=0;
   for(int i = 0; i<xzpos.size();i++){
        if(xpos[i].size()!=0 && xzpos[i].size()!=0 && ypos[i].size()!=0 && yzpos[i].size()!=0){
		     double x_exp = x_fit->Eval(zfinal[counting],0,0);
           double y_exp = y_fit->Eval(zfinal[counting],0,0);

           double xwithout_rot_exp = xwithout_rot_fit->Eval(zfinal[counting],0,0);
		     double ywithout_rot_exp = ywithout_rot_fit->Eval(zfinal[counting],0,0);

           double x_obs = xfinal[counting];
           double y_obs = yfinal[counting];

		     double xwithout_rot_obs = xwithout_rot[counting];
           double ywithout_rot_obs = ywithout_rot[counting];

           counting ++;

           residualx_array[i]->Fill((x_exp-x_obs)*sqrt(12)/0.06);
	        residualy_array[i]->Fill((y_exp-y_obs)*sqrt(12)/0.06);

           baselinex_array[i]->Fill((xwithout_rot_exp-xwithout_rot_obs)*sqrt(12)/0.06);
	        baseliney_array[i]->Fill((ywithout_rot_exp-ywithout_rot_obs)*sqrt(12)/0.06);

           chix_with += pow((x_exp-x_obs)*sqrt(12)/0.06,2);
           chiy_with += pow((y_exp-y_obs)*sqrt(12)/0.06,2);
           chix_without += pow((xwithout_rot_exp-xwithout_rot_obs)*sqrt(12)/0.06,2);
           chiy_without += pow((ywithout_rot_exp-ywithout_rot_obs)*sqrt(12)/0.06,2);

	        }

       }

   }
   
   //All plotting and statistics for x
   TCanvas *ce4 = new TCanvas("ce4","ce4",1400,700);
   ce4->Divide(3,4);
   
   cout<<endl<<"x with rotations"<<endl;
   for(int i=0;i<6;i++){
      cout<<residualx_array[i]->GetMean()<<"   "<<residualx_array[i]->GetStdDev()<<endl;
      ce4->cd(i+1);
      residualx_array[i]->GetYaxis()->SetTitle("counts");
      residualx_array[i]->GetXaxis()->SetTitle("residual");
   	residualx_array[i]->Draw();
	}
   cout<<endl<<"x without rotations"<<endl;
   for(int i=0;i<6;i++){
      cout<<baselinex_array[i]->GetMean()<<"   "<<baselinex_array[i]->GetStdDev()<<endl;
      ce4->cd(i+7);
      baselinex_array[i]->GetYaxis()->SetTitle("counts");
      baselinex_array[i]->GetXaxis()->SetTitle("residual");
      baselinex_array[i]->Draw();
   }


   //All plotting and statistics for y
   TCanvas *c1 = new TCanvas("c1","c1",1400,700);
   c1->Divide(3,4);
   cout<<endl<<"y with rotations"<<endl;
   for(int i=0;i<6;i++){
      cout<<residualy_array[i]->GetMean()<<"   "<<residualy_array[i]->GetStdDev()<<endl;
      c1->cd(i+1);
      residualy_array[i]->GetYaxis()->SetTitle("counts");
      residualy_array[i]->GetXaxis()->SetTitle("residual");
      residualy_array[i]->Draw();
   }
   cout<<endl<<"y without rotations"<<endl;
   for(int i=0;i<6;i++){
      cout<<baseliney_array[i]->GetMean()<<"   "<<baseliney_array[i]->GetStdDev()<<endl;
      c1->cd(i+7);
      baseliney_array[i]->GetYaxis()->SetTitle("counts");
      baseliney_array[i]->GetXaxis()->SetTitle("residual");
      baseliney_array[i]->Draw();
   }

  cout<<endl;
  cout<<" x without rotations: "<<chix_without<<endl;
  cout<<" x with rotations: "<<chix_with<<endl;
  cout<<" y without rotations: "<<chiy_without<<endl;
  cout<<" y with rotations: "<<chiy_with<<endl;


  return {0,0,0,0};//new_residuals; 
}



//.........................................................................


void histograms(){
   
   //Input: vector with SSD offsets
   //Returns a vector with new offsets

   
   //Read in data from SSD.root tree 
   TFile *f = new TFile(infile);
   TTree *t1 = (TTree*)f->Get("SSDtree");
   int fer, module, row;
   t1->SetBranchAddress("fer",&fer);
   t1->SetBranchAddress("module",&module);
   t1->SetBranchAddress("row",&row);


   double size = 40;
   int num_bins = floor(size/0.06);
   TH1F** hits_array = new TH1F*[6]; 

   hits_array[0] = new TH1F("h2","Third Station W-SSD",num_bins,-size,size);
   hits_array[1] = new TH1F("h3","Fourth Station W-SSD",num_bins,-size,size);
   hits_array[2] = new TH1F("h4","Fith Station W-SSD",num_bins,-size,size);
   hits_array[3] = new TH1F("h5","Fith Station W-SSD",num_bins,-size,size);
   //hits_array[4] = new TH1F("h6","Sixth Station Upper W-SSD",num_bins,-size,size);
   //hits_array[5] = new TH1F("h7","Sixth Station Lower W-SSD",num_bins,-size,size);
  
   for(int j=0; j<bounds.size(); j++){
      for(int i=bounds[j]+1;i<=bounds[j+1];i++){
         t1->GetEntry(i);

         //w positions 
         if(fer ==1 && module==4){
            hits_array[0]->Fill((row-320)*ssdPitch);  

         }else if(fer ==1 && module==5){
            hits_array[1]->Fill((row-320)*ssdPitch); 

         }else if(fer ==2 && module==4){
            hits_array[2]->Fill((-row+640)*ssdPitch); 

         }else if(fer ==2 && module==5){
            hits_array[2]->Fill((row-640)*ssdPitch); 

         }else if(fer ==3 && module==4){
            hits_array[3]->Fill((-row+640)*ssdPitch);

         }else if(fer ==3 && module==5){
            hits_array[3]->Fill((row-640)*ssdPitch);  
         }
         
      }
      
   }

   
  TCanvas *ce4 = new TCanvas("ce4","ce4",1300,1300);
   ce4->Divide(2,2);
   for(int i=0;i<4;i++){
        hits_array[i]->GetYaxis()->SetTitle("counts");
        hits_array[i]->GetXaxis()->SetTitle("diagonal position");
        ce4->cd(i+1);
        hits_array[i]->Draw();
   }


}
//.........................................................................

vector<double> beam_profile(){
   
   //Input: vector with SSD offsets
   //Returns a vector with new offsets

   //beam profiles
   TH2F** beam_array = new TH2F*[6];
   double sizes = 30;
   int bin = 100;
   beam_array[0] = new TH2F("h7","First Station",bin,-sizes,sizes,bin,-sizes,sizes);
   beam_array[1] = new TH2F("h8","Second Station",bin,-sizes,sizes,bin,-sizes,sizes);
   beam_array[2] = new TH2F("h9","Third Station",bin,-sizes,sizes,bin,-sizes,sizes);
   beam_array[3] = new TH2F("h10","Fourth Station",bin,-sizes,sizes,bin,-sizes,sizes);
   beam_array[4] = new TH2F("h11","Fith Station",bin,-sizes,sizes,bin,-sizes,sizes);
   beam_array[5] = new TH2F("h12","Fith Station",bin,-sizes,sizes,bin,-sizes,sizes);


   //Read in data from SSD.root tree 
   TFile *f = new TFile(infile);
   TTree *t1 = (TTree*)f->Get("SSDtree");
   int fer, module, row;
   t1->SetBranchAddress("fer",&fer);
   t1->SetBranchAddress("module",&module);
   t1->SetBranchAddress("row",&row);
   vector<double> new_residuals;
   //Histstograms filled with residuals
   TH1F** residual_array = new TH1F*[6]; 
   double size = 0.3;
   int bins = 100;
   residual_array[0] = new TH1F("h0","First Station",bins,-size,size);
   residual_array[1] = new TH1F("h1","Second Station",bins,-size,size);
   residual_array[2] = new TH1F("h2","Third Station",bins,-size,size);
   residual_array[3] = new TH1F("h3","Fourth Station",bins,-size,size);
   residual_array[4] = new TH1F("h4","Fith Station",bins,-size,size);
   residual_array[5] = new TH1F("h5","Fith Station",bins,-size,size);
   //residual_array[6] = new TH1F("h6","Sixth Station Lower X-SSD",bins,-size,size);
   //residual_array[7] = new TH1F("h7","Sixth Station Upper X-SSD",bins,-size,size);

   //loops over single event and fills vector
   for(int j=0; j<bounds.size(); j++){

      vector<vector<double>> ypos = {{},{},{},{},{},{}};
      vector<vector<double>> yzpos = {{},{},{},{},{},{}};
      vector<vector<double>> xpos = {{},{},{},{},{},{}};
      vector<vector<double>> xzpos = {{},{},{},{},{},{}};


      for(int i=bounds[j]+1;i<=bounds[j+1];i++){
      	t1->GetEntry(i);

      	//Filling vectors for xz-plot
         if(fer ==0 && module==1){
            xpos[0].push_back((-row+320)*ssdPitch-xresiduals[0]);
            xzpos[0].push_back(station0_zpos);
         }else if(fer ==0 && module==3){
            xpos[1].push_back((-row+320)*ssdPitch-xresiduals[1]);
            xzpos[1].push_back(station0_zpos+station_thickness);
         }else if(fer ==1 && module==1){
            xpos[2].push_back((-row+320)*ssdPitch-xresiduals[2]);
            xzpos[2].push_back(station1_zpos);
         }else if(fer ==1 && module==3){
            xpos[3].push_back((-row+320)*ssdPitch-xresiduals[3]);
            xzpos[3].push_back(station1_zpos+station_thickness);
         }else if(fer ==2 && module==0){
            xpos[4].push_back((row)*ssdPitch-xresiduals[4]);
            xzpos[4].push_back(station2_zpos);
         }else if(fer ==2 && module==1){
            xpos[4].push_back((-row)*ssdPitch-xresiduals[5]);
            xzpos[4].push_back(station2_zpos);
         }else if(fer ==3 && module==0){
            xpos[5].push_back((-row+640)*ssdPitch-xresiduals[6]);
            xzpos[5].push_back(station2_zpos+station2_thickness);
         }else if(fer ==3 && module==1){
            xpos[5].push_back((row-640)*ssdPitch-xresiduals[7]);
            xzpos[5].push_back(station2_zpos+station2_thickness);  

      	//Filling vectors for yz-plot
      	}else if(fer ==0 && module==0){
      		ypos[0].push_back((row-320)*ssdPitch-yresiduals[0]);
      		yzpos[0].push_back(station0_zpos);
      	}else if(fer ==0 && module==2){
      		ypos[1].push_back((row-320)*ssdPitch-yresiduals[1]);
      		yzpos[1].push_back(station0_zpos+station_thickness);		
      	}else if(fer ==1 && module==0){
      		ypos[2].push_back((row-320)*ssdPitch-yresiduals[2]);
      		yzpos[2].push_back(station1_zpos);		
      	}else if(fer ==1 && module==2){
      		ypos[3].push_back((row-320)*ssdPitch-yresiduals[3]);
      		yzpos[3].push_back(station1_zpos+station_thickness);		
      	}else if(fer ==2 && module==2){
      		ypos[4].push_back((row-640)*ssdPitch-yresiduals[4]);
      		yzpos[4].push_back(station2_zpos);		
      	}else if(fer ==2 && module==3){
      		ypos[4].push_back(-(row-640)*ssdPitch-yresiduals[5]);
      		yzpos[4].push_back(station2_zpos);		
      	}else if(fer ==3 && module==2){
      		ypos[5].push_back((row-640)*ssdPitch-yresiduals[6]);
      		yzpos[5].push_back(station2_zpos+station2_thickness);		
      	}else if(fer ==3 && module==3){
      		ypos[5].push_back(-(row-640)*ssdPitch-yresiduals[7]);
      		yzpos[5].push_back(station2_zpos+station2_thickness);
      		
      	}
      	
	
   }
  
  
   // graph of x_event
   int count = 0;	
   
   //for(int i=0;i<xpos.size()-1;i++)  if(xpos[i].size()!=0 && xzpos[i].size()!=0 && ypos[i].size()!=0 && yzpos[i].size()!=0)  count ++;
   //Double_t ylist[count], zlist[count],zlist[count];
   
   for(int i=0;i<xzpos.size();i++){
      if(xpos[i].size()!=0 && xzpos[i].size()!=0 && ypos[i].size()!=0 && yzpos[i].size()!=0){
   		beam_array[i]->Fill(xpos[i][0],ypos[i][0]);
   	   //num++;
   		}			
      }
   

   /*
	TGraph* x_event = new TGraph(count,zlist,xlist);
   TGraph* y_event = new TGraph(count,zlist,ylist);
	

   TF1 *x_fit = new TF1("x_fit","[1]*x+[0]",-200,1200);
   TF1 *y_fit = new TF1("y_fit","[1]*x+[0]",-200,1200);
   x_event->Fit(x_fit,"Q"); //fit single event
   y_event->Fit(y_fit,"Q"); 
   
   int counting=0;
   for(int i = 0; i<xzpos.size();i++){
        if(xpos[i].size()!=0 && xzpos[i].size()!=0 && ypos[i].size()!=0 && yzpos[i].size()!=0){
		     //Double_t x_exp = x_fit->Eval(zlist[counting],0,0);
		     //Double_t y_exp = y_fit->Eval(zlist[counting],0,0);

	        double x_obs = xlist[counting];
		     double y_obs = ylist[counting];

		     double theta_exp = angle(x_exp,y_exp);
		     double theta_obs =angle(x_obs,y_obs);
	        
	        counting ++;
		     residual_array[i]->Fill(theta_obs-theta_exp);
		     }

      }

      */

   }
   
   
   TCanvas *c1 = new TCanvas("c1","c1",1400,700);
   c1->Divide(3,2);
   for(int i=0;i<6;i++){
        //new_residuals.push_back(residuals[i]+(residual_array[i]->GetMean())*ssdPitch/sqrt(12)); //Fill vector with offsets	
        beam_array[i]->GetYaxis()->SetTitle("y-position mm");
        beam_array[i]->GetXaxis()->SetTitle("x-position mm");
        c1->cd(i+1);
   	  beam_array[i]->Draw("COLZ");
	}
  return {0,0,0,0};//new_residuals; 
}


//.........................................................................


vector<double> ws(int event){
   
   //Input: vector with SSD offsets
   //Returns a vector with new offsets

   
   //Read in data from SSD.root tree 
   TFile *f = new TFile(infile);
   TTree *t1 = (TTree*)f->Get("SSDtree");
   int fer, module, row;
   t1->SetBranchAddress("fer",&fer);
   t1->SetBranchAddress("module",&module);
   t1->SetBranchAddress("row",&row);
  
  
   vector<vector<double>> ypos = {{},{},{},{},{},{}};
   vector<vector<double>> yzpos = {{},{},{},{},{},{}};
   vector<vector<double>> xpos = {{},{},{},{},{},{}};
   vector<vector<double>> xzpos = {{},{},{},{},{},{}};
   vector<vector<double>> wpos = {{},{},{},{},{},{}};
   vector<vector<double>> wzpos = {{},{},{},{},{},{}};

   
   for(int i=bounds[accepted_tracks[event]]+1;i<=bounds[accepted_tracks[event]+1];i++){
      t1->GetEntry(i);

      //Filling vectors for xz-plot
      if(fer ==0 && module==1){
         xpos[0].push_back((-row+320)*ssdPitch-xresiduals[0]);
         xzpos[0].push_back(station0_zpos);
      }else if(fer ==0 && module==3){
         xpos[1].push_back((-row+320)*ssdPitch-xresiduals[1]);
         xzpos[1].push_back(station0_zpos+station_thickness);
      }else if(fer ==1 && module==1){
         xpos[2].push_back((-row+320)*ssdPitch-xresiduals[2]);
         xzpos[2].push_back(station1_zpos);
      }else if(fer ==1 && module==3){
         xpos[3].push_back((-row+320)*ssdPitch-xresiduals[3]);
         xzpos[3].push_back(station1_zpos+station_thickness);
      }else if(fer ==2 && module==0){
         xpos[4].push_back((row)*ssdPitch-xresiduals[4]);
         xzpos[4].push_back(station2_zpos);
      }else if(fer ==2 && module==1){
         xpos[4].push_back((-row)*ssdPitch-xresiduals[5]);
         xzpos[4].push_back(station2_zpos);
      }else if(fer ==3 && module==0){
         xpos[5].push_back((-row+640)*ssdPitch-xresiduals[6]);
         xzpos[5].push_back(station2_zpos+station2_thickness);
      }else if(fer ==3 && module==1){
         xpos[5].push_back((row-640)*ssdPitch-xresiduals[7]);
         xzpos[5].push_back(station2_zpos+station2_thickness);

      //Filling vectors for yz-plot
      }else if(fer ==0 && module==0){
         ypos[0].push_back((row-320)*ssdPitch-yresiduals[0]);
         yzpos[0].push_back(station0_zpos);
      }else if(fer ==0 && module==2){
         ypos[1].push_back((row-320)*ssdPitch-yresiduals[1]);
         yzpos[1].push_back(station0_zpos+station_thickness);     
      }else if(fer ==1 && module==0){
         ypos[2].push_back((row-320)*ssdPitch-yresiduals[2]);
         yzpos[2].push_back(station1_zpos);     
      }else if(fer ==1 && module==2){
         ypos[3].push_back((row-320)*ssdPitch-yresiduals[3]);
         yzpos[3].push_back(station1_zpos+station_thickness);     
      }else if(fer ==2 && module==2){
         ypos[4].push_back((row-640)*ssdPitch-yresiduals[4]);
         yzpos[4].push_back(station2_zpos);     
      }else if(fer ==2 && module==3){
         ypos[4].push_back(-(row-640)*ssdPitch-yresiduals[5]);
         yzpos[4].push_back(station2_zpos);     
      }else if(fer ==3 && module==2){
         ypos[5].push_back((row-640)*ssdPitch-yresiduals[6]);
         yzpos[5].push_back(station2_zpos+station2_thickness);     
      }else if(fer ==3 && module==3){
         ypos[5].push_back(-(row-640)*ssdPitch-yresiduals[7]);
         yzpos[5].push_back(station2_zpos+station2_thickness); 
      }
      
      //w positions 
      else if(fer ==1 && module==4){
         wpos[0].push_back((row-320)*ssdPitch);
         wzpos[0].push_back(station1_zpos);
      }else if(fer ==1 && module==5){
         wpos[1].push_back((row-320)*ssdPitch);
         wzpos[1].push_back(station1_zpos+station_thickness);
      }else if(fer ==2 && module==4){
         wpos[2].push_back((-row+640)*ssdPitch);
         wzpos[2].push_back(station2_zpos);
      }else if(fer ==2 && module==5){
         wpos[3].push_back((row-640)*ssdPitch);
         wzpos[3].push_back(station2_zpos);
      }else if(fer ==3 && module==4){
         wpos[4].push_back((-row+640)*ssdPitch);
         wzpos[4].push_back(station2_zpos+station2_thickness);
      }else if(fer ==3 && module==5){
         wpos[5].push_back((row-640)*ssdPitch);
         wzpos[5].push_back(station2_zpos+station2_thickness);
      }
      
   }
  
  
   
   // get x,y,z positions
   int count = 0; 
   for(int i=0;i<xpos.size();i++)  if(xpos[i].size()!=0 &&  ypos[i].size()!=0 && wpos.size() !=0)  count ++;
   Double_t zfinal[count],xfinal[count],yfinal[count];
   int num = 0;
   

   //Fill lists to be used for plotting
   for(int i=0;i<xpos.size();i++){
      if(xpos[i].size()!=0 && xzpos[i].size()!=0 && ypos[i].size()!=0 && yzpos[i].size()!=0){
         if(yzpos[i][0]!=xzpos[i][0])continue;
         double x = getAverage(xpos[i]);
         double y = getAverage(ypos[i]);
         xfinal[num] = x;
         yfinal[num] = y;
         zfinal[num] = getAverage(yzpos[i]);
         num++;

      }        
   }
   

   
   TF1** w_array = new TF1*[4];
   w_array[0] = new TF1("w_station2","x+[0]",-30,30);
   w_array[1] = new TF1("w_station3","x+[0]",-30,30);
   w_array[2] = new TF1("w_station4","-x+[0]",-30,30);
   w_array[3] = new TF1("w_station5","-x+[0]",-30,30);

   num = 0;
   for(int i=0; i<wpos.size(); i++){
      if(wpos[i].size()==0)continue;
      double intercept = sqrt(2)*getAverage(wpos[i]);
      w_array[num]->SetParameter(0,intercept);
      num++;
   }



   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
   c1->Divide(2,2);
   
   TEllipse** ellipse_array = new TEllipse*[4];
   for(int i=0; i<4; i++){
      c1->cd(i+1);
      ellipse_array[i] = new TEllipse(xfinal[i+2],yfinal[i+2],1.,1.);
      ellipse_array[i]->SetFillColor(6);
      w_array[i]->Draw();
      ellipse_array[i]->Draw();
   }



   /*
   //TGraph* y_event = new TGraph(count,zfinal,yfinal);
   //TGraph* xwithout_rot_event = new TGraph(count,zfinal,xwithout_rot);
   //TGraph* ywithout_rot_event = new TGraph(count,zfinal,ywithout_rot);
   

   //All plotting and statistics for x
   /TCanvas *ce4 = new TCanvas("ce4","ce4",1400,700);
   ce4->Divide(3,4);
   
   cout<<endl<<"x with rotations"<<endl;
   for(int i=0;i<6;i++){
      cout<<residualx_array[i]->GetMean()<<"   "<<residualx_array[i]->GetStdDev()<<endl;
      ce4->cd(i+1);
      residualx_array[i]->GetYaxis()->SetTitle("counts");
      residualx_array[i]->GetXaxis()->SetTitle("residual");
      residualx_array[i]->Draw();
   }
   
   */
  

  return {0,0,0,0};//new_residuals; 
}





//.........................................................................



void diagonal(){
 
 //chi();
  
}




/*
for(int i=0;i<xpos.size();i++){
          print(xpos[i]) ;
          cout<<", ";
       }

       cout<<endl;
       for(int i=0;i<ypos.size();i++){
          print(ypos[i]) ;
          cout<<", ";
       }
      cout<<count;
      cout<<endl;   
*/





