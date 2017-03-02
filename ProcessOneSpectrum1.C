#include "smear_test.C"

const int nbOfHistos = 388;
gStyle->SetOptFit(11111111);

TF1* fitFcn;
TF1* fitFcn1;

fstream foutfit("fitoutput.txt",ios::out);

//my pol2 function
Double_t myPol2(Double_t *x, Double_t *par){
    return par[0]*(x[0]-par[1])*(x[0]-par[1])+par[2];
}

void ProcessOneSpectrum1(){
    TVirtualFitter::SetDefaultFitter("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

    TFile* f =  new TFile("Histograms.root","read");

    TFile* rootfile = new TFile("Histograms_smoothed6-Minuit2.root","recreate");
    TH1F* hspec[nbOfHistos];
    TH1F* hspecSmear[nbOfHistos];
    TH1F* hspecDeri[nbOfHistos];
    TH1F* hspecDeri1[nbOfHistos];
    TCanvas* c[nbOfHistos];
    //TF1* pol2fun[nbOfHistos];
    TGraphErrors* gEdge = new TGraphErrors();
    TGraphErrors* gEdgeRatio = new TGraphErrors();
    TGraph* gChi2 = new TGraph();
    TH1F* hChi2Prob = new TH1F("Chi2Probability","",55*4,0,1.1);

    TCanvas* cderi = new TCanvas();

    TH1F* myError = new TH1F("Error","",4096,0,4096);
    TH1F* ComptonEdgeHist = new TH1F("ComptonEdgeHist","",130,650,780);
    ComptonEdgeHist->SetXTitle("Compton edge (Channel number)");
    ComptonEdgeHist->SetYTitle("Counts (bin size 1 )");

    int stepSize_superpose = 16; // 16 means to take 4 hours step.
    TCanvas* cSpecSuperpose = new TCanvas();
    cSpecSuperpose->SetName("RawDistribution");
    TLegend* leg1 = new TLegend(0.7,0.7,0.85,0.85);
    TCanvas* cSpecSuperpose1 = new TCanvas();
    cSpecSuperpose1->SetName("RawDistribution_smoothed");
    TLegend* leg2 = new TLegend(0.7,0.7,0.85,0.85);

    fstream fout("errors.txt",ios::out);
    fstream fin("MeanSmearedSpectrum.txt",ios::in);
    fstream fin1("CovarianceSmeared.txt",ios::in);
    fstream fin2("RawUncertaintySmeared.txt",ios::in);
    double means[4096], cov[4096], rawerr[4096],a;
    int acnt=0;
    while(fin>>a){means[acnt]=a; acnt++;} fin.close(); 
    acnt=0;
    while(fin1>>a){cov[acnt]=a; acnt++; } fin1.close();
    acnt=0;
    while(fin2>>a){rawerr[acnt]=a; acnt++; } fin2.close();

    string filenames[nbOfHistos];
    fstream fin("filenamelist.txt",ios::in);//filenamelist = open("filenamelist.txt");

    string aname;
    int cnt;
    while(fin>>aname) {
        char tempstr[35]; sprintf(tempstr,"_spec%d",cnt);
        filenames[cnt] = aname + tempstr; 
        //cout<<filenames[cnt]<<endl;
        cnt++;
    }
    fin.close();

gStyle->SetOptFit(11111111);
    int cnt_this = 0;
    for(int i=0; i<nbOfHistos;i++){
        string poststr = ";1";
        string thisname = filenames[i] + poststr;
        cout<<thisname<<endl;
        hspec[i] = (TH1F*)f->Get(thisname.c_str());

        filenames[i].replace(8, 4, "");
        string histname = hspec[i]->GetName(); histname.replace(8,4,"");
        hspec[i]->SetName(histname.c_str());
          
        string smearname = "_smeared";

        string thisnamesmear = filenames[i] + smearname;
        string deriname = "_deri";
        string thisnamederi = filenames[i] + deriname;
        string funname = "_pol2fit";
        string deriname1 = "_deri1";
        string thisnamederi1 = filenames[i] + deriname1;

        string thisfun = filenames[i] + funname;
        cout<<filenames[i]<<endl;foutfit<<filenames[i]<<endl;
        //hspecSmear[i]->Sumw2();
        hspecSmear[i] = smear(hspec[i], 20, thisnamesmear);
        
        if(i%stepSize_superpose==0){
             
             cSpecSuperpose1->cd();hspecSmear[i]->Draw("same"); 
             char leglabel2[15]; sprintf(leglabel2,"%d hours (smoothed)",(int)(stepSize_superpose*15/60)*(cnt_this+1));
             hspecSmear[i]->SetLineColor(cnt_this+1);hspecSmear[i]->SetLineWidth(2); leg2->AddEntry(hspecSmear[i],leglabel2,"l"); leg2->Draw();//cSpecSuperpose1->Close();

             cSpecSuperpose->cd(); hspec[i]->Draw("same");
             char leglabel1[15]; sprintf(leglabel1,"%d hours",(int)(stepSize_superpose*15/60)*(cnt_this+1));
             hspec[i]->SetLineColor(cnt_this+1);hspec[i]->SetLineWidth(2); leg1->AddEntry(hspec[i],leglabel1,"l"); leg1->Draw();//cSpecSuperpose->Close();
             cnt_this++;
        }

//        c[i] = new TCanvas(); 
//        c[i]->SetName(filenames[i].c_str());
//        hspec[i]->SetLineWidth(2); hspec[i]->Draw();
//       hspecSmear[i]->SetLineColor(kRed); hspecSmear[i]->Draw("same");
        //char savefilename[100]; sprintf(savefilename,"./spectrum/%s_%d_%d.png",filenames[i].c_str(),i,15.0*i);
        //c[i]->SaveAs(savefilename); 

        //c[i]->Print("animation_raw_spec.gif+30");
//        c[i]->Write(); c[i]->Close();

        int nBins = hspecSmear[i]->GetNbinsX();
        //cout<<nBins<<endl;

        hspecDeri[i] = new TH1F(thisnamederi.c_str(),"",4096,0,4096);
        hspecDeri1[i] = new TH1F(thisnamederi1.c_str(),"",4096,0,4096);

        double temp=0, temp1=0;

        for (int j=0; j<nBins; j++) {
//            double x = h2->GetBinCenter(i);
            double content = hspecSmear[i]->GetBinContent(j+1);
//
            hspecDeri[i]->SetBinContent(j+1, content-temp);

            double theerr = 0;
            
                if (j==0) {theerr = rawerr[j]**2;}
                else {theerr = rawerr[j]**2+rawerr[j-1]**2-2.0*cov[j];}
            
                fout<<theerr<<endl; 
                //if(j>500 && j<1000) cout<<theerr<<"\t"<<sqrt(TMath::Abs(theerr))<<endl;
            hspecDeri[i]->SetBinError(j+1,sqrt(TMath::Abs(theerr)));
            //myError->SetBinContent(j+1,theerr);
            //int cj=0;
            if(j%2 ==1){
                hspecDeri1[i]->SetBinContent(j+1, hspecSmear[i]->GetBinContent(j+2)-hspecSmear[i]->GetBinContent(j+1));
               // hspecDeri1[i]->SetBinContent((j+2), hspec[i]->GetBinContent(j+2)-hspec[i]->GetBinContent(j+1));
                hspecDeri1[i]->SetBinError(j+1,sqrt(TMath::Abs(theerr)));
                //hspecDeri1[i]->SetBinError((j+2),sqrt(abs(theerr)));
                //cj++;
            }
            temp = content;
        }
  

        hspecDeri[i]->Rebin(8);
        hspecDeri[i]->GetXaxis()->SetRangeUser(500,1000);

        double valey = 1000; double minimum = 0;
        int index ;

        for(int k=62;k<=125;k++){

           if(valey > hspecDeri[i]->GetBinContent(k)) {valey = hspecDeri[i]->GetBinContent(k); minimum = hspecDeri[i]->GetBinCenter(k); index = k;}

          // cout<<hspecDeri[i]->GetBinCenter(k)<<endl;
        }

        //double mean = hspecDeri[i]->GetMean();
        double rms = hspecDeri[i]->GetRMS();
        double slope = hspecDeri[i]->GetBinContent(index+2)-hspecDeri[i]->GetBinContent(index);
        slope /= (hspecDeri[i]->GetBinCenter(index+2)-hspecDeri[i]->GetBinCenter(index));
        double initialcur = TMath::Abs(2.0*slope);
        cout<<minimum<<"\t"<<slope<<endl;foutfit<<minimum<<endl;
        double minimumseed = minimum;
        //prepare for first fit
        fitFcn = new TF1("fitFcn",myPol2,minimum-120,minimum+120,3);
        fitFcn->SetParName(1,"ComptonEdge");
        fitFcn->SetParameter(0,initialcur);
        fitFcn->SetParameter(1,minimum);
        fitFcn->SetParameter(2,valey);
        hspecDeri[i]->Fit("fitFcn","RQE");
        double p0 = fitFcn->GetParameter(0), p2 = fitFcn->GetParameter(2); 
        double width = TMath::Sqrt(-1.0*p2/p0);
        //foutfit<<"Fit once done \t"<<minimum<<"\t"<<width<<"\t"<<p0<<"\t"<<p2<<endl;
        double edge;// = fitFcn1->GetParameter(1);
        double err;// = fitFcn1->GetParError(1); //       double err2 = pol2fun->GetParError(2);
        //double err = edge*sqrt((err1/pol2fun->GetParameter(1))**2+(err2/pol2fun->GetParameter(2))**2);
        double goodness;// = fitFcn1->GetChisquare()/fitFcn1->GetNDF();
        double chi2p;
        //do fit a few more times
        for(int npass = 0; npass<5; npass++){
           if(minimum>=1000 || minimum<=500){ // need to get rid of bad fits
               minimum = minimumseed; width = 240;
           }
           if(width>=400 || width<=0){minimum = minimumseed; width = 200.0;}
           fitFcn1 = new TF1("fitFcn1",myPol2,minimum-width*0.75,minimum+width*0.75,3);
           fitFcn1->SetParName(1,"ComptonEdge");
           fitFcn1->SetParameters(p0,minimum,p2);
           foutfit<<"Looping \t"<<npass<<"\t"<<minimum<<"\t"<<width<<"\t"<<p0<<"\t"<<p2<<endl;
           hspecDeri[i]->Fit("fitFcn1","RQE");
           p0 = fitFcn1->GetParameter(0); p2 = fitFcn1->GetParameter(2); 
           width = TMath::Sqrt(-1.0*p2/p0);
           minimum = fitFcn1->GetParameter(1);
           edge = fitFcn1->GetParameter(1);
           err = fitFcn1->GetParError(1);
           goodness = fitFcn1->GetChisquare()/fitFcn1->GetNDF();
           chi2p = fitFcn1->GetProb();
           //cout<<"-----\t"<<edge<<"\t"<<width<<"\t"<<p0<<"\t"<<p2<<endl;
         }
        //double edge = -0.5*pol2fun->GetParameter(1)/pol2fun->GetParameter(2);
        //double edge = fitFcn1->GetParameter(1);

        //double err = fitFcn1->GetParError(1); //       double err2 = pol2fun->GetParError(2);
        //double err = edge*sqrt((err1/pol2fun->GetParameter(1))**2+(err2/pol2fun->GetParameter(2))**2);
        //double goodness = fitFcn1->GetChisquare()/fitFcn1->GetNDF();
        
        gEdge->SetPoint(i,15.0*(i+1),edge);
        gEdge->SetPointError(i,0,err);
        gEdgeRatio->SetPoint(i,15.0*(i+1),edge/664.2);
        gEdgeRatio->SetPointError(i,0,err/664.2);
        cout<<edge<<"\t"<<err<<"\t"<<goodness<<endl;
        ComptonEdgeHist->Fill(edge);
        gChi2->SetPoint(i,15.0*(i+1),goodness);
        hChi2Prob->Fill(chi2p);
        //cderi->cd();

        //hspecDeri[i]->SetLineColor(kBlue); hspecDeri[i]->Draw("");
      
        //hspecDeri[i]->Draw("");
        hspecDeri[i]->Write();
        //cderi->Print("animation_deri_spec.gif+30");
        //cderi->Close();;
    }

    TCanvas* c2 = new TCanvas();
    gEdge->SetName("ComptonEdgeTime");
    gEdge->GetXaxis()->SetTitle("Time [minute]");
    gEdge->GetYaxis()->SetTitle("Compton edge (Channel number)");
    gEdge->SetMarkerStyle(21);
    gEdge->Draw("ACP");
    c2->Write(); c2->Close();
    TCanvas* c3 = new TCanvas();
    gChi2->GetXaxis()->SetTitle("Time [minute]");
    gChi2->GetYaxis()->SetTitle("Fit goodness (chi2 / ndf)");
    gChi2->SetMarkerStyle(21);
    gChi2->Draw("ACP");
    c3->Write(); c3->Close();

    TCanvas* c4 = new TCanvas();
    gEdgeRatio->SetName("ComptonEdgeRatio");
    gEdgeRatio->GetXaxis()->SetTitle("Time [minute]");
    gEdgeRatio->GetYaxis()->SetTitle("Compton edge ratio (P50 LiLS / LAB)");
    gEdgeRatio->SetMarkerStyle(21);
    gEdgeRatio->Draw("ACP");
    c4->Write(); c4->Close();

    hChi2Prob->GetXaxis()->SetTitle("Fit probability");
    hChi2Prob->GetYaxis()->SetTitle("Counts");
    hChi2Prob->Write();

    ComptonEdgeHist->Write();

    cSpecSuperpose->Write(); cSpecSuperpose1->Write();
}


/*
        //cout<<filenames[i]<<endl;
        string smearname = "_smeared";
        string thisnamesmear = filenames[i] + smearname;
        string deriname = "_deri";
        string thisnamederi = filenames[i] + deriname;
        string funname = "_pol2fit";
        string thisfun = filenames[i] + funname;

        //hspecSmear[i]->Sumw2();
        hspecSmear[i] = smear(hspec[i], 10, thisnamesmear);

        c[i] = new TCanvas(); 
        c[i]->SetName(filenames[i].c_str());
        hspec[i]->SetLineWidth(2); hspec[i]->Draw();
        hspecSmear[i]->SetLineColor(kRed); hspecSmear[i]->Draw("same");
        //char savefilename[100]; sprintf(savefilename,"./spectrum/%s_%d_%d.png",filenames[i].c_str(),i,15.0*i);
        //c[i]->SaveAs(savefilename); 
        //c[i]->Print("animation_raw_spec.gif+30");
        c[i]->Write(); c[i]->Close();

        int nBins = hspecSmear[i]->GetNbinsX();

        hspecDeri[i] = new TH1F(thisnamederi.c_str(),"",4096,0,4096);

        double temp=0;
//        double valey = 100; double minimum = 0;
        for (int j=1; j<=nBins; j++) {
//            double x = h2->GetBinCenter(i);
            double content = hspecSmear[i]->GetBinContent(j);
//
            if(j>=500 && j<=1000){
                if(valey > (content-temp)) {
                    valey = content-temp; minimum = hspecSmear[i]->GetBinCenter(j); 
                    //cout<<valey<<"\t"<<content-temp<<"\t"<<minimum<<endl;

                }
            }
//
            hspecDeri[i]->SetBinContent(j, content-temp);
            //hspecDeri[i]->SetBinContent(j, content);
            temp = content;
        // cout << content << endl;
        }

        hspecDeri[i]->Rebin(8);
        hspecDeri[i]->GetXaxis()->SetRangeUser(500,1000);
        //cout<<hspecDeri[i]->GetNbinsX()<<endl;
        double valey = 1000; double minimum = 0;
        for(int k=62;k<=125;k++){
           if(valey > hspecDeri[i]->GetBinContent(k)) {valey = hspecDeri[i]->GetBinContent(k); minimum = hspecDeri[i]->GetBinCenter(k); }
          // cout<<hspecDeri[i]->GetBinCenter(k)<<endl;
        }
        //double mean = hspecDeri[i]->GetMean();
        double rms = hspecDeri[i]->GetRMS();

        TF1* pol2fun = new TF1("pol2fun","[0]*(x-[1])**2+[2]",minimum-1.2*rms,800);//minimum-2.0*rms,minimum+2.0*rms [0]*(x-[1])**2+[2]
        pol2fun->SetParName(1,"ComptonEdge");
        pol2fun->SetParameter(1,minimum);
        hspecDeri[i]->Fit("pol2fun","RQIE");

        minimum = pol2fun->GetParameter(1);
        double p0 = pol2fun->GetParameter(0), p2 = pol2fun->GetParameter(2); 
        double width = sqrt(-0.5*p2/p0);
        //cout<<minimum<<"\t"<<width<<"\t"<<p0<<"\t"<<p2<<endl;
        TF1* pol2fun1 = new TF1("pol2fun1","[0]*(x-[1])**2+[2]",minimum-width,minimum+width);
        pol2fun1->SetParameters(p0,minimum,p2);
        hspecDeri[i]->Fit("pol2fun1","RQIE");

        //double edge = -0.5*pol2fun->GetParameter(1)/pol2fun->GetParameter(2);
        double edge = pol2fun1->GetParameter(1);
        double err = pol2fun1->GetParError(1); //       double err2 = pol2fun->GetParError(2);
        //double err = edge*sqrt((err1/pol2fun->GetParameter(1))**2+(err2/pol2fun->GetParameter(2))**2);
        double goodness = pol2fun1->GetChisquare()/pol2fun1->GetNDF();
        cout<<edge<<"\t"<<err<<"\t"<<goodness<<endl;
        gEdge->SetPoint(i,15.0*(i+1),edge);
        gEdge->SetPointError(i,0,err);
        gChi2->SetPoint(i,15.0*(i+1),goodness);

        cderi->cd();
        hspecDeri[i]->Draw("pe0");
        hspecDeri[i]->Write();
        cderi->Print("animation_deri_spec.gif+30");
        cderi->Close();;
*/
