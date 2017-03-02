#include "smear_test.C"

const int nbOfHistos = 388;
gStyle->SetOptFit(11111111);

void GetMeanSpectrum(){
    TVirtualFitter::SetDefaultFitter("Minuit2");

    TFile* f =  new TFile("Histograms.root","read");
    fstream fout("MeanSpectrum.txt",ios::out);

    //TFile* rootfile = new TFile("Histograms_smoothed4-test.root","recreate");
    TH1F* hspec[nbOfHistos];
    TH1F* hspecSmear[nbOfHistos];
    TH1F* hspecDeri[nbOfHistos];
    TCanvas* c[nbOfHistos];
    //TF1* pol2fun[nbOfHistos];
    TGraphErrors* gEdge = new TGraphErrors();
    TGraph* gChi2 = new TGraph();

    TCanvas* cderi = new TCanvas();

    string filenames[nbOfHistos];
    fstream fin("filenamelist.txt",ios::in);//filenamelist = open("filenamelist.txt");
    string aname;
    int cnt;
    while(fin>>aname) {
        char tempstr[30]; sprintf(tempstr,"_spec%d",cnt);
        filenames[cnt] = aname + tempstr; 
        //cout<<filenames[cnt]<<endl;
        cnt++;
    }
    fin.close();

    TH1F* hMean = new TH1F("MeanSpectrum","",4096,0,4096); //use SetBinContent()
    TH1F* hMean1 = new TH1F("MeanSpectrum1","",4096,0,4096); //use Fill(x,w)
    double means[4096];

    for(int i=133; i<nbOfHistos;i++){
        string poststr = ";1";
        string thisname = filenames[i] + poststr;
        hspec[i] = (TH1F*)f->Get(thisname.c_str());

        filenames[i].replace(8, 4, "");
        string histname = hspec[i]->GetName(); histname.replace(8,4,"");
        hspec[i]->SetName(histname.c_str());


        for(int j=0;j<4096;j++){ means[j] += hspec[i]->GetBinContent(j+1)/255.0; }

    }

    for(int i=0; i<4096; i++){hMean->SetBinContent(i+1,means[i]); hMean1->Fill(i,means[i]); fout<<means[i]<<endl;}

    TCanvas* c2 = new TCanvas();
    hMean->Draw();
    hMean1->SetLineColor(kRed);
    hMean1->Draw("same");
    
}
