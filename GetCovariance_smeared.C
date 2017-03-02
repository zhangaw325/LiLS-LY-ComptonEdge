#include "smear_test.C"

const int nbOfHistos = 388;
gStyle->SetOptFit(11111111);

void GetCovariance_smeared(){
    TVirtualFitter::SetDefaultFitter("Minuit2");

    TFile* f =  new TFile("Histograms.root","read");

    TFile* rootfile = new TFile("Histograms_smoothed4-test.root","recreate");
    TH1F* hspec[nbOfHistos];
    TH1F* hspecSmear[nbOfHistos];
    TH1F* hspecDeri[nbOfHistos];
    TCanvas* c[nbOfHistos];
    //TF1* pol2fun[nbOfHistos];
    TGraphErrors* gEdge = new TGraphErrors();
    TGraph* gChi2 = new TGraph();

    TCanvas* cderi = new TCanvas();

    TH1F* myCov = new TH1F("Covariance","",4096,0,4096);
    //TH1F* myRawErr = new TH1F("RawError","",4096,0,4096);
    fstream fin("MeanSmearedSpectrum.txt",ios::in);
    fstream fout("CovarianceSmeared.txt",ios::out);
    fstream fsigerr("RawUncertaintySmeared.txt",ios::out);
    double means[4096],a;
    int acnt=0;
    while(fin>>a){means[acnt]=a;acnt++;}
    fin.close();

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

    double cov[4096]; cov[0]=0;
    double temp_cov[388][4096];
    double err[4096];
    //TMatrixD cov(4096,4096);
    //TMatrixD Matrix1(388,4096);
    //TMatrixD Matrix2(4096,4096);
    //TMatrixD vecc(1,4096);

    for(int i=133; i<nbOfHistos;i++){
        string poststr = ";1";
        string thisname = filenames[i] + poststr;
        hspec[i] = (TH1F*)f->Get(thisname.c_str());

        filenames[i].replace(8, 4, "");
        string histname = hspec[i]->GetName(); histname.replace(8,4,"");
        hspec[i]->SetName(histname.c_str());

        string smearname = "_smeared";
        string thisnamesmear = filenames[i] + smearname;
        string deriname = "_deri";
        string thisnamederi = filenames[i] + deriname;
        string funname = "_pol2fit";
        string thisfun = filenames[i] + funname;

        //hspecSmear[i]->Sumw2();
        hspecSmear[i] = smear(hspec[i], 20, thisnamesmear);

        for(int j=0;j<4096;j++){
            temp_cov[i][j] = hspecSmear[i]->GetBinContent(j+1) - means[j];
            //Matrix1(i,j) = hspec[i]->GetBinContent(j+1) - means[j];
            err[j] += (hspecSmear[i]->GetBinContent(j+1)-means[j])**2;
        }
    }
    for(int k=1;k<4096;k++){
        for(int m=133;m<nbOfHistos;m++){
            //cov[k] += Matrix1(m,k-1)*Matrix1(m,k); // nbOfHistos - 1;
            cov[k] += temp_cov[m][k-1]*temp_cov[m][k];
        }
    }
    
    for(int k=0;k<4096;k++) {/*myCov->SetBinContent(k+1,cov[k]/255.0);*/fout<<cov[k]/255.0<<endl;/*myRawErr->SetBinContent(k+1,sqrt(err[k]/255.0));*/fsigerr<<sqrt(err[k]/255.0)<<endl;}

//    TCanvas* c2 = new TCanvas();
//    myCov->Draw();
//    TCanvas* c3 = new TCanvas();
//    myRawErr->Draw();
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
