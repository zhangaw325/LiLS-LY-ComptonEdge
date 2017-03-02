# LiLS-LY-ComptonEdge
My scripts to process data from LS6500 machine for liquid scintillation light yield measurments

1, need python and pyroot to run the .py file, which will then read the raw data files (in excel format) and store histograms (Cs137 spectrum) in a root file.

2, In ROOT, use GetMean, GetCovariance scritps to calculate average and standard deviation of each ADC channel (bin). All stable runs need to be used. By stable, I mean the spectra should not change along time. [in my case, I measure liquid scintillator where Cs137 spectrum shows no peak, and the spectrum changes due to oxygen quenching.]

3, In ROOT, finally run the processonefile to process all runs. 
