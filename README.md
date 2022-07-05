# Probing Heavy Neutral Leptons using Muon Collider

To generate root file from LHE:
go `PROJECT/PATH/data_ISR/script`
and do for example: `./convertLHEtoROOT_v2 ../sig_Maj_E-3_m-1000.lhe ../sig_Maj_E-3_m-1000.root`

To run reconstruction, go: `PROJECT/PATH/reconstruction`
To run single reconstruction on a specific data file, do:
`root -q -l 'allinone_ISR(type, save, num_test)'`
For example: to run and save all samples on the 1000GeV Majorana at sqrt{s}=3TeV, do:
`root -q -l 'allinone_ISR("s_M_3_1000", true, 0)'`

[Not finsihed] To run all signas at sqrt{s}=E TeV, mN=m GeV, do:
`./run_sig_E_m 3 1000`
store the terminal printed output to some file, do: for example
`./run_sig_E_m 3 1000 > ../data_ISR/reco_output.txt`

To run all backgrounds at sqrt{s}=E TeV, do:
`./run_bg_E 3`
append the terminal printed output to the existing file, do: for example
`./run_bg_E 3 >> ../data_ISR/reco_output.txt`


To extract the reconstruction eff. from the printed output, and store to a csv (for later BDT analysis), do:
`cd PROJECT/PATH/data_ISR/script`
`python to_yield_table.py PATH/TO/INPUT PATH/TO/OUTPUT`
for exmaple `python to_yield_table.py ../data_ISR/reco_output.txt ../data_ISR/eff_table.csv`


To run BDT analysis: go `PROJECT/PATH/python`
and use jupyter notebook to read `BDT_ISR.ipynb`
Modify the target sqrt{s}(cm), mass(mN), and the reweighting |V|^2 (can use default)
****Need to modify the cross section accordingly in the .ipynb (cell 6)*



File system: 
There could be multiple final states e.g. N>lW(>jj), N>lW(>lv)
So, make two folders for the two main channels 
Name: N_ljj means the first channels; N_llv means the second channels 
Inside each channel folder:
Having main folder:
1. reconstruction: mainly .C codes, and scripts to run all sig/bg for convenient
2. data:
    a. detector: the raw data after Delphes simulation 
    b. features: the features extracted from the reconstruction, and to be used in the BDT later
3. python: BDT or other related python anaylsis

Some common and useful files:
1. scripts: e.g. the converter for the data, and the eff, and the yields
