import warnings
warnings.filterwarnings("ignore")
import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xgboost as xgb
import sys


cm = int(sys.argv[1])
mN = int(sys.argv[2])
V2 = 1e-2
Lumi = 1000000
V2_sim = 0.2**2

df_eff = pd.read_csv("../data/eff.csv")

# extract the reconstruction eff from the eff table
eff_M = df_eff[(df_eff.Energy==cm) & (df_eff.Mass==mN*1000) & (df_eff.Type=='s_M')]['Efficiency'].values[0]
eff_D = df_eff[(df_eff.Energy==cm) & (df_eff.Mass==mN*1000) & (df_eff.Type=='s_D')]['Efficiency'].values[0]
eff_baMu_qql = df_eff[(df_eff.Energy==cm) & (df_eff.Type=='b_aMu_qql')]['Efficiency'].values[0]
eff_baMu_qqv = df_eff[(df_eff.Energy==cm) & (df_eff.Type=='b_aMu_qqv')]['Efficiency'].values[0]
eff_bMuMu_qqlv = df_eff[(df_eff.Energy==cm) & (df_eff.Type=='b_MuMu_qqlv')]['Efficiency'].values[0]
eff_bMuMu_qqll = df_eff[(df_eff.Energy==cm) & (df_eff.Type=='b_MuMu_qqll')]['Efficiency'].values[0]
eff_baa_qqlv = df_eff[(df_eff.Energy==cm) & (df_eff.Type=='b_aa_qqlv')]['Efficiency'].values[0]

df_Xsec = pd.read_csv("../data/Xsec.csv")
Xsec_M = df_Xsec[(df_Xsec.Energy==cm) & (df_Xsec.Mass==mN*1000) & (df_Xsec.Type=='s_M')]['Efficiency'].values[0] / V2_sim
Xsec_D = df_Xsec[(df_Xsec.Energy==cm) & (df_Xsec.Mass==mN*1000) & (df_Xsec.Type=='s_D')]['Efficiency'].values[0] / V2_sim
Xsec_baMu_qql = df_Xsec[(df_Xsec.Energy==cm) & (df_Xsec.Type=='b_aMu_qql')]['Efficiency'].values[0]
Xsec_baMu_qqv = df_Xsec[(df_Xsec.Energy==cm) & (df_Xsec.Type=='b_aMu_qqv')]['Efficiency'].values[0]
Xsec_bMuMu_qqlv = df_Xsec[(df_Xsec.Energy==cm) & (df_Xsec.Type=='b_MuMu_qqlv')]['Efficiency'].values[0]
Xsec_bMuMu_qqll = df_Xsec[(df_Xsec.Energy==cm) & (df_Xsec.Type=='b_MuMu_qqll')]['Efficiency'].values[0]
Xsec_baa_qqlv = df_Xsec[(df_Xsec.Energy==cm) & (df_Xsec.Type=='b_aa_qqlv')]['Efficiency'].values[0]

# read the files, and calculate the corresponding yields
modes_dt = {1: {'mode': f'../data/features/sig_Maj_E-{cm}_m-{mN*1000}_reco.root', 'yields': Lumi*Xsec_M*eff_M*V2}, 
            2: {'mode': f'../data/features/sig_Dir_E-{cm}_m-{mN*1000}_reco.root', 'yields': Lumi*Xsec_D*eff_D*V2}, 
            3: {'mode': f'../data/features/bg_aMu_qql_E-3_reco.root', 'yields': Lumi*Xsec_baMu_qql*eff_baMu_qql},
            4: {'mode': f'../data/features/bg_aMu_qqv_E-3_reco.root', 'yields': Lumi*Xsec_baMu_qqv*eff_baMu_qqv},
            5: {'mode': f'../data/features/bg_MuMu_qqlv_E-3_reco.root', 'yields': Lumi*Xsec_bMuMu_qqlv*eff_bMuMu_qqlv},
            6: {'mode': f'../data/features/bg_MuMu_qqll_E-3_reco.root', 'yields': Lumi*Xsec_bMuMu_qqll*eff_bMuMu_qqll},
            7: {'mode': f'../data/features/bg_aa_qqlv_E-3_reco.root', 'yields': Lumi*Xsec_baa_qqlv*eff_baa_qqlv}
            }


def load_train_test(modes_dt, size=0.5):
    np.random.seed(9)
    df_train = pd.DataFrame()
    df_test = pd.DataFrame()
    num_lt = []
        
    # loop over different modes
    for i, (k, v) in enumerate(modes_dt.items()):
        file = uproot.open(v['mode'])
        print("reading: ", v['mode'])
        df_i = pd.DataFrame(np.array(file['t']['features'].array()))
        df_i['target'] = k    # add the target label
        df_i['weight'] = v['yields']/len(df_i)
        num_lt.append(len(df_i))

        # shuffle the index for training and testing sets
        idx = df_i.index.tolist()
        np.random.shuffle(idx)
                
        # cut according to the fraction
        cut = int(np.ceil(len(idx) * size))
        df_train_i = df_i.loc[idx[:cut]]
        df_test_i = df_i.loc[idx[cut:]]
                
        # Put to the global dataframs
        df_train = pd.concat([df_train, df_train_i])
        df_test = pd.concat([df_test, df_test_i])
            
    df_train.reset_index(drop=True, inplace=True)
    df_test.reset_index(drop=True, inplace=True)
    
    print('\ntrain size: {} ({:.2f}%); test size: {} ({:.2f}%)'.format(len(df_train), 100*len(df_train)/(len(df_train)+len(df_test)), len(df_test), 100*len(df_test)/(len(df_train)+len(df_test))))
    print('data points per mode: ',num_lt)
    return df_train, df_test

    
tr_te_size = 0.5 
df_train, df_test = load_train_test(modes_dt, tr_te_size)

df_train_o, df_test_o = df_train.copy(), df_test.copy()



def relabel(x):
    if x != 1 and x != 2:    # bkg
        return 0
    elif x == 1:     # Maj signal
        return 1
    elif x == 2:     # Dir signal
        return 2




features = ['ptLep', 'etaLep', 'ELep',      # lepton kinematics info 
            'chargeLep', 'lepisMu',         # lepton type info
            'DeltaPhijjl', 'DeltaRjjl',     # Delta phi & Delta R between jj(from W boson) and lepton
            'ptJJ', 'etaJJ', 'mJJ',         # jj(from W boson) momentum info
            'ptN', 'pzN',                   # reconstructed N 4 momentum info
            'EJet1', 'EJet2',               # energy inbalance of the jets, for W boson classification
            'ptFwMu'                        # Forward muon pT
           ]



# # relabel the, all bkg become one label
df_train['target'] = df_train['target'].apply(relabel)
df_test['target'] = df_test['target'].apply(relabel)


# trainging
X_train = df_train[features]
y_train = df_train['target']

# testing
X_test = df_test[features]
y_test = df_test['target']

print("traning BDT")
# fitting, with reweighting
xgbc1 = xgb.XGBClassifier(seed=0)
xgbc1.fit(X_train, y_train, sample_weight=df_train.weight.values);
# xgbc1.fit(X_train, y_train);
print("Finished training")



# the "preds" means the prob. of the event being HNL (Maj/Dir)
# since the 0 index is the bkg, so the prob. of HNL is 1 - prob. of bkg

# scores for training
df_bdt_train_s = df_train[['target', 'weight']]
df_bdt_train_s.loc[:, 'preds'] = 1 - xgbc1.predict_proba(X_train)[:, 0]

# score for testing 
df_bdt_test_s = df_test[['target', 'weight']]
df_bdt_test_s.loc[:, 'preds'] = 1 - xgbc1.predict_proba(X_test)[:, 0]

def bdt_cut(df, cut):
    df1 = df[(df['preds'] >= cut)]    # pass cut to be classified as HNL
    #    "/ 2" because the samples include M & D. now assume they are just half half
    S1 = df1[df1['target'] != 0]['weight'].values.sum() / 2    # Number of signals
    B1 = df1[df1['target'] == 0]['weight'].values.sum()    # Number of bkg.
    return ((S1+B1)**0.5, (0.1*B1), S1/(np.sqrt(S1+B1)))




# looping to find the optimal cut
def find_opt_cut(df_bdt_train_s, step=0.01):
    # cuts
    cuts = np.arange(0, 1, step)
        
    # store cuts
    c_lt = []
    # store losses
    loss1, loss2 = [], []
    snr = []

    # loop over two cuts
    for i1, c1 in enumerate(cuts):
        print("{}/{}".format(i1, len(cuts)), end='\r')
        res = bdt_cut(df_bdt_train_s, c1)
        c_lt.append(c1)    
        loss1.append(res[0])
        loss2.append(res[1])
        snr.append(res[2])
                
    # store the cuts and corresponding losses
    df_bdt_loss = pd.DataFrame([c_lt, loss1, loss2, snr]).T
    df_bdt_loss.columns = ['c', 'loss1', 'loss2', 'snr']

    df_bdt_loss['loss'] = (df_bdt_loss['loss1']**2 + df_bdt_loss['loss2']**2)**0.5
    df_bdt_loss['tot'] = (df_bdt_loss['loss'] - df_bdt_loss['loss'].min())/(df_bdt_loss['loss'].max() - df_bdt_loss['loss'].min()) - \
                        ((df_bdt_loss['snr'] - df_bdt_loss['snr'].min())/(df_bdt_loss['snr'].max() - df_bdt_loss['snr'].min()))
                            
    return df_bdt_loss.iloc[df_bdt_loss['tot'].argmin()][['c']].values



# the optimal threshold
threshold = find_opt_cut(df_bdt_train_s)[0]



df_bdt_cut = df_bdt_test_s[(df_bdt_test_s['preds'] >= threshold)]


# separate to Maj/Dir/bkg

df_bdt_cut['iEvt'] = df_test_o.loc[df_bdt_cut.index]['iEvt']
df_bdt_cut['mN'] = df_test_o.loc[df_bdt_cut.index]['mN']
df_bdt_cut['lepisEle'] = df_test_o.loc[df_bdt_cut.index]['lepisEle']
df_bdt_cut['lepisMu'] = df_test_o.loc[df_bdt_cut.index]['lepisMu']
df_bdt_cut1 = df_bdt_cut[(df_bdt_cut['target'] == 1)]
df_bdt_cut2 = df_bdt_cut[(df_bdt_cut['target'] == 2)]
df_bdt_cut0 = df_bdt_cut[(df_bdt_cut['target'] != 1) & (df_bdt_cut['target'] != 2)]


partition_u = 0.025
partition_l = 0.05
df_bdt_cut['mN'] = df_test_o.loc[df_bdt_cut.index]['mN']
df_bdt_cut = df_bdt_cut[(df_bdt_cut.mN>=1000*mN*(1-partition_l)) & (df_bdt_cut.mN<=1000*mN*(1+partition_u))]

# limits on |V|^2

# yields after BDT and mass cuts

df_bdt_cut0 = df_bdt_cut[(df_bdt_cut['target'] == 0)]   # Bkg
df_bdt_cut1 = df_bdt_cut[(df_bdt_cut['target'] == 1)]   # M
df_bdt_cut2 = df_bdt_cut[(df_bdt_cut['target'] == 2)]   # D

# yield (without |V|^2 dependence)
nb = df_bdt_cut0.weight.sum() / tr_te_size          # Bkg
nM = df_bdt_cut1.weight.sum() / tr_te_size / V2    # M
nD = df_bdt_cut2.weight.sum() / tr_te_size / V2    #

N = 1.96
M = N**2



print()
print("E=",cm)
print("mN=",mN)

print('Fully Majorana case:')
upperM, lowerM = ((M+(M**2 + 4*M*nb)**0.5)/2)/(nM), ((M-(M**2 + 4*M*nb)**0.5)/2)/(nM)
print(f"upper limit on |V|^2: {upperM:.2e}")
print(f"lower limit on |V|^2: {lowerM:.2e}")
print()

print('Fully Dirac case:')
upperD, lowerD = ((M+(M**2 + 4*M*nb)**0.5)/2)/(nD), ((M-(M**2 + 4*M*nb)**0.5)/2)/(nD)
print(f"upper limit on |V|^2: {upperD:.2e}")
print(f"lower limit on |V|^2: {lowerD:.2e}")



print("STORE: e={}; m={}; V2_M={}; V2_D={} END STORE".format(cm, mN, upperM, upperD))
