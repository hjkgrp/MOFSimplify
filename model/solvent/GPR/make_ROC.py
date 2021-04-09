import numpy as np
import pandas as pd
import sklearn.preprocessing
from scipy.stats import pearsonr, spearmanr
import sys
import os
import json
import pickle
import matplotlib.pyplot as plt

blue = "rgba(0, 0, 255, 1)"
red = "rgba(255, 0, 0, 1)"
green = "rgba(0, 196, 64, 1)"
gray = "rgba(140, 140, 140, 1)"

def normalize_data(df_train, df_test, fnames, lname, unit_trans=1, debug=False):
    _df_train = df_train.copy().dropna(subset=fnames+lname)
    _df_test = df_test.copy().dropna(subset=fnames+lname)
    X_train, X_test = _df_train[fnames].values, _df_test[fnames].values
    y_train, y_test = _df_train[lname].values, _df_test[lname].values
    if debug:
        print("training data reduced from %d -> %d because of nan." % (len(df_train), y_train.shape[0]))
        print("test data reduced from %d -> %d because of nan." % (len(df_test), y_test.shape[0]))
    x_scaler = sklearn.preprocessing.StandardScaler()
    x_scaler.fit(X_train)
    X_train = x_scaler.transform(X_train)
    X_test = x_scaler.transform(X_test)
    y_train = np.array([1 if x == 1 else 0 for x in y_train.reshape(-1, )])
    y_test = np.array([1 if x == 1 else 0 for x in y_test.reshape(-1, )])
    return X_train, X_test, y_train, y_test, x_scaler

def standard_labels(df, key="flag"):
    flags = [1 if row[key] == 1 else 0 for _, row in df.iterrows()]
    df[key] = flags
    return df

RACs = ['D_func-I-0-all','D_func-I-1-all','D_func-I-2-all','D_func-I-3-all',
 'D_func-S-0-all', 'D_func-S-1-all', 'D_func-S-2-all', 'D_func-S-3-all',
 'D_func-T-0-all', 'D_func-T-1-all', 'D_func-T-2-all', 'D_func-T-3-all',
 'D_func-Z-0-all', 'D_func-Z-1-all', 'D_func-Z-2-all', 'D_func-Z-3-all',
 'D_func-chi-0-all', 'D_func-chi-1-all', 'D_func-chi-2-all',
 'D_func-chi-3-all', 'D_lc-I-0-all', 'D_lc-I-1-all', 'D_lc-I-2-all',
 'D_lc-I-3-all', 'D_lc-S-0-all', 'D_lc-S-1-all', 'D_lc-S-2-all',
 'D_lc-S-3-all', 'D_lc-T-0-all', 'D_lc-T-1-all', 'D_lc-T-2-all',
 'D_lc-T-3-all', 'D_lc-Z-0-all', 'D_lc-Z-1-all', 'D_lc-Z-2-all',
 'D_lc-Z-3-all', 'D_lc-chi-0-all', 'D_lc-chi-1-all', 'D_lc-chi-2-all',
 'D_lc-chi-3-all', 'D_mc-I-0-all', 'D_mc-I-1-all', 'D_mc-I-2-all',
 'D_mc-I-3-all', 'D_mc-S-0-all', 'D_mc-S-1-all', 'D_mc-S-2-all',
 'D_mc-S-3-all', 'D_mc-T-0-all', 'D_mc-T-1-all', 'D_mc-T-2-all',
 'D_mc-T-3-all', 'D_mc-Z-0-all', 'D_mc-Z-1-all', 'D_mc-Z-2-all',
 'D_mc-Z-3-all', 'D_mc-chi-0-all', 'D_mc-chi-1-all', 'D_mc-chi-2-all',
 'D_mc-chi-3-all', 'f-I-0-all', 'f-I-1-all', 'f-I-2-all', 'f-I-3-all',
 'f-S-0-all', 'f-S-1-all', 'f-S-2-all', 'f-S-3-all', 'f-T-0-all', 'f-T-1-all',
 'f-T-2-all', 'f-T-3-all', 'f-Z-0-all', 'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all',
 'f-chi-0-all', 'f-chi-1-all', 'f-chi-2-all', 'f-chi-3-all', 'f-lig-I-0',
 'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3', 'f-lig-S-0', 'f-lig-S-1', 'f-lig-S-2',
 'f-lig-S-3', 'f-lig-T-0', 'f-lig-T-1', 'f-lig-T-2', 'f-lig-T-3', 'f-lig-Z-0',
 'f-lig-Z-1', 'f-lig-Z-2', 'f-lig-Z-3', 'f-lig-chi-0', 'f-lig-chi-1',
 'f-lig-chi-2', 'f-lig-chi-3', 'func-I-0-all', 'func-I-1-all',
 'func-I-2-all', 'func-I-3-all', 'func-S-0-all', 'func-S-1-all',
 'func-S-2-all', 'func-S-3-all', 'func-T-0-all', 'func-T-1-all',
 'func-T-2-all', 'func-T-3-all', 'func-Z-0-all', 'func-Z-1-all',
 'func-Z-2-all', 'func-Z-3-all', 'func-chi-0-all', 'func-chi-1-all',
 'func-chi-2-all', 'func-chi-3-all', 'lc-I-0-all', 'lc-I-1-all', 'lc-I-2-all',
 'lc-I-3-all', 'lc-S-0-all', 'lc-S-1-all', 'lc-S-2-all', 'lc-S-3-all',
 'lc-T-0-all', 'lc-T-1-all', 'lc-T-2-all', 'lc-T-3-all', 'lc-Z-0-all',
 'lc-Z-1-all', 'lc-Z-2-all', 'lc-Z-3-all', 'lc-chi-0-all', 'lc-chi-1-all',
 'lc-chi-2-all', 'lc-chi-3-all', 'mc-I-0-all', 'mc-I-1-all', 'mc-I-2-all',
 'mc-I-3-all', 'mc-S-0-all', 'mc-S-1-all', 'mc-S-2-all', 'mc-S-3-all',
 'mc-T-0-all', 'mc-T-1-all', 'mc-T-2-all', 'mc-T-3-all', 'mc-Z-0-all',
 'mc-Z-1-all', 'mc-Z-2-all', 'mc-Z-3-all', 'mc-chi-0-all', 'mc-chi-1-all',
 'mc-chi-2-all', 'mc-chi-3-all']
geo = ['Df','Di', 'Dif','GPOAV','GPONAV','GPOV','GSA','POAV','POAV_vol_frac',
  'PONAV','PONAV_vol_frac','VPOV','VSA','rho']
 
other = ['cif_file','name','filename']
path = os.getcwd()+'/'+str(sys.argv[1])
df_train_all = pd.read_csv(path+"/train.csv").append(pd.read_csv(path+"/val.csv"))
df_train = pd.read_csv(path+"/train.csv")
df_train = df_train.loc[:, (df_train != df_train.iloc[0]).any()]
df_val = pd.read_csv(path+"/val.csv")
df_test = pd.read_csv(path+"/test.csv")
features = [val for val in df_train.columns.values if val in RACs+geo]

df_train = standard_labels(df_train, key="flag")
df_val = standard_labels(df_val, key="flag")
df_test = standard_labels(df_test, key="flag")
df_train_all = standard_labels(df_train_all, key="flag")

X_train, X_test, y_train, y_test, x_scaler = normalize_data(df_train, df_test, features, ["flag"], unit_trans=1, debug=False)
X_train, X_val, y_train, y_val, x_scaler = normalize_data(df_train, df_val, features, ["flag"], unit_trans=1, debug=False)


rfa_features = json.load(open(path+'/rfa_features.json','r'))['features'] # /Users/gianmarcoterrones/Research/MOFs/for_GT/solvent/GPR/dropped_connectivity_dupes_rfa_RBF/rfa_features.json
joint_train_val = pd.concat([df_train,df_val],axis=0)
X_train, X_test, y_train, y_test, x_scaler = normalize_data(joint_train_val, df_test, rfa_features, ["flag"], unit_trans=1, debug=False)
X_train.shape, y_train.reshape(-1, ).shape
import GPy
kernel = GPy.kern.Matern52(input_dim=X_train.shape[1])
gp = GPy.models.GPClassification(X_train, y_train.reshape(-1, 1), kernel)
gp.optimize(messages=True)

from sklearn import metrics
from sklearn.metrics import roc_curve, auc, accuracy_score, roc_auc_score, confusion_matrix
fpr_train, tpr_train, _ = roc_curve(y_train.ravel(), np.squeeze(gp.predict(X_train)[0]).ravel())
fpr_test, tpr_test, _ = roc_curve(y_test.ravel(), np.squeeze(gp.predict(X_test)[0]).ravel())
roc_auc_vals_train = auc(fpr_train, tpr_train)
roc_auc_vals_test = auc(fpr_test, tpr_test)
plt.figure()
lw = 2
plt.plot(fpr_train, tpr_train, color='#8c8c8c',
         lw=lw, label='Train ROC curve (area = %0.2f)' % roc_auc_vals_train)
plt.plot(fpr_test, tpr_test, color='#0000FF',
         lw=lw, label='Test ROC curve (area = %0.2f)' % roc_auc_vals_test)
plt.plot([0, 1], [0, 1], color='k', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.savefig(path+'/final_ROC.png')

train_roc = pd.DataFrame()
train_roc['fpr'] = fpr_train
train_roc['tpr'] = tpr_train
test_roc = pd.DataFrame()
test_roc['fpr'] = fpr_test
test_roc['tpr'] = tpr_test
train_roc.to_csv(path+'/'+'train_roc_curve.csv',index=False)
test_roc.to_csv(path+'/'+'test_roc_curve.csv',index=False)

threshold_3 = []
threshold_4 = []
threshold_6 = []
threshold_7 = []

train_GPC = pd.DataFrame()
train_GPC['CoRE_name'] = joint_train_val['CoRE_name']
train_GPC['actual_y'] = y_train
train_GPC['prediction_y'] = gp.predict(X_train)[0]
for i, row in train_GPC.iterrows():
    if row['prediction_y']>= 0.3:
        threshold_3.append(1)
    else:
        threshold_3.append(0)
for i, row in train_GPC.iterrows():
    if row['prediction_y']>= 0.4:
        threshold_4.append(1)
    else:
        threshold_4.append(0)
for i, row in train_GPC.iterrows():
    if row['prediction_y']>= 0.6:
        threshold_6.append(1)
    else:
        threshold_6.append(0)
for i, row in train_GPC.iterrows():
    if row['prediction_y']>= 0.7:
        threshold_7.append(1)
    else:
        threshold_7.append(0)
train_GPC['threshold_0.3'] = threshold_3
train_GPC['threshold_0.4'] = threshold_4
train_GPC['threshold_0.6'] = threshold_6
train_GPC['threshold_0.7'] = threshold_7

test_GPC = pd.DataFrame()
test_GPC['CoRE_name'] = df_test['CoRE_name']
test_GPC['actual_y'] = y_test
test_GPC['prediction_y'] = gp.predict(X_test)[0]

threshold_3 = []
threshold_4 = []
threshold_6 = []
threshold_7 = []
for i, row in test_GPC.iterrows():
    if row['prediction_y']>= 0.3:
        threshold_3.append(1)
    else:
        threshold_3.append(0)
for i, row in test_GPC.iterrows():
    if row['prediction_y']>= 0.4:
        threshold_4.append(1)
    else:
        threshold_4.append(0)
for i, row in test_GPC.iterrows():
    if row['prediction_y']>= 0.6:
        threshold_6.append(1)
    else:
        threshold_6.append(0)
for i, row in test_GPC.iterrows():
    if row['prediction_y']>= 0.7:
        threshold_7.append(1)
    else:
        threshold_7.append(0)
test_GPC['threshold_0.3'] = threshold_3
test_GPC['threshold_0.4'] = threshold_4
test_GPC['threshold_0.6'] = threshold_6
test_GPC['threshold_0.7'] = threshold_7
train_GPC.to_csv(path+'/train_GPC.csv',index=False) # /Users/gianmarcoterrones/Research/MOFs/for_GT/solvent/GPR/dropped_connectivity_dupes_rfa_RBF/train_GPC.csv
test_GPC.to_csv(path+'/test_GPC.csv',index=False)


from sklearn.metrics import confusion_matrix
print(confusion_matrix(y_test,(np.round(gp.predict(X_test)[0])).astype(int)),'test')
print(confusion_matrix(y_train,(np.round(gp.predict(X_train)[0])).astype(int)),'train')
import numpy as np
np.save(path+'/model_ROC.npy', gp.param_array)
with open(path+'/ROC_model_performance.csv','w') as f:
    f.writelines('train/val acc: '+str(accuracy_score(y_train, np.round(np.squeeze(gp.predict(X_train)[0]))))+', train/val auc: '+str(roc_auc_score(y_train, np.squeeze(gp.predict(X_train)[0])))+'\n')
    f.writelines('test acc: '+str(accuracy_score(y_test, np.round(np.squeeze(gp.predict(X_test)[0]))))+', test auc: '+str(roc_auc_score(y_test, np.squeeze(gp.predict(X_test)[0])))+'\n')

# TODO ask: where to feed in complex for property prediction?