# This script trains a random forest model
# for the prediction of MOF water stability (2-class)
# and MOF acid stability

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold, train_test_split
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    roc_auc_score,
    f1_score,
    confusion_matrix,
    ConfusionMatrixDisplay,
    RocCurveDisplay,
    )
from functools import partial
from hyperopt import (
    fmin, 
    hp,
    space_eval, 
    tpe, 
    Trials
    )
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
import numpy as np
import joblib
import random
import os

df = pd.read_csv('features_and_labels.csv')
column_names = list(df.columns)
column_names.remove('MOF_name')  
column_names.remove('data_set')

# Make a copy of df, that is only MOFs with a water_label of 3 or 4
stable_df = df[df['water_label'].isin([3,4])]

# Getting the feature columns of the dataframe
feature_names = column_names.copy()
feature_names.remove('water_label')
feature_names.remove('acid_label')
feature_names.remove('base_label')
feature_names.remove('boiling_label')

def binarizer(val):
    if val in [1, 2]:
        return 0
    elif val in [3, 4]:
        return 1

def preprocess(X, y, task_name, rand_seed):
    # Accounting for different cases
    # In the case of water, need to make the problem 2-class
    # In other cases, will balance classes
    random.seed(rand_seed)
    if task_name == 'water':
        v_bin = np.vectorize(binarizer)
        y = v_bin(y)
    elif task_name in ['acid', 'base', 'boiling']:
        num_positive_examples = np.count_nonzero(y == 1)
        indices_positive_examples = np.where(y == 1)[0]
        indices_negative_examples = np.where(y == 0)[0]
        indices_negative_examples_kept = np.array(random.sample(list(indices_negative_examples), num_positive_examples))
        indices_kept = np.concatenate((indices_positive_examples, indices_negative_examples_kept))
        X = X[indices_kept, :]
        y = y[indices_kept]

    return X, y  

def preprocess_2(X, y, rand_seed):
    # Train test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, 
        y,
        stratify=y,
        test_size=0.2,
        random_state=rand_seed)

    train_names = X_train[:,0]
    test_names = X_test[:,0]
    X_train, X_test = X_train[:,1:], X_test[:, 1:]

    # Z normalization
    scaler = StandardScaler()
    scaler.fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    return X_train, X_test, y_train, y_test, scaler, train_names, test_names

# Evaluation function 
# args should be a dict
def model_eval(X_train, y_train, rand_seed, task_name, args):

    '''Take suggested arguments and perform model evaluation'''
    clf = RandomForestClassifier(random_state=rand_seed,
    max_depth=int(np.rint(args['max_depth'])),
    )

    if task_name == '4_class_water':
        scoring = 'roc_auc_ovo'
    else:
        scoring = 'roc_auc'
    scores = cross_val_score(clf, X_train, y=y_train, scoring=scoring) 
    cv_score = np.mean(scores)

    # Return the negative of the CV score to ensure we minimize the error by minimizing the negative ROC-AUC
    return -cv_score

def model_train(df, feature_names, label_name, rand_seed, four_class=False):
    task_name = label_name.replace('_label', '')
    if four_class and task_name == 'water':
        task_name = '4_class_' + task_name

    X = df[['MOF_name'] + feature_names].to_numpy() # Use ['MOF_name'] to generate CSV fo train/test splits
    y = df[label_name].to_numpy()
    X, y = preprocess(X, y, task_name, rand_seed)
    X_train, X_test, y_train, y_test, scaler, train_names, test_names = preprocess_2(X, y, rand_seed)

    # Save a CSV file of train/test split.
    train_test_bins = {'train_bin': list(train_names), 'test_bin': list(test_names)}

    # Hyperparameter optimization
    max_features_options = ["sqrt", None]
    parameter_space =  {
                        "max_depth": hp.uniform("max_depth", 3, 30),
                        }
    print("Start trials") 
    trials = Trials()
    best_hyperparams = fmin(
        partial(model_eval, X_train, y_train, rand_seed, task_name), 
        parameter_space, 
        algo=tpe.suggest, 
        max_evals=20, 
        trials=trials,
        rstate= np.random.RandomState(rand_seed) # np.random.default_rng(rand_seed) <-- accounting for older version of numpy
        )
    best_hyperparams['max_depth'] = int(np.rint(best_hyperparams['max_depth']))

    print(f'best_hyperparams is {best_hyperparams}')
    clf = RandomForestClassifier(random_state=rand_seed,
    max_depth=best_hyperparams['max_depth'],
    ).fit(X_train, y_train)

    y_hat_train = clf.predict(X_train)
    y_hat_test = clf.predict(X_test)
    if task_name == '4_class_water':
        y_proba_train = clf.predict_proba(X_train)
        y_proba_test = clf.predict_proba(X_test)
    else:
        y_proba_train = clf.predict_proba(X_train)[:, 1]
        y_proba_test = clf.predict_proba(X_test)[:, 1]

    # Model performance metrics
    train_acc = accuracy_score(y_train, y_hat_train)
    train_bal_acc = balanced_accuracy_score(y_train, y_hat_train)
    test_acc = accuracy_score(y_test, y_hat_test)
    test_bal_acc = balanced_accuracy_score(y_test, y_hat_test)
    if task_name == '4_class_water':
        train_roc_auc = roc_auc_score(y_train, y_proba_train, multi_class='ovo', average='macro')
        test_roc_auc = roc_auc_score(y_test, y_proba_test, multi_class='ovo', average='macro')
    else:
        train_roc_auc = roc_auc_score(y_train, y_proba_train)
        test_roc_auc = roc_auc_score(y_test, y_proba_test)

    train_F1 = f1_score(y_train, y_hat_train, average='macro')
    test_F1 = f1_score(y_test, y_hat_test, average='macro')        

    joblib.dump(clf, f'models/{task_name}_model.joblib')
    joblib.dump(scaler, f'models/{task_name}_scaler.joblib')  
    joblib.dump(X_train, f'models/{task_name}_X_train.joblib')  
    joblib.dump(y_train, f'models/{task_name}_y_train.joblib')  
    joblib.dump(X_test, f'models/{task_name}_X_test.joblib')  
    joblib.dump(y_test, f'models/{task_name}_y_test.joblib')
    with open(f'models/{task_name}_model_details.txt', 'w') as f:
        f.write(f'The features used are {feature_names}\n')
        f.write(f'The best hyperparameters are {best_hyperparams}\n')

    metrics_dict = {
    'train_acc': train_acc,
    'train_bal_acc': train_bal_acc,
    'train_roc_auc': train_roc_auc,
    'train_F1': train_F1,
    'test_acc': test_acc, 
    'test_bal_acc': test_bal_acc, 
    'test_roc_auc': test_roc_auc, 
    'test_F1': test_F1,
    }

    return metrics_dict, train_test_bins

def make_train_test_csv(df, label, all_train_test_bins):
    MOF_names = df['MOF_name'].tolist()
    ttb_df = pd.DataFrame({'MOF_name': MOF_names})
    for _i, ttb in enumerate(all_train_test_bins):
        bin_breakdown = []
        for mn in MOF_names:
            if mn in ttb['train_bin']:
                bin_breakdown.append('train')
            elif mn in ttb['test_bin']:
                bin_breakdown.append('test')
            else:
                # Current MOF is not in the train nor the test set for this split
                bin_breakdown.append(np.nan)
        ttb_df[f'split_{_i}'] = bin_breakdown
    ttb_df.to_csv(f'splits/{label}.csv', index=False)

required_folders = ['models', 'metrics', 'splits']
for i in required_folders:
    if not os.path.isdir(i):
        os.mkdir(i)

rfa_2_class_water_features = ['mc-Z-3-all', 'D_mc-Z-3-all', 'D_mc-Z-2-all', 'D_mc-Z-1-all', 'mc-chi-3-all', 'mc-Z-1-all', 'mc-Z-0-all', 'D_mc-chi-2-all', 'f-lig-Z-2', 'GSA', 'f-lig-I-0', 'func-S-1-all']
rfa_acid_features = ['mc-chi-3-all', 'Dif', 'mc-Z-2-all', 'Di', 'f-T-2-all', 'D_mc-chi-2-all', 'lc-chi-3-all', 'D_mc-S-1-all']
# rfa_base_features = ['D_mc-chi-1-all', 'f-lig-T-1', 'f-lig-S-2', 'Di', 'Df', 'f-lig-S-0', 'GPOV', 'GPOAV', 'f-I-1-all', 'f-lig-chi-0', 'D_lc-Z-3-all', 'D_func-alpha-2-all', 'f-lig-Z-0']
# rfa_boiling_features = ['D_mc-chi-3-all', 'f-S-0-all', 'mc-S-0-all', 'f-S-2-all', 'D_mc-S-2-all', 'Di', 'lc-S-3-all', 'f-Z-3-all']

label_list = ['water_label', 'acid_label'] #, 'base_label', 'boiling_label']
features_lists = [rfa_2_class_water_features, rfa_acid_features] #, rfa_base_features, rfa_boiling_features]
for label, features_to_use in zip(label_list, features_lists):

    print(f'Current label: {label}')
    all_train_test_bins = []
    rand_seed = 0

    if label == 'water_label':
        metrics_dict, train_test_bins = model_train(df, features_to_use, label, rand_seed)
    else:
        # use stable_df for acid, base, and boiling
        metrics_dict, train_test_bins = model_train(stable_df, features_to_use, label, rand_seed)

    all_train_test_bins.append(train_test_bins)

    with open(f'metrics/metrics_{label}.txt', 'w') as f:
        for key in metrics_dict:
            f.write(f'{key}: {metrics_dict[key]}\n')

    make_train_test_csv(df, label, all_train_test_bins)
        