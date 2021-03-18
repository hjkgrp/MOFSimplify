"""Utility for generating PCA maps and models for the training space"""
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.decomposition import PCA # PCA for visualization of variance in space.
from sklearn.preprocessing import StandardScaler # Scaler for PCA
from sklearn.pipeline import Pipeline
import pickle

names20 = []
oxs = []
indices20 = []
# Read in names and indicies of alpha=20 complexes
with open('train_names.csv','r') as file1:
    for i,line in enumerate(file1):
        sline = line.split(',')
        name = sline[0]
        ox = sline[1].split('=')[1]
        if sline[2].split('=')[1].strip('\n') == '0.2':
            names20.append(name) # Name
            oxs.append(ox) # Ox to add to name
            indices20.append(i) # Indices to map to RACs

dent_dict = {
'acac':2,
 'aceticacidbipyridine':2,
 'ammonia':1,
 'benzisc':1,
 'bipy':2,
 'carbonyl':1,
 'cat':2,
 'chloropyridine':1,
 'cl':1,
 'cyanide':1,
 'cyanopyridine':1,
 'en':2,
 'ethbipyridine':2,
 'ethylesteracac':2,
 'furan':1,
 'isothiocyanate':1,
 'mebipyridine':2,
 'mecat':2,
 'methylamine':1,
 'misc':1,
 'ox':1,
 'phen':2,
 'phenacac':2,
 'phenisc':1,
 'phosacidbipyridine':2,
 'pisc':1,
 'porphyrin':4,
 'pyridine':1,
 'pyrrole':1,
 'tbisc':1,
 'tbuc':2,
 'thiopyridine':1,
 'water':1
 }

duplicates_dict = {
  'chloride':'cl',
  'pyr':'pyridine',
  'mec':'mecat'
 }
 
ox_convert_dict = {
    '2': '(II)',
    '3': '(III)'
}

def name_convert(name,ox):
    """ Name converter utility

    Parameters
    ----------
    name : str
        'metal_lig1_lig2_lig3' - assuming lig1 eq, lig2 ax1, lig3 ax2,
        except where all bidentate. Then just listing as "formula"
    ox : str
        oxidation state as either '2' or '3'  
    """
    metal,l1,l2,l3 = name.split('_')
    ss = metal.title() + ox_convert_dict[ox]
    if l1 in duplicates_dict:
        l1 = duplicates_dict[l1]
    if l2 in duplicates_dict:
        l2 = duplicates_dict[l2]
    if l3 in duplicates_dict:
        l3 = duplicates_dict[l3]
    ligs = [l1,l2,l3]
    # print(ligs)
    dents = [dent_dict[x] for x in ligs]
    if dents[0] == 4: # 
        ss += '('+l1+')eq'+'('+l2+')ax1'+'('+l3+')ax2'
    elif (len(set(dents)) == 1) and (dents[0] == 2): # Bidentate
        uligs = list(set(ligs))
        ucount = [ligs.count(x) for x in uligs]
        scount,sligs = zip(*sorted(zip(ucount, uligs)))
        for l,c in zip(sligs, scount):
            ss += '('+l+')'+str(c)
    elif (len(set(dents)) == 1) and (dents[0] == 1): # Monodentate
        ss += '('+l1+')4eq'+'('+l2+')ax1'+'('+l3+')ax2'
    elif (len(set(dents)) == 2): # Bidentate and Monodentate - bidentate eq plane
        ss += '('+l2+')2eq' + '('+l2+')ax1'+'('+l3+')ax2'
    return ss

savenames = []
# Convert the names20 and ox into savenames.
for i,n in enumerate(names20):
    newname = name_convert(n,oxs[i])
    savenames.append(newname)

def rac_name_convert(name): # Model contains old RAC Names. 
    rac_name = name.replace('.','-')
    return rac_name


trainlabels = pd.read_csv('split_y.csv')
reducedlabels = trainlabels.loc[indices20]

traindf = pd.read_csv('split_x.csv')
reduceddf = traindf.loc[indices20]

rac_names = reduceddf.columns.values
# LASSO-28 from https://pubs.acs.org/doi/10.1021/acs.jpca.7b08750
rac_names = ['ox', 'alpha', 'misc.dent.eq', 'f.Z.0.ax', 'f.I.3.ax', 'f.T.3.ax',
            'f.Z.2.eq','lc.chi.0.ax','lc.Z.2.ax','lc.T.3.ax','lc.Z.1.eq','mc.chi.0.all',
            'mc.chi.2.all', 'mc.Z.0.all', 'mc.Z.1.all', 'mc.Z.3.all', 'mc.S.0.all',
            'f.chi.2.all','D_lc.chi.3.ax','D_lc.Z.1.ax','D_lc.Z.3.ax','D_lc.T.2.ax',
            'D_lc.S.3.ax','D_lc.chi.2.eq','D_mc.chi.1.all','D_mc.chi.2.all','D_mc.chi.3.all',
            'D_mc.S.1.all']
save_rac_names = [rac_name_convert(x) for x in rac_names]

X_train = reduceddf[rac_names]

# PCA pipeline with standard scaler
pipeline = Pipeline([('scaling',StandardScaler()),('pca',PCA(n_components=2))])
out = pipeline.fit_transform(X_train)

# Save model
pickle.dump(pipeline,open('PCA_model.pkl','wb'))

# Save RAC names in order
with open('rac_names.txt','w') as file1:
    for rname in save_rac_names:
        file1.write(rname+'\n')


# Save dataset
sses = reducedlabels.values.flatten()
pc1 = out[:,0].flatten()
pc2 = out[:,1].flatten()
newdf = pd.DataFrame({
    'Chemical Name':savenames,
    'Spin Splitting Energy (kcal/mol)':sses,
    'PC1':pc1,
    'PC2':pc2})

newdf.to_csv('pca_for_bokeh.csv')
