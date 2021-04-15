import pandas as pd 
import sys

df = pd.read_csv(sys.argv[1])
original_train = df.sample(frac=0.8,random_state=0)
test = df[~df['refcode'].isin(original_train['refcode'])]
train = original_train.sample(frac=0.8,random_state=0)
val = original_train[~original_train['refcode'].isin(train['refcode'])]
print(df.shape, train.shape, test.shape,val.shape)

train.to_csv('train.csv',index=False)
val.to_csv('val.csv',index=False)
test.to_csv('test.csv',index=False)