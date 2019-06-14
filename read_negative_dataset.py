import pandas as pd 

df=pd.read_hdf("negative_dataset.h5",key="ndf")

print(df.info())
print(df.loc[33333:33433])