import pandas as pd
import glob


csvs = glob.glob('pat/*.csv'))

# generator expression to read files, assig(), to create
# a new colunm and concatenate through concat()
pd.concat((pd.read_csv(file).assign(filename = file) for file in csvs), ignore_index = True)
