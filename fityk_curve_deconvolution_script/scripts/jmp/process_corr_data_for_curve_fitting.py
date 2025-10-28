# Replaces diagonal of ones along corr matrix with avg of closest values

import jmp
import jmputils
import sys

jmputils.jpip('install --upgrade', 'pip setuptools')
jmputils.jpip('install', 'numpy pandas')

import numpy as np
import pandas as pd

# Load current jmp table as df

dt = jmp.current()

# Load df into pandas dataframe

def jmp_to_pandas(jmp_table):
    pandasdf = pd.DataFrame()
    for col in jmp_table[:len(jmp_table)]:
        pandasdf[col.name] = list(col)
    return pandasdf

dataframe = jmp_to_pandas(dt)

# Ensure table is square

print(dataframe)

# assert len(dataframe.columns) == len(dataframe[0]), "Correlation matrix must be square"

# Column processing functions
# Replace [0] in first column with [1]

def process_first_col(col):
    col[0] = col[1]

# Replace [i] in middle columns with np.mean([i-1], [i+1])

def process_middle_col(col):
    for i in np.arange(len(col)):
        if col[i] == 1:
            col[i] = np.mean([col[i-1], col[i+1]])

# Replace last value in last column with second to last value

def process_last_col(col):
    col[len(col) - 1] = col[len(col) - 2]

def process_df(pd_df):
    process_first_col(pd_df.iloc[:, 0])
    for i in np.arange(0, len(pd_df.columns) - 2):
        process_middle_col(pd_df.iloc[:, i+1])
    process_last_col(pd_df.iloc[:, len(pd_df.columns) - 1])

process_df(dataframe)

print(dataframe)

dataframe.to_csv('Z:/User Folders/SRajurkar/Processed CSV/' + str(dt) + '.csv',index=False)


