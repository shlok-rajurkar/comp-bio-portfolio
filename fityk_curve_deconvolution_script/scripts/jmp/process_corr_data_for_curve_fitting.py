# Replaces diagonal of ones along corr matrix with avg of closest values

import jmp
import jmputils

jmputils.jpip('install --upgrade', 'pip setuptools')
jmputils.jpip('install', 'numpy pandas')

import numpy
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

assert len(dataframe.columns) == len(dataframe[0]), "Correlation matrix must be square"

# Column processing functions
# Replace [0] in first column with [1]
# Replace [i] in middle columns with np.mean([i-1], [i+1])
# Replace [len.column] in last column with 







