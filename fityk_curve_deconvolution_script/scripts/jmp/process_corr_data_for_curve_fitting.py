# Replaces diagonal of ones along corr matrix with avg of closest values

import jmp
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

assert len(dataframe.columns) == len(dataframe.rows), "Correlation matrix must be square"







