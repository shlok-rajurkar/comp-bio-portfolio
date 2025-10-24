import jmp
import numpy
import pandas as pd

# Load current jmp table as df

df = jmp.current()

# Load df into pandas dataframe

def jmp_to_pandas(jmp_table):
    pandasdf = pd.DataFrame()
    for col in jmp_table[:len(jmp_table)]:
        pandasdf[col.name] = list(col)
    return pandasdf

jmp_to_pandas(df)





