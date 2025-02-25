import sqlalchemy 
import pandas as pd
import os
from sqlalchemy import create_engine

work_dir = os.getcwd()
print(work_dir)

data_dir = os.path.join(work_dir,'_data')
print(data_dir)

df = pd.read_csv('_data\\airlines.csv')
df
