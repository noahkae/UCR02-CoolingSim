import pandas as pd
import numpy as np

# load the Excel file
file_path = r"data\lapsim_velocity.xlsx"  
sheet_name = "Sim Output Data"

# read the sheet
df = pd.read_excel(file_path, sheet_name=sheet_name)

# find time increment
avg_diff = df.iloc[:250, 0].diff().mean()

# Compute necessary values
num_rows = len(df)

first_col = df.iloc[:, 0] 
other_cols = df.iloc[:, 1:]  

new_num_rows = num_rows * 19  

print(other_cols)

new_first_col = df.iloc[-1, 0] + avg_diff * np.arange(new_num_rows)

new_first_col = np.concatenate([df.iloc[:, 0].to_numpy(), new_first_col]) 

data_repeated = pd.concat([other_cols] * 20, ignore_index=True)

df_expanded = pd.DataFrame(np.column_stack((new_first_col, data_repeated)), columns=df.columns)

df_expanded.to_csv("data\lapsim_velocity_extended-pause.csv", index=False)






