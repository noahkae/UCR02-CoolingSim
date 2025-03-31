import pandas as pd
import numpy as np

def halve_spreadsheet(input_file, output_file, has_header=True):
    # Read the CSV file
    df = pd.read_csv(input_file, header=0 if has_header else None)
    
    # If the number of rows is odd, drop the last row
    if len(df) % 2 != 0:
        df = df.iloc[:-1]
    
    # Reshape the DataFrame so we can take the mean of every two rows
    df_avg = df.groupby(np.arange(len(df)) // 2).mean().reset_index(drop=True)
    
    # Save the new DataFrame
    df_avg.to_csv(output_file, index=False, header=has_header)
    
    print(f"Processed file saved as {output_file}")

# Example usage
halve_spreadsheet(r"data/shortened_velocity.csv", r"data/extra_shortened_velocity.csv")