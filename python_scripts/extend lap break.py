import pandas as pd

# Load the CSV into a DataFrame
df = pd.read_csv("data/shortened_velocity.csv")

# Ensure there are enough rows
if len(df) > 40000:
    start = (len(df) - 40000) // 2  # Calculate start index for middle 80k rows
    end = start + 40000  # Calculate end index

    # Set all values to zero except the first column
    df.iloc[start:end, 1:] = 0

# Save the modified DataFrame if needed
df.to_csv("gapped.csv", index=False)
