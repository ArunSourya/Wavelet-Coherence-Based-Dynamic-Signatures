#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# to determine variability signature 


# In[ ]:


# to extract the discharge data 
import os
import pandas as pd
input_folder = r"E:\camels_de\13837553\camels_de\germany_timeseries"
output_folder = r"E:\camels_de\13837553\camels_de\germany_timeseries_IQR"

os.makedirs(output_folder, exist_ok=True)

for file in os.listdir(input_folder):
    if file.endswith(".csv"):
        input_path = os.path.join(input_folder, file)

        df = pd.read_csv(input_path)

        df = df[["date", "discharge_mm"]]

        df = df.dropna(subset=["discharge_mm"])         
        df = df[df["discharge_mm"].astype(str).str.strip() != ""]  
        output_path = os.path.join(output_folder, file)

        df.to_csv(output_path, index=False)

        print(f"Processed: {file}")

print("All files processed and saved.")



# to determine the inter quartile range values of every month in year
import os
import pandas as pd
import numpy as np

input_folder = r"E:\camels_de\13837553\camels_de\germany_timeseries_IQR"
output_folder = r"E:\phd_results\objective_1\germany\germany_Variability_interpretation"
os.makedirs(output_folder, exist_ok=True)

max_window = 30      
max_use_years = 45   

month_names = ["Jan","Feb","Mar","Apr","May","Jun",
               "Jul","Aug","Sep","Oct","Nov","Dec"]

csv_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]

for file in csv_files:
    file_path = os.path.join(input_folder, file)
    catchment_id = file.split("_")[4]

    df = pd.read_csv(file_path)
    df["date"] = pd.to_datetime(df["date"], format="%Y-%m-%d")
    df.rename(columns={"discharge_mm": "Flow"}, inplace=True)

    df["Year"] = df["date"].dt.year
    df["Month"] = df["date"].dt.month
    df["DOY"] = df["date"].dt.dayofyear

    start_year = df["Year"].min()
    final_year = df["Year"].max()
    available_years = final_year - start_year + 1
    use_years = min(available_years, max_use_years)
    end_year = start_year + use_years - 1

    df_use = df[df["Year"] <= end_year].copy()

    global_mean = df_use["Flow"].mean()

    iqr_matrix = np.full((max_window, 366), np.nan)

    for m in range(1, max_window + 1):

        df_use[f"roll_mean_{m}"] = df_use["Flow"].rolling(
            window=m, center=True, min_periods=1
        ).mean()

        for doy in range(1, 367):
            vals = df_use.loc[df_use["DOY"] == doy, f"roll_mean_{m}"].dropna().values
            if len(vals) > 0:
                q75, q25 = np.percentile(vals, [75, 25])
                iqr_matrix[m-1, doy-1] = (q75 - q25) / global_mean

    
    doy_month_map = df_use.groupby("DOY")["Month"].first().to_dict()

    monthly_iqr = pd.DataFrame(index=range(1, max_window + 1), columns=month_names)

    for m in range(1, max_window + 1):
        for month in range(1, 13):
            
            month_days = [d for d, mo in doy_month_map.items() if mo == month]
            vals = [iqr_matrix[m-1, d-1] for d in month_days if not np.isnan(iqr_matrix[m-1, d-1])]

            if len(vals) > 0:
                monthly_iqr.loc[m, month_names[month-1]] = np.mean(vals)
            else:
                monthly_iqr.loc[m, month_names[month-1]] = np.nan

    out_path = os.path.join(output_folder, f"{catchment_id}_Monthly_IQR.csv")
    monthly_iqr.to_csv(out_path, index_label="m_window_days")

    print(f"Saved Monthly IQR: {out_path}")

print("All Monthly IQR files generated.")



# In[ ]:


# this is to find the mean of variability (IQR) across months
import os
import pandas as pd

input_folder = r"E:\phd_results\objective_1\germany\germany_Variability_interpretation"

output_file = r"E:\phd_results\objective_1\germany\germany_Variability_interpretation.csv"

csv_files = [f for f in os.listdir(input_folder) if f.endswith("_Monthly_IQR.csv")]

final_df = pd.DataFrame()

month_cols = ["Jan","Feb","Mar","Apr","May","Jun",
              "Jul","Aug","Sep","Oct","Nov","Dec"]

rows = []

for file in csv_files:

    file_path = os.path.join(input_folder, file)
    station_id = file.split("_")[0]   

    df = pd.read_csv(file_path)

    monthly_avg = df[month_cols].mean()

    row = {"Station": station_id}
    for month in month_cols:
        row[month] = monthly_avg[month]

    rows.append(row)

final_df = pd.DataFrame(rows)

final_df.to_csv(output_file, index=False)

print("Combined monthly IQR averages saved to:")
print(output_file)


# In[ ]:


# the below code is to determine difference between max variability - minimum variability, second maximum variability - second minimum variability  for every catchment
import pandas as pd
import numpy as np

input_file = r"E:\phd_results\objective_1\germany\germany_Variability_interpretation.csv"
output_file = r"E:\phd_results\objective_1\germany\germany_Variability_interpretation_max_and_second_max_diff.csv"

df = pd.read_csv(input_file)

df.rename(columns={df.columns[0]: 'Station'}, inplace=True)

monthly_data = df.iloc[:, 1:]

values = monthly_data.to_numpy()

sorted_vals = np.sort(values, axis=1)

max_min_diff = sorted_vals[:, -1] - sorted_vals[:, 0]
second_max_min_diff = sorted_vals[:, -2] - sorted_vals[:, 1]

result = pd.DataFrame({
    'Station': df['Station'],
    'max_minus_min': max_min_diff,
    'second_max_minus_second_min': second_max_min_diff
})

result.to_csv(output_file, index=False)

print("File saved successfully.")


# In[ ]:


# to determine apportionment entropy
import os
import pandas as pd
import numpy as np

input_folder = r"E:\camels_de\13837553\camels_de\germany_timeseries"
output_file = r"E:\phd_results\objective_1\germany\germany_AE_values.csv"

years = range(1981, 2015)
months = range(1, 13)

def compute_daily_AE(month_values, annual_total):
    
    if annual_total <= 0:
        return np.nan

    p = month_values / annual_total
    p = p[p > 0]

    if len(p) == 0:
        return np.nan

    return -np.sum(p * np.log2(p))


results = []

for file in os.listdir(input_folder):
    if not file.endswith(".csv"):
        continue

    file_path = os.path.join(input_folder, file)
    catchment = file.split("_")[-1]

    df = pd.read_csv(file_path)

    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df.dropna(subset=["date"])

    df["year"] = df["date"].dt.year
    df["month"] = df["date"].dt.month

    row = {"catchment": catchment}

    precip_monthly_AE = {m: [] for m in months}
    qobs_monthly_AE = {m: [] for m in months}

    for yr in years:
        df_year = df[df["year"] == yr]

        if df_year.empty:
            continue

        annual_precip = df_year["precipitation_mm"].sum(skipna=True)
        annual_qobs = df_year["discharge_mm"].sum(skipna=True)

        for m in months:
            df_month = df_year[df_year["month"] == m]

            precip_vals = df_month["precipitation_mm"].dropna()
            AE_precip = compute_daily_AE(precip_vals, annual_precip)
            precip_monthly_AE[m].append(AE_precip)

            qobs_vals = df_month["discharge_mm"].dropna()
            AE_qobs = compute_daily_AE(qobs_vals, annual_qobs)
            qobs_monthly_AE[m].append(AE_qobs)

    for m in months:
        row[f"AE_precip_mean_{m:02d}"] = np.nanmean(precip_monthly_AE[m])
        row[f"AE_qobs_mean_{m:02d}"] = np.nanmean(qobs_monthly_AE[m])

    results.append(row)

AE_df = pd.DataFrame(results)
AE_df.to_csv(output_file, index=False)

print("Monthly climatological AE (Precip & Qobs) computed successfully.")
print(f"Results saved to:\n{output_file}")


# In[ ]:


# the absolute difference between the AE values for catchments
import pandas as pd

input_file = r"E:\phd_results\objective_1\germany\germany_AE_values.csv"
output_file = r"E:\phd_results\objective_1\germany\germany_mean_diff_AE_values.csv"

df = pd.read_csv(input_file)

df.rename(columns={df.columns[0]: 'catchment'}, inplace=True)

diff_cols = []

for i in range(1, 13):
    m = f"{i:02d}"
    precip_col = f"AE_precip_mean_{m}"
    qobs_col = f"AE_qobs_mean_{m}"
    
    diff_col = f"diff_{m}"
    df[diff_col] = (df[precip_col] - df[qobs_col]).abs()
    diff_cols.append(diff_col)

df['Mean_DIFF_AE_value'] = df[diff_cols].mean(axis=1)

result = df[['catchment', 'Mean_DIFF_AE_value']]

result.to_csv(output_file, index=False)

print("Mean absolute AE difference saved successfully.")


# In[ ]:

# to extract the lag signatures

import os
import pandas as pd
import numpy as np
input_folder = r"E:\phd_results\objective_1\germany\lag_filtered_coh_gt_0p8"
output_csv = r"E:\phd_results\objective_1\germany\mean_lag_overscales1.csv"

rows = []

for filename in os.listdir(input_folder):
    if not filename.endswith("_lag_deg_coh_gt_0p8.csv"):
        continue

    catchment = filename.replace("_lag_deg_coh_gt_0p8.csv", "")
    file_path = os.path.join(input_folder, filename)

    df = pd.read_csv(file_path)

    df_abs = df.abs()

    col_means = df_abs.mean(skipna=True)
    row = col_means.to_dict()
    row["catchment_name"] = catchment

    rows.append(row)
mean_df = pd.DataFrame(rows)
cols = ["catchment_name"] + [c for c in mean_df.columns if c != "catchment_name"]
mean_df = mean_df[cols]
mean_df.to_csv(output_csv, index=False)

print(f"Mean ABSOLUTE lag values saved to:\n{output_csv}")

# this is to extract lag value for each scale
import pandas as pd

input_file = r"E:\phd_results\objective_1\germany\mean_lag_overscales1.csv"
output_file = r"E:\phd_results\objective_1\germany\max_lag_overscales1.csv"

df = pd.read_csv(input_file)

df.rename(columns={df.columns[0]: 'catchment'}, inplace=True)

scale_cols = df.columns[1:]

scale_map = {}
for col in scale_cols:
    try:
        scale_map[col] = float(col)
    except:
        continue

lag1_cols = [col for col, val in scale_map.items() if 2.066087295 <= val <= 6.949461602]
lag2_cols = [col for col, val in scale_map.items() if 7.362698093 <= val <= 13.8989232]
lag3_cols = [col for col, val in scale_map.items() if 14.72539619 <= val <= 29.45079237]
lag4_cols = [col for col, val in scale_map.items() if 31.20202762 <= val <= 176.5053225]
lag5_cols = [col for col, val in scale_map.items() if 187.0008752 <= val <= 353.010645]

result = pd.DataFrame()
result['catchment'] = df['catchment']

result['lag1'] = df[lag1_cols].max(axis=1)
result['lag2'] = df[lag2_cols].max(axis=1)
result['lag3'] = df[lag3_cols].max(axis=1)
result['lag4'] = df[lag4_cols].max(axis=1)
result['lag5'] = df[lag5_cols].max(axis=1)

result.to_csv(output_file, index=False)

print("Max lag values saved successfully.")


# RMI signature for the scale range 

import pandas as pd
import re

input_file = r"E:\phd_results\objective_1\germany\mean_coherence_values.csv"
output_file = r"E:\phd_results\objective_1\germany\mean_coherence_values_rmi.csv"

df = pd.read_csv(input_file)

catchment = df["Catchment"]
data = df.drop(columns=["Catchment"])

periods = data.columns.str.extract(r"P_(\d+\.?\d*)")[0].astype(float)

group1_cols = data.columns[(periods >= 2.07) & (periods <= 6.95)]
group2_cols = data.columns[(periods >= 7.36) & (periods <= 14.73)]
group3_cols = data.columns[(periods >= 15.60) & (periods <= 58.90)]
group4_cols = data.columns[(periods >= 62.40) & (periods <= 176.51)]
group5_cols = data.columns[(periods >= 187.00) & (periods <= 353.01)]

g1 = data[group1_cols].mean(axis=1)
g2 = data[group2_cols].mean(axis=1)
g3 = data[group3_cols].mean(axis=1)
g4 = data[group4_cols].mean(axis=1)
g5 = data[group5_cols].mean(axis=1)
g6 = g1+g2+g3+g4+g5
result = pd.DataFrame({
    "Catchment": catchment,
    "Mean_P_2.07_to_6.95": g1,
    "Mean_P_7.36_to_14.73": g2,
    "Mean_P_15.60_to_58.90": g3,
    "Mean_P_62.40_to_176.51": g4,
    "Mean_P_187.00_to_353.01": g5,
    "G1_by_G6": g1 / g6,
    "G2_by_G6": g2 / g6,
    "G3_by_G6": g3 / g6,
    "G4_by_G6": g4 / g6,
    "G5_by_G6": g5 / g6
})

result.to_csv(output_file, index=False)

print("Correct RMI file saved to:", output_file)


