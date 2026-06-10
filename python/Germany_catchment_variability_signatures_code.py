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

#

# to extract the lag signatures

import os
import pandas as pd
import numpy as np
input_folder = r"E:\phd_results\objective_1\germany\lag_days_all"
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

# to combine all extracted Signatures 

import os
import pandas as pd

input_folder = r"E:\phd_results\objective_1\germany\features"

output_file = r"E:\phd_results\objective_1\germany\germany_ALL_features_combined.csv"

csv_files = [f for f in os.listdir(input_folder) if f.endswith("_features.csv")]

all_rows = []

for file in csv_files:

    file_path = os.path.join(input_folder, file)
    station_id = file.split("_")[4]  
    df = pd.read_csv(file_path)
    if "Station" not in df.columns:
        df.insert(0, "Station", station_id)
    df["Station"] = df["Station"].astype(str)
    all_rows.append(df)

combined_df = pd.concat(all_rows, ignore_index=True)
combined_df.to_csv(output_file, index=False)


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


# to determine flashiness index

import os
import pandas as pd
import numpy as np

input_folder = r"E:\camels_de\13837553\camels_de\germany_timeseries_IQR"
output_folder = r"E:\phd_results\objective_1\germany\germany_FLASHINESSINDEX"

os.makedirs(output_folder, exist_ok=True)

csv_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]

all_summary = []

for file in csv_files:

    file_path = os.path.join(input_folder, file)
    catchment_id = file.split("_")[4]

    df = pd.read_csv(file_path)

    df["date"] = pd.to_datetime(df["date"], format="%Y-%m-%d")
    df.rename(columns={"discharge_mm": "Flow"}, inplace=True)

    df = df.sort_values("date")

    df["Year"] = df["date"].dt.year
    df["Month"] = df["date"].dt.month

    df["water_year"] = df["Year"]
    df.loc[df["Month"] >= 10, "water_year"] += 1

    yearly_results = []

    for year, group in df.groupby("water_year"):

        group = group.sort_values("date")

        diff = group["Flow"].diff().abs().sum()
        total_flow = group["Flow"].sum()

        if total_flow == 0 or np.isnan(total_flow):
            rb = np.nan
        else:
            rb = diff / total_flow

        yearly_results.append({
            "water_year": year,
            "RB_index": rb
        })

    yearly_df = pd.DataFrame(yearly_results)

    yearly_output = os.path.join(output_folder, f"{catchment_id}_RB_yearly.csv")
    yearly_df.to_csv(yearly_output, index=False)

    mean_rb = yearly_df["RB_index"].mean()
    median_rb = yearly_df["RB_index"].median()

    all_summary.append({
        "catchment_id": catchment_id,
        "RB_mean": mean_rb,
        "RB_median": median_rb
    })

    print(f" Processed: {catchment_id}")

summary_df = pd.DataFrame(all_summary)

summary_output = os.path.join(output_folder, "ALL_CATCHMENTS_RB_SUMMARY.csv")
summary_df.to_csv(summary_output, index=False)

print(" DONE!")
print(f"Summary file saved at: {summary_output}")


# the below codes is to draw plots 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import MaxNLocator

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 1.2,
    "lines.linewidth": 1.8,
})

phase_csv = r"E:\phd_results\objective_1\germany\plots_for_attributes.csv"
lag_csv   = r"E:\phd_results\objective_1\germany\plots_for_attributes_lag.csv"

out_phase = r"E:\phd_results\objective_1\germany\plots\phase"
out_lag   = r"E:\phd_results\objective_1\germany\plots\lag"

os.makedirs(out_phase, exist_ok=True)
os.makedirs(out_lag, exist_ok=True)

MAX_PERIOD = 365

phase_df = pd.read_csv(phase_csv)
lag_df   = pd.read_csv(lag_csv)

for df in (phase_df, lag_df):
    df.rename(columns={df.columns[0]: "Catchment"}, inplace=True)
    df.set_index("Catchment", inplace=True)

periods = np.array(phase_df.columns, dtype=float)
valid_cols = periods <= MAX_PERIOD
periods = periods[valid_cols]

phase_df = phase_df.iloc[:, valid_cols]
lag_df   = lag_df.iloc[:, valid_cols]

for catchment in phase_df.index:

    phase = phase_df.loc[catchment].astype(float).values
    lag   = lag_df.loc[catchment].astype(float).values

    fig, ax = plt.subplots(figsize=(7.5, 4.5))

    ax.plot(periods, phase, color="black")
    ax.set_xlabel("Period (days)")
    ax.set_ylabel("Phase (degrees)")
    ax.set_title(f"Catchment {catchment}", loc="left")

    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))

    ax.grid(True, linestyle="--", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.savefig(
        os.path.join(out_phase, f"{catchment}.jpg"),
        format="jpg",
        dpi=600,
        bbox_inches="tight"
    )
    plt.close()

    fig, ax = plt.subplots(figsize=(7.5, 4.5))

    ax.plot(periods, lag, color="black")
    ax.set_xlabel("Period (days)")
    ax.set_ylabel("Lag (days)")
    ax.set_title(f"Catchment {catchment}", loc="left")

    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))

    ax.grid(True, linestyle="--", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.savefig(
        os.path.join(out_lag, f"{catchment}.jpg"),
        format="jpg",
        dpi=600,
        bbox_inches="tight"
    )
    plt.close()

    print(f"Saved paper-quality JPG plots for catchment {catchment}")

print("All catchments processed successfully.")



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import MaxNLocator

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 1.2,
    "lines.linewidth": 1.8,
})
input_csv = r"E:\phd_results\objective_1\germany\plots_for_attributes_mean.csv"
output_folder = r"E:\phd_results\objective_1\germany\plots\mean_coherence"

os.makedirs(output_folder, exist_ok=True)

df = pd.read_csv(input_csv)
df = df.rename(columns={df.columns[0]: "Station"})
df.set_index("Station", inplace=True)

x_labels = df.columns.tolist()
x_positions = np.arange(len(x_labels))  
for station in df.index:
    y_values = df.loc[station].values

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(x_positions, y_values, color="black")

    ax.set_xlabel("Period Range (days)")
    ax.set_ylabel("Coherence")
    ax.set_title(f"Station {station}", loc="left")

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45)

    ax.yaxis.set_major_locator(MaxNLocator(6))

    ax.grid(True, linestyle="--", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    save_path = os.path.join(output_folder, f"{station}.jpg")
    plt.savefig(save_path, dpi=600, bbox_inches="tight")
    plt.close()

    print(f"Saved: {save_path}")

print("All plots created successfully.")


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import MaxNLocator

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 1.2,
    "lines.linewidth": 1.8,
})
input_csv = r"E:\phd_results\objective_1\germany\plots_for_attributes_mean.csv"
output_folder = r"E:\phd_results\objective_1\germany\plots\mean_significance"

os.makedirs(output_folder, exist_ok=True)

df = pd.read_csv(input_csv)
df = df.rename(columns={df.columns[0]: "Station"})
df.set_index("Station", inplace=True)

x_labels = df.columns.tolist()
x_positions = np.arange(len(x_labels))  
for station in df.index:
    y_values = df.loc[station].values

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(x_positions, y_values, color="black")

    ax.set_xlabel("Period Range (days)")
    ax.set_ylabel("Significance (%)")
    ax.set_title(f"Station {station}", loc="left")

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45)

    ax.yaxis.set_major_locator(MaxNLocator(6))

    ax.grid(True, linestyle="--", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    save_path = os.path.join(output_folder, f"{station}.jpg")
    plt.savefig(save_path, dpi=600, bbox_inches="tight")
    plt.close()

    print(f"Saved: {save_path}")

print("All plots created successfully.")


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import MaxNLocator

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 1.2,
    "lines.linewidth": 1.8,
})
input_csv = r"E:\phd_results\objective_1\germany\plots_for_attributes_mean.csv"
output_folder = r"E:\phd_results\objective_1\germany\plots\mean_phase_conc"

os.makedirs(output_folder, exist_ok=True)

df = pd.read_csv(input_csv)
df = df.rename(columns={df.columns[0]: "Station"})
df.set_index("Station", inplace=True)

x_labels = df.columns.tolist()
x_positions = np.arange(len(x_labels)) 
for station in df.index:
    y_values = df.loc[station].values

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(x_positions, y_values, color="black")

    ax.set_xlabel("Period Range (days)")
    ax.set_ylabel("Phase Concentration")
    ax.set_title(f"Station {station}", loc="left")

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45)

    ax.yaxis.set_major_locator(MaxNLocator(6))

    ax.grid(True, linestyle="--", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    save_path = os.path.join(output_folder, f"{station}.jpg")
    plt.savefig(save_path, dpi=600, bbox_inches="tight")
    plt.close()

    print(f"Saved: {save_path}")

print("All plots created successfully.")


import pandas as pd
import matplotlib.pyplot as plt
import os

csv_file = r"E:\camels_de\13837553\camels_de\germany_timeseries\CAMELS_DE_hydromet_timeseries_DEG10010.csv"
out_dir  = r"E:\phd_results\objective_1\germany\plots\timeseries_plot"

os.makedirs(out_dir, exist_ok=True)

catchment_id = "DEG10010"

df = pd.read_csv(csv_file)

df["date"] = pd.to_datetime(df["date"], format="%Y-%m-%d")

fig, ax = plt.subplots(figsize=(14, 5))

ax.vlines(
    df["date"],
    ymin=0,
    ymax=df["precipitation_mm"],
    color="red",
    linewidth=0.8,
    label="Precipitation (mm)"
)

ax.plot(
    df["date"],
    df["discharge_mm"],
    color="blue",
    linewidth=1.2,
    label="Streamflow (mm)"
)

ax.set_xlabel("Date")
ax.set_ylabel("Values")

start_date = df["date"].min().strftime("%d/%m/%Y")
end_date   = df["date"].max().strftime("%d/%m/%Y")

ax.set_title(
    f"Daily Precipitation and Streamflow ({start_date} to {end_date})"
)

ax.grid(True, linestyle="-", alpha=0.6)
ax.legend(loc="upper right")

out_file = os.path.join(out_dir, f"{catchment_id}.jpg")

plt.tight_layout()
plt.savefig(
    out_file,
    format="jpg",
    dpi=600,
    bbox_inches="tight"
)
plt.close()

print(f"Saved 600 dpi JPG for catchment {catchment_id} to:\n{out_file}")



import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    "font.size": 14,
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 15
})

input_file = r"E:\camels_de\13837553\camels_de\germany_timeseries\CAMELS_DE_hydromet_timeseries_DEG10010.csv"
base_output_folder = r"E:\phd_results\objective_1\germany\plots\timeseries_plot\DEG10010"

hydro_folder = os.path.join(base_output_folder, "hydro")
temp_folder = os.path.join(base_output_folder, "temperature")

os.makedirs(hydro_folder, exist_ok=True)
os.makedirs(temp_folder, exist_ok=True)

df = pd.read_csv(input_file)
df.columns = df.columns.str.strip()
df['date'] = pd.to_datetime(df['date'], dayfirst=True, errors='coerce')
df = df.dropna(subset=['date'])

years = range(1981, 2015)

for year in years:
    df_year = df[df['date'].dt.year == year]
    if df_year.empty:
        continue

    fig, ax1 = plt.subplots(figsize=(14, 6))

    p_bar = ax1.bar(
        df_year['date'], df_year['precipitation_mm'],
        color='blue', label='Precipitation (mm)', width=1
    )
    
    ax1.set_ylabel('Precipitation (mm)', fontweight='bold')
    ax1.set_title(f"Daily Precipitation, and Streamflow ({year})", fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.spines['top'].set_visible(False)

    ax2 = ax1.twinx()
    q_line, = ax2.plot(
        df_year['date'], df_year['discharge_mm'],
        color='green', label='Streamflow (mm)', linewidth=2.5
    )
    ax2.set_ylabel('Streamflow (mm)', fontweight='bold', color='green')
    ax2.tick_params(axis='y', labelcolor='green')
    ax2.spines['top'].set_visible(False)

    ax1.set_xlabel('date', fontweight='bold')

    for ax in [ax1, ax2]:
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')

    handles = [p_bar, q_line]
    labels = ['Precipitation (mm)', 'Streamflow (mm)']

    legend = fig.legend(
        handles, labels,
        loc='lower center',
        ncol=3,
        frameon=False
    )

    for text in legend.get_texts():
        text.set_fontweight('bold')

    plt.subplots_adjust(bottom=0.22)

    hydro_path = os.path.join(hydro_folder, f"hydro_{year}.jpg")
    plt.savefig(hydro_path, dpi=1200, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(12, 5))

    ax.plot(df_year['date'], df_year['temp_max'], color='red', label='Tmax (°C)', linewidth=2.5)
    ax.plot(df_year['date'], df_year['temp_min'], color='blue', label='Tmin (°C)', linewidth=2.5)

    ax.set_xlabel("date", fontweight='bold')
    ax.set_ylabel("Temperature (°C)", fontweight='bold')
    ax.set_title(f"Daily Maximum and Minimum Temperature ({year})", fontweight='bold')

    legend = ax.legend(loc='upper right')
    for text in legend.get_texts():
        text.set_fontweight('bold')

    ax.grid(True)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    plt.tight_layout()

    temp_path = os.path.join(temp_folder, f"temperature_{year}.jpg")
    plt.savefig(temp_path, dpi=1200)
    plt.close()

print("All yearly plots created successfully.")


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 16,
    "axes.labelsize": 18,
    "axes.titlesize": 20,
    "axes.titleweight": "bold",
    "axes.labelweight": "bold",
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 17,
    "axes.linewidth": 1.4,
    "lines.linewidth": 2.8
})

input_file = r"E:\camels_de\13837553\camels_de\germany_timeseries\CAMELS_DE_hydromet_timeseries_DEG10010.csv"
out_dir = r"E:\phd_results\objective_1\germany\plots\time_series_plot_events\DEG10010"
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv(input_file)
df.columns = df.columns.str.strip()
df['date'] = pd.to_datetime(df['date'], dayfirst=True)

event_start = pd.to_datetime("1997-01-01")
event_end   = pd.to_datetime("2005-01-01")
event_label = f"{event_start.year}_{event_end.year}"

ma_windows = [256, 128, 64, 32, 16, 8]

for ma_window in ma_windows:

    half_win = ma_window // 2
    max_lag  = ma_window // 2

    df_ext = df[
        (df['date'] >= event_start - pd.Timedelta(days=half_win)) &
        (df['date'] <= event_end   + pd.Timedelta(days=half_win + max_lag))
    ].copy()

    df_ext['P_ma'] = df_ext['precipitation_mm'].rolling(
        window=ma_window, center=True, min_periods=ma_window
    ).mean()

    df_ext['Q_ma'] = df_ext['discharge_mm'].rolling(
        window=ma_window, center=True, min_periods=ma_window
    ).mean()

    P_event = df_ext[
        (df_ext['date'] >= event_start) &
        (df_ext['date'] <= event_end)
    ][['date', 'P_ma']].dropna().reset_index(drop=True)

    correlations = []

    for lag in range(max_lag + 1):
        Q_event = df_ext[
            (df_ext['date'] >= event_start + pd.Timedelta(days=lag)) &
            (df_ext['date'] <= event_end   + pd.Timedelta(days=lag))
        ][['date', 'Q_ma']].dropna().reset_index(drop=True)

        n = min(len(P_event), len(Q_event))
        correlations.append(
            P_event.iloc[:n]['P_ma'].corr(Q_event.iloc[:n]['Q_ma'])
        )

    best_lag  = int(np.nanargmax(correlations))
    best_corr = correlations[best_lag]

    Q_best = df_ext[
        (df_ext['date'] >= event_start + pd.Timedelta(days=best_lag)) &
        (df_ext['date'] <= event_end   + pd.Timedelta(days=best_lag))
    ][['date', 'Q_ma']].dropna().reset_index(drop=True)

    n = min(len(P_event), len(Q_best))
    P_plot = P_event.iloc[:n]
    Q_plot = Q_best.iloc[:n]

    fig, ax1 = plt.subplots(figsize=(14, 6))

    p_line, = ax1.plot(
        P_plot['date'], P_plot['P_ma'],
        color='red',
        label=f'{ma_window}-day Centered MA Precipitation'
    )
    ax1.set_ylabel("Precipitation (mm)", fontweight='bold', color='red')
    ax1.tick_params(axis='y', labelcolor='red')
    ax1.grid(True, linestyle='--', alpha=0.4)
    ax1.spines[['top']].set_visible(False)

    ax2 = ax1.twinx()
    q_line, = ax2.plot(
        Q_plot['date'], Q_plot['Q_ma'],
        color='blue',
        label=f'{ma_window}-day Centered MA Streamflow (Lag = {best_lag} days)'
    )
    ax2.set_ylabel("Streamflow (mm)", fontweight='bold', color='blue')
    ax2.tick_params(axis='y', labelcolor='blue')
    ax2.spines[['top']].set_visible(False)

    ax1.set_xlabel("date", fontweight='bold')
    ax1.set_title(
        f"Event {event_label} | MA = {ma_window} days | Best lag = {best_lag} days",
        fontweight='bold'
    )

    for ax in [ax1, ax2]:
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontweight('bold')

    lines = [p_line, q_line]
    labels = [l.get_label() for l in lines]

    legend = fig.legend(
        lines, labels,
        loc='lower center',
        ncol=2,
        frameon=False
    )

    for txt in legend.get_texts():
        txt.set_fontweight('bold')

    plt.subplots_adjust(bottom=0.22)

    ts_name = f"{event_label}_MA{ma_window}_bestlag{best_lag}_timeseries.jpg"
    plt.savefig(
        os.path.join(out_dir, ts_name),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    plt.figure(figsize=(7, 4))
    plt.plot(range(max_lag + 1), correlations, marker='o', color='black')
    plt.axvline(best_lag, color='red', linestyle='--',
                label=f'Best lag = {best_lag}')
    plt.xlabel("Lag (days)", fontweight='bold')
    plt.ylabel("Correlation", fontweight='bold')
    plt.title(f"Lagged Correlation | MA = {ma_window} days", fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.4)

    legend = plt.legend(frameon=False)
    for txt in legend.get_texts():
        txt.set_fontweight('bold')

    for tick in plt.gca().get_xticklabels() + plt.gca().get_yticklabels():
        tick.set_fontweight('bold')

    plt.tight_layout()
    corr_name = f"{event_label}_MA{ma_window}_correlation_vs_lag.jpg"
    plt.savefig(os.path.join(out_dir, corr_name), dpi=300)
    plt.close()

    print(f"Best lag = {best_lag} days | r = {best_corr:.4f}")

print("All MA window iterations completed successfully.")


import pandas as pd
import matplotlib.pyplot as plt
import os

input_file = r"E:\camels_de\13837553\camels_de\germany_timeseries\CAMELS_DE_hydromet_timeseries_DE110290.csv"
output_folder = r"E:\phd_results\objective_1\germany\plots\time_series_plot_events\DE110290"
os.makedirs(output_folder, exist_ok=True)

df = pd.read_csv(input_file)
df.columns = df.columns.str.strip()
df['date'] = pd.to_datetime(df['date'], format='%Y-%m-%d')

start_date = pd.to_datetime('1986-01-01')
end_date   = pd.to_datetime('1990-01-01')

df_date_range = df[(df['date'] >= start_date) & (df['date'] <= end_date)]

if df_date_range.empty:
    print(f"No data available between {start_date.date()} and {end_date.date()}")

else:

    fig, axs = plt.subplots(2, 1, figsize=(14, 7), sharex=True)

    axs[0].plot(df_date_range['date'], df_date_range['precipitation_mm'],
                color='red', label='Precipitation (mm)', lw=1)

    axs[0].set_ylabel('Precipitation')
    axs[0].set_title(
        f'Daily Hydro-Meteorological Variables\n{start_date.date()} to {end_date.date()}'
    )
    axs[0].grid(True, linestyle='--', alpha=0.4)
    axs[0].spines[['top', 'right']].set_visible(False)
    axs[0].legend(frameon=False)

    axs[1].plot(df_date_range['date'], df_date_range['discharge_mm'],
                color='blue', label='Streamflow (mm)', lw=1)

    axs[1].set_ylabel('Streamflow (mm)')
    axs[1].set_xlabel('Date')
    axs[1].grid(True, linestyle='--', alpha=0.4)
    axs[1].spines[['top', 'right']].set_visible(False)
    axs[1].legend(frameon=False)

    plt.subplots_adjust(hspace=0.15)

    output_path = os.path.join(
        output_folder,
        f"daily_precip_qobs_{start_date.date()}_to_{end_date.date()}.jpg"
    )

    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close()

    print(f"Plot saved successfully:\n{output_path}")


import os
import pandas as pd
import matplotlib.pyplot as plt

csv_file = r"E:\phd_results\objective_1\germany\plots_for_attributes_rmi_AE.csv"

output_dir = r"E:\phd_results\objective_1\germany\plots\rmi"

os.makedirs(output_dir, exist_ok=True)

df = pd.read_csv(csv_file)

rmi_cols = ["RMI1", "RMI2", "RMI3", "RMI4", "RMI5"]

for _, row in df.iterrows():
    gauge_id = str(row["Catchment"])
    rmi_values = row[rmi_cols].values

    plt.figure(figsize=(6, 4))
    plt.plot(rmi_cols, rmi_values, marker='o', linestyle='-')

    plt.xlabel("RMI")
    plt.ylabel("Ratio")
    plt.title("RMI")

    plt.grid(True)

    output_path = os.path.join(output_dir, f"{gauge_id}.jpg")
    plt.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close()

print("All catchment plots have been saved successfully.")

