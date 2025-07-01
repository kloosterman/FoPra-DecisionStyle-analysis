import pandas as pd
import numpy as np
from scipy.stats import norm

# Load data
df = pd.read_csv("Neurodose bereinigte Daten.csv")

# Identify valid FI and KI items (those with both _confidence and _correct)
fi_conf = [col for col in df.columns if col.startswith("FI_") and "_confidence" in col]
fi_corr = [col for col in df.columns if col.startswith("FI_") and "_correct" in col]
ki_conf = [col for col in df.columns if col.startswith("KI_") and "_confidence" in col]
ki_corr = [col for col in df.columns if col.startswith("KI_") and "_correct" in col]

valid_fi = sorted(set(col.replace("_confidence", "").replace("_correct", "") for col in fi_conf) &
                  set(col.replace("_confidence", "").replace("_correct", "") for col in fi_corr))
valid_ki = sorted(set(col.replace("_confidence", "").replace("_correct", "") for col in ki_conf) &
                  set(col.replace("_confidence", "").replace("_correct", "") for col in ki_corr))

# SDT functions
def type2_dprime(correct, conf, threshold=50):
    high_conf = conf >= threshold
    low_conf = conf < threshold
    hit = np.sum(high_conf & correct)
    miss = np.sum(low_conf & correct)
    fa = np.sum(high_conf & ~correct)
    cr = np.sum(low_conf & ~correct)
    hit_rate = (hit + 0.5) / (hit + miss + 1)
    fa_rate = (fa + 0.5) / (fa + cr + 1)
    return norm.ppf(hit_rate) - norm.ppf(fa_rate)

def type1_dprime(correct):
    hits = np.sum(correct == 1)
    misses = np.sum(correct == 0)
    total = hits + misses
    hit_rate = (hits + 0.5) / (total + 1)
    fa_rate = (misses + 0.5) / (total + 1)
    return norm.ppf(hit_rate) - norm.ppf(fa_rate)

# Compute all measures
records = []

for _, row in df.iterrows():
    rec = {'Codes': row['Codes']}
    
    # FI
    fi_correct = row[[f"{item}_correct" for item in valid_fi]].values.astype(float)
    fi_conf = row[[f"{item}_confidence" for item in valid_fi]].values.astype(float)
    rec['FI_accuracy_new'] = np.nanmean(fi_correct)
    rec['FI_mean_conf_new'] = np.nanmean(fi_conf)
    rec['FI_Type2_dprime_new'] = type2_dprime(fi_correct == 1, fi_conf)
    rec['FI_TypeI_dprime_new'] = type1_dprime(fi_correct)
    rec['FI_EmpiricistIndex_new'] = np.nanmean(fi_conf * fi_correct) - np.nanmean(fi_conf)
    rec['FI_EmpiricistIndex_simple_new'] = np.nanmean(fi_conf * fi_correct)

    # KI
    ki_correct = row[[f"{item}_correct" for item in valid_ki]].values.astype(float)
    ki_conf = row[[f"{item}_confidence" for item in valid_ki]].values.astype(float)
    rec['KI_accuracy_new'] = np.nanmean(ki_correct)
    rec['KI_mean_conf_new'] = np.nanmean(ki_conf)
    rec['KI_Type2_dprime_new'] = type2_dprime(ki_correct == 1, ki_conf)
    rec['KI_TypeI_dprime_new'] = type1_dprime(ki_correct)
    rec['KI_EmpiricistIndex_new'] = np.nanmean(ki_conf * ki_correct) - np.nanmean(ki_conf)
    rec['KI_EmpiricistIndex_simple_new'] = np.nanmean(ki_conf * ki_correct)

    records.append(rec)

# Create final DataFrame
recomputed_df = pd.DataFrame(records)

# Merge and export
final_df = df.merge(recomputed_df, on="Codes", how="left")
final_df.to_csv("Neurodose_Daten_with_all_recomputed_measures.csv", index=False)

