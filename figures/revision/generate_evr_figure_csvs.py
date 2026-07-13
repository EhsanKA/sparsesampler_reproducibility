#!/usr/bin/env python

"""
Generate CSV tables for the RF classification panel (middle row) of the
three EVR figures (50k, 100k, 200k).

Each CSV reports, per feature index:
  - Delta F1 (SPS - Random) for macro F1 and each rare cell type
  - Standard deviation of both SPS and Random F1 across the 5 splits
"""

import os
import numpy as np
import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))

DATASET_ORDER = ['mcc', 'mcc_05', 'mcc_01', 'lcmv']

DATASETS = {
    'mcc': {
        'label': 'MCC (1%)',
        'rf_path': os.path.join(project_root, 'jobs/feature_index_classification/results/all_feature_indices_summary.csv'),
        'rare_cell_types': ['osteoblast'],
        'rf_delta_cols': ['osteoblast_f1_improvement'],
        'rf_sps_std_cols': ['sps_osteoblast_f1_std'],
        'rf_random_std_cols': ['random_osteoblast_f1_std'],
    },
    'mcc_05': {
        'label': 'MCC (0.5%)',
        'rf_path': os.path.join(project_root, 'jobs/feature_index_classification/results_mcc_05/all_feature_indices_summary.csv'),
        'rare_cell_types': ['osteoblast'],
        'rf_delta_cols': ['osteoblast_f1_improvement'],
        'rf_sps_std_cols': ['sps_osteoblast_f1_std'],
        'rf_random_std_cols': ['random_osteoblast_f1_std'],
    },
    'mcc_01': {
        'label': 'MCC (0.1%)',
        'rf_path': os.path.join(project_root, 'jobs/feature_index_classification/results_mcc_01/all_feature_indices_summary.csv'),
        'rare_cell_types': ['osteoblast'],
        'rf_delta_cols': ['osteoblast_f1_improvement'],
        'rf_sps_std_cols': ['sps_osteoblast_f1_std'],
        'rf_random_std_cols': ['random_osteoblast_f1_std'],
    },
    'lcmv': {
        'label': 'LCMV',
        'rf_path': os.path.join(project_root, 'jobs/feature_index_classification_lcmv/results/all_feature_indices_summary.csv'),
        'rare_cell_types': ['interacting', 'NK1_1_TCRgd_T', 'cDC2', 'pDCs', 'CD4_LCMV_spec'],
        'rf_delta_cols': [
            'interacting_f1_improvement', 'NK1_1_TCRgd_T_f1_improvement',
            'cDC2_f1_improvement', 'pDCs_f1_improvement', 'CD4_LCMV_spec_f1_improvement'
        ],
        'rf_sps_std_cols': [
            'sps_interacting_f1_std', 'sps_NK1_1_TCRgd_T_f1_std',
            'sps_cDC2_f1_std', 'sps_pDCs_f1_std', 'sps_CD4_LCMV_spec_f1_std'
        ],
        'rf_random_std_cols': [
            'random_interacting_f1_std', 'random_NK1_1_TCRgd_T_f1_std',
            'random_cDC2_f1_std', 'random_pDCs_f1_std', 'random_CD4_LCMV_spec_f1_std'
        ],
    },
}

EVR_RANGE = list(range(1, 31))


def generate_csv(dataset_key):
    cfg = DATASETS[dataset_key]
    if not os.path.exists(cfg['rf_path']):
        print(f"  WARNING: RF data not found: {cfg['rf_path']}")
        return None

    rf_data = pd.read_csv(cfg['rf_path']).sort_values('feature_index').reset_index(drop=True)

    rows = []
    for fi in EVR_RANGE:
        row = {'feature_index': fi}
        if fi in rf_data['feature_index'].values:
            rf_row = rf_data[rf_data['feature_index'] == fi].iloc[0]

            row['delta_macro_f1_pct'] = round(rf_row.get('macro_f1_improvement', np.nan), 4)
            row['sps_macro_f1_std'] = round(rf_row.get('sps_macro_f1_std', np.nan), 6)
            row['random_macro_f1_std'] = round(rf_row.get('random_macro_f1_std', np.nan), 6)

            for rare_ct, delta_col, sps_std_col, rand_std_col in zip(
                cfg['rare_cell_types'],
                cfg['rf_delta_cols'],
                cfg['rf_sps_std_cols'],
                cfg['rf_random_std_cols'],
            ):
                row[f'delta_{rare_ct}_f1_pct'] = round(rf_row.get(delta_col, np.nan), 4)
                row[f'sps_{rare_ct}_f1_std'] = round(rf_row.get(sps_std_col, np.nan), 6)
                row[f'random_{rare_ct}_f1_std'] = round(rf_row.get(rand_std_col, np.nan), 6)

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    output_dir = os.path.join(script_dir, 'evr_figure_csvs')
    os.makedirs(output_dir, exist_ok=True)

    for dataset_key in DATASET_ORDER:
        cfg = DATASETS[dataset_key]
        print(f"  {cfg['label']}...", end=' ')

        df = generate_csv(dataset_key)
        if df is not None:
            filename = f'{dataset_key}_rf_classification.csv'
            path = os.path.join(output_dir, filename)
            df.to_csv(path, index=False)
            print(f"saved ({len(df.columns)} cols)")
        else:
            print("skipped")

    print(f"\nDone. Output: {output_dir}")


if __name__ == "__main__":
    main()
