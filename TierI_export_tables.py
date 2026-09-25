from pathlib import Path
import argparse
import shutil
import pandas as pd

SUBTYPE_COLUMNS = [
    'ARG_subtype', 'ARG_type', 'Hospital_core_candidate', 'Nonhospital_additional_candidate',
    'Cities_detected_any_setting', 'Cities_sample_recurrent_any_setting', 'Cities_site_recurrent_any_setting',
    'Cities_sample_recurrent_nonhospital', 'Cities_site_recurrent_nonhospital',
    'Months_detected_any_setting', 'Months_sample_recurrent_nonhospital', 'Months_site_recurrent_nonhospital',
    'All_sites_sample_rec_both_same_setting', 'All_sites_site_rec_both_same_setting',
    'Common_sites_sample_rec_both_same_setting', 'Common_sites_site_rec_both_same_setting',
    'Common_sites_site_rec_both_same_nonhospital_setting',
    'Balanced_both_same_setting_count', 'Balanced_both_same_nonhospital_setting_count', 'Balanced_repetitions',
]


def main(argv=None):
    parser = argparse.ArgumentParser(description='Export the numerical tables for Supplementary Data S5-2 to S5-6 from saved analysis outputs; no model fitting or resampling.')
    parser.add_argument('--threshold-dir', type=Path, required=True)
    parser.add_argument('--analysis-dir', type=Path, required=True)
    parser.add_argument('--input-dir', type=Path, required=True,
                        help='Directory containing analysis_metadata.csv')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    sub = pd.read_csv(args.analysis_dir / 'subtype_consistency_summary.csv')
    meta = pd.read_csv(args.input_dir / 'analysis_metadata.csv')
    reps = pd.read_csv(args.analysis_dir / 'balanced_subsampling_repetitions.csv')
    if sub.ARG_subtype.duplicated().any() or meta.Sample.duplicated().any():
        raise ValueError('Duplicated subtype or sample IDs')
    n_panel = len(sub)
    n_hospital = int(sub.Hospital_core_candidate.sum())
    n_additional = int(sub.Nonhospital_additional_candidate.sum())
    if not sub.Balanced_repetitions.eq(len(reps)).all():
        raise ValueError('Subtype and repetition files have different run lengths')
    if sub.Balanced_both_same_setting_count.sum() != reps.Both_periods_same_setting_fixed102.sum():
        raise ValueError('Per-subtype and aggregate persistence counts disagree')
    nh = sub.Nonhospital_additional_candidate.eq(1)
    if sub.loc[nh, 'Balanced_both_same_nonhospital_setting_count'].sum() != reps.Both_periods_same_nonhospital_setting_additional10.sum():
        raise ValueError('Non-hospital persistence counts disagree')
    chosen = sub[SUBTYPE_COLUMNS].copy()
    chosen['Balanced_same_setting_fraction'] = chosen.Balanced_both_same_setting_count / chosen.Balanced_repetitions
    chosen['Balanced_same_nonhospital_setting_fraction'] = chosen.Balanced_both_same_nonhospital_setting_count / chosen.Balanced_repetitions
    chosen = chosen.sort_values(['Nonhospital_additional_candidate', 'ARG_subtype'], ascending=[False, True])
    chosen.to_csv(args.output / 'S5-4_subtype_consistency.csv', index=False, encoding='utf-8-sig')
    common = set(meta.groupby('PhysicalSite').Period.nunique().loc[lambda v: v.eq(2)].index)
    schemes = []
    for basis, prefix in [('All available sites', 'All_sites'), ('Sites sampled in both periods', 'Common_sites')]:
        m = meta if prefix == 'All_sites' else meta.loc[meta.PhysicalSite.isin(common)]
        ns = m.groupby('Period').size().reindex(['2024-11 to 2025-05', '2025-06 to 2025-12'])
        if ns.isna().any():
            raise ValueError('Both specified periods must be present')
        for weighting, metric in [('Sample-weighted', 'sample_rec'), ('Equal-site-weighted', 'site_rec')]:
            yes = sub[f'{prefix}_{metric}_both_same_setting'].eq(1)
            outside = sub[f'{prefix}_{metric}_both_same_nonhospital_setting'].eq(1)
            schemes.append([basis, weighting, int(m.PhysicalSite.nunique()), int(ns.iloc[0]), int(ns.iloc[1]),
                            int(yes.sum()), int((yes & sub.Hospital_core_candidate.eq(1)).sum()),
                            int((outside & nh).sum())])
    columns = ['Site basis', 'Weighting', 'Sites in either period', 'Early samples', 'Late samples',
               f'Persistent candidates / {n_panel}', f'Persistent hospital-core candidates / {n_hospital}',
               f'Persistent non-hospital additions / {n_additional}']
    pd.DataFrame(schemes, columns=columns).to_csv(args.output / 'S5-6_temporal_schemes.csv', index=False, encoding='utf-8-sig')
    sources = {
        'S5-2_threshold_summary.csv': args.threshold_dir / 'TierI_threshold_summary.csv',
        'S5-3_subtype_membership.csv': args.threshold_dir / 'TierI_subtype_membership.csv',
        'S5-5_stratum_profiles.csv': args.analysis_dir / 'subtype_stratum_profiles.csv',
        'S5-6_balanced_subsampling.csv': args.analysis_dir / 'balanced_subsampling_summary.csv',
        'S5-6_period_concordance.csv': args.analysis_dir / 'period_concordance.csv',
        'S5-6_matched_sample_coverage.csv': args.analysis_dir / 'stratum_coverage.csv',
    }
    for name, source in sources.items():
        shutil.copyfile(source, args.output / name)
    print(f'Exported S5-2 to S5-6 numerical tables: {args.output}')


if __name__ == '__main__':
    main()
