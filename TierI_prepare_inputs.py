from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

SETTINGS = ['Hospital', 'WWTP', 'Community', 'Wet market']
CORE_COLUMNS = ['hospital_core', 'wwtp_core', 'community_core', 'wet_market_core']


def require(frame, columns, label):
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f'{label}: missing columns {sorted(missing)}')


def check_workbooks(paths, matrix):
    import openpyxl
    seen, blanks, sources = set(), 0, []
    for path in paths:
        book = openpyxl.load_workbook(path, read_only=True, data_only=True)
        try:
            for sheet in book:
                sheet.reset_dimensions()
                rows = sheet.iter_rows(values_only=True)
                header = next(rows)
                if str(header[0]).strip() != 'Subtype':
                    raise ValueError(f'{path.name}/{sheet.title}: first column must be Subtype')
                samples = [str(v).strip() for v in header[1:]]
                if len(samples) != len(set(samples)) or seen.intersection(samples):
                    raise ValueError('Duplicated sample columns across workbooks or sheets')
                seen.update(samples)
                data = {}
                for row in rows:
                    if not any(v is not None for v in row):
                        continue
                    name = str(row[0]).strip()
                    if name in data or name not in matrix.columns:
                        raise ValueError(f'Duplicate or unexpected subtype: {name}')
                    cells = list(row[1:])
                    if len(cells) > len(samples):
                        raise ValueError('More abundance cells than sample columns')
                    cells += [None] * (len(samples) - len(cells))
                    blanks += sum(v is None or str(v).strip() == '' for v in cells)
                    data[name] = [0.0 if v is None or str(v).strip() == '' else float(v) for v in cells]
                if set(data) != set(matrix.columns):
                    raise ValueError('Workbook subtype roster differs from the fixed panel')
                frame = pd.DataFrame(data, index=samples).reindex(columns=matrix.columns)
                np.testing.assert_allclose(frame.to_numpy(), matrix.loc[samples].to_numpy(),
                                           atol=1e-14, rtol=1e-12)
                sources.append({'file': path.name, 'sheet': sheet.title,
                                'samples': len(samples), 'values_checked': frame.size})
        finally:
            book.close()
    if seen != set(matrix.index):
        raise ValueError('Workbook sample IDs differ from the canonical cohort')
    return blanks, sources


def main(argv=None):
    parser = argparse.ArgumentParser(description='Prepare fixed-panel inputs from canonical metadata, a long abundance table and the original Tier I roster.')
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--abundance-long', type=Path, required=True)
    parser.add_argument('--original-tier1', type=Path, required=True)
    parser.add_argument('--membership', type=Path, required=True,
                        help='TierI_subtype_membership.csv from the threshold script')
    parser.add_argument('--workbooks', type=Path, nargs='+', help='Optional original XLSX matrices for reconciliation; blank abundance cells must represent nondetection')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    meta = pd.read_csv(args.metadata)
    original = pd.read_csv(args.original_tier1)
    long = pd.read_csv(args.abundance_long)
    membership = pd.read_csv(args.membership)
    require(meta, ['Sample', 'City', 'Setting', 'PhysicalSite', 'Sample_Date', 'SamplingMonth'], 'metadata')
    require(original, ['arg_subtype'], 'original Tier I roster')
    require(long, ['Sample', 'Subtype', 'abd', 'Site', 'City'], 'long abundance table')
    require(membership, ['arg_subtype', 'arg_type', 'original_TierI', 'nonhospital_additional_core', *CORE_COLUMNS], 'membership')
    for frame, key in [(meta, 'Sample'), (original, 'arg_subtype'), (membership, 'arg_subtype')]:
        if frame[key].isna().any() or frame[key].astype(str).str.strip().eq('').any() or frame[key].duplicated().any():
            raise ValueError(f'Missing or duplicated {key}')
    if meta[['City', 'Setting', 'PhysicalSite']].isna().any().any():
        raise ValueError('Missing sample metadata')
    if set(meta.Setting) != set(SETTINGS):
        raise ValueError(f'Expected setting labels: {SETTINGS}')
    if (meta.groupby('PhysicalSite')[['City', 'Setting']].nunique() > 1).any().any():
        raise ValueError('A physical site maps to more than one city or setting')
    dates = pd.to_datetime(meta.Sample_Date, format='%Y-%m-%d', errors='raise')
    if not dates.between('2024-11-01', '2025-12-31').all():
        raise ValueError('Sampling dates fall outside the specified analysis window')
    if not dates.dt.strftime('%Y-%m').eq(meta.SamplingMonth).all():
        raise ValueError('SamplingMonth does not match Sample_Date')
    names = original.arg_subtype.tolist()
    if set(long.Sample) != set(meta.Sample) or set(long.Subtype) != set(names):
        raise ValueError('Long-table sample or subtype IDs differ from the supplied cohort and fixed roster')
    if long.duplicated(['Sample', 'Subtype']).any():
        raise ValueError('Duplicate sample-subtype records')
    matrix = long.pivot(index='Sample', columns='Subtype', values='abd').reindex(index=meta.Sample, columns=names)
    values = matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Abundances must be complete, finite and non-negative; missing values are not zeros')
    matched = long[['Sample', 'Site', 'City']].drop_duplicates().merge(
        meta, on='Sample', validate='one_to_one', suffixes=('_long', '_official'))
    if not (matched.Site.eq(matched.Setting) & matched.City_long.eq(matched.City_official)).all():
        raise ValueError('Long-table setting or city disagrees with canonical metadata')
    traits = membership.set_index('arg_subtype').loc[names].copy()
    if not traits.original_TierI.eq(1).all() or set(membership.loc[membership.original_TierI.eq(1), 'arg_subtype']) != set(names):
        raise ValueError('The supplied roster differs from baseline Tier I membership')
    baseline, mismatch = {}, []
    for setting, key in zip(SETTINGS, CORE_COLUMNS):
        a = matrix.loc[meta.loc[meta.Setting.eq(setting), 'Sample']]
        n = len(a); positive = a.gt(0).sum(); mean = a.mean()
        recurrent = (positive * 10 >= n * 7) & mean.ge(1e-5)
        baseline[setting] = {'samples': n, 'sites': int(meta.loc[meta.Setting.eq(setting), 'PhysicalSite'].nunique()),
                             'original_core_candidates': int(traits[key].sum()), 'recomputed_recurrent': int(recurrent.sum())}
        for name in names:
            if int(recurrent[name]) != traits.loc[name, key]:
                mismatch.append({'Subtype': name, 'Setting': setting, 'original': int(traits.loc[name, key]),
                                 'recomputed': int(recurrent[name]), 'prevalence': positive[name] / n, 'mean': mean[name]})
    blanks, sheets = check_workbooks(args.workbooks, matrix) if args.workbooks else (None, [])
    meta['Period'] = np.where(dates.lt('2025-06-01'), '2024-11 to 2025-05', '2025-06 to 2025-12')
    inputs = {'metadata': args.metadata, 'abundance_long': args.abundance_long,
              'original_tier1': args.original_tier1, 'membership': args.membership}
    for i, path in enumerate(args.workbooks or []):
        inputs[f'workbook_{i+1}'] = path
    audit = {'samples': len(meta), 'sites': int(meta.PhysicalSite.nunique()), 'candidates': len(names),
             'complete_values': matrix.size, 'workbook_reconciliation_performed': bool(args.workbooks),
             'xlsx_blank_nondetections': blanks, 'source_sheets': sheets, 'baseline': baseline,
             'baseline_flag_mismatches': mismatch,
             'inputs': {key: {'file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for key, path in inputs.items()}}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'matrix_audit.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
    if mismatch:
        raise ValueError('Full-period core flags do not reproduce the fixed roster; see matrix_audit.json')
    matrix.reset_index().to_csv(args.output / 'analysis_abundance_matrix.csv', index=False, encoding='utf-8-sig')
    meta.to_csv(args.output / 'analysis_metadata.csv', index=False, encoding='utf-8-sig')
    traits.reset_index().to_csv(args.output / 'fixed_panel_traits.csv', index=False, encoding='utf-8-sig')
    print(f'Prepared inputs: {args.output}')


if __name__ == '__main__':
    main()
