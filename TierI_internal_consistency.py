from pathlib import Path
import argparse
import hashlib
import json
import platform
import numpy as np
import pandas as pd

SETTINGS = ['Hospital', 'WWTP', 'Community', 'Wet market']
PERIODS = ['2024-11 to 2025-05', '2025-06 to 2025-12']

def recurrent(prevalence, abundance):
    return (prevalence >= 0.70 - 1e-12) & (abundance >= 1e-5)

def main(argv=None):
    parser = argparse.ArgumentParser(description='Fixed Tier I panel: city, time, common-site and balanced-subsampling analyses; Supplementary Data S5-4 to S5-6.')
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repetitions', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=20260925)
    args = parser.parse_args(argv)
    if not __debug__:
        raise RuntimeError('Run Python without -O; input checks must remain enabled.')
    if args.repetitions < 1:
        parser.error('--repetitions must be positive')
    ROOT = args.output
    ROOT.mkdir(parents=True, exist_ok=True)
    REPETITIONS = args.repetitions
    RNG = np.random.default_rng(args.seed)
    paths = {name: args.input_dir / name for name in
             ['analysis_metadata.csv', 'fixed_panel_traits.csv', 'analysis_abundance_matrix.csv']}
    meta = pd.read_csv(paths['analysis_metadata.csv']).set_index('Sample')
    traits = pd.read_csv(paths['fixed_panel_traits.csv']).set_index('arg_subtype')
    names = list(traits.index)
    matrix = pd.read_csv(paths['analysis_abundance_matrix.csv']).set_index('Sample')
    required_meta = {'City', 'Setting', 'PhysicalSite', 'Sample_Date', 'SamplingMonth', 'Period'}
    core_columns = ['hospital_core', 'wwtp_core', 'community_core', 'wet_market_core']
    required_traits = {'arg_type', 'nonhospital_additional_core', *core_columns}
    if required_meta - set(meta.columns) or required_traits - set(traits.columns):
        raise ValueError('Missing metadata or trait columns; see INPUTS.md')
    for label, frame in [('metadata', meta), ('traits', traits), ('abundance', matrix)]:
        if frame.index.has_duplicates or frame.index.isna().any() or not len(frame):
            raise ValueError(f'{label}: missing or duplicate IDs, or empty input')
    if set(matrix.index) != set(meta.index) or set(matrix.columns) != set(names):
        raise ValueError('Abundance IDs differ from the supplied cohort and fixed panel')
    if meta[list(required_meta)].isna().any().any() or traits[list(required_traits)].isna().any().any():
        raise ValueError('Missing metadata or panel traits')
    if set(meta.Setting) != set(SETTINGS) or set(meta.Period) != set(PERIODS):
        raise ValueError('Unexpected setting or period labels')
    dates = pd.to_datetime(meta.Sample_Date, format='%Y-%m-%d', errors='raise')
    if not dates.between('2024-11-01', '2025-12-31').all():
        raise ValueError('Sampling dates fall outside the specified analysis window')
    expected_period = np.where(dates.lt('2025-06-01'), PERIODS[0], PERIODS[1])
    if not dates.dt.strftime('%Y-%m').eq(meta.SamplingMonth).all() or not (expected_period == meta.Period).all():
        raise ValueError('Period or month does not match Sample_Date')
    if (meta.groupby('PhysicalSite')[['City', 'Setting']].nunique() > 1).any().any():
        raise ValueError('A physical site maps to multiple cities or settings')
    if not traits[[*core_columns, 'nonhospital_additional_core']].isin([0, 1]).all().all():
        raise ValueError('Core-membership flags must be 0 or 1')
    expected_outside = traits.hospital_core.eq(0) & traits[core_columns[1:]].sum(axis=1).gt(0)
    if not expected_outside.eq(traits.nonhospital_additional_core.astype(bool)).all():
        raise ValueError('Non-hospital addition flags conflict with core membership')
    matrix = matrix.loc[meta.index, names]
    values = matrix.to_numpy(); detection = values > 0
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Abundances must be finite and non-negative; missing values are not zeros')
    outside = traits.nonhospital_additional_core.to_numpy().astype(bool)
    hosp = traits.hospital_core.to_numpy().astype(bool)
    hospital_only = hosp & traits[['wwtp_core','community_core','wet_market_core']].sum(axis=1).eq(0).to_numpy()
    shared = traits[['hospital_core','wwtp_core','community_core','wet_market_core']].sum(axis=1).eq(4).to_numpy()
    assert (hosp ^ outside).all()
    if not hosp.any() or not outside.any():
        raise ValueError('The panel comparison requires both hospital-core candidates and non-hospital additions')
    site_periods = meta.groupby('PhysicalSite').Period.nunique()
    common_sites = site_periods.index[site_periods.eq(2)]
    common_mask = meta.PhysicalSite.isin(common_sites).to_numpy()
    if set(meta.loc[common_mask, 'Setting']) != set(SETTINGS):
        raise ValueError('Each setting must have a physical site sampled in both periods')
    samples = []
    for i, (sample, row) in enumerate(meta.iterrows()):
        total = values[i].sum(); extra = values[i,outside].sum()
        samples.append({'Sample':sample,'City':row.City,'Setting':row.Setting,'PhysicalSite':row.PhysicalSite,
                        'SamplingMonth':row.SamplingMonth,'Period':row.Period,
                        'Detected_fixed102':int(detection[i].sum()),'Detected_hospital_core92':int(detection[i,hosp].sum()),
                        'Detected_nonhospital10':int(detection[i,outside].sum()),
                        'Abundance_fixed102':float(total),'Abundance_nonhospital10':float(extra),
                        'Nonhospital10_abundance_fraction':float(extra/total) if total else None})
    pd.DataFrame(samples).to_csv(ROOT/'sample_coverage.csv',index=False,encoding='utf-8-sig')
    details, summaries, strata_arrays = [], [], {}

    def add_stratum(scope, label, setting, idx, basis='All available sites'):
        if not len(idx):
            return
        m = meta.iloc[idx]; v = values[idx]; d = detection[idx]
        site_prev = pd.DataFrame(d, index=m.PhysicalSite).groupby(level=0).mean()
        site_mean = pd.DataFrame(v, index=m.PhysicalSite).groupby(level=0).mean()
        count = d.sum(axis=0); prev = count / len(idx); avg = v.mean(axis=0)
        eqprev = site_prev.mean(axis=0).to_numpy(); eqmean = site_mean.mean(axis=0).to_numpy()
        rec = recurrent(prev, avg); eqrec = recurrent(eqprev, eqmean)
        detected = count > 0
        key = (scope,label,setting,basis)
        strata_arrays[key] = {'detected':detected, 'sample_rec':rec, 'site_rec':eqrec,
                             'sample_prev':prev,'site_prev':eqprev,'sample_mean':avg,'site_mean':eqmean}
        for j, name in enumerate(names):
            details.append({'Scope':scope,'Stratum':label,'Setting':setting,'Site_basis':basis,'ARG_subtype':name,
                            'Hospital_core_candidate':int(hosp[j]),'Nonhospital_additional_candidate':int(outside[j]),
                            'Samples':len(idx),'Sites':len(site_prev),'Positive_samples':int(count[j]),
                            'Sample_prevalence':float(prev[j]),'Sample_mean_copies_per_cell':float(avg[j]),
                            'Equal_site_prevalence':float(eqprev[j]),'Equal_site_mean_copies_per_cell':float(eqmean[j]),
                            'Meets_sample_recurrence':int(rec[j]),'Meets_equal_site_recurrence':int(eqrec[j])})
        extra_counts = d[:,outside].sum(axis=1)
        site_extra = pd.DataFrame({'extra_count':extra_counts, 'any_extra':extra_counts>0,
                                   'extra_abundance':v[:,outside].sum(axis=1)},index=m.PhysicalSite).groupby(level=0).mean()
        summaries.append({'Scope':scope,'Stratum':label,'Setting':setting,'Site_basis':basis,
                          'Samples':len(idx),'Sites':len(site_prev),
                          'Detected_fixed102':int(detected.sum()),'Detected_hospital_core92':int(detected[hosp].sum()),
                          'Detected_nonhospital10':int(detected[outside].sum()),
                          'Sample_recurrent_fixed102':int(rec.sum()),'Sample_recurrent_hospital_core92':int(rec[hosp].sum()),
                          'Sample_recurrent_nonhospital10':int(rec[outside].sum()),
                          'Equal_site_recurrent_fixed102':int(eqrec.sum()),'Equal_site_recurrent_hospital_core92':int(eqrec[hosp].sum()),
                          'Equal_site_recurrent_nonhospital10':int(eqrec[outside].sum()),
                          'Samples_with_nonhospital_candidate':int((extra_counts>0).sum()),
                          'Sample_fraction_with_nonhospital_candidate':float((extra_counts>0).mean()),
                          'Mean_detected_hospital_core_candidates':float(d[:,hosp].sum(axis=1).mean()),
                          'Mean_detected_nonhospital_candidates':float(extra_counts.mean()),
                          'Mean_detected_fixed_candidates':float(d.sum(axis=1).mean()),
                          'Equal_site_mean_detected_nonhospital_candidates':float(site_extra.extra_count.mean()),
                          'Mean_nonhospital_candidate_abundance':float(v[:,outside].sum(axis=1).mean()),
                          'Equal_site_mean_nonhospital_candidate_abundance':float(site_extra.extra_abundance.mean())})

    all_idx = np.arange(len(meta))
    for setting in SETTINGS:
        add_stratum('Network','All months',setting,all_idx[meta.Setting.eq(setting)])
    for city in sorted(meta.City.unique()):
        for setting in SETTINGS:
            add_stratum('City',city,setting,all_idx[meta.City.eq(city)&meta.Setting.eq(setting)])
    for month in sorted(meta.SamplingMonth.unique()):
        for setting in SETTINGS:
            add_stratum('Month',month,setting,all_idx[meta.SamplingMonth.eq(month)&meta.Setting.eq(setting)])
    for period in PERIODS:
        for setting in SETTINGS:
            mask = (meta.Period.eq(period)&meta.Setting.eq(setting)).to_numpy()
            add_stratum('Period',period,setting,all_idx[mask])
            add_stratum('Period',period,setting,all_idx[mask&common_mask], 'Sites sampled in both periods')
    detail = pd.DataFrame(details); coverage = pd.DataFrame(summaries)
    detail.to_csv(ROOT/'subtype_stratum_profiles.csv',index=False,encoding='utf-8-sig')
    coverage.to_csv(ROOT/'stratum_coverage.csv',index=False,encoding='utf-8-sig')

    def across(scope, label, key, subset=None, basis='All available sites'):
        take = [s for k,s in strata_arrays.items() if k[0]==scope and k[1]==label and k[3]==basis and (subset is None or k[2] in subset)]
        assert take
        return np.any([s[key] for s in take],axis=0)

    sub_summary = []
    for j,name in enumerate(names):
        row={'ARG_subtype':name,'ARG_type':traits.loc[name,'arg_type'],'Hospital_core_candidate':int(hosp[j]),
             'Hospital_only_core_candidate':int(hospital_only[j]),'Four_setting_shared_candidate':int(shared[j]),
             'Nonhospital_additional_candidate':int(outside[j]),'Cities_observed':meta.City.nunique(),
             'Months_observed':meta.SamplingMonth.nunique()}
        for scope,plural,labels in [('City','Cities',sorted(meta.City.unique())),('Month','Months',sorted(meta.SamplingMonth.unique()))]:
            row[plural+'_detected_any_setting']=sum(bool(across(scope,x,'detected')[j]) for x in labels)
            row[plural+'_sample_recurrent_any_setting']=sum(bool(across(scope,x,'sample_rec')[j]) for x in labels)
            row[plural+'_site_recurrent_any_setting']=sum(bool(across(scope,x,'site_rec')[j]) for x in labels)
            row[plural+'_sample_recurrent_nonhospital']=sum(bool(across(scope,x,'sample_rec',SETTINGS[1:])[j]) for x in labels)
            row[plural+'_site_recurrent_nonhospital']=sum(bool(across(scope,x,'site_rec',SETTINGS[1:])[j]) for x in labels)
        for label,prefix in [('All available sites','All_sites'),('Sites sampled in both periods','Common_sites')]:
            for metric in ['sample_rec','site_rec']:
                early = across('Period',PERIODS[0],metric,basis=label)
                late = across('Period',PERIODS[1],metric,basis=label)
                row[f'{prefix}_{metric}_early_any_setting']=int(early[j])
                row[f'{prefix}_{metric}_late_any_setting']=int(late[j])
                both_same = [strata_arrays[('Period',PERIODS[0],s,label)][metric] & strata_arrays[('Period',PERIODS[1],s,label)][metric] for s in SETTINGS]
                row[f'{prefix}_{metric}_both_same_setting']=int(np.any(both_same,axis=0)[j])
                row[f'{prefix}_{metric}_both_same_nonhospital_setting']=int(np.any(both_same[1:],axis=0)[j])
        sub_summary.append(row)
    sub = pd.DataFrame(sub_summary)
    sub.to_csv(ROOT/'subtype_consistency_summary.csv',index=False,encoding='utf-8-sig')
    concordance=[]
    for basis in ['All available sites','Sites sampled in both periods']:
        for setting in SETTINGS:
            early=strata_arrays[('Period',PERIODS[0],setting,basis)]
            late=strata_arrays[('Period',PERIODS[1],setting,basis)]
            for weight in ['sample','site']:
                a,b = early[weight+'_rec'], late[weight+'_rec']
                both=a&b; union=a|b
                ranks_early = pd.Series(early[weight+'_prev']).rank(method='average')
                ranks_late = pd.Series(late[weight+'_prev']).rank(method='average')
                rho = ranks_early.corr(ranks_late) if ranks_early.nunique()>1 and ranks_late.nunique()>1 else None
                concordance.append({'Site_basis':basis,'Setting':setting,'Weighting':weight,
                                     'Early_recurrent':int(a.sum()),'Late_recurrent':int(b.sum()),
                                     'Both_periods_recurrent':int(both.sum()),'Either_period_recurrent':int(union.sum()),
                                     'Jaccard':float(both.sum()/union.sum()) if union.any() else None,
                                     'Prevalence_rank_correlation':float(rho) if rho is not None else None,
                                     'Both_periods_nonhospital_candidates':int(both[outside].sum())})
    concordance=pd.DataFrame(concordance)
    concordance.to_csv(ROOT/'period_concordance.csv',index=False,encoding='utf-8-sig')

    site_blocks = []
    site_manifest = []
    for site in sorted(common_sites):
        blocks = [all_idx[(meta.PhysicalSite.eq(site)&meta.Period.eq(p)).to_numpy()] for p in PERIODS]
        n = min(map(len,blocks)); setting = meta.loc[meta.PhysicalSite.eq(site),'Setting'].iloc[0]
        site_blocks.append((setting,blocks,n))
        site_manifest.append({'PhysicalSite':site,'City':meta.loc[meta.PhysicalSite.eq(site),'City'].iloc[0],
                              'Setting':setting,'Early_samples':len(blocks[0]),'Late_samples':len(blocks[1]),'Samples_selected_per_period':n})
    pd.DataFrame(site_manifest).to_csv(ROOT/'balanced_site_manifest.csv',index=False,encoding='utf-8-sig')
    reps=[]; selected_counts=np.zeros((len(names),2),dtype=int)
    for rep in range(REPETITIONS):
        prev={s:[] for s in SETTINGS}; avg={s:[] for s in SETTINGS}
        for setting,blocks,n in site_blocks:
            chosen=[RNG.choice(b,size=n,replace=False) for b in blocks]
            prev[setting].append(np.stack([detection[x].mean(axis=0) for x in chosen]))
            avg[setting].append(np.stack([values[x].mean(axis=0) for x in chosen]))
        rules=np.stack([recurrent(np.mean(prev[s],axis=0),np.mean(avg[s],axis=0)) for s in SETTINGS])
        both=rules[:,0,:]&rules[:,1,:]
        union_same=both.any(axis=0); nonhospital_same=both[1:].any(axis=0)
        selected_counts[:,0]+=union_same; selected_counts[:,1]+=nonhospital_same
        reps.append({'Repetition':rep+1,'Both_periods_same_setting_fixed102':int(union_same.sum()),
                     'Both_periods_same_setting_hospital_core92':int((union_same&hosp).sum()),
                     'Both_periods_same_nonhospital_setting_additional10':int((nonhospital_same&outside).sum()),
                     'Both_periods_in_hospital_original92':int((both[0]&hosp).sum()),
                     'Both_periods_in_hospital_original6_exclusive':int((both[0]&hospital_only).sum())})
    repframe=pd.DataFrame(reps)
    repframe.to_csv(ROOT/'balanced_subsampling_repetitions.csv',index=False,encoding='utf-8-sig')
    sub['Balanced_repetitions']=REPETITIONS
    sub['Balanced_both_same_setting_count']=selected_counts[:,0]
    sub['Balanced_both_same_nonhospital_setting_count']=selected_counts[:,1]
    sub.to_csv(ROOT/'subtype_consistency_summary.csv',index=False,encoding='utf-8-sig')
    repsum=[]
    for key in repframe.columns[1:]:
        r=repframe[key]
        repsum.append({'Metric':key,'Repetitions':REPETITIONS,'Minimum':int(r.min()),'P2.5':float(r.quantile(.025)),
                       'Median':float(r.median()),'P97.5':float(r.quantile(.975)),'Maximum':int(r.max())})
    pd.DataFrame(repsum).to_csv(ROOT/'balanced_subsampling_summary.csv',index=False,encoding='utf-8-sig')

    extra=sub[sub.Nonhospital_additional_candidate.eq(1)]
    summary={'analysis':'Fixed-panel internal descriptive consistency, not independent validation',
             'samples':len(meta),'sites':meta.PhysicalSite.nunique(),'cities':meta.City.nunique(),'months':meta.SamplingMonth.nunique(),
             'fixed_candidates':len(names),'hospital_core_candidates':int(hosp.sum()),'nonhospital_candidates':int(outside.sum()),
             'candidates_detected_in_all_cities':int(sub.Cities_detected_any_setting.eq(meta.City.nunique()).sum()),
             'candidates_detected_in_all_months':int(sub.Months_detected_any_setting.eq(meta.SamplingMonth.nunique()).sum()),
             'candidates_recurrent_in_all_cities_sample_weighted':int(sub.Cities_sample_recurrent_any_setting.eq(meta.City.nunique()).sum()),
             'candidates_recurrent_in_all_cities_site_weighted':int(sub.Cities_site_recurrent_any_setting.eq(meta.City.nunique()).sum()),
             'nonhospital_candidates_cities_recurrent_sample_range':[int(extra.Cities_sample_recurrent_nonhospital.min()),int(extra.Cities_sample_recurrent_nonhospital.max())],
             'nonhospital_candidates_cities_recurrent_site_range':[int(extra.Cities_site_recurrent_nonhospital.min()),int(extra.Cities_site_recurrent_nonhospital.max())],
             'nonhospital_candidates_months_recurrent_sample_range':[int(extra.Months_sample_recurrent_nonhospital.min()),int(extra.Months_sample_recurrent_nonhospital.max())],
             'common_sites':len(common_sites),'balanced_samples_per_period':sum(s['Samples_selected_per_period'] for s in site_manifest),
             'period_concordance':concordance.to_dict('records'),'balanced_subsampling':repsum,
             'nonhospital_candidate_detail':extra.to_dict('records')}
    for basis in ['All_sites','Common_sites']:
        for metric in ['sample_rec','site_rec']:
            col=f'{basis}_{metric}_both_same_setting'; nh=f'{basis}_{metric}_both_same_nonhospital_setting'
            summary[f'{col}_fixed102']=int(sub[col].sum())
            summary[f'{nh}_additional10']=int(extra[nh].sum())
    (ROOT/'results.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    run_settings = {'seed': args.seed, 'repetitions': REPETITIONS, 'periods': PERIODS,
                    'prevalence_threshold': 0.70, 'mean_abundance_threshold': 1e-5,
                    'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
                    'inputs': {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()}}
    (ROOT/'run_settings.json').write_text(json.dumps(run_settings, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['nonhospital_candidate_detail','period_concordance']},ensure_ascii=False,indent=2))
    print(extra[['ARG_subtype','Cities_sample_recurrent_nonhospital','Cities_site_recurrent_nonhospital','Months_sample_recurrent_nonhospital','Common_sites_site_rec_both_same_nonhospital_setting','Balanced_both_same_nonhospital_setting_count']].to_string(index=False))

if __name__ == '__main__':
    main()
