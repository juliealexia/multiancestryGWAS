def mrmega_analysis(dfafr, dfamr, dfeur, dfeas, dfsas, dfmid):

    AFR_mrmega = pd.DataFrame({
        'MARKERNAME': dfafr['ID'],
        'EA': dfafr['ALT'],
        'NEA': dfafr['REF'],
        'BETA': dfafr['BETA'],
        'SE': dfafr['SE'],
        'EAF': dfafr['A1_FREQ'],
        'N': dfafr['OBS_CT'],
        'CHROMOSOME': dfafr['#CHROM'],
        'POSITION': dfafr['POS'],
        'P': dfafr['P']
    })

    AMR_mrmega = pd.DataFrame({
        'MARKERNAME': dfamr['ID'],
        'EA': dfamr['ALT'],
        'NEA': dfamr['REF'],
        'BETA': dfamr['BETA'],
        'SE': dfamr['SE'],
        'EAF': dfamr['A1_FREQ'],
        'N': dfamr['OBS_CT'],
        'CHROMOSOME': dfamr['#CHROM'],
        'POSITION': dfamr['POS'],
        'P': dfamr['P']
    })

    EUR_mrmega = pd.DataFrame({
        'MARKERNAME': dfeur['ID'],
        'EA': dfeur['ALT'],
        'NEA': dfeur['REF'],
        'BETA': dfeur['BETA'],
        'SE': dfeur['SE'],
        'EAF': dfeur['A1_FREQ'],
        'N': dfeur['OBS_CT'],
        'CHROMOSOME': dfeur['#CHROM'],
        'POSITION': dfeur['POS'],
        'P': dfeur['P']
    })
    
    EAS_mrmega = pd.DataFrame({
        'MARKERNAME': dfeas['ID'],
        'EA': dfeas['ALT'],
        'NEA': dfeas['REF'],
        'BETA': dfeas['BETA'],
        'SE': dfeas['SE'],
        'EAF': dfeas['A1_FREQ'],
        'N': dfeas['OBS_CT'],
        'CHROMOSOME': dfeas['#CHROM'],
        'POSITION': dfeas['POS'],
        'P': dfeas['P']
    })
    
    SAS_mrmega = pd.DataFrame({
        'MARKERNAME': dfsas['ID'],
        'EA': dfsas['ALT'],
        'NEA': dfsas['REF'],
        'BETA': dfsas['BETA'],
        'SE': dfsas['SE'],
        'EAF': dfsas['A1_FREQ'],
        'N': dfsas['OBS_CT'],
        'CHROMOSOME': dfsas['#CHROM'],
        'POSITION': dfsas['POS'],
        'P': dfsas['P']
    })
    
    MID_mrmega = pd.DataFrame({
        'MARKERNAME': dfmid['ID'],
        'EA': dfmid['ALT'],
        'NEA': dfmid['REF'],
        'BETA': dfmid['BETA'],
        'SE': dfmid['SE'],
        'EAF': dfmid['A1_FREQ'],
        'N': dfmid['OBS_CT'],
        'CHROMOSOME': dfmid['#CHROM'],
        'POSITION': dfmid['POS'],
        'P': dfmid['P']
    })
    MID_mrmega = MID_mrmega.dropna()
    SAS_mrmega = SAS_mrmega.dropna()

    AFR_mrmega.to_csv('AFR_MRMEGA.txt', sep='\t', index=False)
    AMR_mrmega.to_csv('AMR_MRMEGA.txt', sep='\t', index=False)
    EUR_mrmega.to_csv('EUR_MRMEGA.txt', sep='\t', index=False)
    EAS_mrmega.to_csv('EAS_MRMEGA.txt', sep='\t', index=False)
    SAS_mrmega.to_csv('SAS_MRMEGA.txt', sep='\t', index=False)
    MID_mrmega.to_csv('MID_MRMEGA.txt', sep='\t', index=False)
    
    command = f"./MR-MEGA -i mrmega.in --qt --pc 3 --out mrmegares"

    subprocess.run(command, shell=True, check=True)
    
    mrmegasres = pd.read_csv('mrmegares.result',sep='\t')
    mrmegasres = mrmegasres.dropna(subset=['beta_0'])
    
    return mrmegasres


def meta_analysis(dfafr, dfamr, dfeur, dfeas, dfsas, dfmid):
    ids_df1 = dfafr['ID']
    ids_df2 = dfamr['ID']
    ids_df3 = dfeur['ID']
    ids_df4 = dfeas['ID']
    ids_df5 = dfsas['ID']
    ids_df6 = dfmid['ID']
    unique_ids = pd.concat([ids_df1, ids_df2, ids_df3, ids_df4, ids_df5, ids_df6]).unique()

    unique_ids_df = pd.DataFrame({'ID': unique_ids})
    
    combined_df = pd.merge(unique_ids_df,dfafr[['ID','OBS_CT','BETA','SE']],on='ID', how='left')
    combined_df.rename(columns={'OBS_CT': 'OBS_CT_AFR'},inplace=True)
    combined_df.rename(columns={'BETA': 'BETA_AFR'},inplace=True)
    combined_df.rename(columns={'SE': 'SE_AFR'},inplace=True)
    combined_df = pd.merge(combined_df,dfamr[['ID','OBS_CT','BETA','SE']],on='ID', how='left')
    combined_df.rename(columns={'OBS_CT': 'OBS_CT_AMR'},inplace=True)
    combined_df.rename(columns={'BETA': 'BETA_AMR'},inplace=True)
    combined_df.rename(columns={'SE': 'SE_AMR'},inplace=True)
    combined_df = pd.merge(combined_df,dfeur[['ID','OBS_CT','BETA','SE']],on='ID', how='left')
    combined_df.rename(columns={'OBS_CT': 'OBS_CT_EUR'},inplace=True)
    combined_df.rename(columns={'BETA': 'BETA_EUR'},inplace=True)
    combined_df.rename(columns={'SE': 'SE_EUR'},inplace=True)
    combined_df = pd.merge(combined_df,dfeas[['ID','OBS_CT','BETA','SE']],on='ID', how='left')
    combined_df.rename(columns={'OBS_CT': 'OBS_CT_EAS'},inplace=True)
    combined_df.rename(columns={'BETA': 'BETA_EAS'},inplace=True)
    combined_df.rename(columns={'SE': 'SE_EAS'},inplace=True)
    combined_df = pd.merge(combined_df,dfsas[['ID','OBS_CT','BETA','SE']],on='ID', how='left')
    combined_df.rename(columns={'OBS_CT': 'OBS_CT_SAS'},inplace=True)
    combined_df.rename(columns={'BETA': 'BETA_SAS'},inplace=True)
    combined_df.rename(columns={'SE': 'SE_SAS'},inplace=True)
    combined_df = pd.merge(combined_df,dfmid[['ID','OBS_CT','BETA','SE']],on='ID', how='left')
    combined_df.rename(columns={'OBS_CT': 'OBS_CT_MID'},inplace=True)
    combined_df.rename(columns={'BETA': 'BETA_MID'},inplace=True)
    combined_df.rename(columns={'SE': 'SE_MID'},inplace=True)
    
    combined_df['temp_afr'] = combined_df['BETA_AFR']/(combined_df['SE_AFR']**2)
    combined_df['temp_amr'] = combined_df['BETA_AMR']/(combined_df['SE_AMR']**2)
    combined_df['temp_eur'] = combined_df['BETA_EUR']/(combined_df['SE_EUR']**2)
    combined_df['temp_eas'] = combined_df['BETA_EAS']/(combined_df['SE_EAS']**2)
    combined_df['temp_sas'] = combined_df['BETA_SAS']/(combined_df['SE_SAS']**2)
    combined_df['temp_mid'] = combined_df['BETA_MID']/(combined_df['SE_MID']**2)
    
    combined_df['temp2_afr'] = 1/(combined_df['SE_AFR']**2)
    combined_df['temp2_amr'] = 1/(combined_df['SE_AMR']**2)
    combined_df['temp2_eur'] = 1/(combined_df['SE_EUR']**2)
    combined_df['temp2_eas'] = 1/(combined_df['SE_EAS']**2)
    combined_df['temp2_sas'] = 1/(combined_df['SE_SAS']**2)
    combined_df['temp2_mid'] = 1/(combined_df['SE_MID']**2)
    
    combined_df['temp_all'] = combined_df[['temp_afr', 'temp_amr', 'temp_eur', 'temp_eas', 'temp_sas', 'temp_mid']].sum(axis=1, skipna=True)
    combined_df['temp2_all'] = combined_df[['temp2_afr', 'temp2_amr', 'temp2_eur', 'temp2_eas', 'temp2_sas', 'temp2_mid']].sum(axis=1, skipna=True)
    
    combined_df['BETA'] = combined_df['temp_all']/combined_df['temp2_all']
    combined_df['SE'] = np.sqrt(1/combined_df['temp2_all'])
    combined_df['T_STAT'] = combined_df['BETA']/combined_df['SE']
    combined_df['P'] = 2 * (1 - norm.cdf(np.abs(combined_df['T_STAT'])))
    
    meta_tab = combined_df[['ID','BETA','SE','T_STAT','P']].copy()
    
    return meta_tab

