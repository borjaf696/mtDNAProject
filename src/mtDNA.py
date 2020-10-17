# Borja :)
import sys
import pandas as pd

control_lsu = '../control/DV-gen-LSU.csv'

def preprocess_df(df):
    '''
    Method to adjust all columns to 'sensitive' formats
    Columns:
        * Change - Mutation - Reference
        * GB Freq‡ - change to float
        * GB Seqs
    '''
    mutation, reference = [],[]
    for i in df['Nucleotide Change']:
        changes = i.strip().split('-')
        mutation.append(changes[0])
        reference.append(changes[1])
    df['Mutated Base'] = mutation
    df['Reference Base'] = reference
    frequencies = []
    for i in df['GB Freq‡']:
        frequencies.append(float(i[:-1]))
    df['GB Freq‡'] = frequencies

def add_cv(df_section, df_control):
    '''
    Method to add cv:
        * Column: Genomic
    '''
    df_agg = df_control.groupby('Genomic').first()
    cvtot, cvmit, appearences, positions = [], [], [],[]
    for index, row in df_control.iterrows():
        genomic_pos = row['Genomic']
        df_tmp = df_section.loc[df_section['Position'] == genomic_pos]
        for i2, r2 in df_tmp.iterrows():
            cvtot.append(row['CVTOT'])
            cvmit.append(row['CVMIT'])
            appearences.append(r2['GB Seqs'])
            positions.append(r2['Position'])
        if df_tmp.empty:
            cvtot.append(row['CVTOT'])
            cvmit.append(row['CVMIT'])
            appearences.append(0)
            positions.append(row['position'])
    return pd.DataFrame({'Position':positions,'CVTOT':cvtot, 'CVMIT':cvmit,'GB Seqs':appearences})
    df_tmp = df_control.loc[~df_control.Genomic.isin(df_agg.index)]
    new_pot_df = {'Position':[],'GB Seqs':[],'Mutated':[],'Reference Base':[],'CVTOT':[],'CVMIT':[]}

if __name__ == '__main__':
    mitomap, tsv_read = [], False
    print('python mtDNA.py --mitomap comma_separated_files --somatic comma_separated_files')
    df_control_lsu = pd.read_csv(control_lsu, sep = ';')
    df_control_lsu['Genomic'] = pd.to_numeric(df_control_lsu['Genomic'].fillna(-1))
    print(df_control_lsu)
    for i in sys.argv:
        if tsv_read:
            mitomap = i.strip().split(',')
            tsv_read = False
        if i == '--mitomap':
            tsv_read = True
    ssu_df, lsu_df = pd.DataFrame(), pd.DataFrame()
    if len(mitomap) != 0:
        print('MitoMap pages: ', mitomap)
        for i in mitomap:
            if 'RNR1' in i:
                ssu_df = pd.read_csv(i, sep = '\t')
            elif 'RNR2' in i:
                lsu_df = pd.read_csv(i, sep = '\t')
        # Process ssu
        preprocess_df(ssu_df)
        print('SSU MitoMap information: ',ssu_df)
        # Process lsu
        preprocess_df(lsu_df)
        results_df = add_cv(lsu_df, df_control_lsu)
        results_df.to_csv('../output/lsu_df.csv')
        print('LSU MitoMap information: ', results_df)