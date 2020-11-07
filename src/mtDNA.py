# Borja :)
import sys, re
import pandas as pd
from utils.utils import *
import urllib.request

control_lsu, haplogrupos_control = '../control/DV-gen-LSU.csv', '../control/haplogrupos.tsv'
control_ssu = '../control/DV-gen-SSU.csv'

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
        mutation.append(changes[1])
        reference.append(changes[0])
    df['Mutated Base'] = mutation
    df['Reference Base'] = reference
    frequencies = []
    for i in df['GB Freq‡']:
        frequencies.append(float(i[:-1]))
    df['GB Freq‡'] = frequencies

def add_cv(df_section, df_control, haplogroups = None):
    '''
    Method to add cv:
        * Column: Genomic
    '''
    df_control.groupby('Genomic').first()
    cvtot, cvmit, appearences, positions, genomic, corrected_appearences, cells_haplogroups = [], [], [], [], [], [], []
    for index, row in df_control.iterrows():
        genomic_pos = row['Genomic']
        df_tmp = df_section.loc[df_section['Position'] == genomic_pos]
        for i2, r2 in df_tmp.iterrows():
            cvtot.append(row['CVTOT'])
            appearences.append(max(r2['GB Seqs'], r2['Curated References']))
            corrected_appearences.append(max(r2['GB Seqs'] - r2['Nb.Seqs.correction'], r2['Curated References']))
            genomic.append(genomic_pos)
            cells_haplogroups.append(r2[haplogroups])
        if df_tmp.empty and genomic_pos != -1:
            cvtot.append(row['CVTOT'])
            appearences.append(0)
            corrected_appearences.append(0)
            genomic.append(genomic_pos)
            cells_haplogroups.append(np.zeros(len(haplogroups)))
    df_result = pd.DataFrame({'CVTOT':cvtot,'GB Seqs':appearences,'GB.Seqs.Corrected':corrected_appearences,'Genomic':genomic})
    if haplogroups is not None:
        return pd.concat([df_result, pd.DataFrame(data = cells_haplogroups, columns = haplogroups)], axis = 1)
    else:
        return df_result
'''
Obtención de haplogrupos para cada una de las secuencias:
    * Construcción de la query de MitoMap.
    * wget via subprocess
    * Parseo del html
    * mkdir htmls_haplogrupos
'''
def get_haplogroups(df, suffix = 'lsu', df_haplogroups = None):
    '''
    :param df: dataframe with information to build the query
    :return: same dataframe with new columns counting the number of times each haplogroup appears
    base query = https://www.mitomap.org/cgi-bin/index_mitomap.cgi?title=Coding+Polymorphism+G-A+at+rCRS+position+2701&pos=2701&ref=G&alt=A&purge_type=
    '''
    assert df_haplogroups is not None
    haplogrupos = df_haplogrupos['Top Level Haplogroup'].values
    base_query = 'https://www.mitomap.org/cgi-bin/index_mitomap.cgi?title=Coding+Polymorphism+\-\+at+rCRS+position+\&pos=\&ref=\&alt=\&purge_type='
    base_query_split = base_query.split('\\')
    dir_htmls = '../mitomap_html_dirs_'+suffix+'/'
    if not Utils.exists(dir_htmls):
        Utils.mkdir(dir_htmls)
    files = Utils.get_files(dir_htmls,extension=['html'])
    haplogroups_dict = {i:np.zeros(len(df.index)) for i in haplogrupos}
    columna_new_seqs = {'Nb.Seqs.correction':np.zeros(len(df.index))}
    df = pd.concat([df,pd.DataFrame(haplogroups_dict)], axis = 1)
    df = pd.concat([df, pd.DataFrame(columna_new_seqs)], axis = 1)
    print('Haplogroups dataframe columns: ', df_haplogroups.columns.values)
    print('LSU dataframe columns: ',df.columns.values)
    if not len(files) != 0:
        print('Retrieving haplogroup information from MitoMap')

        cols = ['Position', 'Mutated Base','Reference Base']
        for index, row in df.iterrows():
            data = (str(int(row[cols[0]])), str(row[cols[1]]), str(row[cols[2]]))
            custom_query = ''.join(base_query_split[0]+data[2]+base_query_split[1]+
                                   data[1]+base_query_split[2]+data[0]+base_query_split[3]+data[0]
                                   +base_query_split[4] + data[2]+base_query_split[5]+data[1]+base_query_split[6])
            fp = urllib.request.urlopen(custom_query)
            print("Query: ",custom_query)
            mybytes = fp.read()
            my_web_page = mybytes.decode("utf8")
            fp.close()
            html_file = 'haplogroups_'+str(index)+'.html'
            with open(dir_htmls+html_file,'w+') as f:
                f.write(my_web_page)
            haplogroups = [my_web_page[m.end(0):m.end(0)+2] if my_web_page[m.end(0)] == 'H' and my_web_page[m.end(0)+1] == 'V' else my_web_page[m.end(0):m.end(0)+1]
                if my_web_page[m.end(0)] != 'L' else my_web_page[m.end(0):m.end(0)+2] for m in re.finditer('haplogroup=', my_web_page)]
            key_change = str(data[0]+data[1])
            print(key_change)
            for haplogroup in haplogroups:
                set_positions = set([i.strip() for i in df_haplogroups.loc[haplogroup]['Ancestral Marker Motif †'].strip().split(',')])
                if key_change in set_positions:
                    print('change')
                    df.at[index,'Nb.Seqs.correction'] += 1
                df.at[index,haplogroup] += 1
    else:
        print('Assuming information already available')
        df = pd.read_csv('../output/'+suffix+'_haplogroups.csv', sep = ',')
    df.to_csv('../output/lsu_haplogroups.csv')
    return df, haplogrupos

if __name__ == '__main__':
    #pd.set_option('display.max_rows', None)
    mitomap, tsv_read = [], False
    print('python mtDNA.py --mitomap comma_separated_files --somatic comma_separated_files')
    df_control_lsu = pd.read_csv(control_lsu, sep = ';')
    df_control_lsu['Genomic'] = pd.to_numeric(df_control_lsu['Genomic'].fillna(-1))
    print(df_control_lsu)
    df_control_ssu = pd.read_csv(control_ssu, sep=';')
    df_control_ssu['Genomic'] = pd.to_numeric(df_control_ssu['Genomic'].fillna(-1))
    print(df_control_ssu)
    df_haplogrupos = pd.read_csv(haplogrupos_control, sep = '\t')
    df_haplogrupos_indexed = df_haplogrupos.set_index('Top Level Haplogroup')
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
        ssu_df, haplogroups = get_haplogroups(ssu_df, df_haplogroups = df_haplogrupos_indexed, suffix = 'ssu')
        print(ssu_df)
        results_df = add_cv(ssu_df, df_control_ssu, haplogroups)
        results_df.to_csv('../output/ssu_df.csv')
        print('SSU MitoMap information: ', results_df)
        # Process lsu
        preprocess_df(lsu_df)
        # Getting haplogroups
        lsu_df, haplogroups = get_haplogroups(lsu_df, df_haplogroups = df_haplogrupos_indexed)
        results_df = add_cv(lsu_df, df_control_lsu, haplogroups)
        results_df.to_csv('../output/lsu_df.csv')
        print('LSU MitoMap information: ', results_df)