# Borja :)
import sys, re
import pandas as pd
from utils.utils import *
import urllib.request

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
        mutation.append(changes[1])
        reference.append(changes[0])
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
    df_control.groupby('Genomic').first()
    cvtot, cvmit, appearences, positions, genomic = [], [], [], [], []
    for index, row in df_control.iterrows():
        genomic_pos = row['Genomic']
        df_tmp = df_section.loc[df_section['Position'] == genomic_pos]
        for i2, r2 in df_tmp.iterrows():
            cvtot.append(row['CVTOT'])
            cvmit.append(row['CVMIT'])
            appearences.append(max(r2['GB Seqs'], r2['Curated References']))
            positions.append(row['position'])
            genomic.append(genomic_pos)
        if df_tmp.empty:
            cvtot.append(row['CVTOT'])
            cvmit.append(row['CVMIT'])
            appearences.append(0)
            positions.append(row['position'])
            genomic.append(-1)
    return pd.DataFrame({'Position':positions,'CVTOT':cvtot,
                         'CVMIT':cvmit,'GB Seqs':appearences,'Genomic':genomic})

'''
Obtención de haplogrupos para cada una de las secuencias:
    * Construcción de la query de MitoMap.
    * wget via subprocess
    * Parseo del html
    * mkdir htmls_haplogrupos
'''
def get_haplogroups(df, suffix = 'lsu'):
    '''
    :param df: dataframe with information to build the query
    :return: same dataframe with new columns counting the number of times each haplogroup appears
    base query = https://www.mitomap.org/cgi-bin/index_mitomap.cgi?title=Coding+Polymorphism+G-A+at+rCRS+position+2701&pos=2701&ref=G&alt=A&purge_type=
    '''
    base_query = 'https://www.mitomap.org/cgi-bin/index_mitomap.cgi?title=Coding+Polymorphism+\-\+at+rCRS+position+\&pos=\&ref=\&alt=\&purge_type='
    base_query_split = base_query.split('\\')
    dir_htmls = '../mitomap_html_dirs_'+suffix+'/'
    if not Utils.exists(dir_htmls):
        Utils.mkdir(dir_htmls)
    files = Utils.get_files(dir_htmls,extension=['html'])
    haplogroup_column = []
    if not len(files) != 0:
        print('Retrieving haplogroup information from MitoMap')

        cols = ['Position', 'Mutated Base','Reference Base']
        for index, row in df.iterrows():
            data = (str(row[cols[0]]), str(row[cols[1]]), str(row[cols[2]]))
            custom_query = ''.join(base_query_split[0]+data[2]+base_query_split[1]+
                                   data[1]+base_query_split[2]+data[0]+base_query_split[3]+data[0]
                                   +base_query_split[4] + data[2]+base_query_split[5]+data[1]+base_query_split[6])
            fp = urllib.request.urlopen(custom_query)
            mybytes = fp.read()
            my_web_page = mybytes.decode("utf8")
            fp.close()
            html_file = 'haplogroups_'+str(index)+'.html'
            with open(dir_htmls+html_file,'w+') as f:
                f.write(my_web_page)
            haplogroups = [my_web_page[m.end(0):m.end(0)+1] if my_web_page[m.end(0)] != 'L' else my_web_page[m.end(0):m.end(0)+2]
                   for m in re.finditer('haplogroup=', my_web_page)]
            haplogroup_column.append(haplogroups)
    else:
        print('Assuming information already available')
    df['haplogroups'] = haplogroup_column

if __name__ == '__main__':
    #pd.set_option('display.max_rows', None)
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
        # Getting haplogroups
        get_haplogroups(lsu_df)
        results_df = add_cv(lsu_df, df_control_lsu)
        results_df.to_csv('../output/lsu_df.csv')
        print('LSU MitoMap information: ', results_df)