# columns kmer,Observed,Expected,Z_score,E_values,P_values,Frequency

def get_kmer_data_from_csv(filename, column_interest, value_type=int):
    kmer_data = defaultdict(value_type)
    with open(filename, 'r') as fh:
        csvreader = csv.reader(fh)
        col_idx = next(csvreader).index(column_interest)
        for row in csvreader:
            kmer, values = row[0], row[col_idx]
            kmer_data[kmer] = value_type(values)
    return kmer_data
    
    
def get_all_data_from_csvs(filenames, column_interest, value_type=int):
    data_all = defaultdict(dict)
    for filename in filenames:
        name = filename.split('/')[2]
        data_all[name] = get_kmer_data_from_csv(filename, column_interest, value_type=value_type)
    return data_all
    
    
def make_data_correlation(data_dict, out_file_name=None, method='pearson', as_save=False):
    df = pd.DataFrame(data_dict)
    df_cor = df.corr(method=method)
    if as_save:
        csv_name = out_file_name
        df_cor.to_csv(csv_name)
    return df_cor      
    
    
    
    
    
    
    
        
