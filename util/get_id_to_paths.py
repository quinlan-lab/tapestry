import pandas as pd
from pathlib import Path

def get_id_to_paths(array_dir, hifi_dir):
    sample_ids_filename = array_dir / 'K1463 EPIC Array ID_Age v2.betas.csv'
    sample_ids = pd.read_csv(sample_ids_filename)
    sample_ids = sample_ids.astype({
        'Sentrix ID': 'Int64',  # Use 'Int64' for nullable integer type
        'EPIC ID': 'str',
        'LINK ID': 'Int64', 
        'Coriell ID': 'Int64', 
        'LabID': 'Int64',
        'Gender': 'category',
        'Age': 'float',
        'PaID': 'Int64',
        'MaID': 'Int64',
    })
    sample_ids['Coriell ID'] = sample_ids['Coriell ID'].apply(lambda x: f'NA{int(x)}' if pd.notnull(x) else x)

    id_to_paths = {} # map 'LINK ID' to file paths for Illumina and HiFi methylation data

    for key, value in sample_ids.set_index('LINK ID').to_dict('index').items():
        id_to_paths[key] = {}

        if pd.notnull(value['Illumina-methylation-array betas at CpG sites']):
            id_to_paths[key]['Illumina-methylation-array betas at CpG sites'] = Path(value['Illumina-methylation-array betas at CpG sites'].replace('.csv', '.sorted.bed'))

        prefix = key if pd.isnull(value['Coriell ID']) else value['Coriell ID']
        for hap in ['hap1', 'hap2', 'combined']:
            id_to_paths[key][f'HiFi-methylation betas genomewide {hap}'] = hifi_dir / f'{prefix}.GRCh38.{hap}.sorted.bed.gz'

    return id_to_paths

def get_header_path__epic(): 
    return '/scratch/ucgd/lustre-labs/quinlan/data-shared/Neklason_MethylationEPIC_20241008/betas/bed_header.txt'

def get_header_path__hifi(mode): 
    if mode == 'model':
        return '/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/pb-cpg/full_pedigree/bed.header'
    elif mode == 'count':
        return None
    else:
        raise ValueError("Invalid mode. Use 'model' or 'count'.")

def get_uid_to_path__epic_and_hifi():
    repo_dir = '/scratch/ucgd/lustre-labs/quinlan/u6018199/dna-methylation'
    data_dir = '/scratch/ucgd/lustre-labs/quinlan/data-shared'
    array_dir = Path(f'{data_dir}/Neklason_MethylationEPIC_20241008')
    hifi_dir_model_mode = Path(f'{data_dir}/datasets/Palladium/pb-cpg/full_pedigree')
    hifi_dir_count_mode = Path(f'{data_dir}/dna-methylation/palladium-count-mode')

    epic_paths_filename = array_dir / 'K1463 EPIC Array ID_Age v2.betas.csv'
    epic_paths = pd.read_csv(epic_paths_filename)
    epic_paths = epic_paths.astype({
        'Sentrix ID': 'Int64',  # Use 'Int64' for nullable integer type
        'EPIC ID': 'str',
        'LINK ID': 'Int64', 
        'Coriell ID': 'Int64', 
        'LabID': 'Int64',
        'Gender': 'category',
        'Age': 'float',
        'PaID': 'Int64',
        'MaID': 'Int64',
    })
    epic_paths['Coriell ID'] = epic_paths['Coriell ID'].apply(lambda x: f'NA{int(x)}' if pd.notnull(x) else x)

    # 1. incorporate ages at blood draw for all platinum samples, including those not assayed on EPIC
    # 2. exclude 4 samples in generation 1 when investigating the relationship between HiFi-based DNA-methylation and age, 
    #    as these are cell lines, and DNA methylation in a cell line may be quite far from what it was in the individual at the time the cell line was established
    #    c.f., https://quinlangroup.slack.com/archives/C02KEHXJ274/p1744997614832519
    ages_generations_filename = Path(repo_dir) / '2024 05 13 k1463 DNA Ages_DEIDENTIFIED.xlsx'
    ages_generations = pd.read_excel(ages_generations_filename)
    df = pd.merge(epic_paths, ages_generations, how='outer', left_on='LINK ID', right_on='LINK ID')
    df = df[[
        'LINK ID',
        'Coriell ID',
        'Gender',
        'Age at blood draw',
        'Generation?',
        'Illumina-methylation-array betas at CpG sites',
    ]]
    df = df.rename(columns={  # type: ignore
        'Age at blood draw': 'Age',
        'Generation?': 'Generation',
    })
    df = df[~df['Generation'].str.startswith('1st')]

    uid_to_path = {} # map 'LINK ID' to file paths for Illumina and HiFi methylation data

    for key, value in df.set_index('LINK ID').to_dict('index').items():
        key = str(key)

        uid_to_path[key] = {}

        uid_to_path[key]['age'] = value['Age']
        uid_to_path[key]['generation'] = value['Generation']
        uid_to_path[key]['gender'] = value['Gender']
        
        if pd.notnull(value['Illumina-methylation-array betas at CpG sites']):
            uid_to_path[key]['Illumina-methylation-array betas at CpG sites'] = Path(value['Illumina-methylation-array betas at CpG sites'].replace('.csv', '.unique.sorted.bed'))

        prefix = key if pd.isnull(value['Coriell ID']) else value['Coriell ID']
        uid_to_path[key]['prefix'] = prefix
        for hap in ['hap1', 'hap2', 'combined']:
            uid_to_path[key][f'HiFi-methylation {hap} (model mode)'] = hifi_dir_model_mode / f'{prefix}.GRCh38.{hap}.sorted.bed.gz'
            uid_to_path[key][f'HiFi-methylation {hap} (count mode)'] = hifi_dir_count_mode / f'{prefix}.GRCh38.haplotagged.{hap}.bed.gz'

    return uid_to_path

def get_prefixes(uid_to_path):
    return [uid_to_path[uid]['prefix'] for uid in uid_to_path.keys()]

