from tqdm import tqdm 
from pathlib import Path 

from tile import get_tiles 
from get_palladium_prefixes import get_prefixes_wrapper
from read_data import read_tapestry
from methylation import compute_methylation, compute_delta_methylation
from version_sort import version_sort

def prefix_columns(df, prefix, join_keys):
    # add prefixes to all columns except the join keys
    return df.rename({col: f"{prefix}_{col}" for col in df.columns if col not in join_keys})

def compute_delta_methylation_all_samples(reference_genome, tile_size, meth_read_phased_dir): 
    df_tiles = get_tiles(reference_genome, tile_size)
    prefixes = get_prefixes_wrapper()
    df_all_samples = None
    for prefix in tqdm(prefixes):
        bed_meth = f"{meth_read_phased_dir}/{prefix}.dna-methylation.founder-phased.all_cpgs.bed"
        if Path(bed_meth).exists():
            df_meth = read_tapestry(bed_meth)
        else: 
            print(f"Could not read CpG sites at which count- and model-based methylation levels have been phased to founder haplotypes") 
            print(f"Required file does not exist: '{bed_meth}'")
            print(f"This may be because this sample is a founder and therefore cannot be inheritance-based phased")
            continue

        # df_meth = df_meth.head(10000) # TESTING
        df_tiles_with_meth = compute_methylation(df_tiles, df_meth)
        df_tiles_with_delta_meth = compute_delta_methylation(df_tiles_with_meth)
        join_keys = ['chrom', 'start', 'end']
        df_tiles_with_delta_meth = (
            df_tiles_with_delta_meth
            .select(join_keys + ['delta_of_count_based_meth', 'delta_of_model_based_meth'])
            .rename({
                'delta_of_count_based_meth': 'count',
                'delta_of_model_based_meth': 'model'
            })
        )
        df_tiles_with_delta_meth = prefix_columns(df_tiles_with_delta_meth, prefix=prefix, join_keys=join_keys)

        if df_all_samples is None:
            df_all_samples = df_tiles_with_delta_meth
        else:
            df_all_samples = df_all_samples.join(
                df_tiles_with_delta_meth,
                on=join_keys,
                # capture tiles in which at least of the two dfs has a record: 
                how="full", 
                coalesce=True, 
            )   
    return version_sort(df_all_samples)

def compute_methylation_all_samples_at_given_loci(df_loci, meth_read_phased_dir): 
    df_loci = df_loci.select(['chrom', 'start', 'end'])

    prefixes = get_prefixes_wrapper()
    # prefixes = prefixes[:2] # TESTING
    df_all_samples = None
    for prefix in tqdm(prefixes):
        bed_meth = f"{meth_read_phased_dir}/{prefix}.dna-methylation.founder-phased.all_cpgs.bed"
        if Path(bed_meth).exists():
            df_meth = read_tapestry(bed_meth)
        else: 
            print(f"Could not read CpG sites at which count- and model-based methylation levels have been phased to founder haplotypes") 
            print(f"Required file does not exist: '{bed_meth}'")
            print(f"This may be because this sample is a founder and therefore cannot be inheritance-based phased")
            continue

        df_loci_with_meth = compute_methylation(df_loci, df_meth)
        df_loci_with_meth = df_loci_with_meth.drop([
            'num_cpgs', 
            'count_based_meth', 
            'model_based_meth', 
            'founder_pat', 
            'founder_mat'
        ])
        join_keys = ['chrom', 'start', 'end']
        df_loci_with_meth = prefix_columns(df_loci_with_meth, prefix=prefix, join_keys=join_keys)

        if df_all_samples is None:
            df_all_samples = df_loci_with_meth
        else:
            df_all_samples = df_all_samples.join(
                df_loci_with_meth,
                on=join_keys,
                # capture loci in which at least one of the two dfs has a record: 
                how="full", 
                coalesce=True, 
            )   

    return version_sort(df_all_samples)
