from tqdm import tqdm 
from pathlib import Path 
import polars as pl 

from tile import get_tiles 
from get_palladium_prefixes import get_prefixes_wrapper
from read_data import read_tapestry
from methylation import compute_methylation
from version_sort import version_sort
from prefix_columns import prefix_columns

def compute_delta_methylation(df): 
    # List of the metric types you want to calculate a delta for
    metrics_to_diff = [
        "count_based_meth",
        "model_based_meth"
    ]
    delta_expressions = []
    for metric in metrics_to_diff:
        expr = (pl.col(f"{metric}_pat") - pl.col(f"{metric}_mat")).alias(f"delta_of_{metric}")
        delta_expressions.append(expr)
    df = df.with_columns(delta_expressions)

    df = df.drop([
        'num_cpgs',
        'count_based_meth',
        'model_based_meth'
    ])

    return df 

def compute_delta_methylation_all_samples(reference_genome, tile_size, meth_read_phased_dir): 
    df_tiles = get_tiles(reference_genome, tile_size)
    prefixes = get_prefixes_wrapper()
    prefixes = prefixes[:2] # TESTING
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

        df_meth_free_from_allele_specific_cpgs = df_meth.filter(~pl.col('cpg_is_allele_specific'))
        df_tiles_with_meth = compute_methylation(df_tiles, df_meth_free_from_allele_specific_cpgs)        
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
                # capture tiles in which at least one of the two dfs has a record: 
                how="full", 
                coalesce=True, 
            )   
    return version_sort(df_all_samples)

