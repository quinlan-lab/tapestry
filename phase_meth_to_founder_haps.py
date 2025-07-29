import argparse
from pathlib import Path
from get_all_phasing import (
    get_read_based_phasing, 
    get_phase_blocks, 
    get_parental_phasing, 
    get_iht_blocks, 
    get_all_phasing
)
from get_hap_map import (
    get_hap_map,
    write_hap_map_blocks,
    write_bit_vector_sites_and_mismatches,    
)
from get_meth_hap1_hap2 import get_meth_hap1_hap2
from util.write_data import write_bed

# Constants from the notebook
NUMBER_VARIANTS = 100000  # TODO: testing
READ_BACKED_PHASED_DIR = Path('/scratch/ucgd/lustre-labs/quinlan/data-shared/read-backed-phasing')
HAPLOTYPE_MAPS_DIR = Path('/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38')
PB_CPG_TOOL_MODE = 'model'
METH_READ_BACKED_PHASED_DIR = Path(f'/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.{PB_CPG_TOOL_MODE}.read-backed-phased')
METH_FOUNDER_PHASED_DIR = Path(f'/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.{PB_CPG_TOOL_MODE}.founder-phased') 

def get_all_phasing_from_uid(uid):
    """Wrapper function that creates all required DataFrames and calls get_all_phasing"""

    # Get read-based phasing data
    df_read_based_phasing = get_read_based_phasing(
        uid=uid, 
        read_backed_phased_dir=READ_BACKED_PHASED_DIR, 
        number_variants=NUMBER_VARIANTS
    )
    
    # Get read-basedphase blocks
    df_phase_blocks = get_phase_blocks(
        uid=uid, 
        read_backed_phased_dir=READ_BACKED_PHASED_DIR
    )
    
    # Get inheritance-based phasing data
    df_parental_phasing = get_parental_phasing(
        uid=uid, 
        haplotype_maps_dir=HAPLOTYPE_MAPS_DIR, 
        number_variants=NUMBER_VARIANTS
    )
    
    # Get IHT blocks
    df_iht_blocks = get_iht_blocks(
        uid=uid, 
        haplotype_maps_dir=HAPLOTYPE_MAPS_DIR
    )
    
    # Combine all phasing data
    df_all_phasing = get_all_phasing(
        df_read_based_phasing, 
        df_phase_blocks, 
        df_parental_phasing, 
        df_iht_blocks
    )
    
    return df_all_phasing

def main(uid):
    # Create output directory if it doesn't exist
    METH_FOUNDER_PHASED_DIR.mkdir(parents=True, exist_ok=True)
    
    df_all_phasing = get_all_phasing_from_uid(uid)
    df_hap_map, df_sites, df_sites_mismatch = get_hap_map(df_all_phasing)

    write_hap_map_blocks(df_hap_map, uid, "paternal", METH_FOUNDER_PHASED_DIR)
    write_hap_map_blocks(df_hap_map, uid, "maternal", METH_FOUNDER_PHASED_DIR)
    write_bed(METH_FOUNDER_PHASED_DIR, df_hap_map, f"{uid}.hap-map-blocks")
    write_bit_vector_sites_and_mismatches(df_sites, df_sites_mismatch, uid, METH_FOUNDER_PHASED_DIR)

    df_meth_hap1_hap2 = get_meth_hap1_hap2(uid=uid, pb_cpg_tool_mode=PB_CPG_TOOL_MODE, meth_read_backed_phased_dir=METH_READ_BACKED_PHASED_DIR)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Phase DNA methylation data to founder haplotypes')
    parser.add_argument('-uid', required=True, help='Sample UID to process')
    args = parser.parse_args()

    main(args.uid)