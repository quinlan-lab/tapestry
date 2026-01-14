import re
from pathlib import Path
import polars as pl

def print_fraction_cpg_sites(log_directory_path, file_pattern, methylation_logic):
    """
    Parses log files to create a Polars DataFrame of sample IDs and their
    fraction of phased CpG sites.

    Args:
        log_directory_path (str): The path to the directory containing log files.
        file_pattern (str): The pattern to match for log files (e.g., "*.log").
    """
    log_dir = Path(log_directory_path)
    if not log_dir.is_dir():
        print(f"Error: Directory not found at '{log_directory_path}'")
        return

    id_pattern = re.compile(r'/(\w+)\.dna-methylation.founder-phased.all_cpgs.bed')
    escaped_logic = re.escape(methylation_logic)
    percent_pattern = re.compile(rf'{escaped_logic}: (\d+\.\d+)%')

    # 1. Collect data into a list of dictionaries
    data_for_df = []

    for log_file in log_dir.glob(file_pattern):
        found_id = None
        found_percent = None
        with open(log_file, 'r') as f:
            for line in f:
                id_match = id_pattern.search(line)
                if id_match:
                    found_id = id_match.group(1)
                percent_match = percent_pattern.search(line)
                if percent_match:
                    found_percent = float(percent_match.group(1))
        
        if found_id and found_percent is not None:
            data_for_df.append({"sample_id": found_id, "percentage_of_cpg_sites": found_percent})

    # 2. Check if data was found and create the DataFrame
    if not data_for_df:
        print("No data found to create a DataFrame.")
        return

    df = pl.DataFrame(data_for_df)

    # 3. Use Polars expressions to create the final DataFrame and print it
    final_df = (
        df
        .sort("sample_id")
    )

    pl.Config.set_tbl_rows(50)

    print(f"{methylation_logic}:")
    print(final_df)

if __name__ == '__main__':
    print_fraction_cpg_sites(
        log_directory_path='/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.founder-phased.all-cpgs',
        file_pattern="*.log",
        methylation_logic='Percentage of CpG sites (in reference and sample genomes, and on phasable chroms) at which count-based methylation is phased to at least one parental haplotype',
    )

    print_fraction_cpg_sites(
        log_directory_path='/scratch/ucgd/lustre-labs/quinlan/data-shared/dna-methylation/CEPH1463.GRCh38.hifi.founder-phased.all-cpgs',
        file_pattern="*.log",
        methylation_logic='Percentage of CpG sites (in reference and sample genomes, and on phasable chroms) at which count-based unphased methylation is reported',
    )    