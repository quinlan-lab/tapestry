# https://quinlangroup.slack.com/archives/C02KEHXJ274/p1724332317643529
# /scratch/ucgd/lustre-labs/quinlan/u6018199/cyvcf2
# https://quinlangroup.slack.com/archives/C449KJT3J/p1751389842484399
from cyvcf2 import VCF  # type: ignore
from contextlib import closing

# from tqdm import tqdm # testing
import polars as pl
import bioframe as bf # https://bioframe.readthedocs.io/en/latest/index.html

from phasing import stringify, is_snv_het

READ_PHASING_SCHEMA = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "REF": pl.String,
    "ALT": pl.String,
    "phase_block_id": pl.String,
    "allele_hap1": pl.String,
    "allele_hap2": pl.String,
}

IHT_PHASING_SCHEMA = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "REF": pl.String,
    "ALT": pl.String,
    "allele_pat": pl.String,
    "allele_mat": pl.String,
}


def get_read_phasing(vcf, chrom):
    """Read phased heterozygous SNVs from one indexed VCF chromosome."""
    columns = {name: [] for name in READ_PHASING_SCHEMA}
    with closing(VCF(vcf, strict_gt=True)) as vcf_reader:
        samples = vcf_reader.samples
        assert len(samples) == 1, f"Expected single sample in VCF, found {len(samples)} samples: {samples}"
        sample_index = 0

        for variant in vcf_reader(chrom):
            if not is_snv_het(variant, sample_index):
                continue

            pos = variant.POS # pos is 1-based
            start = pos - 1 
            end = pos 
            REF = variant.REF
            ALT = variant.ALT[0]  

            phase_block_id_all_samples = variant.format('PS')
            if phase_block_id_all_samples is None:   
                phase_block_id = '.' 
            else:
                phase_block_id = phase_block_id_all_samples[sample_index, 0] 
                phase_block_id = stringify(phase_block_id)

            # single sample of 0|1 in vcf becomes [[0, 1, True]]
            # 2 samples of 0/0 and 1|1 would be [[0, 0, False], [1, 1, True]]
            # https://brentp.github.io/cyvcf2/#cyvcf2
            genotype_all_samples = variant.genotypes
            genotype = genotype_all_samples[sample_index]

            # "-1" in "genotype" indicates missing data: 
            # c.f., "test_set_gts" at: https://github.com/brentp/cyvcf2/issues/31#issuecomment-275195917
            # hap1 | hap2 
            #   https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md#haplotagged-bam-files
            #   https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag-algorithm
            allele_hap1 = str(genotype[0]) if genotype[0] != -1 else '.'
            allele_hap2 = str(genotype[1]) if genotype[1] != -1 else '.'

            phased = genotype[2]
            if not phased: 
                continue 

            columns["chrom"].append(variant.CHROM)
            columns["start"].append(start)
            columns["end"].append(end)
            columns["REF"].append(REF)
            columns["ALT"].append(ALT)
            columns["phase_block_id"].append(phase_block_id)
            columns["allele_hap1"].append(allele_hap1)
            columns["allele_hap2"].append(allele_hap2)

    return pl.DataFrame(columns, schema=READ_PHASING_SCHEMA)

def get_read_phase_blocks(tsv): 
    df = (
        pl
        .read_csv(
            tsv,
            separator="\t",
            has_header=True,
            infer_schema_length=1000000,
            # n_rows=100000  # testing
        )
        .cast({
            "phase_block_id": pl.String,
        })
        # HiPhase reports 1-based inclusive block coordinates. Internal overlap
        # operations use zero-based half-open coordinates.
        .with_columns((pl.col("start") - 1).alias("start"))
    )
    return df

def get_iht_phasing(uid, vcf, chrom):
    """Read inheritance-phased heterozygous SNVs for one chromosome."""
    columns = {name: [] for name in IHT_PHASING_SCHEMA}
    with closing(VCF(vcf, strict_gt=True)) as vcf_reader:
        samples = vcf_reader.samples
        if uid not in samples:
            raise ValueError(f"Sample {uid!r} is absent from inheritance VCF: {samples}")
        sample_index = samples.index(uid)

        for variant in vcf_reader(chrom):
            if not is_snv_het(variant, sample_index):
                continue

            pos = variant.POS # pos is 1-based
            start = pos - 1 
            end = pos
            REF = variant.REF
            ALT = variant.ALT[0]

            # single sample of 0|1 in vcf becomes [[0, 1, True]]
            # 2 samples of 0/0 and 1|1 would be [[0, 0, False], [1, 1, True]]
            # https://brentp.github.io/cyvcf2/#cyvcf2
            genotype_all_samples = variant.genotypes
            genotype = genotype_all_samples[sample_index]

            # "-1" in "genotype" indicates missing data: 
            # c.f., "test_set_gts" at: https://github.com/brentp/cyvcf2/issues/31#issuecomment-275195917
            allele_pat = str(genotype[0]) if genotype[0] != -1 else '.'
            allele_mat = str(genotype[1]) if genotype[1] != -1 else '.'

            phased = genotype[2]

            if not phased: 
                raise ValueError(f"Expected phased genotype, but found unphased: {genotype}")

            columns["chrom"].append(variant.CHROM)
            columns["start"].append(start)
            columns["end"].append(end)
            columns["REF"].append(REF)
            columns["ALT"].append(ALT)
            columns["allele_pat"].append(allele_pat)
            columns["allele_mat"].append(allele_mat)

    return pl.DataFrame(columns, schema=IHT_PHASING_SCHEMA)

def get_iht_blocks(uid, txt):
    records = []
    with open(txt, 'r') as f:
        header = f.readline().strip().strip('#')
        assert header.startswith("chrom start end"), f"Unexpected header format: {header}"
        assert header.endswith("marker_count len markers"), f"Unexpected header format: {header}"
        samples = header.split()[3:-3]  # skip first 3 columns (chrom, start, end) and last 3 columns (marker_count len markers)
        if uid not in samples:
            raise ValueError(f"Sample {uid!r} is absent from IHT header: {samples}")
        sample_index = samples.index(uid)

        for line in f:
            line = line.strip()
            fields = line.split()

            chrom, start, end = fields[:3] 

            genotypes_all_samples = fields[3:-3] # skip first 3 columns (chrom, start, end) and last 3 columns (marker_count len markers)
            genotype = genotypes_all_samples[sample_index] # type: ignore
            if '/' in genotype:
                raise ValueError(f"This sample is a founder and therefore cannot be inheritance-based phased")
            founder_label_pat, founder_label_mat = genotype.split('|')

            if founder_label_pat == '?' and founder_label_mat == '?':
                continue

            records.append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "founder_label_pat": founder_label_pat,
                "founder_label_mat": founder_label_mat,
            })
            
    return pl.DataFrame(
        records,
        schema={
            "chrom": pl.String,
            "start": pl.Int64,
            "end": pl.Int64,
            "founder_label_pat": pl.String,
            "founder_label_mat": pl.String,
        },
    )

def get_all_phasing(df_read_phasing, df_read_phase_blocks, df_iht_phasing, df_iht_blocks):
    if any(
        frame.is_empty()
        for frame in (
            df_read_phasing,
            df_read_phase_blocks,
            df_iht_phasing,
            df_iht_blocks,
        )
    ):
        return pl.DataFrame(
            schema={
                "chrom": pl.String,
                "start": pl.Int64,
                "end": pl.Int64,
                "REF": pl.String,
                "ALT": pl.String,
            }
        )
    df = (
        df_read_phasing
        .join(
            df_read_phase_blocks,
            on=["chrom", "phase_block_id"], 
            how="inner"
        )
        .rename({
            "start_right": "start_phase_block",
            "end_right": "end_phase_block",
            "num_variants": "num_variants_phase_block",
        })
        .join(
            df_iht_phasing, 
            # Join on REF and ALT, in addition to joining on chrom, start, end
            # This is important because the parental-phased vcf results from joint-calling, 
            # whereas the read-backed-phased vcf results from single-sample calling,
            # so REF and ALT may differ between the two data sources 
            on=["chrom", "start", "end", "REF", "ALT"], 
            how="inner"
        )
    )

    df = bf.overlap(
        df.to_pandas(), 
        df_iht_blocks.to_pandas(), 
        how='inner', 
        suffixes=('','_iht_block'),
        cols2=["chrom", "start", "end"]
    ).drop([
        "source_block_index", 
        "sample_name", 
        "phase_block_id", 
        "num_variants_phase_block", 
        "chrom_iht_block"
    ], axis=1)  
    
    df = pl.from_pandas(df)
    df = df.sort(["chrom", "start", "end"])

    return df
