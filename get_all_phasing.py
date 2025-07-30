# https://quinlangroup.slack.com/archives/C02KEHXJ274/p1724332317643529
# /scratch/ucgd/lustre-labs/quinlan/u6018199/cyvcf2
# https://quinlangroup.slack.com/archives/C449KJT3J/p1751389842484399
from cyvcf2 import VCF  # type: ignore

# from tqdm import tqdm # testing
import polars as pl
import bioframe as bf # https://bioframe.readthedocs.io/en/latest/index.html

# NUMBER_VARIANTS = 100000 # testing 

def stringify(phase_block_id): 
    if phase_block_id == -2147483648: # https://github.com/brentp/cyvcf2/issues/31#issuecomment-273214346
        return '.'
    else: 
        return str(phase_block_id)

def is_snv_het(variant, sample_index):
    # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
    # https://brentp.github.io/cyvcf2/#cyvcf2
    gt_types = variant.gt_types
    is_het = gt_types[sample_index] == 1 # 1 indicates heterozygous

    # Even without filtering to hets, all sites in bit-vectors in first 40mb of chr1 are observed to be hets.
    # Therefore constructing bitvectors from hets only is not actually restrictive, AND allows us to track only 2 bitvectors, instead of 4.

    is_snp = variant.is_snp

    # is_snp allows multiallelic sites, so multiplicity of ALT alleles must be tested separately 
    # https://github.com/brentp/cyvcf2/blob/541ab16a255a5287c331843d8180ed6b9ef10e00/cyvcf2/cyvcf2.pyx#L1903-L1911
    has_single_ALT_allele = len(variant.ALT) == 1

    return is_het and is_snp and has_single_ALT_allele

def get_read_phasing(vcf):
    # assume vcf follows this phasing format: hap1 | hap2 

    records = []
    with VCF(vcf, strict_gt=True) as vcf_reader: # cyvcf2 handles .vcf.gz directly
        # # assume multi-sample vcf:
        # samples = vcf_reader.samples
        # sample_index = samples.index(uid) if uid in samples else None

        # assume single-sample vcf: 
        samples = vcf_reader.samples
        assert len(samples) == 1, f"Expected single sample in VCF, found {len(samples)} samples: {samples}"
        sample_index = 0

        for variant in vcf_reader:
        # for variant in tqdm(vcf_reader, total=vcf_reader.num_records): # testing 
        # for i, variant in enumerate(vcf_reader): # testing
        #     if i >= NUMBER_VARIANTS: # testing
        #         break # testing

            if not is_snv_het(variant, sample_index):
                continue
        
            chrom = variant.CHROM
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

            records.append({ 
                "chrom": chrom,
                "start": start,
                "end": end,
                "REF": REF,
                "ALT": ALT,
                "phase_block_id": phase_block_id,
                "allele_hap1": allele_hap1,
                "allele_hap2": allele_hap2,
            })

    df = pl.DataFrame(records)
    return df 

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
    )
    return df

def get_iht_phasing(uid, vcf): 
    # assume vcf is phased as: "paternal | maternal" 
    # https://quinlangroup.slack.com/archives/C08U7NLC9PZ/p1748885496941579

    records = []
    with VCF(vcf, strict_gt=True) as vcf_reader: # cyvcf2 handles .vcf.gz directly
        samples = vcf_reader.samples
        sample_index = samples.index(uid) if uid in samples else None

        # for variant in tqdm(vcf_reader, total=vcf_reader.num_records): # testing 
        for variant in vcf_reader:
        # for i, variant in enumerate(vcf_reader): # testing
        #     if i >= NUMBER_VARIANTS: # testing
        #         break # testing

            if not is_snv_het(variant, sample_index):
                continue

            chrom = variant.CHROM
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

            records.append({ 
                "chrom": chrom,
                "start": start,
                "end": end,
                "REF": REF,
                "ALT": ALT,
                "allele_pat": allele_pat,
                "allele_mat": allele_mat,
            })
 
    df = pl.DataFrame(records)
    return df 

def get_iht_blocks(uid, txt):
    records = []
    with open(txt, 'r') as f:
        header = f.readline().strip().strip('#')
        assert header.startswith("chrom start end"), f"Unexpected header format: {header}"
        assert header.endswith("marker_count len markers"), f"Unexpected header format: {header}"
        samples = header.split()[3:-3]  # skip first 3 columns (chrom, start, end) and last 3 columns (marker_count len markers)
        sample_index = samples.index(uid) if uid in samples else None

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
            
    df = pl.DataFrame(records)
    return df 

def get_all_phasing(df_read_phasing, df_read_phase_blocks, df_iht_phasing, df_iht_blocks):
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
