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
