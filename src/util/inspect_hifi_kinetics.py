import pysam
import os 

BAM_PATH = "/scratch/ucgd/lustre-core/UCGD_Staging/keck_longread/unaligned_bams__ped__notes/r84286_20251101_013201_2_D01/m84286_251101_183156_s2.hifi_reads.bc2057.2008268_2033.bam"

def inspect_tags(bam_path, num_reads=5):
    """
    Prints all tags and their values for the first N reads.
    """
    # check_sq=False is required for uBAMs 
    # as otherwise pysam checks if SQ entries are present in header
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for i, read in enumerate(bam):
            if i >= num_reads:
                break
                
            print(f"--- Read: {read.query_name} ---")
            
            # read.tags returns a list of tuples: [('RG', 'A'), ('HP', 1), ...]
            # It automatically handles the values (int, string, float)
            tags = read.tags 
            
            if not tags:
                print("  (No tags found)")
                continue

            for tag_name, tag_value in tags:
                print(f"  {tag_name}: {tag_value}")

def inspect_hifi_kinetics(bam_path, num_reads=2):
    """
    Scans a HiFi BAM for specific methylation-relevant kinetic tags (fi, fp, ri, rp).
    """
    # The standard kinetic tags required by Jasmine/Primrose
    kinetic_tags = {
        'fi': 'Forward IPD',
        'fp': 'Forward PulseWidth',
        'ri': 'Reverse IPD',
        'rp': 'Reverse PulseWidth'
    }

    print(f"Checking {bam_path} for HiFi kinetics...\n")

    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for i, read in enumerate(bam):
            if i >= num_reads:
                break

            print(f"=== Read: {read.query_name} ===")
            found_any = False

            for tag, desc in kinetic_tags.items():
                if read.has_tag(tag):
                    found_any = True

                    # These return array('B') - unsigned chars (0-255)
                    values = read.get_tag(tag)
                    
                    # Basic stats
                    length = len(values)
                    
                    # Preview the data 
                    preview = list(values[:10])
                    
                    print(f"  [{tag}] {desc}")
                    print(f"    ├─ Length:  {length} (matches seq len: {length == read.query_length})")
                    print(f"    └─ Values:  {preview} ...")
            
            if not found_any:
                print("  [!] No HiFi kinetic tags found.")
                print("      (Ensure CCS was run with --hifi-kinetics)")
            
            print("")

def compute_bam_size(bam_path): 
    if not os.path.exists(bam_path):
        print(f"Error: File not found at {bam_path}")
        return

    size_bytes = os.path.getsize(bam_path)
    size_gb = size_bytes / (1024 ** 3)
    
    print(f"{'='*60}")
    print(f"Bam Size:      {size_gb:.2f} GB")    

# inspect_tags(BAM_PATH)
inspect_hifi_kinetics(BAM_PATH)
compute_bam_size(BAM_PATH)