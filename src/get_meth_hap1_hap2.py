import gzip
from pathlib import Path
import polars as pl
import pysam
from remove_funky_chromosomes import remove_funky_chromosomes
from version_sort import version_sort


METHYLATION_SCHEMA = {
    "chromosome": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "total_read_count": pl.Int64,
    "methylation_level": pl.Float64,
}


class IndexedMethylationBed:
    """Read one chromosome at a time from an indexed pb-CpG-tools BED."""

    def __init__(self, bed: str | Path, pb_cpg_tool_mode: str):
        self.path = Path(bed)
        self.index = Path(f"{self.path}.tbi")
        if not self.path.is_file():
            raise FileNotFoundError(f"Methylation BED does not exist: {self.path}")
        if not self.index.is_file():
            raise FileNotFoundError(f"Methylation BED index does not exist: {self.index}")

        header: list[str] | None = None
        found_pileup_mode = False
        with gzip.open(self.path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##pileup-mode="):
                    observed = line.rstrip("\n").split("=", 1)[1]
                    if observed != pb_cpg_tool_mode:
                        raise ValueError(
                            f"Expected pileup-mode {pb_cpg_tool_mode!r} but found "
                            f"{observed!r} in {self.path}"
                        )
                    found_pileup_mode = True
                elif line.startswith("##"):
                    continue
                elif line.startswith("#"):
                    header = line.rstrip("\n")[1:].split("\t")
                    break
                elif line.strip():
                    break
        if header is None:
            raise ValueError(f"File {self.path} is missing a # header line")
        if not found_pileup_mode:
            raise ValueError(f"File {self.path} is missing ##pileup-mode metadata")
        required = {"chrom", "begin", "end", "mod_score", "cov"}
        missing = sorted(required - set(header))
        if missing:
            raise ValueError(f"File {self.path} is missing columns: {missing}")

        self._column = {name: header.index(name) for name in required}
        self._tabix = pysam.TabixFile(str(self.path), index=str(self.index))
        self._contigs = set(self._tabix.contigs)

    def close(self) -> None:
        self._tabix.close()

    def __enter__(self) -> "IndexedMethylationBed":
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()

    def read_chromosome(self, chrom: str) -> pl.DataFrame:
        columns: dict[str, list] = {
            "chromosome": [],
            "start": [],
            "end": [],
            "total_read_count": [],
            "methylation_level": [],
        }
        if chrom not in self._contigs:
            return pl.DataFrame(columns, schema=METHYLATION_SCHEMA)

        for line in self._tabix.fetch(chrom):
            values = line.split("\t")
            columns["chromosome"].append(values[self._column["chrom"]])
            columns["start"].append(int(values[self._column["begin"]]))
            columns["end"].append(int(values[self._column["end"]]))
            columns["total_read_count"].append(int(values[self._column["cov"]]))
            columns["methylation_level"].append(
                float(values[self._column["mod_score"]]) / 100
            )
        return pl.DataFrame(columns, schema=METHYLATION_SCHEMA)

def read_meth_level(bed: Path, pb_cpg_tool_mode: str) -> pl.DataFrame:
    """
    Reads a methylation level bed file from pb-cpg-tools into a Polars DataFrame.
    
    The function asserts that the file contains header lines starting with '##' 
    and a single header line starting with '#'. It also asserts that the
    pileup-mode in the file matches the expected mode.
    """
    is_gzipped = str(bed).endswith('.gz')
    _open = gzip.open if is_gzipped else open

    has_comment_lines = False
    has_header_line = False
    found_pileup_mode = False

    with _open(bed, 'rt') as f:
        for line in f:
            if line.startswith('##pileup-mode='):
                pileup_mode_in_file = line.strip().split('=')[1]
                assert pileup_mode_in_file == pb_cpg_tool_mode, \
                    f"Expected pileup-mode '{pb_cpg_tool_mode}' but found '{pileup_mode_in_file}' in {bed}"
                found_pileup_mode = True
            
            if line.startswith('##'):
                has_comment_lines = True
            elif line.startswith('#'):
                has_header_line = True
                break
            else:
                # Data line reached before header
                break
    
    assert has_comment_lines, f"File {bed} is missing comment lines starting with '##'"
    assert has_header_line, f"File {bed} is missing a header line starting with '#'"
    assert found_pileup_mode, f"File {bed} is missing '##pileup-mode' comment line."

    df = (
        pl
        .read_csv(
            bed,
            separator='\t',
            comment_prefix='##',
            has_header=True,
            # The header line starts with '#', which polars handles automatically.
        )
        # "The bed file columns will differ between the model and count pileup methods, but both share the first six columns"
        # https://github.com/PacificBiosciences/pb-CpG-tools?tab=readme-ov-file#bed-file-format
        .rename({
            '#chrom': 'chromosome',
            'begin': 'start',
            'end': 'end',
            'mod_score': 'methylation_level_percent',
            'type': 'type',
            'cov': 'total_read_count',
        })
        .select([
            'chromosome', 
            'start', 
            'end', 
            'methylation_level_percent', 
            'total_read_count'
        ]) 
        .with_columns(
            (pl.col('methylation_level_percent') / 100)
            .alias('methylation_level')
        )
        .drop('methylation_level_percent')
    )

    df = remove_funky_chromosomes(df, chrom_column='chromosome')

    return df

def read_meth_hap1_hap2(pb_cpg_tool_mode, bed_hap1, bed_hap2):
    df_meth_hap1 = read_meth_level(
        bed_hap1,
        pb_cpg_tool_mode
    ).select(pl.all().name.suffix("_hap1"))

    df_meth_hap2 = read_meth_level(
        bed_hap2,
        pb_cpg_tool_mode
    ).select(pl.all().name.suffix("_hap2"))

    df_meth = (
        df_meth_hap1        
        .join(
            df_meth_hap2,
            left_on=["chromosome_hap1", "start_hap1", "end_hap1"],
            right_on=["chromosome_hap2", "start_hap2", "end_hap2"],
            # "how=full" captures CpG sites for which methylation level is reported for one haplotype, but not the other, 
            # e.g., because CpG sites are created or destroyed in an individual 
            # https://quinlangroup.slack.com/archives/C0803TM7X0X/p1759354796808349?thread_ts=1759349045.861589&cid=C0803TM7X0X 
            how="full", 
            coalesce=True
        )
        .rename({
            "chromosome_hap1": "chrom", 
            "start_hap1": "start", 
            "end_hap1": "end",
        })
    )
    return version_sort(df_meth)


def read_meth_hap1_hap2_chromosome(
    hap1_reader: IndexedMethylationBed,
    hap2_reader: IndexedMethylationBed,
    chrom: str,
) -> pl.DataFrame:
    """Full-join one chromosome of hap1/hap2 CpGs, preserving unique sites."""
    df_meth_hap1 = hap1_reader.read_chromosome(chrom).select(
        pl.all().name.suffix("_hap1")
    )
    df_meth_hap2 = hap2_reader.read_chromosome(chrom).select(
        pl.all().name.suffix("_hap2")
    )
    return version_sort(
        df_meth_hap1
        .join(
            df_meth_hap2,
            left_on=["chromosome_hap1", "start_hap1", "end_hap1"],
            right_on=["chromosome_hap2", "start_hap2", "end_hap2"],
            how="full",
            coalesce=True,
        )
        .rename(
            {
                "chromosome_hap1": "chrom",
                "start_hap1": "start",
                "end_hap1": "end",
            }
        )
    )
