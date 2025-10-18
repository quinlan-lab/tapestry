import bioframe as bf # https://bioframe.readthedocs.io/en/latest/index.html
import polars as pl 

def get_tiles(reference_genome, tile_size):
    chromsizes = bf.fetch_chromsizes(db=reference_genome)
    windows = []
    for chrom, length in chromsizes.items():
        for start in range(0, length, tile_size): # type: ignore
            end = start + tile_size
            if end > length:
                continue
            windows.append({"chrom": chrom, "start": start, "end": end})

    return pl.DataFrame(windows)

