"""Regenerate the second deck's four scientific figures from chr22 run records.

This is a fixture analysis, intentionally restricted to the deck's 1 Mb region.
All intervals plotted on the Mb axis are zero-based, half-open.
"""

import argparse
import csv
import json
import logging
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import polars as pl

LEFT, RIGHT = 20499999, 21500000
SIDES = ('paternal', 'maternal')
COLORS = {'A': '#2878b5', 'B': '#66b7e2', 'C': '#af408c', 'D': '#e193bd'}


def setup_axis(ax, title, limits=(LEFT, RIGHT)):
    ax.set_title(title, loc='left')
    ax.set_xlim(limits[0] / 1e6, limits[1] / 1e6)
    ax.set_xlabel('GRCh38 chr22 position (Mb; zero-based)')
    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.spines[['top', 'right']].set_visible(False)


def save(fig, path):
    fig.savefig(path, dpi=160, bbox_inches='tight')
    plt.close(fig)
    logging.info('Wrote %s', path)


def map_lanes(ax, maps):
    for y, side in enumerate(SIDES):
        for row in maps[side].iter_rows(named=True):
            ax.broken_barh([(row['start'] / 1e6, (row['end'] - row['start']) / 1e6)],
                           (y - .25, .5), facecolors=COLORS[row[f'{side}_haplotype']])
            # Only label long intervals; the zoom displays the short left-edge blocks.
            visible_start = max(row['start'] / 1e6, ax.get_xlim()[0])
            visible_end = min(row['end'] / 1e6, ax.get_xlim()[1])
            if visible_end - visible_start > .14 * (ax.get_xlim()[1] - ax.get_xlim()[0]):
                n_snvs = row['num_het_SNVs_in_parent']
                ax.text((visible_start + visible_end) / 2, y,
                        f"{row[f'{side}_haplotype']} · {n_snvs} {'SNV' if n_snvs == 1 else 'SNVs'}\n"
                        f"{row['haplotype_concordance']:.0%} concordance",
                        ha='center', va='center', color='white', fontsize=9)
    ax.set_yticks([0, 1], ['Paternal', 'Maternal'])
    ax.set_ylim(1.5, -.5)


def shade_unresolved(ax, mapping):
    cursor = LEFT
    for row in mapping.sort('start').iter_rows(named=True):
        if cursor < row['start']:
            ax.axvspan(cursor / 1e6, row['start'] / 1e6, color='#dfe2e7', zorder=0)
        cursor = row['end']
    if cursor < RIGHT:
        ax.axvspan(cursor / 1e6, RIGHT / 1e6, color='#dfe2e7', zorder=0)


def run(run_dir, out):
    logging.info('Reading chr22 model experiment: %s', run_dir)
    map_dir = run_dir / 'per_chrom_hap_maps/chr22'
    maps = {side: pl.read_csv(map_dir / f'{side}.hap-map.tsv', separator='\t') for side in SIDES}
    calls = pl.read_csv(run_dir / 'per_chrom_phased/chr22.crossovers.tsv', separator='\t')
    frame = pl.read_csv(run_dir / 'per_chrom_methylation/chr22/trio.dna-methylation.all_cpgs.sorted.bed.gz',
                        separator='\t', comment_prefix='##', infer_schema_length=None)
    if set(frame['#chrom']) != {'chr22'} or not frame['start_cpg'].is_sorted():
        raise ValueError('Expected sorted chr22 CpGs')
    if frame.is_empty() or frame.filter(~pl.col('start_cpg').is_between(LEFT, RIGHT - 1)).height:
        raise ValueError('CpGs outside the slide-deck fixture')
    intervals = []
    for side, mapping in maps.items():
        mismatch_path = map_dir / f'NA12879.bit-vector-sites-mismatches.{side}.bed'
        mismatch_positions = []
        with mismatch_path.open() as handle:
            for line in handle:
                if line.strip() and not line.startswith('#'):
                    chrom, start, *_ = line.split()
                    if chrom == 'chr22':
                        mismatch_positions.append(int(start))
        previous_end = LEFT
        for row in mapping.sort('start').iter_rows(named=True):
            if row['chrom'] != 'chr22' or not previous_end <= row['start'] < row['end'] <= RIGHT:
                raise ValueError(f'Invalid or overlapping {side} map interval: {row}')
            previous_end = row['end']
            intervals.append(dict(side=side, chrom=row['chrom'], start=row['start'], end=row['end'],
                                  haplotype=row[f'{side}_haplotype'],
                                  concordance=row['haplotype_concordance'],
                                  informative_snvs=row['num_het_SNVs_in_parent'],
                                  mismatches=sum(row['start'] <= pos < row['end'] for pos in mismatch_positions)))
        labels = ('kid_pat', 'dad_A', 'dad_B') if side == 'paternal' else ('kid_mat', 'mom_C', 'mom_D')
        unresolved = frame.filter(pl.col(f'{side}_haplotype').is_null())
        if any(unresolved[f'methylation_level_{label}_model'].drop_nulls().len() for label in labels):
            raise ValueError(f'Measurements assigned in unresolved {side} regions')
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'resolved-intervals.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(intervals[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(intervals)
    labels = ('kid_pat', 'kid_mat', 'dad_A', 'dad_B', 'mom_C', 'mom_D', 'kid', 'dad', 'mom')
    summary = {'rows': frame.height, 'measurements': {
        label: frame[f'methylation_level_{label}_model'].drop_nulls().len() for label in labels},
        'crossovers': {'total': calls.height, 'clustered': calls['clustered'].sum()},
        'coordinate_system': 'zero-based, half-open',
        'run_dir': str(run_dir.resolve())}
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')

    # Slide 2: phase blocks use each sample's own PS membership, via whatshap stats.
    fig, ax = plt.subplots(figsize=(12, 3.5), layout='constrained')
    setup_axis(ax, 'WhatsHap phase blocks · chr22 trio fixture')
    for y, sample in enumerate(('NA12879', 'NA12877', 'NA12878')):
        blocks = pl.read_csv(run_dir / f'phase-blocks/{sample}.blocks.tsv', separator='\t')
        for row in blocks.iter_rows(named=True):
            ax.broken_barh([((row['from'] - 1) / 1e6, (row['to'] - row['from'] + 1) / 1e6)],
                           (y - .25, .5), facecolors='#2878b5')
            ax.text((row['from'] + row['to'] - 1) / 2e6, y,
                    f"PS {row['phase_set']} · {row['variants']} variants", ha='center', va='center', color='white')
    ax.set_yticks(range(3), ['Child NA12879', 'Father NA12877', 'Mother NA12878'])
    ax.set_ylim(2.5, -.5)
    (out / 'phase-blocks').mkdir(exist_ok=True)
    save(fig, out / 'phase-blocks/phase-blocks.png')

    # Slide 3: retain every candidate, including both parental events in one bracket.
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), layout='constrained')
    for ax, limits, title in zip(axes, ((20500000, 20512000), (20950000, 20970000)),
                                 ('Left-edge cluster · review flagged', 'Isolated paternal candidate · read supported')):
        setup_axis(ax, title, limits)
        for call in calls.iter_rows(named=True):
            y = 0 if call['parent'] == 'father' else 1
            left, right = call['position1'] - 1, call['position2']
            ax.broken_barh([(left / 1e6, (right - left) / 1e6)], (y - .22, .44),
                           facecolors='#d96850' if call['clustered'] else '#d8a326')
        ax.set_yticks([0, 1], ['Paternal', 'Maternal'])
        ax.set_ylim(1.5, -.5)
    save(fig, out / 'crossover-brackets.png')

    # Slide 4: show the small edge blocks rather than implying genome-wide support.
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), layout='constrained')
    for ax, limits, title in zip(axes, ((LEFT, RIGHT), (20500000, 20512000)),
                                 ('Resolved parental intervals', 'Left-edge zoom · few informative SNVs')):
        setup_axis(ax, title, limits)
        map_lanes(ax, maps)
    save(fig, out / 'resolved-intervals.png')

    # Slide 5: never replace nulls by zero or interpolate across missing calls.
    fig, axes = plt.subplots(4, 1, figsize=(12, 11), layout='constrained')
    setup_axis(axes[0], 'Parental haplotype assignments')
    map_lanes(axes[0], maps)
    for ax, label, title, color in zip(axes[1:], ('kid_pat', 'kid_mat', 'kid'),
                                     ('Child paternal', 'Child maternal', 'Child combined'),
                                     ('#2878b5', '#af408c', '#37695e')):
        observed = frame.filter(pl.col(f'methylation_level_{label}_model').is_not_null())
        setup_axis(ax, f'{title} · {observed.height:,} measured CpG rows')
        if label != 'kid':
            shade_unresolved(ax, maps['paternal' if label == 'kid_pat' else 'maternal'])
        ax.scatter(observed['start_cpg'].to_numpy() / 1e6,
                   observed[f'methylation_level_{label}_model'].to_numpy() * 100,
                   s=2, alpha=.35, color=color, rasterized=True)
        ax.set_ylim(-3, 103)
        ax.set_ylabel('Methylation (%)')
    fig.suptitle(f"Chr22 trio fixture · {frame.height:,} all-CpG rows · model methylation\n"
                 'Gray: unresolved assignment; no point: no reported call; combined includes untagged reads')
    save(fig, out / 'parent_phased_methylation.png')
    logging.info('Completed summary: %s', summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, default=Path('trio_dev_data/output/nextflow-trio-model'))
    parser.add_argument('--output-dir', type=Path, help='Write figures and summaries elsewhere; default: run directory')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    run(args.run_dir, args.output_dir or args.run_dir)


if __name__ == '__main__':
    main()
