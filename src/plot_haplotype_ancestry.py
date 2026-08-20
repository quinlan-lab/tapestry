#!/usr/bin/env python3
"""Write an offline, lazy-loaded Plotly bundle of pedigree haplotype blocks."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import tempfile
from pathlib import Path
from typing import Any

from plotly.offline import get_plotlyjs

from plot_pedigree_haplotypes import layout as pedigree_layout
from plot_transmission_qc import TransmissionPlotError, write_qc_page


class HaplotypePlotError(ValueError):
    """Raised when a PED or inheritance map cannot be visualized safely."""


MISSING_PARENTS = {"0", "NA", "N/A", ".", "-"}
FIXED_IHT_COLUMNS = {"#chrom", "chrom", "start", "end", "marker_count", "len", "markers"}


def _parent(value: str) -> str | None:
    return None if value.upper() in MISSING_PARENTS else value


def read_pedigree(path: Path) -> dict[str, dict[str, str | None]]:
    """Read a validated six-column PED into a sample-indexed mapping."""
    people: dict[str, dict[str, str | None]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) != 6:
                raise HaplotypePlotError(
                    f"{path} line {line_number}: expected six PED columns"
                )
            family, sample, father, mother, sex, phenotype = fields
            if sample in people:
                raise HaplotypePlotError(f"{path}: duplicate sample {sample!r}")
            people[sample] = {
                "family": family,
                "father": _parent(father),
                "mother": _parent(mother),
                "sex": sex,
                "phenotype": phenotype,
            }
    if not people:
        raise HaplotypePlotError(f"{path}: no pedigree records")
    return people


def _chromosome_key(chromosome: str) -> tuple[int, int | str]:
    match = re.fullmatch(r"chr(\d+)", chromosome, flags=re.IGNORECASE)
    if match:
        return (0, int(match.group(1)))
    return (1, chromosome)


def _split_haplotypes(cell: str) -> tuple[list[str], bool]:
    if "|" in cell:
        parts = cell.split("|")
        phased = True
    elif "/" in cell:
        parts = cell.split("/")
        phased = False
    else:
        return [], False
    if len(parts) != 2:
        return [], phased
    return parts, phased


def read_inheritance_blocks(
    path: Path, people: dict[str, dict[str, str | None]]
) -> tuple[list[dict[str, Any]], list[str], list[str]]:
    """Read and merge adjacent equal-label runs from a gtg ``.iht.txt`` file."""
    header: list[str] | None = None
    rows: list[list[str]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split()
            if header is None:
                header = fields
                continue
            if len(fields) != len(header):
                raise HaplotypePlotError(
                    f"{path} line {line_number}: expected {len(header)} fields, "
                    f"found {len(fields)}"
                )
            rows.append(fields)
    if header is None or not rows:
        raise HaplotypePlotError(f"{path}: inheritance map has no data rows")

    names = {name: index for index, name in enumerate(header)}
    chrom_name = "#chrom" if "#chrom" in names else "chrom"
    for required in (chrom_name, "start", "end"):
        if required not in names:
            raise HaplotypePlotError(f"{path}: missing required column {required!r}")
    samples = [
        name
        for name in header
        if name not in FIXED_IHT_COLUMNS and name in people
    ]
    if not samples:
        raise HaplotypePlotError(f"{path}: no PED samples occur in the inheritance map")

    blocks: list[dict[str, Any]] = []
    last_by_track: dict[tuple[str, str, int], dict[str, Any]] = {}
    chromosomes: set[str] = set()
    for fields in rows:
        chromosome = fields[names[chrom_name]]
        start = int(fields[names["start"]])
        end = int(fields[names["end"]])
        if end < start:
            raise HaplotypePlotError(
                f"{path}: invalid interval {chromosome}:{start}-{end}"
            )
        # gtg represents an inheritance run supported by a single marker as a
        # point (start == end, len == 0). Preserve it in the visualization as
        # the smallest BED-like interval instead of rejecting or dropping it.
        if end == start:
            end = start + 1
        chromosomes.add(chromosome)
        for sample in samples:
            haplotypes, phased = _split_haplotypes(fields[names[sample]])
            for hap_index, label in enumerate(haplotypes, 1):
                key = (chromosome, sample, hap_index)
                previous = last_by_track.get(key)
                if (
                    previous is not None
                    and previous["end"] == start
                    and previous["label"] == label
                    and previous["phased"] == phased
                ):
                    previous["end"] = end
                    continue
                block: dict[str, Any] = {
                    "chrom": chromosome,
                    "start": start,
                    "end": end,
                    "sample": sample,
                    "hap": hap_index,
                    "label": label,
                    "phased": phased,
                }
                blocks.append(block)
                last_by_track[key] = block
    chromosomes_sorted = sorted(chromosomes, key=_chromosome_key)
    return blocks, chromosomes_sorted, samples


def read_selected_samples(path: Path | None) -> list[str]:
    if path is None:
        return []
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
    if not lines:
        return []
    return [line.split("\t", 1)[0] for line in lines[1:]]


def _generation(sample: str, people: dict[str, dict[str, str | None]]) -> int:
    memo: dict[str, int] = {}

    def visit(person: str, active: set[str]) -> int:
        if person in memo:
            return memo[person]
        if person in active:
            raise HaplotypePlotError(f"pedigree ancestry cycle at {person!r}")
        active.add(person)
        parents = [
            parent
            for parent in (people[person]["father"], people[person]["mother"])
            if parent in people
        ]
        value = 0 if not parents else max(visit(parent, active) for parent in parents) + 1
        active.remove(person)
        memo[person] = value
        return value

    return visit(sample, set())


def _founder_label_names(
    blocks: list[dict[str, Any]], people: dict[str, dict[str, str | None]]
) -> dict[str, str]:
    """Map gtg's compact founder codes to PED sample haplotype names."""
    labels_by_track: dict[tuple[str, int], set[str]] = {}
    for block in blocks:
        key = (str(block["sample"]), int(block["hap"]))
        label = str(block["label"])
        if label != "?":
            labels_by_track.setdefault(key, set()).add(label)

    names: dict[str, str] = {}
    for sample, person in people.items():
        if person["father"] is not None or person["mother"] is not None:
            continue
        for hap in (1, 2):
            labels = labels_by_track.get((sample, hap), set())
            if len(labels) != 1:
                continue
            label = next(iter(labels))
            display = f"{sample} hap{hap}"
            previous = names.get(label)
            if previous is not None and previous != display:
                raise HaplotypePlotError(
                    f"founder label {label!r} maps to both {previous!r} and {display!r}"
                )
            names[label] = display
    return names


def build_payload(
    ped_path: Path, iht_path: Path, selected_samples_path: Path | None = None
) -> dict[str, Any]:
    people = read_pedigree(ped_path)
    blocks, chromosomes, samples = read_inheritance_blocks(iht_path, people)
    # Retain the explicit ancestry-cycle check for standalone visualizer use.
    for sample in samples:
        _generation(sample, people)
    plotted_people: dict[str, dict[str, Any]] = {
        sample: {**people[sample]}
        for sample in samples
    }
    layout_people = {
        sample: {
            "sex": plotted_people[sample]["sex"],
            "father": (
                plotted_people[sample]["father"]
                if plotted_people[sample]["father"] in plotted_people
                else None
            ),
            "mother": (
                plotted_people[sample]["mother"]
                if plotted_people[sample]["mother"] in plotted_people
                else None
            ),
        }
        for sample in samples
    }
    positions, generations = pedigree_layout(layout_people)
    for sample in samples:
        plotted_people[sample]["generation"] = generations[sample]
    requested = [
        sample
        for sample in read_selected_samples(selected_samples_path)
        if sample in plotted_people
    ]
    fallback = max(samples, key=lambda sample: (plotted_people[sample]["generation"], sample))
    labels = sorted({str(block["label"]) for block in blocks if block["label"] != "?"})
    return {
        "people": plotted_people,
        "samples": samples,
        "selectedSamples": requested,
        "initialSample": requested[0] if requested else fallback,
        "chromosomes": chromosomes,
        "blocks": blocks,
        "labels": labels,
        "labelNames": _founder_label_names(blocks, people),
        "positions": {
            sample: {"x": position[0], "y": position[1]}
            for sample, position in positions.items()
        },
        "sources": {"ped": ped_path.name, "inheritanceMap": iht_path.name},
    }


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Tapestry haplotype ancestry</title>
<script src="plotly.min.js"></script>
<script src="metadata.js"></script>
<style>
:root { color-scheme: light; font-family: Inter, ui-sans-serif, system-ui, sans-serif; }
body { margin: 0; color: #172033; background: #f5f7fb; }
.page { max-width: 1700px; margin: 0 auto; padding: 24px; }
h1 { margin: 0 0 6px; font-size: 25px; }
h2 { margin: 22px 0 6px; font-size: 18px; }
.lede,.subtle { color: #5d687b; }
.lede { margin: 0 0 18px; }
.subtle { margin: 0 0 10px; font-size: 13px; }
.controls { display: flex; flex-wrap: wrap; gap: 14px; align-items: end; margin-bottom: 14px; }
label { display: grid; gap: 5px; color: #445067; font-size: 13px; font-weight: 650; }
select { min-width: 210px; padding: 7px 30px 7px 9px; border: 1px solid #cbd2df; border-radius: 6px; background: white; }
.card { background: white; border: 1px solid #dde2eb; border-radius: 10px; box-shadow: 0 1px 4px rgb(30 45 70 / 8%); }
#overview { min-height: 300px; }
#pedigree { min-height: 520px; }
.detail { margin: 12px 0 0; padding: 12px 16px; line-height: 1.45; }
.detail strong { color: #111827; }
.badge { display: inline-block; padding: 2px 8px; border-radius: 999px; background: #edf1f8; font-size: 12px; }
.loading { color: #536174; }
.nav { display: flex; gap: 8px; margin: 0 0 18px; }
.nav a { color: #334e91; background: white; border: 1px solid #ccd5e5; border-radius: 6px; padding: 6px 10px; text-decoration: none; font-size: 13px; }
</style>
</head>
<body>
<main class="page">
  <h1>Haplotype ancestry across the pedigree</h1>
  <p class="lede">Start with one person's whole-genome chromosome painting, then click a block to paint that chromosome and inherited founder haplotype across the family.</p>
  <nav class="nav"><a href="index.html" aria-current="page">Haplotype ancestry</a><a href="transmission-qc.html">Transmission QC</a></nav>
  <div class="controls">
    <label>Overview sample<select id="sample"></select></label>
    <span id="source" class="badge"></span>
    <span id="loading" class="badge loading">Loading data…</span>
  </div>
  <h2>Whole-genome overview</h2>
  <p class="subtle">Each chromosome has two stripes: hap1 above hap2 (paternal above maternal when phased). Shorter chromosomes end earlier on the shared base-pair scale. Blank spans have no inheritance-map block; they do not indicate sample-specific missing sequencing. Click any block for the pedigree view.</p>
  <div id="overview" class="card"></div>
  <div id="detail" class="card detail"></div>
  <h2>Selected chromosome painted onto the pedigree</h2>
  <p id="pedigree-note" class="subtle">The clicked founder haplotype remains saturated and outlined wherever it occurs; other inherited haplotypes are muted.</p>
  <div id="pedigree" class="card"></div>
</main>
<script>
const DATA = window.TAPESTRY_METADATA;
const overview = document.getElementById('overview');
const pedigree = document.getElementById('pedigree');
const detail = document.getElementById('detail');
const sampleSelect = document.getElementById('sample');
const loading = document.getElementById('loading');
const palette = ['#3f2f84','#4185e1','#2cc4bb','#67e65c','#c7eb22','#f5a02c','#e8480f','#aa1b08','#8057b8','#d2609e','#798230','#9b6b43','#596a93','#43a6a1','#b14d63','#6f9951','#8a70c2','#cf9a2c','#5574a6','#bb5d9b'];
const colors = Object.fromEntries(DATA.labels.map((label, index) => [label, palette[index % palette.length]]));
colors['?'] = '#c7ccd5';
let selectedBlock = null;
let activeSampleBlocks = [];
let activeChromosomeBlocks = [];
let activeMethylationRows = [];
const sampleCache = new Map();
const chromosomeCache = new Map();
const methylationCache = new Map();
const CACHE_LIMIT = 3;
const query = new URLSearchParams(location.search);
const requestedSample = query.get('sample');
const requestedChromosome = query.get('chrom');

function naturalCompare(a, b) { return a.localeCompare(b, undefined, {numeric: true, sensitivity: 'base'}); }
function labelName(label) { return DATA.labelNames[label] || `founder ${label}`; }
function sameBlock(a, b) { return a && b && a.chrom === b.chrom && a.sample === b.sample && a.hap === b.hap && a.start === b.start && a.end === b.end && a.label === b.label; }
function cacheGet(cache, key) {
  if (!cache.has(key)) return null;
  const value = cache.get(key);
  cache.delete(key); cache.set(key, value);
  return value;
}
function cachePut(cache, key, value) {
  cache.set(key, value);
  while (cache.size > CACHE_LIMIT) cache.delete(cache.keys().next().value);
}
function loadShard(kind, key) {
  const isSample = kind === 'sample';
  const cache = isSample ? sampleCache : chromosomeCache;
  const cached = cacheGet(cache, key);
  if (cached) return Promise.resolve(cached);
  const files = isSample ? DATA.sampleFiles : DATA.chromosomeFiles;
  const globalName = isSample ? 'TAPESTRY_SAMPLE_SHARDS' : 'TAPESTRY_CHROMOSOME_SHARDS';
  loading.textContent = `Loading ${isSample ? 'sample' : 'chromosome'} ${key}…`;
  return new Promise((resolve, reject) => {
    const script = document.createElement('script');
    script.src = files[key]; script.async = true;
    script.onload = () => {
      const rows = window[globalName]?.[key];
      if (!rows) { reject(new Error(`Data shard did not define ${key}`)); return; }
      const blocks = rows.map(row => isSample ? {
        chrom: DATA.chromosomes[row[0]], start: row[1], end: row[2], sample: key,
        hap: row[3], label: DATA.blockLabels[row[4]], phased: Boolean(row[5])
      } : {
        chrom: key, start: row[1], end: row[2], sample: DATA.samples[row[0]],
        hap: row[3], label: DATA.blockLabels[row[4]], phased: Boolean(row[5])
      });
      delete window[globalName][key]; script.remove();
      cachePut(cache, key, blocks); loading.textContent = 'Ready';
      resolve(blocks);
    };
    script.onerror = () => {
      script.remove(); loading.textContent = 'Load failed';
      reject(new Error(`Could not load ${files[key]}`));
    };
    document.head.appendChild(script);
  });
}
function loadMethylationShard(chromosome) {
  const cached = cacheGet(methylationCache, chromosome);
  if (cached) return Promise.resolve(cached);
  const file = DATA.methylationChromosomeFiles[chromosome];
  if (!file) return Promise.resolve([]);
  loading.textContent = `Loading ${chromosome} methylation…`;
  return new Promise((resolve, reject) => {
    const script = document.createElement('script');
    script.src = file; script.async = true;
    script.onload = () => {
      const rows = window.TAPESTRY_METHYLATION_SHARDS?.[chromosome];
      if (!rows) { reject(new Error(`Methylation shard did not define ${chromosome}`)); return; }
      const decoded = rows.map(row => ({
        sample: DATA.samples[row[0]], start: row[1], end: row[2],
        patLabel: DATA.blockLabels[row[3]], matLabel: DATA.blockLabels[row[4]],
        patSum: row[5], patN: row[6], matSum: row[7], matN: row[8],
        unphasedSum: row[9], unphasedN: row[10], deltaSum: row[11], deltaN: row[12],
        mismatchN: row[13], alleleSpecificN: row[14], cpgN: row[15]
      }));
      delete window.TAPESTRY_METHYLATION_SHARDS[chromosome]; script.remove();
      cachePut(methylationCache, chromosome, decoded); loading.textContent = 'Ready';
      resolve(decoded);
    };
    script.onerror = () => {
      script.remove(); loading.textContent = 'Load failed';
      reject(new Error(`Could not load ${file}`));
    };
    document.head.appendChild(script);
  });
}
function chooseInitialBlock() {
  const resolved = activeSampleBlocks.filter(block => block.label !== '?');
  const requested = requestedChromosome
    ? resolved.find(block => block.chrom === requestedChromosome)
      || activeSampleBlocks.find(block => block.chrom === requestedChromosome)
    : null;
  selectedBlock = requested || resolved[0] || activeSampleBlocks[0] || null;
}

function renderOverview() {
  const sample = sampleSelect.value;
  const blocks = activeSampleBlocks;
  const traces = [];
  const chromosomeY = Object.fromEntries(DATA.chromosomes.map((chromosome, index) => [chromosome, DATA.chromosomes.length - index]));
  for (const label of [...DATA.labels, '?']) {
    const labelled = blocks.filter(block => block.label === label);
    if (!labelled.length) continue;
    traces.push({
      type: 'bar', orientation: 'h', name: label === '?' ? 'unresolved' : labelName(label),
      x: labelled.map(block => block.end - block.start), base: labelled.map(block => block.start),
      y: labelled.map(block => chromosomeY[block.chrom] + (block.hap === 1 ? 0.16 : -0.16)), width: 0.29,
      marker: {color: colors[label], line: {color: 'rgba(255,255,255,.55)', width: 0.5}}, opacity: 0.94,
      customdata: labelled.map(block => ({...block, labelName: labelName(block.label)})),
      hovertemplate: '<b>%{customdata.sample} · %{customdata.chrom} · hap%{customdata.hap}</b><br>Founder haplotype: %{customdata.labelName}<br>%{customdata.start:,}–%{customdata.end:,} bp<extra></extra>'
    });
  }
  if (selectedBlock && selectedBlock.sample === sample) {
    traces.push({
      type: 'bar', orientation: 'h', x: [selectedBlock.end - selectedBlock.start], base: [selectedBlock.start],
      y: [chromosomeY[selectedBlock.chrom] + (selectedBlock.hap === 1 ? 0.16 : -0.16)], width: 0.38,
      marker: {color: colors[selectedBlock.label], line: {color: '#111827', width: 3}},
      hoverinfo: 'skip', showlegend: false
    });
  }
  const maximumEnd = Math.max(...blocks.map(block => block.end), 1);
  const height = Math.max(300, DATA.chromosomes.length * 29 + 150);
  overview.style.height = `${height}px`;
  Plotly.newPlot(overview, traces, {
    barmode: 'overlay', template: 'plotly_white', height,
    margin: {l: 72, r: 30, t: 52, b: 105},
    title: {text: `${sample} · all configured chromosomes`, x: 0.01, xanchor: 'left'},
    xaxis: {title: 'Genomic position (bp)', range: [0, maximumEnd], tickformat: '~s', gridcolor: '#edf0f5'},
    yaxis: {tickmode: 'array', tickvals: DATA.chromosomes.map(chromosome => chromosomeY[chromosome]), ticktext: DATA.chromosomes, range: [0.4, DATA.chromosomes.length + 0.6], fixedrange: true},
    legend: {orientation: 'h', y: -0.13, x: 0, title: {text: 'Founder haplotype'}},
    hoverlabel: {align: 'left'}, uirevision: `overview|${sample}`
  }, {responsive: true, displaylogo: false, scrollZoom: true}).then(() => {
    overview.on('plotly_click', async event => {
      const block = event.points[0]?.customdata;
      if (!block) return;
      selectedBlock = block;
      renderOverview();
      [activeChromosomeBlocks, activeMethylationRows] = await Promise.all([
        loadShard('chromosome', block.chrom), loadMethylationShard(block.chrom)
      ]);
      renderPedigree();
    });
  });
}

function pedigreeConnectors() {
  const x = [], y = [];
  const couples = new Map();
  for (const [child, person] of Object.entries(DATA.people)) {
    if (!person.father || !person.mother || !DATA.positions[person.father] || !DATA.positions[person.mother]) continue;
    const key = `${person.father}|${person.mother}`;
    if (!couples.has(key)) couples.set(key, {father: person.father, mother: person.mother, children: []});
    couples.get(key).children.push(child);
  }
  for (const family of couples.values()) {
    const father = DATA.positions[family.father], mother = DATA.positions[family.mother];
    const middleX = (father.x + mother.x) / 2, middleY = (father.y + mother.y) / 2;
    const childY = Math.max(...family.children.map(child => DATA.positions[child].y));
    const sibY = (middleY + childY) / 2;
    x.push(father.x, mother.x, null, middleX, middleX, null);
    y.push(father.y, mother.y, null, middleY, sibY, null);
    const childXs = family.children.map(child => DATA.positions[child].x);
    x.push(Math.min(...childXs, middleX), Math.max(...childXs, middleX), null);
    y.push(sibY, sibY, null);
    for (const child of family.children) {
      x.push(DATA.positions[child].x, DATA.positions[child].x, null);
      y.push(sibY, DATA.positions[child].y, null);
    }
  }
  return {type: 'scatter', mode: 'lines', x, y, line: {color: '#7c8491', width: 1.2}, hoverinfo: 'skip', showlegend: false};
}

function overlapsLocus(item, start, end) { return item.start < end && item.end > start; }
function methylationAtLocus(sample, label, start, end) {
  const result = {
    carrier: false, selectedSum: 0, selectedN: 0, otherSum: 0, otherN: 0,
    deltaSum: 0, deltaN: 0, mismatchN: 0, alleleSpecificN: 0, cpgN: 0
  };
  for (const row of activeMethylationRows) {
    if (row.sample !== sample || !overlapsLocus(row, start, end)) continue;
    if (row.patLabel === label) {
      result.carrier = true;
      result.selectedSum += row.patSum; result.selectedN += row.patN;
      result.otherSum += row.matSum; result.otherN += row.matN;
      result.deltaSum += row.deltaSum; result.deltaN += row.deltaN;
    } else if (row.matLabel === label) {
      result.carrier = true;
      result.selectedSum += row.matSum; result.selectedN += row.matN;
      result.otherSum += row.patSum; result.otherN += row.patN;
      result.deltaSum -= row.deltaSum; result.deltaN += row.deltaN;
    } else {
      continue;
    }
    result.mismatchN += row.mismatchN;
    result.alleleSpecificN += row.alleleSpecificN;
    result.cpgN += row.cpgN;
  }
  result.selectedMean = result.selectedN ? result.selectedSum / result.selectedN : null;
  result.otherMean = result.otherN ? result.otherSum / result.otherN : null;
  result.deltaMean = result.deltaN ? result.deltaSum / result.deltaN : null;
  return result;
}

function renderPedigree() {
  if (!selectedBlock) return;
  const chromosome = selectedBlock.chrom;
  const selectedLabel = selectedBlock.label;
  const chromosomeBlocks = activeChromosomeBlocks;
  const chromosomeStart = Math.min(...chromosomeBlocks.map(block => block.start));
  const chromosomeEnd = Math.max(...chromosomeBlocks.map(block => block.end));
  const chromosomeSpan = Math.max(chromosomeEnd - chromosomeStart, 1);
  const paintWidth = 1.9, trackHeight = 0.24, paintOffset = 0.9;
  const shapes = [];
  for (const block of chromosomeBlocks) {
    const position = DATA.positions[block.sample];
    if (!position) continue;
    const left = position.x - paintWidth / 2 + (block.start - chromosomeStart) / chromosomeSpan * paintWidth;
    const right = position.x - paintWidth / 2 + (block.end - chromosomeStart) / chromosomeSpan * paintWidth;
    const top = position.y - paintOffset - (block.hap - 1) * trackHeight;
    const selected = sameBlock(block, selectedBlock);
    const matching = block.label === selectedLabel;
    shapes.push({
      type: 'rect', x0: left, x1: right, y0: top - trackHeight, y1: top,
      fillcolor: colors[block.label], opacity: matching ? 1 : 0.17,
      line: {color: selected ? '#111827' : (matching ? colors[block.label] : 'rgba(255,255,255,.15)'), width: selected ? 3 : (matching ? 1.1 : 0.2)},
      layer: 'above'
    });
  }
  const people = Object.keys(DATA.positions);
  const nodeTrace = {
    type: 'scatter', mode: 'markers+text',
    x: people.map(person => DATA.positions[person].x), y: people.map(person => DATA.positions[person].y),
    text: people, textposition: 'middle center', textfont: {size: 9, color: '#111827'},
    marker: {size: 52, color: 'white', line: {color: '#20242b', width: 1.5}, symbol: people.map(person => DATA.people[person].sex === '1' ? 'square' : 'circle')},
    customdata: people, hovertemplate: '<b>%{customdata}</b><extra></extra>', showlegend: false
  };
  const carrierSet = new Set(chromosomeBlocks.filter(block =>
    block.label === selectedLabel && overlapsLocus(block, selectedBlock.start, selectedBlock.end)
  ).map(block => block.sample));
  const carriers = [...carrierSet].sort(naturalCompare);
  const glyphLineX = [], glyphLineY = [];
  const selectedX = [], selectedY = [], selectedCustom = [], selectedText = [];
  const otherX = [], otherY = [], otherCustom = [];
  const statusX = [], statusY = [], statusText = [];
  const glyphWidth = 1.55;
  if (DATA.methylationSamples.length) {
    for (const sample of carriers) {
      const position = DATA.positions[sample];
      if (!position) continue;
      const glyphY = position.y - 1.7;
      if (!DATA.methylationSamples.includes(sample)) {
        statusX.push(position.x); statusY.push(glyphY); statusText.push('no methylation output');
        continue;
      }
      const summary = methylationAtLocus(sample, selectedLabel, selectedBlock.start, selectedBlock.end);
      if (!summary.carrier || summary.selectedMean === null) {
        statusX.push(position.x); statusY.push(glyphY); statusText.push('no phased CpGs');
        continue;
      }
      const left = position.x - glyphWidth / 2;
      const right = position.x + glyphWidth / 2;
      const selectedPosition = left + summary.selectedMean * glyphWidth;
      glyphLineX.push(left, right, null); glyphLineY.push(glyphY, glyphY, null);
      if (summary.otherMean !== null) {
        const otherPosition = left + summary.otherMean * glyphWidth;
        glyphLineX.push(otherPosition, selectedPosition, null);
        glyphLineY.push(glyphY, glyphY, null);
        otherX.push(otherPosition); otherY.push(glyphY);
        otherCustom.push([sample, summary.otherMean, summary.otherN]);
      }
      selectedX.push(selectedPosition); selectedY.push(glyphY);
      selectedText.push(`${Math.round(summary.selectedMean * 100)}%`);
      selectedCustom.push([
        sample, summary.selectedMean, summary.selectedN, summary.otherMean,
        summary.otherN, summary.deltaMean, summary.deltaN, summary.mismatchN,
        summary.alleleSpecificN, summary.cpgN
      ]);
    }
  }
  const glyphLines = {
    type: 'scatter', mode: 'lines', x: glyphLineX, y: glyphLineY,
    line: {color: '#9aa3b2', width: 1.4}, hoverinfo: 'skip', showlegend: false
  };
  const otherDots = {
    type: 'scatter', mode: 'markers', x: otherX, y: otherY,
    marker: {size: 8, color: 'white', line: {color: '#667085', width: 1.8}},
    customdata: otherCustom,
    hovertemplate: '<b>%{customdata[0]} · other haplotype</b><br>Mean model methylation: %{customdata[1]:.3f}<br>Phased CpGs: %{customdata[2]:,}<extra></extra>',
    showlegend: false
  };
  const selectedDots = {
    type: 'scatter', mode: 'markers+text', x: selectedX, y: selectedY,
    text: selectedText, textposition: 'bottom center', textfont: {size: 8, color: '#344054'},
    marker: {size: 9, color: colors[selectedLabel], line: {color: '#111827', width: 1}},
    customdata: selectedCustom,
    hovertemplate: `<b>%{customdata[0]} · ${labelName(selectedLabel)}</b><br>Mean model methylation: %{customdata[1]:.3f}<br>Phased CpGs: %{customdata[2]:,}<br>Other haplotype: %{customdata[3]:.3f} (%{customdata[4]:,} CpGs)<br>Mean paired difference: %{customdata[5]:+.3f} (%{customdata[6]:,} CpGs)<br>Mismatch-proximal CpGs: %{customdata[7]:,}<br>Allele-specific CpGs: %{customdata[8]:,}<extra></extra>`,
    showlegend: false
  };
  const statusTrace = {
    type: 'scatter', mode: 'text', x: statusX, y: statusY, text: statusText,
    textfont: {size: 8, color: '#7c8491'}, hoverinfo: 'skip', showlegend: false
  };
  const generations = Math.max(...people.map(person => DATA.people[person].generation), 0) + 1;
  const height = Math.max(520, generations * 185 + 150);
  pedigree.style.height = `${height}px`;
  const xs = people.map(person => DATA.positions[person].x);
  const ys = people.map(person => DATA.positions[person].y);
  Plotly.newPlot(pedigree, [pedigreeConnectors(), nodeTrace, glyphLines, otherDots, selectedDots, statusTrace], {
    template: 'plotly_white', height, shapes,
    margin: {l: 35, r: 35, t: 65, b: 35},
    title: {text: `${chromosome} · ${labelName(selectedLabel)}`, x: 0.01, xanchor: 'left'},
    xaxis: {visible: false, range: [Math.min(...xs) - 1.5, Math.max(...xs) + 1.5], fixedrange: true},
    yaxis: {visible: false, range: [Math.min(...ys) - 2.05, Math.max(...ys) + 0.8], fixedrange: true, scaleanchor: 'x', scaleratio: 1},
    hovermode: 'closest', uirevision: `${chromosome}|${selectedLabel}`
  }, {responsive: true, displaylogo: false});
  detail.innerHTML = `<strong>${selectedBlock.sample} · ${chromosome} · ${selectedBlock.phased ? (selectedBlock.hap === 1 ? 'paternal' : 'maternal') : `hap${selectedBlock.hap}`} · ${labelName(selectedLabel)}</strong> &nbsp; ${selectedBlock.start.toLocaleString()}–${selectedBlock.end.toLocaleString()} bp<br><span class="subtle">Carried at this locus by ${carriers.length} sample${carriers.length === 1 ? '' : 's'}: ${carriers.join(', ')}</span>`;
}

for (const sample of [...DATA.samples].sort(naturalCompare)) {
  const option = document.createElement('option'); option.value = sample;
  const generation = `F${DATA.people[sample].generation}`;
  option.textContent = DATA.selectedSamples.includes(sample)
    ? `${generation} ${sample}`
    : `${generation} ${sample} (inheritance map only)`;
  sampleSelect.appendChild(option);
}
sampleSelect.value = requestedSample && DATA.samples.includes(requestedSample)
  ? requestedSample : DATA.initialSample;
document.getElementById('source').textContent = `${DATA.sources.inheritanceMap} + ${DATA.sources.ped}`;
document.getElementById('pedigree-note').textContent = DATA.methylationSamples.length
  ? 'The clicked founder haplotype remains saturated. Beneath each carrier, the solid dot is its mean model methylation and the hollow dot is the other haplotype on a shared 0–1 scale; percentages label the selected founder.'
  : 'The clicked founder haplotype remains saturated and outlined wherever it occurs; this bundle contains no methylation summaries.';
async function selectSample(sample) {
  try {
    activeSampleBlocks = await loadShard('sample', sample);
    chooseInitialBlock();
    if (selectedBlock) {
      [activeChromosomeBlocks, activeMethylationRows] = await Promise.all([
        loadShard('chromosome', selectedBlock.chrom), loadMethylationShard(selectedBlock.chrom)
      ]);
    } else {
      activeChromosomeBlocks = []; activeMethylationRows = [];
    }
    renderOverview(); renderPedigree();
  } catch (error) {
    detail.innerHTML = `<strong>Could not load visualization data.</strong> ${error.message}`;
  }
}
sampleSelect.addEventListener('change', () => selectSample(sampleSelect.value));
selectSample(sampleSelect.value);
</script>
</body>
</html>
"""


def _filename(prefix: str, value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._") or "item"
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()[:10]
    return f"{prefix}-{slug[:48]}-{digest}.js"


def _javascript_assignment(global_name: str, key: str, value: Any) -> str:
    encoded_key = json.dumps(key).replace("</", "<\\/")
    encoded_value = json.dumps(value, separators=(",", ":")).replace("</", "<\\/")
    return (
        f"window.{global_name}=window.{global_name}||{{}};"
        f"window.{global_name}[{encoded_key}]={encoded_value};\n"
    )


def _methylation_rows(
    paths: list[Path],
    sample_indexes: dict[str, int],
    chromosomes: list[str],
    label_indexes: dict[str, int],
) -> tuple[dict[str, list[list[Any]]], list[str]]:
    rows: dict[str, list[list[Any]]] = {chromosome: [] for chromosome in chromosomes}
    summary_samples: list[str] = []
    seen: set[str] = set()
    for path in paths:
        document = json.loads(path.read_text(encoding="utf-8"))
        if document.get("schema_version") != 1:
            raise HaplotypePlotError(f"{path}: unsupported methylation summary schema")
        sample = document.get("sample_id")
        if not isinstance(sample, str) or sample not in sample_indexes:
            raise HaplotypePlotError(f"{path}: invalid methylation summary sample {sample!r}")
        if sample in seen:
            raise HaplotypePlotError(f"duplicate methylation summary for {sample!r}")
        seen.add(sample)
        summary_samples.append(sample)
        chromosome_values = document.get("chromosomes")
        if not isinstance(chromosome_values, dict):
            raise HaplotypePlotError(f"{path}: missing methylation chromosomes")
        for chromosome, values in chromosome_values.items():
            if chromosome not in rows:
                continue
            if not isinstance(values, list):
                raise HaplotypePlotError(f"{path}: invalid rows for {chromosome}")
            for value in values:
                if not isinstance(value, list) or len(value) != 15:
                    raise HaplotypePlotError(
                        f"{path}: invalid methylation row for {chromosome}"
                    )
                pat_label = "?" if value[2] is None else str(value[2])
                mat_label = "?" if value[3] is None else str(value[3])
                if pat_label not in label_indexes or mat_label not in label_indexes:
                    raise HaplotypePlotError(
                        f"{path}: unknown founder label in {chromosome} row"
                    )
                rows[chromosome].append(
                    [
                        sample_indexes[sample],
                        int(value[0]),
                        int(value[1]),
                        label_indexes[pat_label],
                        label_indexes[mat_label],
                        *value[4:],
                    ]
                )
    summary_samples.sort(key=lambda sample: sample_indexes[sample])
    return rows, summary_samples


def write_visualization(
    payload: dict[str, Any],
    output: Path,
    methylation_summaries: list[Path] | None = None,
    transmission_qc_summary: Path | None = None,
) -> None:
    """Write an offline bundle with lazily loaded sample/chromosome shards."""
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise HaplotypePlotError(f"output already exists: {output}")
    samples = list(payload["samples"])
    chromosomes = list(payload["chromosomes"])
    display_labels = list(payload["labels"])
    block_labels = [*display_labels, "?"]
    sample_indexes = {sample: index for index, sample in enumerate(samples)}
    chromosome_indexes = {
        chromosome: index for index, chromosome in enumerate(chromosomes)
    }
    label_indexes = {label: index for index, label in enumerate(block_labels)}
    methylation_rows, methylation_samples = _methylation_rows(
        methylation_summaries or [], sample_indexes, chromosomes, label_indexes
    )
    sample_rows: dict[str, list[list[int]]] = {sample: [] for sample in samples}
    chromosome_rows: dict[str, list[list[int]]] = {
        chromosome: [] for chromosome in chromosomes
    }
    for block in payload["blocks"]:
        common = [
            int(block["start"]),
            int(block["end"]),
            int(block["hap"]),
            label_indexes[str(block["label"])],
            int(bool(block["phased"])),
        ]
        sample_rows[str(block["sample"])].append(
            [chromosome_indexes[str(block["chrom"])], *common]
        )
        chromosome_rows[str(block["chrom"])].append(
            [sample_indexes[str(block["sample"])], *common]
        )

    sample_files = {
        sample: f"data/samples/{_filename('sample', sample)}" for sample in samples
    }
    chromosome_files = {
        chromosome: f"data/chromosomes/{_filename('chromosome', chromosome)}"
        for chromosome in chromosomes
    }
    methylation_files = {
        chromosome: f"data/methylation/chromosomes/{_filename('methylation', chromosome)}"
        for chromosome in chromosomes
        if methylation_rows[chromosome]
    }
    metadata = {
        key: value
        for key, value in payload.items()
        if key not in {"blocks", "labels"}
    }
    metadata.update(
        {
            "labels": display_labels,
            "blockLabels": block_labels,
            "sampleFiles": sample_files,
            "chromosomeFiles": chromosome_files,
            "methylationChromosomeFiles": methylation_files,
            "methylationSamples": methylation_samples,
        }
    )

    with tempfile.TemporaryDirectory(prefix=f".{output.name}.", dir=output.parent) as tmp:
        temporary = Path(tmp)
        (temporary / "data" / "samples").mkdir(parents=True)
        (temporary / "data" / "chromosomes").mkdir(parents=True)
        if methylation_files:
            (temporary / "data" / "methylation" / "chromosomes").mkdir(parents=True)
        (temporary / "index.html").write_text(HTML_TEMPLATE, encoding="utf-8")
        (temporary / "plotly.min.js").write_text(get_plotlyjs(), encoding="utf-8")
        (temporary / "metadata.js").write_text(
            "window.TAPESTRY_METADATA="
            + json.dumps(metadata, separators=(",", ":")).replace("</", "<\\/")
            + ";\n",
            encoding="utf-8",
        )
        transmission_qc = write_qc_page(
            temporary, metadata, transmission_qc_summary
        )
        for sample, relative in sample_files.items():
            (temporary / relative).write_text(
                _javascript_assignment(
                    "TAPESTRY_SAMPLE_SHARDS", sample, sample_rows[sample]
                ),
                encoding="utf-8",
            )
        for chromosome, relative in chromosome_files.items():
            (temporary / relative).write_text(
                _javascript_assignment(
                    "TAPESTRY_CHROMOSOME_SHARDS",
                    chromosome,
                    chromosome_rows[chromosome],
                ),
                encoding="utf-8",
            )
        for chromosome, relative in methylation_files.items():
            (temporary / relative).write_text(
                _javascript_assignment(
                    "TAPESTRY_METHYLATION_SHARDS",
                    chromosome,
                    methylation_rows[chromosome],
                ),
                encoding="utf-8",
            )
        files = sorted(
            str(path.relative_to(temporary))
            for path in temporary.rglob("*")
            if path.is_file()
        )
        bundle_manifest = {
            "schema_version": 1,
            "entrypoint": "index.html",
            "files": [*files, "bundle-manifest.json"],
            "samples": samples,
            "chromosomes": chromosomes,
            "blocks": len(payload["blocks"]),
            "methylation_samples": methylation_samples,
            "transmission_qc_comparisons": len(transmission_qc["comparisons"]),
        }
        (temporary / "bundle-manifest.json").write_text(
            json.dumps(bundle_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        # Nextflow's container process may run as root while its host-side
        # publisher and cleanup run as the invoking user. Make every generated
        # directory writable across that boundary; unlinking the root-owned
        # files only requires write permission on their parent directories.
        directories = [temporary, *(path for path in temporary.rglob("*") if path.is_dir())]
        for directory in directories:
            directory.chmod(0o777)
        temporary.replace(output)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ped", required=True, type=Path, help="normalized six-column PED")
    parser.add_argument("--iht", required=True, type=Path, help="sorted gtg inheritance map")
    parser.add_argument(
        "--selected-samples",
        type=Path,
        help="optional selected-samples.tsv used to choose the initial sample",
    )
    parser.add_argument(
        "--methylation-summaries",
        nargs="*",
        type=Path,
        default=[],
        help="optional per-sample methylation summary JSON files",
    )
    parser.add_argument(
        "--transmission-qc-summary",
        type=Path,
        help="optional chromosome-level parent-child methylation QC JSON",
    )
    parser.add_argument("--output", required=True, type=Path, help="output bundle directory")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        payload = build_payload(args.ped, args.iht, args.selected_samples)
        write_visualization(
            payload,
            args.output,
            args.methylation_summaries,
            args.transmission_qc_summary,
        )
    except (HaplotypePlotError, TransmissionPlotError, OSError, ValueError) as exc:
        print(f"haplotype ancestry visualization failed: {exc}", file=sys.stderr)
        return 2
    print(
        f"wrote {args.output} ({len(payload['samples'])} samples, "
        f"{len(payload['chromosomes'])} chromosomes, {len(payload['blocks'])} blocks)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
