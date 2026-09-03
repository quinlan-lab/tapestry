#!/usr/bin/env python3
"""Write the chromosome-level transmitted-haplotype methylation QC page."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


class TransmissionPlotError(ValueError):
    """Raised when a transmission QC summary cannot be plotted safely."""


COMPLETENESS_METRICS = (
    "missing_chromosomes",
    "callable_fraction",
    "mismatch_excluded_fraction",
    "paired_cpgs",
)
CONCORDANCE_METRICS = (
    "agreement",
    "discordant_fraction",
    "inherited_specificity",
    "mean_difference",
)
OVERVIEW_METRICS = COMPLETENESS_METRICS + CONCORDANCE_METRICS
HIGH_WHEN_CONCERNING = {
    "discordant_fraction",
    "missing_chromosomes",
    "mismatch_excluded_fraction",
}
LOW_WHEN_CONCERNING = {
    "agreement",
    "callable_fraction",
    "inherited_specificity",
    "paired_cpgs",
}


def edge_key(parent_id: str, child_id: str, relationship: str) -> str:
    return f"{parent_id}\0{child_id}\0{relationship}"


def median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[middle])
    return (ordered[middle - 1] + ordered[middle]) / 2.0


def scoring_value(metric: str, value: float) -> float:
    if metric == "paired_cpgs":
        return math.log10(value + 1.0)
    return float(value)


def concerning_deviation(value: float, center: float, metric: str) -> float:
    if metric in HIGH_WHEN_CONCERNING:
        return max(0.0, value - center)
    if metric in LOW_WHEN_CONCERNING:
        return max(0.0, center - value)
    return abs(value - center)


def robust_outliers(
    entries: list[tuple[str, float]], metric: str
) -> dict[str, dict[str, float | int | None]]:
    result: dict[str, dict[str, float | int | None]] = {}
    if not entries:
        return result
    raw_center = median([value for _key, value in entries])
    scored = [scoring_value(metric, value) for _key, value in entries]
    score_center = median(scored)
    if raw_center is None or score_center is None:
        return result
    absolute_deviations = [abs(value - score_center) for value in scored]
    scale = 1.4826 * (median(absolute_deviations) or 0.0)
    if scale <= 1e-12:
        mean = sum(scored) / len(scored)
        scale = math.sqrt(sum((value - mean) ** 2 for value in scored) / len(scored))
    countable = len(entries)
    for key, value in entries:
        scored_value = scoring_value(metric, value)
        score = None
        if countable >= 4:
            score = (
                concerning_deviation(scored_value, score_center, metric) / scale
                if scale > 1e-12
                else 0.0
            )
        result[key] = {
            "value": value,
            "center": raw_center,
            "score": score,
            "count": countable,
        }
    return result


def meets_metric_minimum(
    row: dict[str, Any], metric: str, minimum_paired_cpgs: int
) -> bool:
    if metric == "inherited_specificity":
        return int(row.get("specificity_cpgs") or 0) >= minimum_paired_cpgs
    if metric in {
        "missing_chromosomes",
        "mismatch_excluded_fraction",
        "paired_cpgs",
    }:
        return int(row.get("eligible_cpgs") or 0) > 0 or int(row.get("paired_cpgs") or 0) > 0
    return int(row.get("paired_cpgs") or 0) >= minimum_paired_cpgs


def _chromosome_value(
    row: dict[str, Any], metric: str
) -> tuple[float | None, int]:
    if metric == "mismatch_excluded_fraction":
        eligible = int(row.get("eligible_cpgs") or 0)
        excluded = int(row.get("mismatch_excluded_cpgs") or 0)
        return (excluded / eligible if eligible else None, eligible)
    if metric == "paired_cpgs":
        paired = int(row.get("paired_cpgs") or 0)
        return (float(paired), paired)
    if metric == "callable_fraction":
        return (
            None if row.get("callable_fraction") is None else float(row["callable_fraction"]),
            int(row.get("evaluated_cpgs") or 0),
        )
    if metric == "inherited_specificity":
        return (
            None
            if row.get("inherited_specificity") is None
            else float(row["inherited_specificity"]),
            int(row.get("specificity_cpgs") or 0),
        )
    if row.get(metric) is None:
        return None, int(row.get("paired_cpgs") or 0)
    return float(row[metric]), int(row.get("paired_cpgs") or 0)


def global_metric(
    rows: list[dict[str, Any]],
    chromosomes: list[str],
    metric: str,
    minimum_paired_cpgs: int,
) -> float | None:
    if metric == "missing_chromosomes":
        present = {
            str(row["chromosome"])
            for row in rows
            if meets_metric_minimum(row, "agreement", minimum_paired_cpgs)
            and row.get("agreement") is not None
        }
        return float(sum(1 for chromosome in chromosomes if chromosome not in present))
    if metric == "mismatch_excluded_fraction":
        eligible = sum(int(row.get("eligible_cpgs") or 0) for row in rows)
        excluded = sum(int(row.get("mismatch_excluded_cpgs") or 0) for row in rows)
        return None if eligible <= 0 else excluded / eligible
    if metric == "paired_cpgs":
        return float(sum(int(row.get("paired_cpgs") or 0) for row in rows))
    if metric == "callable_fraction":
        evaluated = sum(int(row.get("evaluated_cpgs") or 0) for row in rows)
        paired = sum(int(row.get("paired_cpgs") or 0) for row in rows)
        return None if evaluated <= 0 else paired / evaluated
    weight_key = (
        "specificity_cpgs" if metric == "inherited_specificity" else "paired_cpgs"
    )
    usable = [
        row
        for row in rows
        if meets_metric_minimum(row, metric, minimum_paired_cpgs)
        and row.get(metric) is not None
    ]
    weight = sum(int(row.get(weight_key) or 0) for row in usable)
    if weight < minimum_paired_cpgs:
        return None
    return sum(float(row[metric]) * int(row[weight_key]) for row in usable) / weight


def driver_chromosomes(
    rows: list[dict[str, Any]],
    chromosomes: list[str],
    metric: str,
    minimum_paired_cpgs: int,
    cohort_centers: dict[str, float],
    limit: int = 5,
) -> list[dict[str, float | int | str | None]]:
    if metric == "missing_chromosomes":
        present = {
            str(row["chromosome"])
            for row in rows
            if meets_metric_minimum(row, "agreement", minimum_paired_cpgs)
            and row.get("agreement") is not None
        }
        return [
            {
                "chromosome": chromosome,
                "value": None,
                "weight": 0,
                "contribution": None,
                "center": None,
            }
            for chromosome in chromosomes
            if chromosome not in present
        ][:limit]
    items: list[dict[str, float | int | str | None]] = []
    for row in rows:
        if metric in CONCORDANCE_METRICS and not meets_metric_minimum(
            row, metric, minimum_paired_cpgs
        ):
            continue
        value, weight = _chromosome_value(row, metric)
        if value is None:
            continue
        chromosome = str(row["chromosome"])
        center = cohort_centers.get(chromosome)
        contribution = float(weight)
        if center is not None:
            if metric == "paired_cpgs":
                contribution = max(0.0, center - value)
            else:
                scored_center = scoring_value(metric, center)
                contribution = weight * concerning_deviation(
                    scoring_value(metric, value), scored_center, metric
                )
        items.append(
            {
                "chromosome": chromosome,
                "value": value,
                "weight": weight,
                "contribution": contribution,
                "center": center,
            }
        )
    items.sort(
        key=lambda item: (
            -(float(item["contribution"] or 0.0)),
            -(int(item["weight"] or 0)),
            str(item["chromosome"]),
        )
    )
    return items[:limit]


def chromosome_centers(
    rows_by_edge: dict[str, list[dict[str, Any]]],
    metric: str,
    minimum_paired_cpgs: int,
) -> dict[str, float]:
    values: dict[str, list[float]] = {}
    for rows in rows_by_edge.values():
        for row in rows:
            if metric in CONCORDANCE_METRICS and not meets_metric_minimum(
                row, metric, minimum_paired_cpgs
            ):
                continue
            value, _weight = _chromosome_value(row, metric)
            if value is not None:
                values.setdefault(str(row["chromosome"]), []).append(value)
    return {
        chromosome: center
        for chromosome, chromosome_values in values.items()
        if (center := median(chromosome_values)) is not None
    }


def compute_overview(
    summary: dict[str, Any], chromosomes: list[str]
) -> dict[str, Any]:
    minimum = int(summary.get("minimum_paired_cpgs") or 100)
    rows_by_edge: dict[str, list[dict[str, Any]]] = {}
    for row in summary.get("comparisons") or []:
        key = edge_key(
            str(row["parent_id"]), str(row["child_id"]), str(row["relationship"])
        )
        rows_by_edge.setdefault(key, []).append(row)
    for edge in summary.get("edges") or []:
        key = edge_key(
            str(edge["parent_id"]), str(edge["child_id"]), str(edge["relationship"])
        )
        rows_by_edge.setdefault(key, [])
    metrics: dict[str, Any] = {}
    for metric in OVERVIEW_METRICS:
        entries: list[tuple[str, float]] = []
        values: dict[str, float | None] = {}
        for key, rows in rows_by_edge.items():
            value = global_metric(rows, chromosomes, metric, minimum)
            values[key] = value
            if value is not None:
                entries.append((key, value))
        stats = robust_outliers(entries, metric)
        raw_center = median([value for _key, value in entries])
        per_chromosome_centers = (
            {}
            if metric == "missing_chromosomes"
            else chromosome_centers(rows_by_edge, metric, minimum)
        )
        metrics[metric] = {
            "center": raw_center,
            "count": len(entries),
            "edges": {
                key: {
                    "value": values[key],
                    "center": stats.get(key, {}).get("center", raw_center),
                    "score": stats.get(key, {}).get("score"),
                    "count": stats.get(key, {}).get("count", len(entries)),
                    "drivers": driver_chromosomes(
                        rows_by_edge[key],
                        chromosomes,
                        metric,
                        minimum,
                        per_chromosome_centers,
                    ),
                }
                for key in rows_by_edge
            },
        }
    return {"metrics": metrics}


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Tapestry transmission methylation QC</title>
<script src="plotly.min.js"></script>
<script src="metadata.js"></script>
<script src="transmission-qc.js"></script>
<style>
:root { color-scheme: light; font-family: Inter, ui-sans-serif, system-ui, sans-serif; }
body { margin: 0; color: #172033; background: #f5f7fb; }
.page { max-width: 1700px; margin: 0 auto; padding: 24px; }
h1 { margin: 0 0 6px; font-size: 25px; }
.lede,.subtle { color: #5d687b; }
.lede { margin: 0 0 14px; max-width: 1100px; line-height: 1.45; }
.subtle { margin: 10px 0; font-size: 13px; line-height: 1.45; }
.nav { display: flex; gap: 8px; margin: 0 0 18px; }
.nav a { color: #334e91; background: white; border: 1px solid #ccd5e5; border-radius: 6px; padding: 6px 10px; text-decoration: none; font-size: 13px; }
.controls { display: flex; flex-wrap: wrap; gap: 14px; align-items: end; margin-bottom: 14px; }
label { display: grid; gap: 5px; color: #445067; font-size: 13px; font-weight: 650; }
select { min-width: 175px; padding: 7px 30px 7px 9px; border: 1px solid #cbd2df; border-radius: 6px; background: white; }
#pair { min-width: 300px; }
.badge { display: inline-block; padding: 2px 8px; border-radius: 999px; background: #edf1f8; font-size: 12px; }
.card { background: white; border: 1px solid #dde2eb; border-radius: 10px; box-shadow: 0 1px 4px rgb(30 45 70 / 8%); }
.split { display: grid; grid-template-columns: minmax(0, 1.45fr) minmax(280px, 0.85fr); gap: 12px; margin: 12px 0; }
.strip-row { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; margin: 12px 0; }
@media (max-width: 1100px) { .split, .strip-row { grid-template-columns: 1fr; } }
#overview { min-height: 430px; }
#strip-agreement, #strip-discordance { min-height: 290px; }
#heatmap { min-height: 260px; }
#missing { min-height: 250px; margin-top: 18px; }
#empty { display: none; padding: 24px; line-height: 1.5; }
#rank-table, #drivers { padding: 12px 14px 16px; overflow: auto; }
#rank-table h2, #drivers h2 { margin: 0 0 8px; font-size: 16px; }
table.rank { width: 100%; border-collapse: collapse; font-size: 12.5px; }
table.rank th, table.rank td { padding: 6px 7px; border-bottom: 1px solid #e6ebf2; text-align: left; vertical-align: top; }
table.rank th { color: #5d687b; font-weight: 650; }
table.rank tr { cursor: pointer; }
table.rank tr:hover { background: #f3f6fb; }
table.rank tr.selected { background: #e8eef9; }
.mono { font-variant-numeric: tabular-nums; }
</style>
</head>
<body>
<main class="page">
  <h1>Parent–child methylation concordance</h1>
  <p class="lede">Each comparison uses model-based methylation at the same callable CpGs on the haplotype transmitted from parent to child. Mismatch-proximal CpGs are excluded. Completeness scores ask whether the pair is well measured; concordance scores ask whether the transmitted homolog matches. Neither establishes epigenetic inheritance by itself.</p>
  <nav class="nav"><a href="index.html">Haplotype ancestry</a><a href="transmission-qc.html" aria-current="page">Transmission QC</a></nav>
  <div class="controls">
    <label>Parent–child pair<select id="pair"></select></label>
    <label>Cell metric<select id="metric"></select></label>
    <label>Chromosome cell colors<select id="color-mode"></select></label>
    <label>Pair order<select id="sort-order"></select></label>
    <label>Compare pairs<select id="row-limit"></select></label>
    <span id="threshold" class="badge"></span>
    <span id="minimum" class="badge"></span>
    <span id="coverage" class="badge"></span>
  </div>
  <div id="overview" class="card"></div>
  <p class="subtle">The genome-wide overview keeps every pedigree transmission pair visible. Completeness rows are missing chromosomes, callable fraction, mismatch-excluded fraction, and paired-CpG yield. Concordance rows are agreement, large-discordance fraction, inherited specificity, and signed difference. Color is robust outlier severity versus the cohort median; cells stay blank when fewer than four pairs can be scored. Hover reports the raw value and cohort median. Click a column to inspect that pair.</p>
  <div class="split">
    <div id="rank-table" class="card"></div>
    <div id="drivers" class="card"></div>
  </div>
  <div class="strip-row">
    <div id="strip-agreement" class="card"></div>
    <div id="strip-discordance" class="card"></div>
  </div>
  <p class="subtle">Strip plots show genome-wide CpG-weighted agreement and large-discordance fraction, split by paternal versus maternal transmission. The selected pair is outlined. Click a point to inspect that pair.</p>
  <div id="heatmap" class="card"></div>
  <div id="empty" class="card"></div>
  <p class="subtle">Selecting one pair shows only that pair’s chromosomes. Use Compare pairs to overlay the strongest outliers. In chromosome-relative mode, each cell is compared with all eligible pairs on that chromosome, even when one pair is selected. Concordance genome-wide means use only chromosomes that meet the configured minimum paired CpGs. Agreement is 1 − mean absolute parent–child difference. Inherited specificity is the non-transmitted parental haplotype MAE minus transmitted-haplotype MAE at CpGs callable in all three measurements; positive values mean the child is closer to the transmitted haplotype. Click a populated chromosome cell to open that child in the ancestry view.</p>
  <div id="missing" class="card"></div>
  <p class="subtle">Missing comparisons count pedigree transmission edges with no chromosome-level estimate for the selected metric. The chart separates absent sample output from edges whose samples have output but fewer than the configured minimum contributing CpGs on that chromosome.</p>
</main>
<script>
const DATA = window.TAPESTRY_METADATA;
const QC = window.TAPESTRY_TRANSMISSION_QC;
const overview = document.getElementById('overview');
const heatmap = document.getElementById('heatmap');
const missingPlot = document.getElementById('missing');
const empty = document.getElementById('empty');
const rankTable = document.getElementById('rank-table');
const drivers = document.getElementById('drivers');
const stripAgreement = document.getElementById('strip-agreement');
const stripDiscordance = document.getElementById('strip-discordance');
const pairSelect = document.getElementById('pair');
const metricSelect = document.getElementById('metric');
const colorModeSelect = document.getElementById('color-mode');
const sortOrderSelect = document.getElementById('sort-order');
const rowLimitSelect = document.getElementById('row-limit');
const metrics = {
  agreement: {label: 'Inherited-haplotype agreement', range: [0, 1], colorscale: [[0,'#b2182b'],[0.5,'#f7f7f7'],[1,'#1a9850']], family: 'concordance', signed: false, integer: false},
  mean_difference: {label: 'Signed child − parent difference', range: [-1, 1], colorscale: [[0,'#2166ac'],[0.5,'#f7f7f7'],[1,'#b2182b']], family: 'concordance', signed: true, integer: false},
  discordant_fraction: {label: 'Large-discordance fraction', range: [0, 1], colorscale: [[0,'#fff7ec'],[1,'#b30000']], family: 'concordance', signed: false, integer: false},
  inherited_specificity: {label: 'Inherited specificity', range: [-1, 1], colorscale: [[0,'#b2182b'],[0.5,'#f7f7f7'],[1,'#2166ac']], family: 'concordance', signed: true, integer: false},
  callable_fraction: {label: 'Callable fraction', range: [0, 1], colorscale: [[0,'#f7fbff'],[1,'#08519c']], family: 'completeness', signed: false, integer: false},
  missing_chromosomes: {label: 'Missing chromosomes', family: 'completeness', signed: false, integer: true},
  mismatch_excluded_fraction: {label: 'Mismatch-excluded fraction', family: 'completeness', signed: false, integer: false},
  paired_cpgs: {label: 'Paired CpG yield', family: 'completeness', signed: false, integer: true}
};
const chromosomeMetrics = ['agreement', 'mean_difference', 'discordant_fraction', 'callable_fraction', 'inherited_specificity'];
const completenessMetrics = ['missing_chromosomes', 'callable_fraction', 'mismatch_excluded_fraction', 'paired_cpgs'];
const concordanceMetrics = ['agreement', 'discordant_fraction', 'inherited_specificity', 'mean_difference'];
const overviewMetrics = [...completenessMetrics, ...concordanceMetrics];
for (const key of chromosomeMetrics) {
  const option = document.createElement('option'); option.value = key; option.textContent = metrics[key].label;
  metricSelect.appendChild(option);
}
for (const [value, label] of [['outlier', 'Chromosome-relative outliers'], ['absolute', 'Absolute metric values']]) {
  const option = document.createElement('option'); option.value = value; option.textContent = label; colorModeSelect.appendChild(option);
}
for (const [value, label] of [['concordance', 'Concordance outliers first'], ['completeness', 'Completeness outliers first'], ['pedigree', 'Pedigree label']]) {
  const option = document.createElement('option'); option.value = value; option.textContent = label; sortOrderSelect.appendChild(option);
}
for (const [value, label] of [['selected', 'Selected pair only'], ['25', 'Top 25 pairs'], ['50', 'Top 50 pairs'], ['all', 'All pairs']]) {
  const option = document.createElement('option'); option.value = value; option.textContent = label; rowLimitSelect.appendChild(option);
}
document.getElementById('threshold').textContent = `large difference ≥ ${QC.discordance_threshold.toFixed(2)}`;
document.getElementById('minimum').textContent = `minimum ${QC.minimum_paired_cpgs.toLocaleString()} contributing CpGs`;
document.getElementById('coverage').textContent = `${QC.samples.length} sample${QC.samples.length === 1 ? '' : 's'} with methylation output`;
function generation(sample) { return DATA.people[sample] ? `F${DATA.people[sample].generation}` : ''; }
function edgeKey(row) { return `${row.parent_id}\u0000${row.child_id}\u0000${row.relationship}`; }
function recordKey(row) { return `${edgeKey(row)}\u0000${row.chromosome}`; }
function edgeLabel(row) { return `${generation(row.parent_id)} ${row.parent_id} → ${generation(row.child_id)} ${row.child_id} · ${row.relationship}`; }
function edgeShortLabel(row) { return `${row.parent_id}→${row.child_id} · ${row.relationship === 'paternal' ? 'pat' : 'mat'}`; }
function fmt(value, signed=false, integer=false) {
  if (value === null || value === undefined) return 'not available';
  if (integer) return Number(value).toLocaleString();
  const text = Number(value).toFixed(3);
  return signed && value > 0 ? `+${text}` : text;
}
function fmtScore(score) {
  return score === null || score === undefined ? 'n<4' : Number(score).toFixed(2);
}
const edgeDefinitions = new Map(QC.edges.map(edge => [edgeKey(edge), edge]));
for (const row of QC.comparisons) if (!edgeDefinitions.has(edgeKey(row))) edgeDefinitions.set(edgeKey(row), row);
const allEdges = [...edgeDefinitions.keys()].sort((left, right) => edgeLabel(edgeDefinitions.get(left)).localeCompare(edgeLabel(edgeDefinitions.get(right)), undefined, {numeric: true}));
const edgeIndexes = new Map(allEdges.map((edge, index) => [edge, index]));
const allOption = document.createElement('option'); allOption.value = ''; allOption.textContent = `All pairs (${allEdges.length})`; pairSelect.appendChild(allOption);
for (const [index, edge] of allEdges.entries()) {
  const option = document.createElement('option'); option.value = String(index); option.textContent = edgeLabel(edgeDefinitions.get(edge)); pairSelect.appendChild(option);
}
const records = new Map(QC.comparisons.map(row => [recordKey(row), row]));
const rowsByEdge = new Map(allEdges.map(edge => [edge, []]));
for (const row of QC.comparisons) {
  const edge = edgeKey(row);
  if (!rowsByEdge.has(edge)) rowsByEdge.set(edge, []);
  rowsByEdge.get(edge).push(row);
}
function selectedEdgeKey() { return pairSelect.value === '' ? null : allEdges[Number(pairSelect.value)]; }
function selectedEdges() { return selectedEdgeKey() === null ? allEdges : [selectedEdgeKey()]; }
function rankFamily() { return sortOrderSelect.value === 'completeness' ? completenessMetrics : concordanceMetrics; }
function meetsMetricMinimum(row, metric) {
  if (metric === 'inherited_specificity') return row.specificity_cpgs >= QC.minimum_paired_cpgs;
  if (metric === 'missing_chromosomes' || metric === 'mismatch_excluded_fraction' || metric === 'paired_cpgs') {
    return row.eligible_cpgs > 0 || row.paired_cpgs > 0;
  }
  return row.paired_cpgs >= QC.minimum_paired_cpgs;
}
function median(values) {
  if (!values.length) return null;
  const ordered = [...values].sort((left, right) => left - right);
  const middle = Math.floor(ordered.length / 2);
  return ordered.length % 2 ? ordered[middle] : (ordered[middle - 1] + ordered[middle]) / 2;
}
function concerningDeviation(value, center, metric) {
  if (metric === 'discordant_fraction' || metric === 'missing_chromosomes' || metric === 'mismatch_excluded_fraction') return Math.max(0, value - center);
  if (metric === 'agreement' || metric === 'callable_fraction' || metric === 'inherited_specificity' || metric === 'paired_cpgs') return Math.max(0, center - value);
  return Math.abs(value - center);
}
function robustOutliers(entries, metric) {
  const result = new Map();
  const values = entries.map(entry => entry.value);
  const center = median(values);
  if (center === null) return result;
  const scored = metric === 'paired_cpgs' ? values.map(value => Math.log10(value + 1)) : values;
  const scoreCenter = metric === 'paired_cpgs' ? median(scored) : center;
  const absoluteDeviations = scored.map(value => Math.abs(value - scoreCenter));
  let scale = 1.4826 * median(absoluteDeviations);
  if (!(scale > 1e-12)) {
    const mean = scored.reduce((total, value) => total + value, 0) / scored.length;
    scale = Math.sqrt(scored.reduce((total, value) => total + (value - mean) ** 2, 0) / scored.length);
  }
  for (const [index, entry] of entries.entries()) {
    const score = entries.length < 4 ? null : (scale > 1e-12 ? concerningDeviation(scored[index], scoreCenter, metric) / scale : 0);
    result.set(entry.key, {center, count: entries.length, score});
  }
  return result;
}
function overviewRow(edge, metric) {
  return QC.overview?.metrics?.[metric]?.edges?.[edge] || {};
}
function metricValue(edge, metric) {
  const value = overviewRow(edge, metric).value;
  return value === undefined ? null : value;
}
function metricStat(edge, metric) {
  const row = overviewRow(edge, metric);
  if (row.value === null || row.value === undefined) return undefined;
  return {center: row.center, count: row.count, score: row.score};
}
function edgeOutlierRank(edge, family) {
  return Math.max(0, ...family.map(metric => metricStat(edge, metric)?.score || 0));
}
function strongestMetric(edge, family) {
  let best = family[0], bestScore = -1;
  for (const metric of family) {
    const score = metricStat(edge, metric)?.score;
    if (score !== null && score !== undefined && score > bestScore) {
      best = metric; bestScore = score;
    }
  }
  if (bestScore >= 0) return {metric: best, score: bestScore};
  for (const metric of family) {
    if (metricValue(edge, metric) !== null) return {metric, score: null};
  }
  return {metric: family[0], score: null};
}
function orderedEdges() {
  if (sortOrderSelect.value === 'pedigree') return [...allEdges];
  return rankedEdges(rankFamily());
}
function rankedEdges(family) {
  return [...allEdges].sort((left, right) => edgeOutlierRank(right, family) - edgeOutlierRank(left, family) || edgeIndexes.get(left) - edgeIndexes.get(right));
}
function selectPair(index) {
  pairSelect.value = index === null || index === undefined || index === '' ? '' : String(index);
  if (pairSelect.value !== '') rowLimitSelect.value = 'selected';
  renderAll();
}
function renderOverview() {
  const edges = orderedEdges();
  const z = [], customdata = [], y = [];
  for (const metric of overviewMetrics) {
    const definition = metrics[metric];
    y.push(`${definition.family === 'completeness' ? 'Completeness' : 'Concordance'} · ${definition.label}`);
    const zRow = [], customRow = [];
    for (const edge of edges) {
      const value = metricValue(edge, metric);
      const stat = metricStat(edge, metric);
      zRow.push(value === null || stat?.score === null || stat?.score === undefined ? null : Math.min(stat.score, 4));
      customRow.push([
        edgeLabel(edgeDefinitions.get(edge)), definition.label,
        fmt(value, definition.signed, definition.integer),
        stat ? fmt(stat.center, definition.signed, definition.integer) : 'not available',
        fmtScore(stat?.score),
        stat?.count || 0, edgeIndexes.get(edge)
      ]);
    }
    z.push(zRow); customdata.push(customRow);
  }
  const showLabels = edges.length <= 40;
  const height = showLabels ? 520 : 430; overview.style.height = `${height}px`;
  Plotly.newPlot(overview, [{
    type: 'heatmap', x: edges.map(edge => edgeShortLabel(edgeDefinitions.get(edge))), y,
    z, customdata, zmin: 0, zmax: 4, colorscale: [[0,'#fff7ec'],[0.35,'#fdd49e'],[0.7,'#fc8d59'],[1,'#b30000']], hoverongaps: false,
    colorbar: {title: {text: 'Robust outlier<br>score', side: 'right'}},
    hovertemplate: '<b>%{customdata[0]}</b><br>%{customdata[1]}: %{customdata[2]}<br>Cohort median: %{customdata[3]}<br>Robust outlier score: %{customdata[4]}<br>Eligible pairs: %{customdata[5]}<extra></extra>'
  }], {
    template: 'plotly_white', height, margin: {l: 320, r: 100, t: 65, b: showLabels ? 115 : 35},
    title: {text: `Genome-wide overview · all ${edges.length} parent–child pairs · click a column to inspect`, x: 0.01, xanchor: 'left'},
    xaxis: {title: showLabels ? 'Parent–child pair' : '', showticklabels: showLabels, tickangle: -45, fixedrange: false},
    yaxis: {fixedrange: true, autorange: 'reversed'}, hoverlabel: {align: 'left'}, uirevision: `overview-${sortOrderSelect.value}`
  }, {responsive: true, displaylogo: false}).then(() => {
    overview.removeAllListeners?.('plotly_click');
    overview.on('plotly_click', event => {
      const index = event.points[0]?.customdata?.[6]; if (index === undefined) return;
      selectPair(index);
      drivers.scrollIntoView({behavior: 'smooth', block: 'nearest'});
    });
  });
}
function renderRankTable() {
  const family = rankFamily();
  const familyName = sortOrderSelect.value === 'completeness' ? 'completeness' : 'concordance';
  const rows = rankedEdges(family).slice(0, 10);
  const selected = selectedEdgeKey();
  const body = rows.map(edge => {
    const definition = edgeDefinitions.get(edge);
    const strongest = strongestMetric(edge, family);
    const otherFamily = familyName === 'completeness' ? concordanceMetrics : completenessMetrics;
    const other = strongestMetric(edge, otherFamily);
    const value = metricValue(edge, strongest.metric);
    const stat = metricStat(edge, strongest.metric);
    const paired = metricValue(edge, 'paired_cpgs');
    const driversText = (overviewRow(edge, strongest.metric).drivers || []).slice(0, 3).map(item => String(item.chromosome).replace(/^chr/, '')).join(', ') || '—';
    return `<tr data-index="${edgeIndexes.get(edge)}" class="${edge === selected ? 'selected' : ''}">
      <td>${edgeLabel(definition)}</td>
      <td>${metrics[strongest.metric].label}<br><span class="subtle">${fmt(value, metrics[strongest.metric].signed, metrics[strongest.metric].integer)} vs ${stat ? fmt(stat.center, metrics[strongest.metric].signed, metrics[strongest.metric].integer) : 'n/a'}</span></td>
      <td class="mono">${fmtScore(strongest.score)}</td>
      <td class="mono">${fmtScore(other.score)} <span class="subtle">${metrics[other.metric].label}</span></td>
      <td class="mono">${fmt(paired, false, true)}</td>
      <td>${driversText}</td>
    </tr>`;
  }).join('');
  rankTable.innerHTML = `<h2>Ranked ${familyName} outliers</h2>
    <p class="subtle">Top ${rows.length} pair${rows.length === 1 ? '' : 's'} by ${familyName} score. Click a row to inspect chromosomes.</p>
    <table class="rank"><thead><tr>
      <th>Pair</th><th>Strongest ${familyName} metric</th><th>Score</th><th>Other family</th><th>Paired CpGs</th><th>Driver chr</th>
    </tr></thead><tbody>${body || '<tr><td colspan="6">No pairs available.</td></tr>'}</tbody></table>`;
  rankTable.querySelectorAll('tr[data-index]').forEach(row => {
    row.addEventListener('click', () => {
      selectPair(row.getAttribute('data-index'));
    });
  });
}
function renderDrivers() {
  const edge = selectedEdgeKey();
  if (edge === null) {
    drivers.innerHTML = '<h2>Chromosomes that drive the score</h2><p class="subtle">Select one pair to see the chromosomes that dominate its genome-wide completeness or concordance number.</p>';
    return;
  }
  const family = rankFamily();
  const strongest = strongestMetric(edge, family);
  const items = overviewRow(edge, strongest.metric).drivers || [];
  const definition = metrics[strongest.metric];
  const rows = items.map(item => `<tr>
      <td>${item.chromosome}</td>
      <td class="mono">${item.value === null || item.value === undefined ? 'missing' : fmt(item.value, definition.signed, definition.integer)}</td>
      <td class="mono">${item.center === null || item.center === undefined ? '—' : fmt(item.center, definition.signed, definition.integer)}</td>
      <td class="mono">${item.contribution === null || item.contribution === undefined ? '—' : Number(item.contribution).toLocaleString(undefined, {maximumFractionDigits: 2})}</td>
      <td class="mono">${item.weight ? Number(item.weight).toLocaleString() : '—'}</td>
    </tr>`).join('');
  drivers.innerHTML = `<h2>Chromosomes that drive the score</h2>
    <p class="subtle">${edgeLabel(edgeDefinitions.get(edge))}<br>Strongest ${sortOrderSelect.value === 'completeness' ? 'completeness' : 'concordance'} metric: ${definition.label} (${fmtScore(strongest.score)}). Drivers are ranked by their metric-specific deficit from the same chromosome’s cohort median.</p>
    <table class="rank"><thead><tr><th>Chromosome</th><th>Value</th><th>Chromosome median</th><th>Driver contribution</th><th>Weight</th></tr></thead>
    <tbody>${rows || '<tr><td colspan="5">No chromosome-level values for this metric.</td></tr>'}</tbody></table>`;
}
function renderStrip(node, metric, title) {
  const definition = metrics[metric];
  const selected = selectedEdgeKey();
  const traces = ['paternal', 'maternal'].map((relationship, index) => {
    const edges = allEdges.filter(edge => edgeDefinitions.get(edge).relationship === relationship);
    const x = [], y = [], customdata = [], selectedFlags = [];
    edges.forEach((edge, order) => {
      const value = metricValue(edge, metric);
      if (value === null) return;
      const jitter = ((order * 17) % 10) / 10;
      x.push(index + (jitter - 0.45) * 0.28);
      y.push(value);
      selectedFlags.push(edge === selected);
      customdata.push([edgeLabel(edgeDefinitions.get(edge)), fmt(value, definition.signed, definition.integer), edgeIndexes.get(edge)]);
    });
    return {
      type: 'scatter', mode: 'markers', name: relationship,
      x, y, customdata,
      marker: {
        size: selectedFlags.map(flag => flag ? 12 : 8),
        color: relationship === 'paternal' ? '#2b6cb0' : '#c05621',
        line: {width: selectedFlags.map(flag => flag ? 2 : 0), color: '#111827'}
      },
      hovertemplate: `<b>%{customdata[0]}</b><br>${definition.label}: %{customdata[1]}<extra></extra>`
    };
  });
  Plotly.newPlot(node, traces, {
    template: 'plotly_white', height: 290, margin: {l: 55, r: 20, t: 45, b: 45},
    title: {text: title, x: 0.01, xanchor: 'left'},
    xaxis: {tickvals: [0, 1], ticktext: ['paternal', 'maternal'], range: [-0.5, 1.5], fixedrange: true, zeroline: false},
    yaxis: {title: definition.label, autorange: true, fixedrange: true},
    showlegend: false, hoverlabel: {align: 'left'}, uirevision: `strip-${metric}-${selected || 'all'}`
  }, {responsive: true, displaylogo: false}).then(() => {
    node.removeAllListeners?.('plotly_click');
    node.on('plotly_click', event => {
      const index = event.points[0]?.customdata?.[2]; if (index === undefined) return;
      selectPair(index);
    });
  });
}
function chromosomeOutliers(metric) {
  const result = new Map();
  for (const chromosome of DATA.chromosomes) {
    const entries = QC.comparisons
      .filter(row => row.chromosome === chromosome && meetsMetricMinimum(row, metric) && (
        metric === 'mismatch_excluded_fraction'
          ? row.eligible_cpgs > 0
          : metric === 'paired_cpgs'
            ? row.paired_cpgs !== null && row.paired_cpgs !== undefined
            : row[metric] !== null && row[metric] !== undefined
      ))
      .map(row => ({
        key: recordKey(row),
        value: metric === 'mismatch_excluded_fraction'
          ? row.mismatch_excluded_cpgs / row.eligible_cpgs
          : metric === 'paired_cpgs' ? row.paired_cpgs : row[metric]
      }));
    for (const [key, value] of robustOutliers(entries, metric)) result.set(key, value);
  }
  return result;
}
function renderMissingEdges() {
  const selected = new Set(selectedEdges());
  const metric = metricSelect.value || 'agreement';
  const scopedEdges = QC.edges.filter(edge => selected.has(edgeKey(edge)));
  const edgesWithOutput = scopedEdges.filter(edge => edge.has_methylation_outputs);
  const absentOutput = scopedEdges.length - edgesWithOutput.length;
  const absentByChromosome = DATA.chromosomes.map(() => absentOutput);
  const insufficientByChromosome = DATA.chromosomes.map(chromosome => {
    const comparable = new Set(QC.comparisons.filter(row => selected.has(edgeKey(row)) && row.chromosome === chromosome && meetsMetricMinimum(row, metric)).map(edgeKey));
    return Math.max(edgesWithOutput.length - comparable.size, 0);
  });
  Plotly.newPlot(missingPlot, [
    {type: 'bar', name: 'sample output absent', x: DATA.chromosomes, y: absentByChromosome, marker: {color: '#9da7b8'}, hovertemplate: '<b>%{x}</b><br>Missing sample output: %{y} edges<extra></extra>'},
    {type: 'bar', name: 'insufficient contributing CpGs', x: DATA.chromosomes, y: insufficientByChromosome, marker: {color: '#dd6b55'}, hovertemplate: `<b>%{x}</b><br>Fewer than ${QC.minimum_paired_cpgs.toLocaleString()} CpGs for ${metrics[metric].label}: %{y} edges<extra></extra>`}
  ], {
    barmode: 'stack', template: 'plotly_white', height: 250, margin: {l: 62, r: 30, t: 55, b: 55},
    title: {text: `Missing transmission comparisons · ${scopedEdges.length} parent–child pair${scopedEdges.length === 1 ? '' : 's'} selected`, x: 0.01, xanchor: 'left'},
    xaxis: {title: 'Chromosome', fixedrange: true}, yaxis: {title: 'Missing edges', rangemode: 'tozero', dtick: 1, fixedrange: true},
    legend: {orientation: 'h', y: 1.16, x: 1, xanchor: 'right'}
  }, {responsive: true, displaylogo: false});
}
function heatmapEdges() {
  const selected = selectedEdgeKey();
  const metric = metricSelect.value || 'agreement';
  const comparable = orderedEdges().filter(edge => (rowsByEdge.get(edge) || []).some(row => meetsMetricMinimum(row, metric) && row[metric] !== null && row[metric] !== undefined));
  if (selected !== null || rowLimitSelect.value === 'selected') {
    return {edges: selected !== null && comparable.includes(selected) ? [selected] : [], available: comparable.length};
  }
  if (rowLimitSelect.value === 'all') return {edges: comparable, available: comparable.length};
  return {edges: comparable.slice(0, Number(rowLimitSelect.value)), available: comparable.length};
}
function render() {
  const metric = metricSelect.value || 'agreement';
  const colorMode = colorModeSelect.value || 'absolute';
  const outliers = chromosomeOutliers(metric);
  const {edges, available} = heatmapEdges();
  rowLimitSelect.disabled = selectedEdgeKey() !== null;
  if (!edges.length) {
    heatmap.style.display = 'none'; empty.style.display = 'block';
    empty.innerHTML = selectedEdgeKey() === null
      ? '<strong>Select a parent–child pair to inspect chromosomes.</strong><br>Use the overview, ranked table, or strip plot. Compare pairs overlays multiple transmissions. The missing-comparison chart below reports edges without enough CpGs.'
      : '<strong>No chromosome-level comparison is available for the selected parent–child pair and metric.</strong><br>Both samples need downstream founder-phased methylation output and at least the configured minimum shared evaluable CpGs. The missing-comparison chart below identifies the reason.';
    return;
  }
  heatmap.style.display = 'block'; empty.style.display = 'none';
  const definition = metrics[metric];
  const z = [], customdata = [];
  for (const edge of edges) {
    const zRow = [], customRow = [];
    for (const chromosome of DATA.chromosomes) {
      const row = records.get(`${edge}\u0000${chromosome}`);
      const usable = row && meetsMetricMinimum(row, metric) && row[metric] !== null && row[metric] !== undefined;
      const stat = usable ? outliers.get(recordKey(row)) : null;
      zRow.push(!usable || (colorMode === 'outlier' && (stat?.score === null || stat?.score === undefined)) ? null : colorMode === 'outlier' ? Math.min(stat.score, 4) : row[metric]);
      customRow.push(row ? [
        edgeLabel(row), chromosome, fmt(row.agreement), fmt(row.mean_difference, true),
        fmt(row.discordant_fraction), row.paired_cpgs.toLocaleString(),
        row.evaluated_cpgs.toLocaleString(), fmt(row.callable_fraction),
        row.mismatch_excluded_cpgs.toLocaleString(), fmt(row.inherited_specificity, true),
        row.specificity_cpgs.toLocaleString(), row.child_id,
        stat ? fmt(stat.center, metric === 'mean_difference') : 'not available',
        fmtScore(stat?.score),
        stat?.count || 0, definition.label, fmt(row[metric], metric === 'mean_difference'),
        Number(row.ambiguous_cpgs || 0).toLocaleString()
      ] : [edgeLabel(edgeDefinitions.get(edge)), chromosome, 'not available', 'not available',
           'not available', '0', '0', 'not available', '0', 'not available', '0', edgeDefinitions.get(edge).child_id,
           'not available', 'not available', '0', definition.label, 'not available', '0']);
    }
    z.push(zRow); customdata.push(customRow);
  }
  const height = Math.max(280, Math.min(900, edges.length * 42 + 190)); heatmap.style.height = `${height}px`;
  const outlierMode = colorMode === 'outlier';
  const titleSuffix = selectedEdgeKey() !== null
    ? edgeLabel(edgeDefinitions.get(selectedEdgeKey()))
    : available === edges.length ? `${edges.length} pairs` : `showing ${edges.length} of ${available} comparable pairs`;
  Plotly.newPlot(heatmap, [{
    type: 'heatmap', x: DATA.chromosomes, y: edges.map(edge => edgeLabel(edgeDefinitions.get(edge))),
    z, customdata, zmin: outlierMode ? 0 : definition.range[0], zmax: outlierMode ? 4 : definition.range[1], zmid: !outlierMode && definition.range[0] < 0 ? 0 : undefined,
    colorscale: outlierMode ? [[0,'#fff7ec'],[0.35,'#fdd49e'],[0.7,'#fc8d59'],[1,'#b30000']] : definition.colorscale, hoverongaps: false,
    colorbar: {title: {text: outlierMode ? 'Robust outlier<br>score' : definition.label, side: 'right'}},
    hovertemplate: '<b>%{customdata[0]}</b><br>%{customdata[1]}<br>%{customdata[15]}: %{customdata[16]}<br>Chromosome median: %{customdata[12]}<br>Robust outlier score: %{customdata[13]} (%{customdata[14]} eligible pairs)<br><br>Agreement: %{customdata[2]}<br>Signed child − parent: %{customdata[3]}<br>Large-discordance fraction: %{customdata[4]}<br>Paired CpGs: %{customdata[5]} / %{customdata[6]} evaluated<br>Callable fraction: %{customdata[7]}<br>Ambiguous duplicate CpGs excluded: %{customdata[17]}<br>Mismatch-proximal excluded: %{customdata[8]}<br>Inherited specificity: %{customdata[9]} (%{customdata[10]} CpGs)<extra></extra>'
  }], {
    template: 'plotly_white', height, margin: {l: 300, r: 105, t: 65, b: 65},
    title: {text: `${outlierMode ? 'Chromosome-relative outliers' : definition.label} · ${titleSuffix}`, x: 0.01, xanchor: 'left'},
    xaxis: {title: 'Chromosome', side: 'bottom', fixedrange: true},
    yaxis: {title: 'Transmitted haplotype', automargin: true, fixedrange: true, autorange: 'reversed', tickfont: {size: edges.length > 20 ? 10 : 12}},
    hoverlabel: {align: 'left'}, uirevision: `${metric}-${colorMode}-${sortOrderSelect.value}-${pairSelect.value}-${rowLimitSelect.value}`
  }, {responsive: true, displaylogo: false}).then(() => {
    heatmap.removeAllListeners?.('plotly_click');
    heatmap.on('plotly_click', event => {
      const point = event.points[0]; if (point?.z === null || !point?.customdata) return;
      const chromosome = point.customdata[1], child = point.customdata[11];
      location.href = `index.html?sample=${encodeURIComponent(child)}&chrom=${encodeURIComponent(chromosome)}`;
    });
  });
}
function renderAll() {
  renderOverview();
  renderRankTable();
  renderDrivers();
  renderStrip(stripAgreement, 'agreement', 'Genome-wide agreement');
  renderStrip(stripDiscordance, 'discordant_fraction', 'Genome-wide large-discordance fraction');
  renderMissingEdges();
  render();
}
metricSelect.value = 'agreement';
colorModeSelect.value = allEdges.length >= 4 ? 'outlier' : 'absolute';
sortOrderSelect.value = 'concordance';
rowLimitSelect.value = 'selected';
const initialRanked = orderedEdges();
if (initialRanked.length) pairSelect.value = String(edgeIndexes.get(initialRanked[0]));
metricSelect.addEventListener('change', () => { renderMissingEdges(); render(); });
pairSelect.addEventListener('change', () => {
  if (pairSelect.value !== '') rowLimitSelect.value = 'selected';
  renderAll();
});
colorModeSelect.addEventListener('change', render);
sortOrderSelect.addEventListener('change', renderAll);
rowLimitSelect.addEventListener('change', render);
renderAll();
</script>
</body>
</html>
"""


def read_summary(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {
            "schema_version": 1,
            "measurement": "model-based methylation",
            "discordance_threshold": 0.4,
            "minimum_paired_cpgs": 100,
            "samples": [],
            "edges": [],
            "comparisons": [],
            "unavailable_edges": [],
            "sources": {},
        }
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema_version") != 1:
        raise TransmissionPlotError(f"{path}: unsupported transmission QC schema")
    for key in ("samples", "edges", "comparisons", "unavailable_edges"):
        if not isinstance(document.get(key), list):
            raise TransmissionPlotError(f"{path}: invalid {key}")
    minimum = document.get("minimum_paired_cpgs", 100)
    if not isinstance(minimum, int) or minimum < 1:
        raise TransmissionPlotError(f"{path}: invalid minimum_paired_cpgs")
    document["minimum_paired_cpgs"] = minimum
    return document


def write_qc_page(
    output: Path, payload: dict[str, Any], summary_path: Path | None
) -> dict[str, Any]:
    summary = read_summary(summary_path)
    samples = set(payload["samples"])
    chromosomes = set(payload["chromosomes"])
    for row in summary["comparisons"]:
        if not isinstance(row, dict):
            raise TransmissionPlotError("transmission QC comparison must be an object")
        if row.get("parent_id") not in samples or row.get("child_id") not in samples:
            raise TransmissionPlotError("transmission QC comparison contains unknown sample")
        if row.get("chromosome") not in chromosomes:
            raise TransmissionPlotError("transmission QC comparison contains unknown chromosome")
    document = {
        **summary,
        "overview": compute_overview(summary, list(payload["chromosomes"])),
    }
    (output / "transmission-qc.html").write_text(HTML_TEMPLATE, encoding="utf-8")
    (output / "transmission-qc.js").write_text(
        "window.TAPESTRY_TRANSMISSION_QC="
        + json.dumps(document, separators=(",", ":")).replace("</", "<\\/")
        + ";\n",
        encoding="utf-8",
    )
    return document
