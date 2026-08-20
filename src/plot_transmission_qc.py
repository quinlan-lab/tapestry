#!/usr/bin/env python3
"""Write the chromosome-level transmitted-haplotype methylation QC page."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class TransmissionPlotError(ValueError):
    """Raised when a transmission QC summary cannot be plotted safely."""


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
.lede { margin: 0 0 14px; max-width: 1050px; line-height: 1.45; }
.subtle { margin: 10px 0; font-size: 13px; line-height: 1.45; }
.nav { display: flex; gap: 8px; margin: 0 0 18px; }
.nav a { color: #334e91; background: white; border: 1px solid #ccd5e5; border-radius: 6px; padding: 6px 10px; text-decoration: none; font-size: 13px; }
.controls { display: flex; flex-wrap: wrap; gap: 14px; align-items: end; margin-bottom: 14px; }
label { display: grid; gap: 5px; color: #445067; font-size: 13px; font-weight: 650; }
select { min-width: 245px; padding: 7px 30px 7px 9px; border: 1px solid #cbd2df; border-radius: 6px; background: white; }
.badge { display: inline-block; padding: 2px 8px; border-radius: 999px; background: #edf1f8; font-size: 12px; }
.card { background: white; border: 1px solid #dde2eb; border-radius: 10px; box-shadow: 0 1px 4px rgb(30 45 70 / 8%); }
#heatmap { min-height: 420px; }
#missing { min-height: 250px; margin-bottom: 18px; }
#empty { display: none; padding: 24px; line-height: 1.5; }
</style>
</head>
<body>
<main class="page">
  <h1>Parent–child methylation concordance</h1>
  <p class="lede">Each cell compares model-based methylation at the same callable CpGs on the haplotype transmitted from the named parent to the child. Mismatch-proximal CpGs are excluded. This measures methylation observed on an inherited haplotype; it does not by itself establish epigenetic inheritance.</p>
  <nav class="nav"><a href="index.html">Haplotype ancestry</a><a href="transmission-qc.html" aria-current="page">Transmission QC</a></nav>
  <div class="controls">
    <label>Parent–child pair<select id="pair"></select></label>
    <label>Cell metric<select id="metric"></select></label>
    <span id="threshold" class="badge"></span>
    <span id="minimum" class="badge"></span>
    <span id="coverage" class="badge"></span>
  </div>
  <div id="missing" class="card"></div>
  <p class="subtle">Missing comparisons count pedigree transmission edges with no chromosome-level estimate for the selected metric. The chart separates absent sample output from edges whose samples have output but fewer than the configured minimum contributing CpGs on that chromosome.</p>
  <div id="heatmap" class="card"></div>
  <div id="empty" class="card"></div>
  <p class="subtle">Agreement is 1 − mean absolute parent–child difference. Inherited specificity is the non-transmitted parental haplotype MAE minus transmitted-haplotype MAE at CpGs callable in all three measurements; positive values mean the child is closer to the transmitted haplotype. Click a populated cell to open that child and chromosome in the ancestry view.</p>
</main>
<script>
const DATA = window.TAPESTRY_METADATA;
const QC = window.TAPESTRY_TRANSMISSION_QC;
const heatmap = document.getElementById('heatmap');
const missingPlot = document.getElementById('missing');
const empty = document.getElementById('empty');
const pairSelect = document.getElementById('pair');
const metricSelect = document.getElementById('metric');
const metrics = {
  agreement: {label: 'Inherited-haplotype agreement', range: [0, 1], colorscale: [[0,'#b2182b'],[0.5,'#f7f7f7'],[1,'#1a9850']]},
  mean_difference: {label: 'Signed child − parent difference', range: [-1, 1], colorscale: [[0,'#2166ac'],[0.5,'#f7f7f7'],[1,'#b2182b']]},
  discordant_fraction: {label: 'Large-discordance fraction', range: [0, 1], colorscale: [[0,'#fff7ec'],[1,'#b30000']]},
  callable_fraction: {label: 'Callable fraction', range: [0, 1], colorscale: [[0,'#f7fbff'],[1,'#08519c']]},
  inherited_specificity: {label: 'Inherited specificity', range: [-1, 1], colorscale: [[0,'#b2182b'],[0.5,'#f7f7f7'],[1,'#2166ac']]}
};
for (const [key, value] of Object.entries(metrics)) {
  const option = document.createElement('option'); option.value = key; option.textContent = value.label;
  metricSelect.appendChild(option);
}
document.getElementById('threshold').textContent = `large difference ≥ ${QC.discordance_threshold.toFixed(2)}`;
document.getElementById('minimum').textContent = `minimum ${QC.minimum_paired_cpgs.toLocaleString()} contributing CpGs`;
document.getElementById('coverage').textContent = `${QC.samples.length} sample${QC.samples.length === 1 ? '' : 's'} with methylation output`;
function generation(sample) { return DATA.people[sample] ? `F${DATA.people[sample].generation}` : '' ; }
function edgeKey(row) { return `${row.parent_id}\u0000${row.child_id}\u0000${row.relationship}`; }
function edgeLabel(row) { return `${generation(row.parent_id)} ${row.parent_id} → ${generation(row.child_id)} ${row.child_id} · ${row.relationship}`; }
function fmt(value, signed=false) {
  if (value === null || value === undefined) return 'not available';
  const text = Number(value).toFixed(3); return signed && value > 0 ? `+${text}` : text;
}
const edgeDefinitions = new Map(QC.edges.map(edge => [edgeKey(edge), edge]));
for (const row of QC.comparisons) if (!edgeDefinitions.has(edgeKey(row))) edgeDefinitions.set(edgeKey(row), row);
const allEdges = [...edgeDefinitions.keys()].sort((left, right) => edgeLabel(edgeDefinitions.get(left)).localeCompare(edgeLabel(edgeDefinitions.get(right)), undefined, {numeric: true}));
const allOption = document.createElement('option'); allOption.value = ''; allOption.textContent = `All pairs (${allEdges.length})`; pairSelect.appendChild(allOption);
for (const [index, edge] of allEdges.entries()) {
  const option = document.createElement('option'); option.value = String(index); option.textContent = edgeLabel(edgeDefinitions.get(edge)); pairSelect.appendChild(option);
}
const records = new Map(QC.comparisons.map(row => [`${edgeKey(row)}\u0000${row.chromosome}`, row]));
function selectedEdges() { return pairSelect.value === '' ? allEdges : [allEdges[Number(pairSelect.value)]]; }
function meetsMetricMinimum(row, metric) {
  return metric === 'inherited_specificity'
    ? row.specificity_cpgs >= QC.minimum_paired_cpgs
    : row.paired_cpgs >= QC.minimum_paired_cpgs;
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
function render() {
  const selected = new Set(selectedEdges());
  const metric = metricSelect.value || 'agreement';
  const visibleComparisons = QC.comparisons.filter(row => selected.has(edgeKey(row)));
  const usableComparisons = visibleComparisons.filter(row => meetsMetricMinimum(row, metric));
  const edges = selectedEdges().filter(edge => usableComparisons.some(row => edgeKey(row) === edge));
  if (!usableComparisons.length) {
    heatmap.style.display = 'none'; empty.style.display = 'block';
    empty.innerHTML = `<strong>No chromosome-level comparison is available for the selected parent–child pair and metric.</strong><br>Both samples need downstream founder-phased methylation output and at least the configured minimum shared evaluable CpGs. The missing-comparison chart above identifies the reason.`;
    return;
  }
  heatmap.style.display = 'block'; empty.style.display = 'none';
  const definition = metrics[metric];
  const z = [], customdata = [];
  for (const edge of edges) {
    const zRow = [], customRow = [];
    for (const chromosome of DATA.chromosomes) {
      const row = records.get(`${edge}\u0000${chromosome}`);
      zRow.push(row && meetsMetricMinimum(row, metric) && row[metric] !== null ? row[metric] : null);
      customRow.push(row ? [
        edgeLabel(row), chromosome, fmt(row.agreement), fmt(row.mean_difference, true),
        fmt(row.discordant_fraction), row.paired_cpgs.toLocaleString(),
        row.evaluated_cpgs.toLocaleString(), fmt(row.callable_fraction),
        row.mismatch_excluded_cpgs.toLocaleString(), fmt(row.inherited_specificity, true),
        row.specificity_cpgs.toLocaleString(), row.child_id
      ] : [edgeLabel(edgeDefinitions.get(edge)), chromosome, 'not available', 'not available',
           'not available', '0', '0', 'not available', '0', 'not available', '0', edgeDefinitions.get(edge).child_id]);
    }
    z.push(zRow); customdata.push(customRow);
  }
  const height = Math.max(420, edges.length * 42 + 190); heatmap.style.height = `${height}px`;
  Plotly.newPlot(heatmap, [{
    type: 'heatmap', x: DATA.chromosomes, y: edges.map(edge => edgeLabel(edgeDefinitions.get(edge))),
    z, customdata, zmin: definition.range[0], zmax: definition.range[1], zmid: definition.range[0] < 0 ? 0 : undefined,
    colorscale: definition.colorscale, hoverongaps: false,
    colorbar: {title: {text: definition.label, side: 'right'}},
    hovertemplate: '<b>%{customdata[0]}</b><br>%{customdata[1]}<br>Agreement: %{customdata[2]}<br>Signed child − parent: %{customdata[3]}<br>Large-discordance fraction: %{customdata[4]}<br>Paired CpGs: %{customdata[5]} / %{customdata[6]} evaluated<br>Callable fraction: %{customdata[7]}<br>Mismatch-proximal excluded: %{customdata[8]}<br>Inherited specificity: %{customdata[9]} (%{customdata[10]} CpGs)<extra></extra>'
  }], {
    template: 'plotly_white', height, margin: {l: 300, r: 95, t: 60, b: 65},
    title: {text: definition.label, x: 0.01, xanchor: 'left'},
    xaxis: {title: 'Chromosome', side: 'bottom', fixedrange: true},
    yaxis: {title: 'Transmitted haplotype', automargin: true, fixedrange: true},
    hoverlabel: {align: 'left'}, uirevision: metric
  }, {responsive: true, displaylogo: false}).then(() => {
    heatmap.removeAllListeners?.('plotly_click');
    heatmap.on('plotly_click', event => {
      const point = event.points[0]; if (point?.z === null || !point?.customdata) return;
      const chromosome = point.customdata[1], child = point.customdata[11];
      location.href = `index.html?sample=${encodeURIComponent(child)}&chrom=${encodeURIComponent(chromosome)}`;
    });
  });
}
metricSelect.value = 'agreement';
metricSelect.addEventListener('change', () => { renderMissingEdges(); render(); });
pairSelect.addEventListener('change', () => { renderMissingEdges(); render(); });
renderMissingEdges();
render();
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
    (output / "transmission-qc.html").write_text(HTML_TEMPLATE, encoding="utf-8")
    (output / "transmission-qc.js").write_text(
        "window.TAPESTRY_TRANSMISSION_QC="
        + json.dumps(summary, separators=(",", ":")).replace("</", "<\\/")
        + ";\n",
        encoding="utf-8",
    )
    return summary
