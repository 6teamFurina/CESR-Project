#!/usr/bin/env python3
"""Build the inline interactive pair of sparse affinity heatmaps."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
LATEST_RESULTS = HERE / "results" / "scibmad_latest"


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--scores",
        type=Path,
        default=LATEST_RESULTS / "affinity" / "quadrupole_affinity_scores.csv",
    )
    result.add_argument(
        "--screen",
        type=Path,
        default=LATEST_RESULTS / "responses" / "quadrupole_optics_screen.csv",
    )
    result.add_argument("--output", type=Path, required=True)
    return result


def quadrupole_key(name: str) -> tuple[int, str]:
    match = re.fullmatch(r"Q(\d+)([EW])", name.upper())
    return (int(match.group(1)), match.group(2)) if match else (10_000, name)


FRAGMENT = r'''
<div id="sext-quad-affinity-20260816">
  <h1>Sextupole–quadrupole affinity</h1>
  <div class="card aff-selection" aria-live="polite"></div>
  <section class="aff-panel" aria-labelledby="aff-info-heading">
    <h2 id="aff-info-heading">Nuisance-marginalized information gain</h2>
    <div class="aff-chart" data-metric="info"></div>
    <div class="aff-legend text-small" data-legend="info"></div>
  </section>
  <section class="aff-panel" aria-labelledby="aff-precision-heading">
    <h2 id="aff-precision-heading">Worst-axis precision improvement</h2>
    <div class="aff-chart" data-metric="precision"></div>
    <div class="aff-legend text-small" data-legend="precision"></div>
  </section>
  <div class="tooltip aff-tooltip" role="tooltip" hidden></div>
</div>
<script src="https://cdn.jsdelivr.net/npm/d3@7.9.0/dist/d3.min.js"></script>
<style>
  #sext-quad-affinity-20260816 { position: relative; color: var(--foreground); width: 100%; }
  #sext-quad-affinity-20260816 .aff-selection { margin-block: 12px 18px; }
  #sext-quad-affinity-20260816 .aff-selection .viz-row { justify-content: flex-start; gap: 18px; }
  #sext-quad-affinity-20260816 .aff-panel { margin-block: 18px 28px; }
  #sext-quad-affinity-20260816 .aff-chart { width: 100%; }
  #sext-quad-affinity-20260816 .aff-chart svg { display: block; width: 100%; height: auto; }
  #sext-quad-affinity-20260816 .aff-frame { fill: var(--muted); fill-opacity: 0.18; stroke: var(--border); }
  #sext-quad-affinity-20260816 .aff-grid-line { stroke: var(--border); stroke-opacity: 0.45; shape-rendering: crispEdges; }
  #sext-quad-affinity-20260816 .aff-cell { cursor: pointer; }
  #sext-quad-affinity-20260816 .aff-column-best { stroke: var(--foreground); stroke-width: 1; }
  #sext-quad-affinity-20260816 .aff-current { stroke: var(--primary); stroke-width: 2; }
  #sext-quad-affinity-20260816 .aff-axis text { fill: var(--foreground); font-size: 12px; }
  #sext-quad-affinity-20260816 .aff-axis path,
  #sext-quad-affinity-20260816 .aff-axis line { stroke: var(--border); }
  #sext-quad-affinity-20260816 text.axis-title { fill: var(--foreground); font-size: 12px; }
  #sext-quad-affinity-20260816 .aff-legend { display: flex; align-items: center; gap: 8px; color: var(--muted-foreground); }
  #sext-quad-affinity-20260816 .aff-gradient { width: min(240px, 46%); height: 10px; }
  #sext-quad-affinity-20260816 .aff-tooltip { position: absolute; z-index: 4; pointer-events: none; background: var(--popover); color: var(--popover-foreground); border: 1px solid var(--border); padding: 8px 10px; max-width: 280px; }
</style>
<script>
(() => {
  const root = document.getElementById("sext-quad-affinity-20260816");
  const payload = __DATA_JSON__;
  const pairs = payload.pairs.map(d => ({
    s: payload.sextupoles[d[0]], q: payload.quadrupoles[d[1]], si: d[0], qi: d[1],
    info: d[2], precision: d[3], rank: d[4], leverage: d[5], sx: d[6], sy: d[7]
  }));
  const lookup = new Map(pairs.map(d => [`${d.si}:${d.qi}`, d]));
  const tooltip = root.querySelector(".aff-tooltip");
  const selection = root.querySelector(".aff-selection");
  let selected = pairs.reduce((best, d) => d.info > best.info ? d : best, pairs[0]);

  const formatInfo = d3.format(".3f");
  const formatPrecision = d3.format(".3f");
  const formatLeverage = d3.format(".3g");

  function updateSelection(d) {
    selected = d;
    selection.innerHTML = `<div class="viz-row"><strong>${d.s} × ${d.q}</strong><span>Δ log det F: ${formatInfo(d.info)}</span><span>precision: ${formatPrecision(d.precision)}×</span><span>screen rank: ${d.rank}/${payload.topK}</span><span>σ: (${d.sx.toFixed(1)}, ${d.sy.toFixed(1)}) µm</span></div>`;
    root.querySelectorAll(".aff-cell").forEach(node => {
      const same = +node.dataset.si === d.si && +node.dataset.qi === d.qi;
      node.classList.toggle("aff-current", same);
    });
  }

  function showTooltip(event, d) {
    tooltip.hidden = false;
    tooltip.innerHTML = `<strong>${d.s} × ${d.q}</strong><br>Δ log det F: ${formatInfo(d.info)}<br>precision: ${formatPrecision(d.precision)}×<br>optics leverage: ${formatLeverage(d.leverage)}<br>screen rank: ${d.rank}/${payload.topK}`;
    const box = root.getBoundingClientRect();
    const left = Math.min(event.clientX - box.left + 12, box.width - tooltip.offsetWidth - 8);
    const top = event.clientY - box.top + 12;
    tooltip.style.left = `${Math.max(8, left)}px`;
    tooltip.style.top = `${Math.max(8, top)}px`;
  }

  function hideTooltip() { tooltip.hidden = true; }

  function draw(container) {
    const metric = container.dataset.metric;
    const width = Math.max(320, Math.floor(container.getBoundingClientRect().width));
    const height = width < 500 ? 720 : 800;
    const margin = {top: 8, right: 18, bottom: 92, left: width < 500 ? 88 : 86};
    const innerWidth = width - margin.left - margin.right;
    const innerHeight = height - margin.top - margin.bottom;
    const x = d3.scaleBand(payload.sextupoles, [0, innerWidth]);
    const y = d3.scaleBand(payload.quadrupoles, [0, innerHeight]);
    const extent = d3.extent(pairs, d => d[metric]);
    const color = d3.scaleSequential(extent, d3.interpolateViridis);
    container.replaceChildren();
    const svg = d3.select(container).append("svg")
      .attr("viewBox", `0 0 ${width} ${height}`)
      .attr("role", "img")
      .attr("aria-label", metric === "info" ? "Sparse heatmap of nuisance-marginalized information gain" : "Sparse heatmap of worst-axis precision improvement");
    svg.append("title").text(metric === "info" ? "Information gain heatmap" : "Precision improvement heatmap");
    svg.append("desc").text(`Columns are ${payload.sextupoles.length} target sextupoles in ring order. Rows are the ${payload.quadrupoles.length} quadrupoles retained for at least one target. Only ${payload.topK} pre-screened cells per sextupole are colored.`);
    const plot = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
    plot.append("rect").attr("data-chart-frame", "").attr("class", "aff-frame")
      .attr("width", innerWidth).attr("height", innerHeight);
    d3.range(10, payload.sextupoles.length, 10).forEach(i => plot.append("line")
      .attr("class", "aff-grid-line").attr("x1", i * x.bandwidth()).attr("x2", i * x.bandwidth())
      .attr("y1", 0).attr("y2", innerHeight));
    d3.range(10, payload.quadrupoles.length, 10).forEach(i => plot.append("line")
      .attr("class", "aff-grid-line").attr("x1", 0).attr("x2", innerWidth)
      .attr("y1", i * y.bandwidth()).attr("y2", i * y.bandwidth()));
    const bestBySextupole = d3.rollup(pairs, group => d3.max(group, d => d[metric]), d => d.si);
    plot.selectAll("rect.aff-cell").data(pairs).join("rect")
      .attr("class", d => `aff-cell${d[metric] === bestBySextupole.get(d.si) ? " aff-column-best" : ""}${d.si === selected.si && d.qi === selected.qi ? " aff-current" : ""}`)
      .attr("data-si", d => d.si).attr("data-qi", d => d.qi)
      .attr("x", d => x(d.s)).attr("y", d => y(d.q))
      .attr("width", Math.max(0.6, x.bandwidth())).attr("height", Math.max(0.6, y.bandwidth()))
      .attr("fill", d => color(d[metric]))
      .attr("aria-label", d => `${d.s}, ${d.q}, information gain ${formatInfo(d.info)}, precision ${formatPrecision(d.precision)} times`)
      .on("pointermove", showTooltip).on("pointerleave", hideTooltip)
      .on("click", (event, d) => updateSelection(d));

    const xEvery = width < 760 ? 10 : 6;
    const yEvery = width < 500 ? 14 : width < 760 ? 10 : 7;
    const finalSextupole = payload.sextupoles.length - 1;
    let xTickIndices;
    if (width < 500) {
      xTickIndices = d3.range(4).map(i => Math.round(i * finalSextupole / 3));
    } else {
      xTickIndices = d3.range(0, payload.sextupoles.length, xEvery);
      if (finalSextupole - xTickIndices[xTickIndices.length - 1] >= 0.75 * xEvery) {
        xTickIndices.push(finalSextupole);
      }
    }
    const xTicks = xTickIndices.map(i => payload.sextupoles[i]);
    const yTicks = payload.quadrupoles.filter((_, i) => i % yEvery === 0 || i === payload.quadrupoles.length - 1);
    svg.append("g").attr("class", "aff-axis")
      .attr("transform", `translate(${margin.left},${margin.top + innerHeight})`)
      .call(d3.axisBottom(x).tickValues(xTicks).tickSizeOuter(0))
      .selectAll("text").attr("transform", "rotate(-55)").attr("text-anchor", "end");
    svg.append("g").attr("class", "aff-axis")
      .attr("transform", `translate(${margin.left},${margin.top})`)
      .call(d3.axisLeft(y).tickValues(yTicks).tickSizeOuter(0));
    svg.append("text").attr("class", "axis-title").attr("data-axis", "x")
      .attr("x", margin.left + innerWidth / 2).attr("y", height - 6).attr("text-anchor", "middle")
      .text(width < 500 ? "Target sextupoles" : "Target sextupole — blank means not evaluated");
    svg.append("text").attr("class", "axis-title").attr("data-axis", "y")
      .attr("transform", `translate(13,${margin.top + innerHeight / 2}) rotate(-90)`)
      .attr("text-anchor", "middle").text("Quadrupole");

    const legend = root.querySelector(`[data-legend="${metric}"]`);
    const label = metric === "info" ? "Δ log det F" : "precision ratio";
    legend.innerHTML = `<span>${label}</span><span>${metric === "info" ? formatInfo(extent[0]) : formatPrecision(extent[0])}</span><span class="aff-gradient" aria-hidden="true"></span><span>${metric === "info" ? formatInfo(extent[1]) : formatPrecision(extent[1])}</span>`;
    const gradientStops = d3.range(7).map(i => d3.interpolateViridis(i / 6));
    legend.querySelector(".aff-gradient").style.background = `linear-gradient(90deg, ${gradientStops.join(", ")})`;
  }

  function drawAll() { root.querySelectorAll(".aff-chart").forEach(draw); updateSelection(selected); }
  const observer = new ResizeObserver(() => drawAll());
  root.querySelectorAll(".aff-chart").forEach(node => observer.observe(node));
  new MutationObserver(drawAll).observe(document.documentElement, {attributes: true, attributeFilter: ["class", "style", "data-theme"]});
  drawAll();
})();
</script>
'''.strip()


def main() -> int:
    args = parser().parse_args()
    with args.scores.expanduser().resolve().open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    with args.screen.expanduser().resolve().open(encoding="utf-8", newline="") as stream:
        screen_rows = list(csv.DictReader(stream))
    sextupole_s = {row["sextupole"]: float(row["sextupole_s_m"]) for row in rows}
    sextupoles = sorted(sextupole_s, key=lambda name: sextupole_s[name])
    retained_quadrupoles = {row["quadrupole"] for row in rows}
    if screen_rows and all(row.get("quadrupole_s_m", "").strip() for row in screen_rows):
        quadrupole_s = {
            row["quadrupole"]: float(row["quadrupole_s_m"])
            for row in screen_rows
            if row["quadrupole"] in retained_quadrupoles
        }
        if not all(math.isfinite(value) for value in quadrupole_s.values()):
            raise ValueError("Non-finite quadrupole s coordinate")
        quadrupoles = sorted(quadrupole_s, key=lambda name: quadrupole_s[name])
    else:
        quadrupoles = sorted(retained_quadrupoles, key=quadrupole_key)
    s_index = {name: index for index, name in enumerate(sextupoles)}
    q_index = {name: index for index, name in enumerate(quadrupoles)}
    pairs = [
        [
            s_index[row["sextupole"]],
            q_index[row["quadrupole"]],
            round(float(row["information_gain_logdet"]), 6),
            round(float(row["precision_improvement_worst_axis"]), 6),
            int(row["screen_rank"]),
            round(float(row["optics_leverage"]), 7),
            round(float(row["candidate_sigma_x_um"]), 2),
            round(float(row["candidate_sigma_y_um"]), 2),
        ]
        for row in rows
    ]
    payload = {
        "sextupoles": sextupoles,
        "quadrupoles": quadrupoles,
        "topK": max(int(row["screen_rank"]) for row in rows),
        "pairs": pairs,
    }
    fragment = FRAGMENT.replace(
        "__DATA_JSON__", json.dumps(payload, ensure_ascii=True, separators=(",", ":"))
    )
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(fragment + "\n", encoding="utf-8")
    print(f"Wrote {output} ({output.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
