// Therapeutic targets: what distinguishes these tumour cells, what can be reached,
// and which mechanism the biology supports. Research hypotheses, never advice.
(function () {
  const $ = (s) => document.querySelector(s);
  const pct = (v) => (v === null || v === undefined ? 'n/a' : Math.round(v * 100) + '%');
  const num = (v) => (v === null || v === undefined ? 'unknown' : Number(v).toFixed(2));
  const tri = (v) => (v === null || v === undefined ? 'unknown' : v ? 'yes' : 'no');
  const esc = (s) => String(s == null ? '' : s).replace(/[<>&]/g, (c) => ({ '<': '&lt;', '>': '&gt;', '&': '&amp;' }[c]));

  const CLASS_LABEL = {
    direct_surface: 'Direct surface target',
    neoantigen_hla: 'Neoantigen / HLA target',
    pathway_induced_surface: 'Indirect / pathway-induced',
    secreted: 'Secreted',
    intracellular_only: 'Intracellular only',
    unsuitable: 'Unsuitable',
    unknown: 'Unknown',
  };
  const RISK_COLOUR = { low: 'var(--ok)', moderate: 'var(--c3)', high: 'var(--bad)', unknown: 'var(--muted)' };

  let last = null;

  function bar(v, label) {
    if (v === null || v === undefined) return `<div class="tbar"><span class="tbar-l">${esc(label)}</span><span class="muted">unknown</span></div>`;
    return `<div class="tbar"><span class="tbar-l">${esc(label)}</span><span class="tbar-t"><i style="width:${Math.round(v * 100)}%"></i></span><span class="tbar-v">${num(v)}</span></div>`;
  }

  function overview(r) {
    const by = {};
    r.therapeutic_candidates.forEach((c) => (by[c.target_class] = (by[c.target_class] || 0) + 1));
    const alt = r.therapeutic_candidates.reduce((n, c) => n + c.origin_variants.length, 0);
    const neo = r.therapeutic_candidates.filter((c) => c.neoantigen && c.neoantigen.peptides > 0).length;
    return `<div class="kpi">
      <div><b>${alt}</b><span>somatic alterations</span></div>
      <div><b>${by.direct_surface || 0}</b><span>direct surface targets</span></div>
      <div><b>${neo}</b><span>neoantigen candidates</span></div>
      <div><b>${by.pathway_induced_surface || 0}</b><span>indirect candidates</span></div>
      <div><b>${(by.intracellular_only || 0) + (by.unknown || 0) + (by.secreted || 0)}</b><span>unsuitable or unknown</span></div>
      <div><b>${r.data_level}/9</b><span>patient data level</span></div>
    </div>`;
  }

  function mechanismPanel(mechs) {
    if (!mechs || !mechs.length) return '<div class="muted">no mechanism passed its hard requirements</div>';
    return mechs
      .map((m) => {
        const flag = m.status === 'experimental' ? ' <span class="pill experimental">EXPERIMENTAL</span>' : '';
        const cargo = m.payload_required ? `cargo ${esc(m.cargo)}` : 'no payload';
        return `<div class="tbar"><span class="tbar-l">${esc(m.mechanism)}</span><span class="tbar-t"><i style="width:${Math.round(m.compatibility * 100)}%"></i></span><span class="tbar-v">${pct(m.compatibility)}</span>${flag}</div>
                <div class="muted" style="font-size:11.5px;margin:-4px 0 6px 0">${cargo} · effector ${esc(m.effector)}</div>`;
      })
      .join('');
  }

  function card(c, design) {
    const d = design || {};
    const risk = (d.normal_tissue || {}).on_target_off_tumour_risk || 'unknown';
    const loc = d.localization || {};
    const traf = d.trafficking || {};
    const ready = d.readiness || {};
    const s = c.scores || {};
    return `<div class="tcard">
      <div class="thead">
        <div><b class="mono">${esc(c.gene)}</b> <span class="muted">${esc(c.protein || '')}</span></div>
        <div><span class="pill">${esc(CLASS_LABEL[c.target_class] || c.target_class)}</span>
             <span class="pill" style="color:${RISK_COLOUR[risk]};border-color:${RISK_COLOUR[risk]}">normal-tissue risk ${esc(risk)}</span>
             <b style="font-variant-numeric:tabular-nums">${num(s.overall)}</b></div>
      </div>
      <div class="muted" style="font-size:12.5px">${esc(c.why)}</div>
      <div class="muted" style="font-size:12px;margin-top:2px">${esc(d.class_reason || '')}</div>
      <div class="grid cols-2" style="margin-top:8px;gap:14px">
        <div>
          ${bar(s.surface_accessibility, 'surface accessibility')}
          ${bar(s.tumour_selectivity, 'tumour specificity')}
          ${bar(s.normal_tissue_safety, 'normal-tissue safety')}
          ${bar(s.internalisation, 'internalisation')}
          ${bar(s.evidence_strength, 'evidence strength')}
          <div class="muted" style="font-size:11.5px">scored on ${pct(s.coverage)} of the dimensions · confidence ${num(s.confidence)}</div>
        </div>
        <div>${mechanismPanel(c.recommended_mechanisms)}
          <div class="muted" style="font-size:11.5px;margin-top:4px">Computational compatibility score, not a clinical response probability.</div>
        </div>
      </div>
      <details style="margin-top:8px"><summary>Localisation and trafficking</summary>
        <div style="font-size:12.5px;margin-top:6px">
          <div>compartment <b>${esc((loc.primary || 'unknown').replace(/_/g, ' '))}</b> · topology ${esc(loc.topology || 'unknown')} · orientation ${esc(loc.orientation || 'unknown')}</div>
          ${(loc.extracellular_regions || []).map((r) => `<div class="mono">extracellular ${r.start}-${r.end}</div>`).join('')}
          ${(loc.transmembrane_regions || []).map((r) => `<div class="mono">transmembrane ${r.start}-${r.end}</div>`).join('')}
          <div style="margin-top:4px">internalises ${tri(traf.internalises)} · endosomal ${tri(traf.endosomal)} · lysosomal ${tri(traf.lysosomal)} · shedding ${tri(traf.shedding)}</div>
          ${traf.reason ? `<div class="muted">${esc(traf.reason)}</div>` : ''}
        </div>
      </details>
      <details><summary>Tumour DNA origin</summary>
        <div style="font-size:12.5px;margin-top:6px">${(c.origin_variants || []).map((o) => `<div class="mono">${esc(o.chromosome)}:${o.position} ${esc(o.reference)}&gt;${esc(o.alternate)} ${esc(o.variant_type)} ${esc(o.protein_change || '')} ${o.driver_frequency ? '· driver ' + (o.driver_frequency * 100).toFixed(1) + '%' : ''} · clonality ${esc(o.clonality)}</div>`).join('') || '<div class="muted">not altered in this tumour</div>'}</div>
      </details>
      ${c.neoantigen ? `<details><summary>Neoantigen / HLA route</summary>
        <div style="font-size:12.5px;margin-top:6px">
          <div>novel peptide sequence: <b>${tri(c.neoantigen.novel_peptide_sequence)}</b></div>
          <div class="muted">${esc(c.neoantigen.reason)}</div>
          <div>peptides: ${c.neoantigen.peptides} · HLA: ${esc((c.neoantigen.hla_alleles || []).join(', ') || 'not supplied')} · presentation: ${esc(c.neoantigen.presentation)}</div>
        </div></details>` : ''}
      <details><summary>What a binder must NOT recognise (${(d.negative_targets || []).length})</summary>
        <div style="font-size:12.5px;margin-top:6px">${(d.negative_targets || []).map((n) => `<div>• <b>${esc(n.label)}</b> <span class="muted">${esc(n.reason)}</span></div>`).join('')}</div>
      </details>
      <details><summary>Evidence (${(c.evidence || []).length})</summary>
        <div style="font-size:12px;margin-top:6px">${(c.evidence || []).map((e) => `<div>[<span class="pill">${esc(e.level)}</span> ${esc(e.source)}] ${esc(e.claim)}</div>`).join('')}</div>
      </details>
      <details><summary>Design readiness ${ready.overall !== undefined ? num(ready.overall) : ''}</summary>
        <div style="font-size:12.5px;margin-top:6px">
          ${Object.entries(ready.components || {}).map(([k, v]) => `<div>${esc(k.replace(/_/g, ' '))}: ${v === null ? '<span class="muted">unknown</span>' : num(v)}</div>`).join('')}
          ${(ready.blocking_unknowns || []).map((b) => `<div class="muted">blocked by: ${esc(b)}</div>`).join('')}
        </div>
      </details>
      <details><summary>Limitations (${(c.limitations || []).length})</summary>
        <div style="font-size:12.5px;margin-top:6px">${(c.limitations || []).map((l) => `<div>• ${esc(l)}</div>`).join('')}</div>
      </details>
    </div>`;
  }

  async function run() {
    const status = $('#th-status');
    status.textContent = 'annotating variants, compiling proteins, querying expression and precedent…';
    $('#th-out').innerHTML = '';
    try {
      const r = await window.api('/api/therapeutics', {
        vcf: $('#th-vcf').value,
        hla: $('#th-hla').value.trim(),
        top: +$('#th-top').value || 8,
      });
      last = r;
      status.textContent = `${r.therapeutic_candidates.length} candidates · data level ${r.data_level}/9`;
      const cards = r.therapeutic_candidates.map((c) => card(c, (r.design || {})[c.gene])).join('');
      $('#th-out').innerHTML = `
        <div class="hint" style="margin-bottom:8px">${esc(r.disclaimer)}</div>
        ${overview(r)}
        ${cards}
        <div class="card" style="margin-top:10px">
          <h2><i class="ic">🧩</i>Combination target logic</h2>
          ${(r.combinations || []).length ? r.combinations.map((x) => `<div style="font-size:12.5px"><b class="mono">${esc(x.targets.join(' ' + x.operator + ' '))}</b> ${esc(x.rationale)}${x.normal_tissues_excluded.length ? `<div class="muted">healthy tissues excluded: ${esc(x.normal_tissues_excluded.join(', '))}</div>` : ''}</div>`).join('') : '<div class="muted">no combination assessed</div>'}
          <div class="hint">Bulk tissue RNA cannot show two genes are in the same cell; single-cell data is required before a combination is treated as real.</div>
        </div>
        <div class="card" style="margin-top:10px">
          <h2><i class="ic">📋</i>Required missing data</h2>
          ${(r.missing_data || []).map((m) => `<div style="font-size:12.5px">• <b>${esc(m.input)}</b> <span class="muted">${esc(m.why)}</span></div>`).join('')}
        </div>`;
    } catch (e) {
      status.textContent = e.message;
    }
  }

  window.therapeuticsInit = function (files) {
    const vopts = (files.vcfs || []).map((v) => `<option value="${v.path}">${v.path}</option>`).join('');
    $('#th-vcf').innerHTML = vopts;
    const t = (files.vcfs || []).find((v) => v.path.includes('cancer_tumour'));
    if (t) $('#th-vcf').value = t.path;
    $('#th-run').onclick = run;
    $('#th-report').onclick = () => {
      $('#th-reportout').textContent = last ? last.report : 'run an analysis first';
    };
    $('#th-json').onclick = () => {
      if (!last) return;
      const compact = { sample: last.sample, data_level: last.data_level, therapeutic_candidates: last.therapeutic_candidates, combinations: last.combinations, missing_data: last.missing_data };
      $('#th-reportout').textContent = JSON.stringify(compact, null, 1);
    };
  };
})();
