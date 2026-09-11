// GenomeOS molecules view: a gene's protein from UniProt and its AlphaFold structure,
// drawn as a rotating C-alpha backbone coloured by per-residue confidence (pLDDT).
(function () {
  const $ = s => document.querySelector(s);
  const canvas = $('#m-canvas'); if (!canvas) return;
  const ctx = canvas.getContext('2d');
  let ca = [], plddt = [], angle = 0, spin = true, highlight = null, seq = '';
  const colour = p => p >= 90 ? '#1f6feb' : p >= 70 ? '#39c5c5' : p >= 50 ? '#e3b341' : '#ff6b66';
  function draw() {
    const r = canvas.parentElement.getBoundingClientRect(); if (canvas.width !== Math.floor(r.width)) { canvas.width = Math.floor(r.width); canvas.height = 420; }
    const W = canvas.width, H = canvas.height; ctx.clearRect(0, 0, W, H);
    if (!ca.length) { ctx.fillStyle = '#888'; ctx.font = '13px system-ui'; ctx.fillText('no structure loaded', 16, 24); return; }
    const cx = ca.reduce((a, p) => a + p[0], 0) / ca.length, cy = ca.reduce((a, p) => a + p[1], 0) / ca.length, cz = ca.reduce((a, p) => a + p[2], 0) / ca.length;
    const cos = Math.cos(angle), sin = Math.sin(angle), tilt = 0.35;
    const pts = ca.map(([x, y, z], i) => { const dx = x - cx, dy = y - cy, dz = z - cz; const rx = dx * cos - dz * sin, rz = dx * sin + dz * cos; const ry = dy * Math.cos(tilt) - rz * Math.sin(tilt), rz2 = dy * Math.sin(tilt) + rz * Math.cos(tilt); return {x: rx, y: ry, z: rz2, i}; });
    const ext = Math.max(...pts.map(p => Math.max(Math.abs(p.x), Math.abs(p.y)))) || 1;
    const scale = Math.min(W, H) * 0.42 / ext;
    const X = p => W / 2 + p.x * scale, Y = p => H / 2 - p.y * scale;
    for (let i = 1; i < pts.length; i++) { const a = pts[i - 1], b = pts[i]; const depth = (b.z / ext + 1) / 2; ctx.strokeStyle = colour(plddt[i]); ctx.globalAlpha = 0.35 + 0.65 * depth; ctx.lineWidth = 1.2 + 2.2 * depth; ctx.beginPath(); ctx.moveTo(X(a), Y(a)); ctx.lineTo(X(b), Y(b)); ctx.stroke(); }
    ctx.globalAlpha = 1;
    if (highlight != null && pts[highlight - 1]) { const p = pts[highlight - 1]; ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(X(p), Y(p), 7, 0, Math.PI * 2); ctx.stroke(); ctx.fillStyle = '#fff'; ctx.font = 'bold 12px system-ui'; ctx.fillText(`${seq[highlight - 1] || ''}${highlight} (pLDDT ${plddt[highlight - 1]})`, X(p) + 10, Y(p) - 8); }
    ctx.fillStyle = '#8b949e'; ctx.font = '11px system-ui'; ctx.fillText('pLDDT: ', 10, H - 10);
    [['#1f6feb', '>90 high'], ['#39c5c5', '70–90 confident'], ['#e3b341', '50–70 low'], ['#ff6b66', '<50 very low']].forEach(([c, l], k) => { ctx.fillStyle = c; ctx.fillRect(60 + k * 120, H - 18, 10, 10); ctx.fillStyle = '#8b949e'; ctx.fillText(l, 74 + k * 120, H - 10); });
  }
  function tick() { if (spin) angle += 0.012; draw(); if ($('#molecules').classList.contains('active')) requestAnimationFrame(tick); else setTimeout(() => requestAnimationFrame(tick), 500); }
  requestAnimationFrame(tick);
  canvas.onclick = () => { spin = !spin; };
  function pill(sec) { const k = sec && sec.evidence || 'none'; return `<span class="pill ${k === 'curated' ? 'curated' : k === 'predicted' ? 'predicted' : k === 'experimental' ? 'measured' : 'none'}" data-tip="${(sec && sec.source || '').replace(/"/g, '')} · confidence ${sec ? sec.confidence : 0}">${k}</span>`; }
  function renderDefinition(d) {
    const s = d.sections || {}; const id = (s.identity || {}).items;
    if (!id) { $('#m-def').innerHTML = `<span class="muted">no reviewed UniProt entry${s.identity && s.identity.error ? ': ' + s.identity.error : ''}</span>`; return; }
    const go = (s.genomic_origin || {}).items; const iso = (s.isoforms || {}).items || []; const fn = (s.function || {}).items || {};
    const dom = (s.domains || {}).items || {}; const exp = s.structures_experimental || {}; const pred = s.structures_predicted || {};
    const pw = (s.pathways || {}).items || []; const it = (s.interactions || {}); const its = it.items || []; const ex = (s.expression || {}).items || {};
    const dis = (s.diseases || {}).items || []; const mods = (s.modifications || {}).items || [];
    const methods = {}; for (const x of exp.items || []) methods[x.method] = (methods[x.method] || 0) + 1;
    const cov = d.coverage || {}; const covRow = Object.entries(cov).map(([k, v]) => `<span class="chip" style="opacity:${v ? 1 : .45}" data-tip="${k.replace(/_/g, ' ')}">${v ? '✓' : '·'} ${k.replace(/_/g, ' ')}</span>`).join(' ');
    const canon = go ? (go.transcripts.find(t => t.canonical) || {}) : {};
    const ntpm = ex.tissue_ntpm ? Object.entries(ex.tissue_ntpm).sort((a, b) => b[1] - a[1]).slice(0, 8) : [];
    const sec = (title, body, p) => `<div style="margin-top:8px"><b>${title}</b> ${p}<div style="font-size:12.5px;margin-top:2px">${body}</div></div>`;
    $('#m-def').innerHTML = `
      <div style="font-size:12px;margin-bottom:6px">${covRow}</div>
      <div><b>${d.id}</b> ${id.name} · ${id.length} aa · ${id.existence || ''} ${pill(s.identity)}</div>
      ${go ? sec('Genomic origin', `${go.gene_id} ${go.locus} · ${go.transcripts.length} transcripts → ${go.protein_products} protein products; canonical ${canon.name || '?'} (${canon.protein_length || '?'} aa). Gene → transcripts → isoforms, never gene → protein.`, pill(s.genomic_origin)) : sec('Genomic origin', `<span class="muted">${(s.genomic_origin || {}).error || 'not fetched'}</span>`, pill(s.genomic_origin))}
      ${sec('Isoforms (UniProt)', iso.length ? iso.map(i => `<span class="mono">${i.id}</span>${i.name ? ' ' + i.name : ''}${i.status === 'Displayed' ? ' <span class="muted">(canonical sequence)</span>' : ''}`).join(' · ') : '<span class="muted">one form recorded</span>', pill(s.isoforms))}
      ${sec('Function', `${(fn.summary || [])[0] ? fn.summary[0].slice(0, 400) + (fn.summary[0].length > 400 ? '…' : '') : '<span class="muted">not characterised</span>'}${fn.location && fn.location.length ? `<div class="muted">location: ${fn.location.slice(0, 5).join('; ')}</div>` : ''}`, pill(s.function))}
      ${sec('Domains', `${(dom.interpro || []).map(x => `<span class="chip" data-tip="${x.id}">${x.name || x.id}</span>`).join(' ') || '<span class="muted">none</span>'}${dom.features && dom.features.length ? `<div class="muted">${dom.features.length} annotated features (domains, regions, sites)</div>` : ''}`, pill(s.domains))}
      ${sec('Modifications', mods.length ? `${mods.length} sites: ${Object.entries(mods.reduce((a, f) => { a[f.type] = (a[f.type] || 0) + 1; return a; }, {})).map(([k, v]) => `${k} ${v}`).join(', ')}` : '<span class="muted">none recorded</span>', pill(s.modifications))}
      ${sec('Structures', `experimental <b>${exp.count || 0}</b> ${Object.keys(methods).length ? '(' + Object.entries(methods).map(([k, v]) => `${k} ${v}`).join(', ') + ')' : ''} ${pill(exp)} · predicted <b>${(pred.items || []).length}</b> ${pill(pred)} <span class="muted">— predicted is never treated as observed</span>`, '')}
      ${sec('Pathways (Reactome)', pw.length ? pw.slice(0, 12).map(x => `<span class="chip" data-pw="${x.id}" data-tip="${x.id} · click to run with ${d.gene} knocked out">${x.name}</span>`).join(' ') + (pw.length > 12 ? ` <span class="muted">+${pw.length - 12} more</span>` : '') : '<span class="muted">none</span>', pill(s.pathways))}
      ${sec('Associations (STRING ≥ 0.7)', its.length ? its.slice(0, 20).map(x => `<span class="chip" data-tip="combined ${x.score}; experimental ${x.experimental}, database ${x.database}, text-mining ${x.textmining}${x.physical_evidence ? ' · experimental support' : ''}" style="${x.physical_evidence ? '' : 'opacity:.6'}">${x.partner}</span>`).join(' ') + ` <span class="muted">${its.length} total, ${its.filter(x => x.physical_evidence).length} with experimental support; dimmed = association without physical evidence</span>` : `<span class="muted">${it.error || 'none'}</span>`, pill(it))}
      ${sec('Expression (Human Protein Atlas)', ex.tissue_specificity ? `${ex.tissue_specificity}; ${ex.tissue_distribution || ''}; main location ${(ex.subcellular_main || []).join(', ') || '?'}; cell types: ${ex.cell_type_specificity || '?'}${ntpm.length ? `<table style="margin-top:4px">${ntpm.map(([k, v]) => `<tr><td>${k}</td><td class="num mono">${v.toFixed(0)} nTPM</td></tr>`).join('')}</table>` : ''}` : `<span class="muted">${(s.expression || {}).error || 'not fetched'}</span>`, pill(s.expression))}
      ${dis.length ? sec('Diseases (UniProt)', dis.map(x => `<span class="chip" data-tip="${(x.description || '').replace(/"/g, '')}">${x.name}${x.mim ? ' · MIM ' + x.mim : ''}</span>`).join(' '), pill(s.diseases)) : ''}
      ${d.states && d.states.length ? `<div class="muted" style="margin-top:8px;font-size:12px">${d.states.length} ProteinState records derived (protein × tissue/cell type × level × location); the definition above is what the protein <i>is</i>, a state is where it is and how much.</div>` : ''}`;
  }
  window.moleculesPathway = async function (id, ko) {
    $('#m-pw-status').textContent = 'fetching Reactome export and running…';
    try {
      const d = await window.api(`/api/pathway?id=${encodeURIComponent(id)}${ko ? '&knockout=' + encodeURIComponent(ko) : ''}`);
      const s = d.summary, k = d.knockout;
      $('#m-pw-status').textContent = `${s.name} · Reactome v${s.version}`;
      const lostIds = new Set(k ? k.reactions_lost.map(r => r.id) : []);
      $('#m-pw-out').innerHTML = `
        <div><b>${s.pathway}</b> ${s.name} <span class="pill curated">curated</span> · ${s.species} entities (${s.proteins} proteins), ${s.reactions} reactions, ${s.reactions_reachable_from_sources} reachable from the pathway's inputs, ${s.catalysed} catalysed, ${s.inhibited} inhibited</div>
        ${k ? `<div style="margin-top:6px"><b>knockout ${k.symbol || k.knockout}</b> <span class="pill inferred">inferred</span> <span class="muted">${k.logic} · confidence ${k.confidence}</span><br>${k.entities_containing.length} entities contain it · <b>${k.reactions_lost.length} of ${k.reactions_reachable_baseline} reachable reactions lost (${(k.fraction_lost * 100).toFixed(0)}%)</b> · ${k.products_unreachable.length} products can no longer be made${k.products_unreachable.length ? ': <span class="muted">' + k.products_unreachable.slice(0, 8).join('; ') + (k.products_unreachable.length > 8 ? '…' : '') + '</span>' : ''}</div>` : ''}
        <table style="margin-top:8px;font-size:12px"><tr><th></th><th>reaction</th><th>inputs</th><th>outputs</th><th>catalyst / inhibitor</th></tr>
        ${d.reactions.map(r => `<tr style="${lostIds.has(r.id) ? 'color:var(--bad)' : ''}"><td>${lostIds.has(r.id) ? '✗' : '✓'}</td><td>${r.name}</td><td class="muted">${r.inputs.join(' + ')}</td><td class="muted">${r.outputs.join(' + ')}</td><td class="muted">${r.catalysts.map(c => '⚙ ' + c).concat(r.inhibitors.map(i => '⊣ ' + i)).join('; ')}</td></tr>`).join('')}</table>`;
    } catch (e) { $('#m-pw-status').textContent = e.message; }
  };
  // ---- knowledge-graph neighbourhood: a small force layout on a canvas
  const G = {nodes: [], edges: [], drag: null, raf: 0, ticks: 0};
  const KIND_COL = {protein: '#1f6feb', pathway: '#1a7f37', domain: '#8250df', tissue: '#bf8700'};
  function gCanvas() { return $('#m-graph'); }
  function gDraw() {
    const c = gCanvas(); if (!c) return; const dpr = window.devicePixelRatio || 1; const W = c.clientWidth, H = c.clientHeight; if (!W) return;
    c.width = W * dpr; c.height = H * dpr; const ctx = c.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
    const byId = Object.fromEntries(G.nodes.map(n => [n.id, n]));
    for (const e of G.edges) { const a = byId[e.a], b = byId[e.b]; if (!a || !b) continue; ctx.strokeStyle = e.rel === 'associates' ? (e.physical ? '#1f6feb' : '#8b949e') : KIND_COL[b.kind] || '#999'; ctx.globalAlpha = 0.25 + e.confidence * 0.5; ctx.lineWidth = e.rel === 'associates' ? 0.6 + (e.score || 0.7) : 1; ctx.setLineDash(e.rel === 'associates' && !e.physical ? [3, 3] : []); ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke(); }
    ctx.setLineDash([]); ctx.globalAlpha = 1; ctx.font = '11px system-ui'; ctx.textAlign = 'center';
    const ink = getComputedStyle(document.documentElement).getPropertyValue('--ink').trim() || '#222';
    for (const n of G.nodes) { const r = n.centre ? 9 : n.kind === 'protein' ? 6 : 5; ctx.fillStyle = KIND_COL[n.kind] || '#999'; ctx.globalAlpha = n.compiled === false ? 0.55 : 1; ctx.beginPath(); ctx.arc(n.x, n.y, r, 0, Math.PI * 2); ctx.fill(); ctx.globalAlpha = 1; if (n.centre || n.kind === 'protein' || G.nodes.length <= 25) { ctx.fillStyle = ink; ctx.fillText(n.kind === 'tissue' ? n.name : (n.kind === 'pathway' || n.kind === 'domain') ? (n.name || n.id).slice(0, 22) : n.id, n.x, n.y - r - 3); } }
  }
  function gStep() {
    const c = gCanvas(); const W = c.clientWidth, H = c.clientHeight; const byId = Object.fromEntries(G.nodes.map(n => [n.id, n]));
    for (const n of G.nodes) { n.vx = (n.vx || 0) * 0.85; n.vy = (n.vy || 0) * 0.85; }
    for (let i = 0; i < G.nodes.length; i++) for (let k = i + 1; k < G.nodes.length; k++) { const a = G.nodes[i], b = G.nodes[k]; let dx = b.x - a.x, dy = b.y - a.y; const d2 = dx * dx + dy * dy + 0.01; const f = 900 / d2; const d = Math.sqrt(d2); dx /= d; dy /= d; a.vx -= dx * f; a.vy -= dy * f; b.vx += dx * f; b.vy += dy * f; }
    for (const e of G.edges) { const a = byId[e.a], b = byId[e.b]; if (!a || !b) continue; const dx = b.x - a.x, dy = b.y - a.y; const d = Math.sqrt(dx * dx + dy * dy) + 0.01; const want = e.rel === 'associates' ? 90 : 70; const f = (d - want) * 0.02; a.vx += dx / d * f; a.vy += dy / d * f; b.vx -= dx / d * f; b.vy -= dy / d * f; }
    for (const n of G.nodes) { if (n === G.drag) continue; n.vx += (W / 2 - n.x) * 0.002; n.vy += (H / 2 - n.y) * 0.002; n.x = Math.max(10, Math.min(W - 10, n.x + n.vx)); n.y = Math.max(14, Math.min(H - 10, n.y + n.vy)); }
    gDraw(); if (G.ticks++ < 300) G.raf = requestAnimationFrame(gStep);
  }
  window.moleculesGraph = async function (gene) {
    $('#m-graph-status').textContent = 'building from the local definitions…';
    try {
      const n = await window.api(`/api/graph?gene=${encodeURIComponent(gene)}&max=40`);
      if (!n.nodes.length) { $('#m-graph-status').textContent = 'not compiled locally yet'; G.nodes = []; G.edges = []; gDraw(); return; }
      const c = gCanvas(); const W = c.clientWidth || 800, H = c.clientHeight || 380;
      G.nodes = n.nodes.map((x, i) => ({...x, centre: x.id === n.centre, x: W / 2 + (x.id === n.centre ? 0 : Math.cos(i) * 120), y: H / 2 + (x.id === n.centre ? 0 : Math.sin(i) * 120)}));
      G.edges = n.edges; G.ticks = 0; cancelAnimationFrame(G.raf); gStep();
      const d = n.degree; $('#m-graph-status').textContent = `${n.centre}: ${d.associates} associations, ${d.member_of} pathways, ${d.has_domain} domains, ${d.expressed_in} tissues with nTPM`;
      const uncompiled = n.nodes.filter(x => x.kind === 'protein' && x.compiled === false).map(x => x.id);
      $('#m-graph-info').innerHTML = `showing ${n.nodes.length - 1} of its neighbours (highest confidence first) · blue protein (faded = not compiled locally), green pathway, purple domain, amber tissue · solid blue association = STRING experimental channel, dashed = other channels${uncompiled.length ? ` · <button class="ghost" id="m-graph-grow" style="padding:1px 8px;font-size:11px" data-tip="Compile the ${Math.min(8, uncompiled.length)} nearest uncompiled partners from the public databases (a few seconds each) so the graph grows around this protein.">🌱 compile ${Math.min(8, uncompiled.length)} partners</button>` : ''}`;
      const grow = $('#m-graph-grow');
      if (grow) grow.onclick = async () => {
        const todo = uncompiled.slice(0, 8);
        for (let i = 0; i < todo.length; i++) {
          $('#m-graph-status').textContent = `compiling ${todo[i]} (${i + 1}/${todo.length})…`;
          try { await window.api(`/api/protein_definition?gene=${encodeURIComponent(todo[i])}`); } catch (e) { /* a partner without a reviewed entry is skipped */ }
        }
        window.moleculesGraph(gene);
      };
    } catch (e) { $('#m-graph-status').textContent = e.message; }
  };
  document.addEventListener('DOMContentLoaded', () => {
    const c = gCanvas(); if (!c) return;
    const pick = ev => { const r = c.getBoundingClientRect(); const x = ev.clientX - r.left, y = ev.clientY - r.top; return G.nodes.find(n => (n.x - x) ** 2 + (n.y - y) ** 2 < 100); };
    c.addEventListener('mousedown', ev => { G.drag = pick(ev); });
    c.addEventListener('mousemove', ev => { if (G.drag) { const r = c.getBoundingClientRect(); G.drag.x = ev.clientX - r.left; G.drag.y = ev.clientY - r.top; gDraw(); } else { const n = pick(ev); c.title = n ? `${n.kind}: ${n.name || n.id}` : ''; } });
    window.addEventListener('mouseup', () => { G.drag = null; });
  });
  document.addEventListener('DOMContentLoaded', () => {
    const b = $('#m-report'); if (!b) return;
    b.onclick = async () => {
      const out = $('#m-report-out'); out.style.display = 'block'; out.textContent = 'assembling the dossier…';
      try { const r = await window.api(`/api/report?gene=${encodeURIComponent($('#m-gene').value.trim())}&chrom=${encodeURIComponent($('#m-chrom').value.trim())}`); out.textContent = r.markdown; }
      catch (e) { out.textContent = e.message; }
    };
  });
  window.moleculesRna = async function (gene, chrom) {
    $('#m-rna-status').textContent = 'listing transcripts, fetching GTEx…';
    try {
      const r = await window.api(`/api/rna?gene=${encodeURIComponent(gene)}&chrom=${encodeURIComponent(chrom || '')}`);
      const tx = r.transcripts, e = r.expression;
      let html = '';
      if (tx) {
        html += `<div><b>${tx.count} transcripts</b>, ${tx.coding_isoforms} coding · ${Object.entries(tx.by_biotype).map(([k, v]) => `${k.replace(/_/g, ' ')} ${v}`).join(', ')} <span class="pill curated">curated: GENCODE</span></div>
          <div class="scroll" style="max-height:220px"><table style="margin-top:6px;font-size:12px"><tr><th>transcript</th><th>biotype</th><th class="num">exons</th><th class="num">spliced nt</th><th class="num">CDS nt</th><th class="num">aa</th><th>tags</th></tr>
          ${tx.transcripts.map(x => `<tr><td class="mono">${x.name}</td><td>${x.biotype_label}</td><td class="num">${x.exons}</td><td class="num">${x.spliced_nt ?? '-'}</td><td class="num">${x.cds_nt ?? '-'}</td><td class="num">${x.protein_aa ?? '-'}</td><td class="muted">${x.tags.join(', ')}</td></tr>`).join('')}</table></div>`;
      } else if (chrom) html += `<div class="muted">no local gene models for ${chrom}; transcripts need chr21 or chrM</div>`;
      if (e && e.tissues && Object.keys(e.tissues).length) {
        const mx = e.max_tpm || 1;
        html += `<div style="margin-top:8px"><b>expression</b> ${e.pattern}; median ${e.median_tpm} TPM over ${e.tissues_measured} tissues, ${e.tissues_expressed} with ≥ 1 TPM <span class="pill experimental">measured: GTEx v8</span></div>
          <div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:2px 14px;font-size:11.5px;margin-top:4px">${Object.entries(e.tissues).map(([t, v]) => `<div style="display:flex;align-items:center;gap:6px"><span style="width:120px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="${t}">${t.replace(/_/g, ' ')}</span><i style="display:inline-block;height:8px;width:${Math.max(1, v / mx * 90)}px;background:var(--c1);border-radius:2px"></i><span class="muted num">${v.toFixed(0)}</span></div>`).join('')}</div>`;
      } else if (e) html += `<div class="muted" style="margin-top:6px">expression: ${e.error || 'unavailable'}</div>`;
      $('#m-rna').innerHTML = html || '<span class="muted">nothing found</span>';
      $('#m-rna-status').textContent = tx ? `${tx.gene} · ${tx.gene_type}` : '';
    } catch (err) { $('#m-rna-status').textContent = err.message; }
  };
  window.moleculesDefinition = async function (gene) {
    $('#m-def-status').textContent = 'compiling from Ensembl, UniProt, STRING, HPA…';
    try { const d = await window.api(`/api/protein_definition?gene=${encodeURIComponent(gene)}`); renderDefinition(d); $('#m-def-status').textContent = 'from local knowledge cache after the first compile'; }
    catch (e) { $('#m-def-status').textContent = e.message; }
    $('#m-def').querySelectorAll('[data-pw]').forEach(c => c.onclick = () => { $('#m-pw').value = c.dataset.pw; $('#m-ko').value = gene; window.moleculesPathway(c.dataset.pw, gene); });
  };
  document.addEventListener('DOMContentLoaded', () => { const b = $('#m-pw-run'); if (b) b.onclick = () => window.moleculesPathway($('#m-pw').value.trim(), $('#m-ko').value.trim()); });
  window.moleculesLoad = async function (gene, chrom) {
    window.moleculesDefinition(gene);
    window.moleculesRna(gene, chrom);
    window.moleculesGraph(gene);
    $('#m-status').textContent = 'fetching UniProt and AlphaFold…';
    try {
      const r = await window.api(`/api/protein?gene=${encodeURIComponent(gene)}&chrom=${encodeURIComponent(chrom || '')}`);
      const u = r.uniprot;
      if (!u) { $('#m-status').textContent = 'no reviewed UniProt entry'; ca = []; return; }
      $('#m-status').textContent = `${u.accession} · ${u.name} · ${u.length} aa`;
      $('#m-info').innerHTML = `<div><b>${u.name}</b> <span class="pill curated">curated</span> <span class="muted">${u.accession}, ${u.length} aa</span></div>
        ${r.ours ? `<div class="muted">our translation of ${r.ours.transcript}: ${r.ours.length} aa; identity with UniProt ${(r.agreement.identity * 100).toFixed(1)}%${r.agreement.same_length ? '' : ' (different isoform length)'}</div>` : ''}
        ${u.location ? `<div class="muted" style="font-size:12px">location: ${u.location}</div>` : ''}
        ${u.function ? `<div style="font-size:12.5px;margin-top:6px">${u.function}</div>` : ''}
        <table style="margin-top:8px">${u.features.slice(0, 20).map(f => `<tr><td>${f.type}</td><td class="mono num">${f.start}–${f.end}</td><td class="muted">${f.description}</td><td><button class="ghost" data-res="${f.start}" style="padding:1px 7px;font-size:11px">show</button></td></tr>`).join('')}</table>`;
      $('#m-info').querySelectorAll('button[data-res]').forEach(b => b.onclick = () => { highlight = +b.dataset.res; });
      const s = r.structure;
      if (s && s.ca) { ca = s.ca.map(p => [p[0], p[1], p[2]]); plddt = s.ca.map(p => p[3]); seq = s.residues || ''; $('#m-struct').innerHTML = `<b>AlphaFold ${s.entry}</b> <span class="pill predicted">predicted</span> · ${s.length} residues · mean pLDDT <b>${s.mean_plddt}</b> · ${(s.confident_fraction * 100).toFixed(0)}% confident · Rg ${s.radius_of_gyration_A} Å${s.sequence_matches_uniprot ? ' · sequence matches UniProt' : ''}<div class="muted" style="font-size:11.5px">click the structure to pause/spin; "show" on a feature marks its first residue</div>`; }
      else { ca = []; $('#m-struct').innerHTML = '<span class="muted">no AlphaFold model</span>'; }
      const res = +($('#m-residue').value || 0); highlight = res > 0 ? res : null;
    } catch (e) { $('#m-status').textContent = e.message; }
  };
})();
