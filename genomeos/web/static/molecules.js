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
      ${sec('Pathways (Reactome)', pw.length ? pw.slice(0, 12).map(x => `<span class="chip" data-tip="${x.id}">${x.name}</span>`).join(' ') + (pw.length > 12 ? ` <span class="muted">+${pw.length - 12} more</span>` : '') : '<span class="muted">none</span>', pill(s.pathways))}
      ${sec('Associations (STRING ≥ 0.7)', its.length ? its.slice(0, 20).map(x => `<span class="chip" data-tip="combined ${x.score}; experimental ${x.experimental}, database ${x.database}, text-mining ${x.textmining}${x.physical_evidence ? ' · experimental support' : ''}" style="${x.physical_evidence ? '' : 'opacity:.6'}">${x.partner}</span>`).join(' ') + ` <span class="muted">${its.length} total, ${its.filter(x => x.physical_evidence).length} with experimental support; dimmed = association without physical evidence</span>` : `<span class="muted">${it.error || 'none'}</span>`, pill(it))}
      ${sec('Expression (Human Protein Atlas)', ex.tissue_specificity ? `${ex.tissue_specificity}; ${ex.tissue_distribution || ''}; main location ${(ex.subcellular_main || []).join(', ') || '?'}; cell types: ${ex.cell_type_specificity || '?'}${ntpm.length ? `<table style="margin-top:4px">${ntpm.map(([k, v]) => `<tr><td>${k}</td><td class="num mono">${v.toFixed(0)} nTPM</td></tr>`).join('')}</table>` : ''}` : `<span class="muted">${(s.expression || {}).error || 'not fetched'}</span>`, pill(s.expression))}
      ${dis.length ? sec('Diseases (UniProt)', dis.map(x => `<span class="chip" data-tip="${(x.description || '').replace(/"/g, '')}">${x.name}${x.mim ? ' · MIM ' + x.mim : ''}</span>`).join(' '), pill(s.diseases)) : ''}
      ${d.states && d.states.length ? `<div class="muted" style="margin-top:8px;font-size:12px">${d.states.length} ProteinState records derived (protein × tissue/cell type × level × location); the definition above is what the protein <i>is</i>, a state is where it is and how much.</div>` : ''}`;
  }
  window.moleculesDefinition = async function (gene) {
    $('#m-def-status').textContent = 'compiling from Ensembl, UniProt, STRING, HPA…';
    try { const d = await window.api(`/api/protein_definition?gene=${encodeURIComponent(gene)}`); renderDefinition(d); $('#m-def-status').textContent = 'from local knowledge cache after the first compile'; }
    catch (e) { $('#m-def-status').textContent = e.message; }
  };
  window.moleculesLoad = async function (gene, chrom) {
    window.moleculesDefinition(gene);
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
