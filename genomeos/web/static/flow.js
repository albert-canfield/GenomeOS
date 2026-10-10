// Flow: the DNA → RNA → protein relationship of one gene, drawn as three linked lanes.
// Hover any lane to see where that base / nucleotide / residue sits in the other two.
(function () {
  const $ = s => document.querySelector(s);
  const st = {data: null, hover: null, variant: null};
  const AA_COL = {A: '#8b949e', V: '#8b949e', L: '#8b949e', I: '#8b949e', M: '#8b949e', F: '#6e7681', W: '#6e7681', Y: '#6e7681',
                  S: '#1a7f37', T: '#1a7f37', N: '#1a7f37', Q: '#1a7f37', C: '#bf8700', G: '#bf8700', P: '#bf8700',
                  D: '#d1242f', E: '#d1242f', K: '#1f6feb', R: '#1f6feb', H: '#1f6feb', '*': '#000'};
  const AA3 = {A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Ter'};

  function canvas() { return $('#f-canvas'); }
  function geom() {
    const c = canvas(); const W = c.clientWidth, H = c.clientHeight;
    return {W, H, pad: 70, dnaY: 56, rnaY: H / 2 + 6, protY: H - 40, laneH: 22};
  }
  function xDNA(pos) { const d = st.data, g = geom(); return g.pad + (pos - d.gene_start) / (d.gene_end - d.gene_start) * (g.W - 2 * g.pad); }
  function xRNA(i) { const d = st.data, g = geom(); return g.pad + i / d.mrna_length * (g.W - 2 * g.pad); }
  function xProt(r) { const d = st.data, g = geom(); return g.pad + (r - 1) / Math.max(1, d.protein_length) * (g.W - 2 * g.pad); }

  function draw() {
    const c = canvas(); if (!c) return;
    const dpr = window.devicePixelRatio || 1;
    const g = geom(); if (!g.W) return;
    c.width = g.W * dpr; c.height = g.H * dpr;
    const ctx = c.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, g.W, g.H);
    const d = st.data; if (!d) return;
    const ink = getComputedStyle(document.documentElement).getPropertyValue('--ink').trim() || '#222';
    const muted = getComputedStyle(document.documentElement).getPropertyValue('--muted').trim() || '#777';
    ctx.font = '12px system-ui'; ctx.fillStyle = muted; ctx.textAlign = 'left';
    ctx.fillText(`1 DNA  ${d.chrom} ${d.strand}`, 6, g.dnaY - 16);
    if (d.regulation) ctx.fillText('regulation', 6, g.dnaY - 30);
    ctx.fillText('2 RNA  mRNA', 6, g.rnaY - 16);
    ctx.fillText('3 protein', 6, g.protY - 16);
    // DNA lane: gene span, exons as bars
    ctx.strokeStyle = '#1f6feb'; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(xDNA(d.gene_start), g.dnaY); ctx.lineTo(xDNA(d.gene_end), g.dnaY); ctx.stroke();
    for (const e of d.exons) {
      const x0 = xDNA(e.start), x1 = xDNA(e.end);
      ctx.fillStyle = '#8250df'; ctx.fillRect(x0, g.dnaY - 9, Math.max(1.5, x1 - x0), 18);
    }
    // CDS on DNA (map mRNA cds range back to genomic via exons)
    if (d.cds_start >= 0) {
      for (const e of d.exons) {
        const a = Math.max(e.mrna_start, d.cds_start), b = Math.min(e.mrna_end, d.cds_end);
        if (a >= b) continue;
        const ga = d.mrna_to_genomic[a], gb = d.mrna_to_genomic[b - 1];
        const x0 = xDNA(Math.min(ga, gb)), x1 = xDNA(Math.max(ga, gb) + 1);
        ctx.fillStyle = '#d1242f'; ctx.fillRect(x0, g.dnaY - 9, Math.max(1.5, x1 - x0), 18);
      }
    }
    // regulation: promoter (green) and enhancers (amber) that reach this gene, above the DNA lane
    const reg = d.regulation;
    if (reg) {
      for (const e of [...reg.promoters, ...reg.enhancers]) {
        if (e.end < d.gene_start - 2000 || e.start > d.gene_end + 2000) continue;
        const x0 = xDNA(Math.max(d.gene_start, e.start)), x1 = xDNA(Math.min(d.gene_end, e.end));
        ctx.fillStyle = e.class === 'enhancer' ? '#e3b341' : '#2ea043';
        ctx.fillRect(x0 - 1, g.dnaY - 20, Math.max(2, x1 - x0), 7);
      }
    }
    // splice lines exon → mRNA
    ctx.strokeStyle = muted; ctx.lineWidth = 0.6; ctx.globalAlpha = 0.5;
    for (const e of d.exons) {
      const gx0 = d.strand === '-' ? xDNA(e.end) : xDNA(e.start), gx1 = d.strand === '-' ? xDNA(e.start) : xDNA(e.end);
      ctx.beginPath(); ctx.moveTo(gx0, g.dnaY + 9); ctx.lineTo(xRNA(e.mrna_start), g.rnaY - 9); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(gx1, g.dnaY + 9); ctx.lineTo(xRNA(e.mrna_end), g.rnaY - 9); ctx.stroke();
    }
    ctx.globalAlpha = 1;
    // RNA lane: UTR5 | CDS | UTR3 with exon junction ticks
    ctx.fillStyle = '#8250df'; ctx.fillRect(xRNA(0), g.rnaY - 9, xRNA(d.mrna_length) - xRNA(0), 18);
    if (d.cds_start >= 0) { ctx.fillStyle = '#d1242f'; ctx.fillRect(xRNA(d.cds_start), g.rnaY - 9, xRNA(d.cds_end) - xRNA(d.cds_start), 18); }
    ctx.strokeStyle = '#fff'; ctx.lineWidth = 1;
    for (const e of d.exons.slice(1)) { ctx.beginPath(); ctx.moveTo(xRNA(e.mrna_start), g.rnaY - 9); ctx.lineTo(xRNA(e.mrna_start), g.rnaY + 9); ctx.stroke(); }
    // translation lines CDS → protein
    if (d.cds_start >= 0) {
      ctx.strokeStyle = muted; ctx.lineWidth = 0.6; ctx.globalAlpha = 0.5;
      ctx.beginPath(); ctx.moveTo(xRNA(d.cds_start), g.rnaY + 9); ctx.lineTo(xProt(1), g.protY - 9); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(xRNA(d.cds_end), g.rnaY + 9); ctx.lineTo(xProt(d.protein_length + 1), g.protY - 9); ctx.stroke();
      ctx.globalAlpha = 1;
      // protein lane: residues coloured by chemistry
      const w = Math.max(0.5, (g.W - 2 * g.pad) / Math.max(1, d.protein_length));
      for (let i = 0; i < d.protein.length; i++) { ctx.fillStyle = AA_COL[d.protein[i]] || '#999'; ctx.fillRect(xProt(i + 1), g.protY - 9, w + 0.3, 18); }
    }
    // hover link
    const h = st.hover;
    if (h && h.mrna !== undefined) {
      ctx.strokeStyle = ink; ctx.lineWidth = 1.5;
      const gx = xDNA(h.genomic) , rx = xRNA(h.mrna);
      ctx.beginPath(); ctx.moveTo(gx, g.dnaY - 14); ctx.lineTo(gx, g.dnaY + 14); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(rx, g.rnaY - 14); ctx.lineTo(rx, g.rnaY + 14); ctx.stroke();
      ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(gx, g.dnaY + 14); ctx.lineTo(rx, g.rnaY - 14); ctx.stroke();
      if (h.residue) { const px = xProt(h.residue) + (g.W - 2 * g.pad) / d.protein_length / 2; ctx.beginPath(); ctx.moveTo(rx, g.rnaY + 14); ctx.lineTo(px, g.protY - 14); ctx.stroke(); ctx.setLineDash([]); ctx.beginPath(); ctx.moveTo(px, g.protY - 14); ctx.lineTo(px, g.protY + 14); ctx.stroke(); }
      ctx.setLineDash([]);
    }
    // individuals: one row of ticks per person under the DNA lane (measured genotypes; coding SNVs in red)
    const people = d.individuals || [];
    people.forEach((p, i) => {
      const y = g.dnaY + 16 + i * 9;
      if (y > g.rnaY - 30) return;
      ctx.fillStyle = muted; ctx.font = '10px system-ui'; ctx.textAlign = 'right'; ctx.fillText(p.name, g.pad - 6, y + 4); ctx.textAlign = 'left';
      ctx.fillStyle = '#8b949e'; ctx.globalAlpha = 0.7;
      for (const [pos] of p.positions) { const x = xDNA(pos); ctx.fillRect(x, y, 1, 6); }
      ctx.globalAlpha = 1; ctx.fillStyle = '#d1242f';
      for (const c of p.coding_snvs) { if (c.consequence === 'synonymous') continue; const x = xDNA(c.pos - 1); ctx.fillRect(x - 1, y - 1, 3, 8); }
    });
    // variant marker
    const v = st.variant;
    if (v && v.mrna !== undefined) {
      const gx = xDNA(v.genomic); ctx.fillStyle = '#d1242f';
      ctx.beginPath(); ctx.moveTo(gx, g.dnaY - 22); ctx.lineTo(gx - 5, g.dnaY - 30); ctx.lineTo(gx + 5, g.dnaY - 30); ctx.closePath(); ctx.fill();
      if (v.residue) { const px = xProt(v.residue); ctx.beginPath(); ctx.moveTo(px, g.protY + 22); ctx.lineTo(px - 5, g.protY + 30); ctx.lineTo(px + 5, g.protY + 30); ctx.closePath(); ctx.fill(); }
    }
  }

  function whereMrna(m) {
    const d = st.data; const out = {mrna: m, genomic: d.mrna_to_genomic[m], base: d.mrna[m]};
    if (d.cds_start < 0) out.region = 'non-coding';
    else if (m < d.cds_start) out.region = "5'UTR";
    else if (m >= d.cds_end) out.region = "3'UTR";
    else { const c = m - d.cds_start; out.region = 'CDS'; out.cds = c; out.residue = Math.floor(c / 3) + 1; out.codon = d.mrna.slice(d.cds_start + Math.floor(c / 3) * 3, d.cds_start + Math.floor(c / 3) * 3 + 3); out.pos = c % 3 + 1; out.aa = d.protein[out.residue - 1] || '*'; }
    const ex = d.exons.find(e => m >= e.mrna_start && m < e.mrna_end); out.exon = ex ? ex.index : null;
    return out;
  }
  function describe(h) {
    if (!h) return '';
    const d = st.data;
    let s = `<b>${d.chrom}:${(h.genomic + 1).toLocaleString()}</b> ${h.base.replace('U', 'T')} <span class="muted">(exon ${h.exon} of ${d.exons.length})</span> → mRNA nt <b>${h.mrna + 1}</b> ${h.base} <span class="muted">${h.region}</span>`;
    if (h.residue) s += ` → codon <b>${h.codon}</b> position ${h.pos} → residue <b>${AA3[h.aa] || h.aa}${h.residue}</b>`;
    return s;
  }

  function onMove(ev) {
    const d = st.data; if (!d) return;
    const r = canvas().getBoundingClientRect(); const x = ev.clientX - r.left, y = ev.clientY - r.top; const g = geom();
    const t = (x - g.pad) / (g.W - 2 * g.pad); let h = null;
    if (Math.abs(y - g.dnaY) < 16) {
      const pos = Math.floor(d.gene_start + t * (d.gene_end - d.gene_start));
      const m = d.mrna_to_genomic.indexOf(pos); // exact base; intronic → nearest exon edge
      if (m >= 0) h = whereMrna(m); else { h = {intron: true, genomic: pos}; }
    } else if (Math.abs(y - g.rnaY) < 16) { const m = Math.min(d.mrna_length - 1, Math.max(0, Math.floor(t * d.mrna_length))); h = whereMrna(m); }
    else if (Math.abs(y - g.protY) < 16 && d.cds_start >= 0) { const res = Math.min(d.protein_length, Math.max(1, Math.floor(t * d.protein_length) + 1)); h = whereMrna(d.cds_start + (res - 1) * 3); }
    st.hover = h; draw();
    $('#f-hover').innerHTML = h ? (h.intron ? `<b>${d.chrom}:${(h.genomic + 1).toLocaleString()}</b> <span class="muted">intron: removed by splicing, not in the mRNA</span>` : describe(h)) : '<span class="muted">hover a lane</span>';
  }

  // ---- the tissue's own isoform beside the canonical (/api/flow?tissue=): GTEx v8 dominant isoform when
  // the rule establishes it, else UNKNOWN with the API's reason shown verbatim; the canonical is never drawn
  // or described as the tissue's.
  const esc = t => String(t ?? '').replace(/[&<>"]/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c]));
  const REASON_GLOSS = {
    no_isoform_data: () => 'no GTEx isoform table is cached for this gene',
    tissue_not_measured: () => 'GTEx did not measure this tissue for the gene (names match exactly or after folding case and punctuation)',
    low_expression: r => `the top isoform's median is below ${r.min_tpm} TPM`,
    tie: r => `the next isoform's median is at least ${r.tie_ratio} of the top's`,
    dominant_not_in_local_models: () => 'the top isoform is not in the local GENCODE models',
    dominant_non_coding: () => 'the top isoform codes no protein in the local GENCODE models',
  };
  function cdsSpans(t) { // genomic [lo, hi) pieces of a trace's CDS, read back from its mRNA coordinates
    const out = []; if (!t || t.cds_start < 0) return out;
    for (const e of t.exons) {
      const a = Math.max(e.mrna_start, t.cds_start), b = Math.min(e.mrna_end, t.cds_end); if (a >= b) continue;
      const ga = t.mrna_to_genomic[a], gb = t.mrna_to_genomic[b - 1]; out.push([Math.min(ga, gb), Math.max(ga, gb) + 1]);
    }
    return out;
  }
  const exonsLacking = (a, b) => b.exons.filter(e => !a.exons.some(f => f.start < e.end && e.start < f.end)); // exons of b that a lacks
  function callText(c) {
    const nxt = c.next_transcript ? `, next ${esc(c.next_transcript)} ${c.next_tpm} TPM` : '';
    return `<span class="mono">${esc(c.transcript || '-')}</span> ${c.tpm != null ? c.tpm + ' TPM' : 'no TPM'}${c.share != null ? `, share ${c.share}` : ''}${nxt}`;
  }
  function consText(s) { return s ? `${s.hgvs_c ? esc(s.hgvs_c) + ' → ' : ''}${s.hgvs_p ? esc(s.hgvs_p) + ' ' : ''}<b>${esc(s.consequence)}</b>` : ''; }
  function tissueHtml(d) {
    const t = d && d.tissue; if (!t) return '';
    const c = t.call || {}, rule = c.rule || {}, name = `<span class="mono">${esc(t.tissue)}</span>`;
    const canon = `canonical ${esc(t.canonical_transcript)}, ${t.canonical_protein_length} aa`;
    let h;
    if (t.status === 'ESTABLISHED') {
      const iso = t.trace, lacks = iso ? exonsLacking(iso, d) : [], extra = iso ? exonsLacking(d, iso) : [];
      const diff = lacks.length || extra.length ? ` · exons: ${lacks.length ? `lacks canonical exon${lacks.length === 1 ? '' : 's'} ${lacks.map(e => e.index).join(', ')}` : ''}${lacks.length && extra.length ? '; ' : ''}${extra.length ? `${extra.length} exon${extra.length === 1 ? '' : 's'} the canonical lacks` : ''} <span class="muted">(dashed in the drawing)</span>` : '';
      h = `<div><b>tissue</b> ${name}: dominant isoform <b>${esc(t.transcript)}</b> <span class="mono">${esc(t.transcript_id)}</span>${t.same_as_canonical ? ' (the canonical itself)' : ''}, ${t.protein_length} aa; ${canon}${diff} · ${callText(c)} <span class="pill experimental">measured: GTEx v8 median, not this person's expression</span></div>`;
    } else {
      const gloss = REASON_GLOSS[t.reason];
      h = `<div class="f-tissue-unknown"><b>tissue</b> ${name}: <span class="pill none">UNKNOWN</span> <span class="mono f-tissue-reason">${esc(t.reason)}</span>${gloss ? ` <span class="muted">(${esc(gloss(rule))})</span>` : ''} · ${callText(c)}; ${canon}
        <div class="muted" style="font-size:11.5px">The tissue's isoform is not established, so it is not drawn and its consequence is UNKNOWN, never the canonical's.</div></div>`;
    }
    const v = d.variant_in_tissue;
    if (v) {
      const there = v.status === 'ESTABLISHED' ? consText({...v, hgvs_c: v.isoform ? v.isoform.hgvs_c : null}) : `<b>UNKNOWN</b> <span class="mono f-tissue-reason">${esc(v.reason)}</span>`;
      h += `<div><b>variant</b> in ${esc(v.tissue)}${v.transcript ? ` on ${esc(v.transcript)}` : ''}: ${there} · on the canonical ${esc(v.canonical_transcript)}: ${consText(v.canonical)}${v.agreement ? ` · <span class="mono">${esc(v.agreement)}</span>` : ''}${v.status === 'ESTABLISHED' ? '' : ` <span class="muted">(the canonical's consequence, not the tissue's)</span>`} <span class="pill predicted">predicted by the local trace</span></div>`;
    }
    return h + `<div class="muted" style="font-size:11.5px">${esc(t.note)}; rule: the top isoform at least ${rule.min_tpm} TPM and the next below ${rule.tie_ratio} × the top's.</div>`;
  }
  window.flowTissueHtml = tissueHtml;
  function drawTissue() {
    const c = $('#f-tcanvas'), d = st.data; if (!c) return;
    const t = d && d.tissue; c.style.display = t ? 'block' : 'none'; if (!t) return;
    const W = c.clientWidth, H = c.clientHeight; if (!W) return;
    const dpr = window.devicePixelRatio || 1; c.width = W * dpr; c.height = H * dpr;
    const ctx = c.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
    const css = getComputedStyle(document.documentElement);
    const muted = css.getPropertyValue('--muted').trim() || '#777', bad = css.getPropertyValue('--bad').trim() || '#d1242f';
    const pad = 70, span = W - 2 * pad, xg = p => pad + (p - d.gene_start) / (d.gene_end - d.gene_start) * span;
    const iso = t.status === 'ESTABLISHED' ? t.trace : null, v = d.variant_in_tissue;
    const text = (s, x, y, col, size, align) => { ctx.font = `${size || 12}px system-ui`; ctx.fillStyle = col || muted; ctx.textAlign = align || 'left'; ctx.fillText(s, x, y); ctx.textAlign = 'left'; };
    text('exons on the DNA, same scale as above', 6, 12); text('protein, one residue scale', 6, 78);
    text('canonical', pad - 6, 30, muted, 10, 'right'); text('canonical', pad - 6, 102, muted, 10, 'right');
    text('tissue', pad - 6, 48, muted, 10, 'right'); text('tissue', pad - 6, 134, muted, 10, 'right');
    const exRow = (tr, other, y) => {
      ctx.strokeStyle = '#1f6feb'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(xg(d.gene_start), y + 6); ctx.lineTo(xg(d.gene_end), y + 6); ctx.stroke();
      for (const e of tr.exons) { ctx.fillStyle = '#8250df'; ctx.fillRect(xg(e.start), y, Math.max(1.5, xg(e.end) - xg(e.start)), 12); }
      for (const [a, b] of cdsSpans(tr)) { ctx.fillStyle = '#d1242f'; ctx.fillRect(xg(a), y, Math.max(1.5, xg(b) - xg(a)), 12); }
      ctx.strokeStyle = muted; ctx.setLineDash([2, 2]);
      for (const e of exonsLacking(tr, other)) ctx.strokeRect(xg(e.start) - 1, y - 1, Math.max(2, xg(e.end) - xg(e.start)) + 2, 14);
      ctx.setLineDash([]);
    };
    const L = Math.max(d.protein_length || 0, iso ? iso.protein_length : 0, 1), w = span / L;
    const protRow = (tr, y, residue, label) => {
      for (let i = 0; i < (tr.protein || '').length; i++) { ctx.fillStyle = AA_COL[tr.protein[i]] || '#999'; ctx.fillRect(pad + i * w, y, w + 0.3, 12); }
      if (!residue) return;
      const x = pad + (residue - 0.5) * w; ctx.fillStyle = '#d1242f';
      ctx.beginPath(); ctx.moveTo(x, y - 1); ctx.lineTo(x - 4, y - 8); ctx.lineTo(x + 4, y - 8); ctx.closePath(); ctx.fill();
      ctx.font = '11px system-ui'; const right = x + 7 + ctx.measureText(label).width > W;
      text(label, right ? x - 7 : x + 7, y - 2, '#d1242f', 11, right ? 'right' : 'left');
    };
    exRow(d, iso || d, 20); protRow(d, 92, v ? v.canonical.residue : null, v ? (v.canonical.hgvs_p || v.canonical.consequence) : '');
    if (iso) { exRow(iso, d, 38); protRow(iso, 124, v ? v.residue : null, v ? (v.hgvs_p || v.consequence) : ''); }
    else for (const y of [48, 134]) text(`UNKNOWN: ${t.reason} (not established, so the canonical is not drawn in its place)`, pad, y, bad, 11);
    const gpos = [d.variant && d.variant.genomic, v && v.isoform && v.isoform.genomic].find(x => x != null);
    if (gpos != null && gpos >= d.gene_start && gpos < d.gene_end) { ctx.strokeStyle = '#d1242f'; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(xg(gpos), 17); ctx.lineTo(xg(gpos), 53); ctx.stroke(); }
  }
  function fillTissues(d, asked) { // the selector: the tissues this gene's cached GTEx isoform table measures
    const sel = $('#f-tissue'); if (!sel) return;
    const list = d.tissues_measured || [], chosen = d.tissue ? (list.includes(d.tissue.tissue) ? d.tissue.tissue : asked) : '';
    const opt = (val, label, off) => `<option value="${esc(val)}"${off ? ' disabled' : ''}>${esc(label)}</option>`;
    sel.innerHTML = opt('', 'canonical only') + (list.length ? '' : opt('', `no GTEx isoform table cached for ${d.gene}`, true)) +
      (chosen && !list.includes(chosen) ? opt(chosen, `${chosen} (not measured for ${d.gene})`) : '') + list.map(x => opt(x, x.replace(/_/g, ' '))).join('');
    sel.value = chosen;
  }
  function setTissue(val) { const sel = $('#f-tissue'); if (!sel || !val) return; if (![...sel.options].some(o => o.value === val)) { const o = document.createElement('option'); o.value = val; o.textContent = val; sel.appendChild(o); } sel.value = val; }
  function drawAll() { draw(); drawTissue(); }

  window.flowLoad = async function () {
    const gene = $('#f-gene').value.trim(), chrom = $('#f-chrom').value.trim(), variant = $('#f-variant').value.trim();
    const tissue = $('#f-tissue') ? $('#f-tissue').value : '';
    $('#f-status').textContent = 'tracing…';
    try {
      const d = await window.api(`/api/flow?gene=${encodeURIComponent(gene)}&chrom=${encodeURIComponent(chrom)}${variant ? '&variant=' + encodeURIComponent(variant) : ''}${tissue ? '&tissue=' + encodeURIComponent(tissue) : ''}`);
      st.data = d; st.variant = d.variant || null; st.hover = null;
      fillTissues(d, tissue);
      $('#f-tissue-out').innerHTML = tissueHtml(d);
      setTimeout(() => document.querySelectorAll('.f-ind-var').forEach(a => a.onclick = ev => { ev.preventDefault(); $('#f-variant').value = a.dataset.v; $('#f-load').click(); }), 0);
      $('#f-status').textContent = `${d.gene} · ${d.transcript} (canonical of ${d.transcripts})`;
      const ev = d.evidence;
      $('#f-summary').innerHTML = `
        <div><b>1 DNA</b> ${d.chrom}:${(d.gene_start + 1).toLocaleString()}–${d.gene_end.toLocaleString()} (${d.strand}) · ${(d.gene_end - d.gene_start).toLocaleString()} bp · ${d.exons.length} exons <span class="pill curated">${ev.exons}</span></div>
        ${d.regulation ? `<div><b>1b regulation</b> node ${d.regulation.domain ? d.regulation.domain.id : '?'}${d.regulation.domain ? ` (${(d.regulation.domain.length / 1000).toFixed(0)} kb, ${d.regulation.domain.coding_genes} coding genes)` : ''}: ${d.regulation.promoters.length} promoter element${d.regulation.promoters.length === 1 ? '' : 's'}, ${d.regulation.enhancers_in_domain} enhancers can reach it (${d.regulation.enhancers_nearest_to_this_gene} nearest to this gene, ${d.regulation.enhancers_inside_gene} inside the gene), ${d.regulation.insulators_bounding.length} bounding insulators <span class="pill curated">curated: ENCODE elements</span> <span class="pill inferred">inferred: reach bounded by the CTCF domain</span>${d.regulation.enhancers_predicted ? ` · ${d.regulation.enhancers_predicted} scored by AlphaGenome, ${d.regulation.enhancers_predicted_this_gene} name this gene <span class="pill predicted">predicted</span>` : ''}<div class="muted" style="font-size:11.5px">green tick = promoter, amber = enhancer, drawn where they fall in the gene span; the full list is in <span class="mono">genomeos regulation --gene ${d.gene}</span></div></div>` : ''}
        ${(d.individuals || []).length ? `<div><b>1c individuals</b> ${d.individuals.map(p => `<span class="mono">${p.name}</span> ${p.variants_in_gene.toLocaleString()} variants inside the gene${p.coding_snvs.length ? ', coding: ' + p.coding_snvs.filter(c => c.consequence !== 'synonymous').slice(0, 8).map(c => `<a href="#" class="f-ind-var" data-v="${d.chrom}:${c.pos} ${c.ref}>${c.alt}">${c.hgvs_p || c.consequence}</a> <span class="muted">(${c.genotype})</span>`).join('; ') + (p.coding_snvs.some(c => c.consequence === 'synonymous') ? ` + ${p.coding_snvs.filter(c => c.consequence === 'synonymous').length} synonymous` : '') : ''}`).join(' · ')} <span class="pill experimental">measured</span> <span class="muted">genotypes from each person's file; grey ticks under the DNA lane, red = protein-changing; click a change to trace it</span></div>` : ''}
        <div><b>2 RNA</b> spliced mRNA ${d.mrna_length.toLocaleString()} nt = 5'UTR ${d.utr5} + CDS ${d.cds_length} + 3'UTR ${d.utr3}; introns removed: ${((d.gene_end - d.gene_start) - d.mrna_length).toLocaleString()} bp (${(100 - d.mrna_length / (d.gene_end - d.gene_start) * 100).toFixed(1)}% of the gene) <span class="pill derived">${ev.mrna}</span></div>
        <div><b>3 protein</b> ${d.protein_length} aa, ${d.codon_table} code${d.tags.length ? ' · ' + d.tags.join(', ') : ''} <span class="pill derived">${ev.protein}</span>
          <div class="mono" style="font-size:11px;word-break:break-all;max-height:64px;overflow:auto;margin-top:4px">${d.protein}</div></div>
        ${d.variant ? `<div style="margin-top:6px"><b>variant</b> ${d.variant.consequence === 'reference mismatch' ? `reference base at that position is ${d.variant.expected}` : d.variant.hgvs_p ? `${d.variant.hgvs_c} → ${d.variant.hgvs_p} <b>${d.variant.consequence}</b> (codon ${d.variant.codon_before}→${d.variant.codon_after})${d.variant.residue ? ` <button class="ghost" id="f-to-structure" data-tip="Open the Molecules tab with residue ${d.variant.residue} marked on the AlphaFold model and the UniProt features" style="padding:1px 8px;font-size:11px">🧊 residue ${d.variant.residue} on the structure</button>` : ''}` : `<b>${d.variant.consequence}</b>: no change to the protein`}</div>` : ''}`;
      const toS = $('#f-to-structure');
      if (toS) toS.onclick = () => {
        $('#m-gene').value = d.gene; $('#m-chrom').value = d.chrom; $('#m-residue').value = d.variant.residue;
        const b = document.querySelector('nav button[data-tab="molecules"]'); if (b) b.click();
        setTimeout(() => $('#m-load').click(), 200);
      };
      drawAll();
    } catch (e) { $('#f-status').textContent = e.message; }
  };
  window.lookupVariant = async function () {
    const v = $('#l-variant').value.trim(); if (!v) return;
    $('#l-status').textContent = 'asking VEP, reading the local layers…';
    try {
      const r = await window.api(`/api/lookup?variant=${encodeURIComponent(v)}`);
      const x = r.vep, lt = r.local_trace, p = r.protein, pw = r.pathways;
      const sev = x.consequence.includes('stop') || x.consequence.includes('frameshift') || x.consequence.includes('start_lost') ? 'var(--bad)' : x.consequence.includes('missense') ? 'var(--warn)' : 'var(--ok)';
      $('#l-out').innerHTML = `
        <div><b>${r.variant}</b> → <b style="color:${sev}">${x.consequence.replace(/_/g, ' ')}</b> in <b>${x.gene || '-'}</b> <span class="mono">${x.hgvsc || ''} ${x.hgvsp || ''}</span>${x.sift ? ` · SIFT ${x.sift}` : ''}${x.polyphen ? ` · PolyPhen ${x.polyphen.replace(/_/g, ' ')}` : ''} <span class="pill curated">curated: VEP</span></div>
        <div><b>known</b>: ${r.known.join('; ')}${x.dbsnp && x.dbsnp.length ? ` (${x.dbsnp.join(', ')})` : ''}${x.pubmed && x.pubmed.length ? `<div class="muted" style="font-size:11.5px">PubMed: ${x.pubmed.slice(0, 8).map(id => `<a href="https://pubmed.ncbi.nlm.nih.gov/${id}/" target="_blank" rel="noopener">${id}</a>`).join(', ')}${x.pubmed_count > 8 ? ' …' : ''}</div>` : ''}</div>
        ${(r.individuals || []).map(i => `<div><b>${i.individual}</b>: ${i.carries ? `carries it (genotype ${i.genotype})` : 'does not carry it'} <span class="pill experimental">measured</span></div>`).join('')}
        ${lt ? `<div><b>local trace</b> (${lt.transcript}): ${lt.consequence} <span class="mono">${lt.hgvs_c || ''} ${lt.hgvs_p || ''}</span> · ${lt.agrees_with_vep ? 'agrees with VEP' : 'differs from VEP (isoform numbering or region)'} <span class="pill derived">derived</span></div>` : ''}
        ${p && !p.error ? `<div><b>protein</b> ${p.accession} ${p.name} (${p.length} aa), residue ${p.residue}${p.reference_residue_matches ? '' : ' <span class="muted">(reference residue differs: isoform?)</span>'}: ${p.features_at_residue.length ? p.features_at_residue.map(f => `<span class="chip">${f.type}${f.description ? ' · ' + f.description : ''} ${f.start}–${f.end}</span>`).join(' ') : '<span class="muted">no annotated feature at this residue</span>'} · structures: ${p.structures_experimental} experimental${p.alphafold ? ', AlphaFold' : ''} <span class="pill curated">curated: UniProt</span> <button class="ghost" id="l-to-structure" style="padding:1px 8px;font-size:11px">🧊 show on structure</button></div>` : p ? `<div class="muted">protein: ${p.error}</div>` : ''}
        ${pw ? `<div><b>pathways</b>: ${pw.reactions_lost} reactions lost over ${pw.checked} pathways with the protein ${pw.treated_as}${pw.most_affected ? `; most affected ${pw.most_affected.name} (${(pw.most_affected.fraction_lost * 100).toFixed(0)}%)` : ''} <span class="pill inferred">inferred</span></div>` : ''}`;
      const b = $('#l-to-structure');
      if (b) b.onclick = () => { $('#m-gene').value = x.gene; $('#m-chrom').value = r.chrom; $('#m-residue').value = p.residue; const t = document.querySelector('nav button[data-tab="molecules"]'); if (t) t.click(); setTimeout(() => $('#m-load').click(), 200); };
      $('#l-status').textContent = '';
    } catch (e) { $('#l-status').textContent = e.message; }
  };
  document.addEventListener('DOMContentLoaded', () => { const b = $('#l-run'); if (b) { b.onclick = window.lookupVariant; $('#l-variant').onkeydown = e => { if (e.key === 'Enter') window.lookupVariant(); }; } });
  window.flowResize = drawAll;
  document.addEventListener('DOMContentLoaded', () => {
    const c = canvas(); if (!c) return;
    c.addEventListener('mousemove', onMove); c.addEventListener('mouseleave', () => { st.hover = null; draw(); });
    $('#f-load').onclick = window.flowLoad;
    $('#f-gene').onkeydown = e => { if (e.key === 'Enter') window.flowLoad(); };
    if ($('#f-tissue')) $('#f-tissue').onchange = window.flowLoad;
    window.addEventListener('resize', drawAll);
  });
  document.addEventListener('DOMContentLoaded', () => {
    // deep link: #flow?gene=APP&chrom=chr21[&variant=chr21:25897620 C>T][&tissue=Brain_Cerebellum]
    const q = ((window.initialHash || location.hash).split('?')[1] || '');
    const gm = q.match(/gene=([^&]+)/), cm = q.match(/chrom=([^&]+)/), vm = q.match(/variant=([^&]+)/), tm = q.match(/tissue=([^&]+)/);
    if (gm && $('#f-gene') && (window.initialHash || location.hash).startsWith('#flow')) { $('#f-gene').value = decodeURIComponent(gm[1]); if (cm) $('#f-chrom').value = decodeURIComponent(cm[1]); if (vm) $('#f-variant').value = decodeURIComponent(vm[1]); if (tm) setTissue(decodeURIComponent(tm[1])); setTimeout(() => $('#f-load').click(), 400); }
  });
})();
