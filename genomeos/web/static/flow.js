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

  window.flowLoad = async function () {
    const gene = $('#f-gene').value.trim(), chrom = $('#f-chrom').value.trim(), variant = $('#f-variant').value.trim();
    $('#f-status').textContent = 'tracing…';
    try {
      const d = await window.api(`/api/flow?gene=${encodeURIComponent(gene)}&chrom=${encodeURIComponent(chrom)}${variant ? '&variant=' + encodeURIComponent(variant) : ''}`);
      st.data = d; st.variant = d.variant || null; st.hover = null;
      $('#f-status').textContent = `${d.gene} · ${d.transcript} (canonical of ${d.transcripts})`;
      const ev = d.evidence;
      $('#f-summary').innerHTML = `
        <div><b>1 DNA</b> ${d.chrom}:${(d.gene_start + 1).toLocaleString()}–${d.gene_end.toLocaleString()} (${d.strand}) · ${(d.gene_end - d.gene_start).toLocaleString()} bp · ${d.exons.length} exons <span class="pill curated">${ev.exons}</span></div>
        ${d.regulation ? `<div><b>1b regulation</b> node ${d.regulation.domain ? d.regulation.domain.id : '?'}${d.regulation.domain ? ` (${(d.regulation.domain.length / 1000).toFixed(0)} kb, ${d.regulation.domain.coding_genes} coding genes)` : ''}: ${d.regulation.promoters.length} promoter element${d.regulation.promoters.length === 1 ? '' : 's'}, ${d.regulation.enhancers_in_domain} enhancers can reach it (${d.regulation.enhancers_nearest_to_this_gene} nearest to this gene, ${d.regulation.enhancers_inside_gene} inside the gene), ${d.regulation.insulators_bounding.length} bounding insulators <span class="pill curated">curated: ENCODE elements</span> <span class="pill inferred">inferred: reach bounded by the CTCF domain</span><div class="muted" style="font-size:11.5px">green tick = promoter, amber = enhancer, drawn where they fall in the gene span; the full list is in <span class="mono">genomeos regulation --gene ${d.gene}</span></div></div>` : ''}
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
      draw();
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
  window.flowResize = draw;
  document.addEventListener('DOMContentLoaded', () => {
    const c = canvas(); if (!c) return;
    c.addEventListener('mousemove', onMove); c.addEventListener('mouseleave', () => { st.hover = null; draw(); });
    $('#f-load').onclick = window.flowLoad;
    $('#f-gene').onkeydown = e => { if (e.key === 'Enter') window.flowLoad(); };
    window.addEventListener('resize', draw);
  });
})();
