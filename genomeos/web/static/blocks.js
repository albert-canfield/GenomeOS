// GenomeOS block viewer / editor: a 2-D map of the delimited blocks of a chromosome.
// Genes contain transcripts contain exons; UNKNOWN spans sit between genes; islands,
// repeats, gaps and telomeric blocks come from the sequence. Sizes are to scale.
(function () {
  const $ = s => document.querySelector(s);
  const COLORS = {gene: '#1f6feb', transcript: '#4c8dff', exon: '#8250df', cds: '#d1242f', unknown: '#57606a',
                  cpg_island: '#1a7f37', repeat: '#bf8700', telomere: '#0e8a8a', gap: '#30363d', ccre: '#2ea043', domain: '#6639ba'};
  const REPEAT_COLORS = {LINE: '#bf8700', SINE: '#e3b341', LTR: '#9a6700', DNA: '#7a5901', Satellite: '#0e8a8a', Simple_repeat: '#d29922', Low_complexity: '#b08800'};
  const CCRE_COLORS = {PLS: '#2ea043', pELS: '#e3b341', dELS: '#d29922', 'CTCF-only': '#39c5c5', 'DNase-H3K4me3': '#4c8dff'};
  const st = {path: null, chrom: null, length: 0, view: [0, 1], data: null, sel: null, filter: new Set(), highlight: null,
              mode: 'view', edits: [], drag: null, hover: null, query: ''};
  const canvas = $('#b-canvas'), ctx = canvas.getContext('2d');
  const mini = $('#b-mini'), mctx = mini.getContext('2d');
  const muted = () => getComputedStyle(document.documentElement).getPropertyValue('--muted').trim();

  function resize() {
    const r = canvas.parentElement.getBoundingClientRect();
    canvas.width = Math.floor(r.width); canvas.height = Math.max(440, Math.floor(window.innerHeight - 420)); canvas.style.height = canvas.height + 'px'; mini.width = Math.floor(r.width); mini.height = 26;
    draw();
  }
  window.addEventListener('resize', resize);
  window.blocksResize = resize;

  const X = pos => (pos - st.view[0]) / (st.view[1] - st.view[0]) * canvas.width;
  const P = x => st.view[0] + x / canvas.width * (st.view[1] - st.view[0]);

  async function load(keepView) {
    const span = st.view[1] - st.view[0];
    const fetchWin = (!keepView || span > 3e6) ? [0, 0] : [Math.max(0, Math.floor(st.view[0] - span)), Math.min(st.length || 1e12, Math.ceil(st.view[1] + span))];
    $('#b-status').textContent = 'loading…';
    try {
      const r = await window.api(`/api/blocks?path=${encodeURIComponent(st.path)}&chrom=${encodeURIComponent(st.chrom)}&start=${fetchWin[0]}&end=${fetchWin[1]}`);
      st.data = r; st.length = r.length;
      if (canvas.width < 10) { const rr = canvas.parentElement.getBoundingClientRect(); if (rr.width > 10) resize(); }
      if (!keepView) st.view = [0, r.length];
      else st.view = [Math.max(0, st.view[0]), Math.min(r.length, st.view[1])];
      applyEdits();
      $('#b-status').textContent = `${st.chrom}: ${Object.entries(r.counts).map(([k, v]) => `${v} ${k}`).join(', ')}${r.coarse ? ' (zoom in for transcripts and exons)' : ''}`;
      layout(); draw(); summary();
    } catch (e) { $('#b-status').textContent = e.message; }
  }

  // ---- layout: pack genes into rows per strand; children drawn inside their parent
  let rows = [];
  function layout() {
    const bs = st.data.blocks;
    const byId = Object.fromEntries(bs.map(b => [b.id, b]));
    st.byId = byId;
    const top = bs.filter(b => b.type === 'gene' || (!b.parent && b.type !== 'gene'));
    // each track: list of rows; each row: array of end positions for packing
    const tracks = {domain: [], plus: [], minus: [], other: []};
    for (const b of top.sort((a, c) => a.start - c.start)) {
      const key = b.type === 'domain' ? 'domain' : b.type !== 'gene' ? 'other' : b.strand === '+' ? 'plus' : 'minus';
      const t = tracks[key];
      let placed = false;
      for (let i = 0; i < t.length; i++) { if (t[i].end <= b.start) { t[i].end = b.end; b._row = i; placed = true; break; } }
      if (!placed) { t.push({end: b.end}); b._row = t.length - 1; }
      b._track = key;
    }
    st.tracks = {domain: tracks.domain.length, plus: tracks.plus.length, minus: tracks.minus.length, other: tracks.other.length};
    // transcripts inside gene rows
    const kids = {};
    for (const b of bs) if (b.parent) (kids[b.parent] = kids[b.parent] || []).push(b);
    st.kids = kids;
    for (const g of top) {
      const ts = (kids[g.id] || []).filter(t => t.type === 'transcript').sort((a, b) => (b.attrs.canonical ? 1 : 0) - (a.attrs.canonical ? 1 : 0) || a.start - b.start);
      ts.forEach((t, i) => { t._sub = i; t._n = ts.length; });
    }
  }

  function visible(b) {
    if (st.filter.size && !st.filter.has(b.type)) return false;
    return b.end > st.view[0] && b.start < st.view[1];
  }
  const classOf = b => (b.type === 'unknown' && b.attrs.class) ? 'unknown:' + b.attrs.class : (b.type === 'ccre' ? 'ccre:' + b.attrs.cls : b.type);
  function isMuted(b) {
    if (st.highlight) return classOf(b) !== st.highlight && !(b.parent && st.highlight === 'gene' && b.type !== 'gene' && false);
    if (st.query) { const q = st.query.toLowerCase(); const hit = (b.name || '').toLowerCase().includes(q) || (st.byId[b.parent]?.name || '').toLowerCase().includes(q) || (st.byId[st.byId[b.parent]?.parent]?.name || '').toLowerCase().includes(q); if (!hit) return true; }
    if (st.sel) { const s = st.byId[st.sel]; const chain = new Set([b.id, b.parent, st.byId[b.parent]?.parent]); const related = chain.has(s.id) || s.parent === b.id || st.byId[s.parent]?.parent === b.id || b.parent === s.id || st.byId[b.parent]?.parent === s.id; if (!related) return true; }
    return false;
  }

  function draw() {
    if (!st.data) return;
    const W = canvas.width, H = canvas.height;
    ctx.clearRect(0, 0, W, H);
    const line = getComputedStyle(document.documentElement).getPropertyValue('--line').trim();
    // ruler
    ctx.fillStyle = muted(); ctx.font = '11px system-ui';
    const span = st.view[1] - st.view[0];
    let step = Math.pow(10, Math.floor(Math.log10(span / 6)));
    while (step / span * W < 90) step *= step / span * W < 36 ? 5 : 2;
    for (let p = Math.ceil(st.view[0] / step) * step; p < st.view[1]; p += step) {
      const x = X(p); ctx.fillStyle = line; ctx.fillRect(x, 0, 1, H); ctx.fillStyle = muted(); ctx.fillText(fmtPos(p), x + 3, 12);
    }
    // use the whole height: split it among the packed rows of the three tracks
    const rowsPlus = Math.max(1, st.tracks.plus), rowsMinus = Math.max(1, st.tracks.minus), rowsOther = Math.max(1, st.tracks.other);
    const totalRows = rowsPlus + rowsMinus + rowsOther;
    const gapY = 4, labels = 3 * 18 + 22;
    const geneH = Math.max(10, Math.min(120, (H - labels) / totalRows - gapY));
    const domainH = 12, rowsDomain = st.tracks.domain;
    const yDomain = 22;
    const yPlus = 22 + (rowsDomain ? rowsDomain * (domainH + gapY) + 18 : 0);
    const yMinus = yPlus + rowsPlus * (geneH + gapY) + 18;
    const yOther = yMinus + rowsMinus * (geneH + gapY) + 18;
    ctx.fillStyle = muted(); if (rowsDomain) ctx.fillText('nodes: domains between CTCF boundaries (inferred)', 4, yDomain - 4); ctx.fillText('+ strand', 4, yPlus - 4); ctx.fillText('− strand', 4, yMinus - 4); ctx.fillText('sequence elements & UNKNOWN', 4, yOther - 4);
    st.hitboxes = [];
    const bs = st.data.blocks;
    for (const b of bs) {
      if (b.parent || !visible(b)) continue;
      const y = b._track === 'domain' ? yDomain + b._row * (domainH + gapY) : b._track === 'plus' ? yPlus + b._row * (geneH + gapY) : b._track === 'minus' ? yMinus + b._row * (geneH + gapY) : yOther + b._row * (geneH + gapY);
      drawBlock(b, y, b._track === 'domain' ? domainH : geneH, span);
    }
    if (st.drag) { ctx.fillStyle = 'rgba(255,255,255,.08)'; ctx.fillRect(0, 0, W, H); const b = st.drag.block; const x = X(st.drag.newStart), w = Math.max(2, (b.end - b.start) / span * W); ctx.strokeStyle = st.drag.ok ? '#3fb950' : '#ff6b66'; ctx.lineWidth = 2; ctx.strokeRect(x, st.drag.y, w, geneH); ctx.fillStyle = st.drag.ok ? '#3fb950' : '#ff6b66'; ctx.fillText(`${b.name} → ${fmtPos(st.drag.newStart)} ${st.drag.ok ? '' : '(' + st.drag.reasons.join('; ') + ')'}`, Math.min(x, W - 300), st.drag.y - 4); }
    drawMini();
  }

  function drawBlock(b, y, h, span) {
    const W = canvas.width;
    const x0 = Math.max(-2, X(b.start)), x1 = Math.min(W + 2, X(b.end)), w = Math.max(1.5, x1 - x0);
    const m = isMuted(b);
    ctx.globalAlpha = m ? 0.18 : 1;
    let col = COLORS[b.type] || '#999';
    if (b.type === 'gene') {
      ctx.fillStyle = col; ctx.globalAlpha = m ? 0.12 : 0.22; ctx.fillRect(x0, y, w, h); ctx.globalAlpha = m ? 0.18 : 1;
      ctx.strokeStyle = b.id === st.sel ? '#fff' : col; ctx.lineWidth = b.id === st.sel ? 2 : 1; ctx.strokeRect(x0 + 0.5, y + 0.5, w - 1, h - 1);
      if (w > 30) { ctx.fillStyle = m ? muted() : '#fff'; ctx.font = 'bold 11px system-ui'; ctx.fillText(b.name + (b.attrs.gene_type !== 'protein_coding' ? ' · ' + b.attrs.gene_type : ''), x0 + 4, y + 11); }
      // transcripts inside when zoomed
      const ts = (st.kids[b.id] || []).filter(t => t.type === 'transcript');
      if (span < 2e6 && ts.length && h >= 22) {
        const showN = Math.min(ts.length, Math.max(1, Math.floor((h - 14) / 5)));
        const subH = (h - 14) / showN;
        for (const t of ts) {
          if (t._sub >= showN) break;
          const ty = y + 13 + t._sub * subH;
          if (!visible(t)) continue;
          const tm = isMuted(t);
          ctx.globalAlpha = tm ? 0.15 : 0.8;
          ctx.strokeStyle = COLORS.transcript; ctx.lineWidth = 1;
          ctx.beginPath(); ctx.moveTo(Math.max(0, X(t.start)), ty + subH / 2); ctx.lineTo(Math.min(W, X(t.end)), ty + subH / 2); ctx.stroke();
          for (const e of (st.kids[t.id] || [])) {
            if (!visible(e)) continue;
            const ex0 = X(e.start), ew = Math.max(1, X(e.end) - ex0);
            ctx.fillStyle = COLORS[e.type]; ctx.globalAlpha = isMuted(e) ? 0.15 : (e.type === 'cds' ? 1 : 0.6);
            ctx.fillRect(ex0, ty + (e.type === 'cds' ? 0 : subH * 0.25), ew, e.type === 'cds' ? subH : subH * 0.5);
            st.hitboxes.push({b: e, x: ex0, y: ty, w: ew, h: subH});
          }
          st.hitboxes.push({b: t, x: Math.max(0, X(t.start)), y: ty, w: Math.min(W, X(t.end)) - Math.max(0, X(t.start)), h: subH});
        }
      }
    } else {
      if (b.type === 'ccre') col = CCRE_COLORS[b.attrs.cls] || col;
      if (b.type === 'repeat' && b.attrs.cls) col = REPEAT_COLORS[b.attrs.cls] || '#bf8700';
      const ucol = {interspersed_repeat: '#bf8700', interspersed_repeat_SINE: '#bf8700', tandem_repeat: '#e3b341', low_complexity: '#9a6700', long_orf: '#d1242f', satellite_array: '#7a5901', mixed_intergenic: '#6e7681', unique_intergenic: '#8b949e', promoter_like: '#1a7f37', centromere: '#0e8a8a', telomere: '#0e8a8a', gene_desert: '#3d444d'};
      ctx.fillStyle = (b.type === 'unknown' && b.attrs.class && ucol[b.attrs.class]) || col; ctx.globalAlpha = m ? 0.15 : (b.type === 'unknown' ? (b.attrs.class ? 0.6 : 0.35) : 0.9);
      ctx.fillRect(x0, y + (b.type === 'unknown' ? h * 0.3 : 0), w, b.type === 'unknown' ? h * 0.4 : h);
      if (b.id === st.sel) { ctx.strokeStyle = '#fff'; ctx.lineWidth = 2; ctx.strokeRect(x0, y, w, h); }
      if (w > 40 && h >= 20) { ctx.fillStyle = m ? muted() : '#fff'; ctx.font = '10px system-ui'; ctx.fillText(b.name, x0 + 3, y + 13); }
    }
    ctx.globalAlpha = 1;
    st.hitboxes.push({b, x: x0, y, w, h});
  }

  // coverage of the loaded window by block type (union of intervals, so overlapping genes count once)
  function summary() {
    const d = st.data; if (!d) return;
    classSummary();
    const span = d.end - d.start;
    const types = ['domain', 'gene', 'transcript', 'exon', 'cds', 'unknown', 'ccre', 'cpg_island', 'repeat', 'telomere', 'gap'];
    const rows = [];
    for (const t of types) {
      const iv = d.blocks.filter(b => b.type === t).map(b => [Math.max(d.start, b.start), Math.min(d.end, b.end)]).filter(([a, b]) => b > a).sort((a, b) => a[0] - b[0]);
      let bp = 0, cur = null;
      for (const [a, b] of iv) { if (!cur || a > cur[1]) { if (cur) bp += cur[1] - cur[0]; cur = [a, b]; } else cur[1] = Math.max(cur[1], b); }
      if (cur) bp += cur[1] - cur[0];
      if (iv.length) rows.push({t, n: iv.length, bp, pct: bp / span * 100});
    }
    const note = d.coarse ? ' · zoom under 3 Mb to count transcripts, exons, CDS, repeats and islands' : '';
    $('#b-summary').innerHTML = `<div class="muted" style="font-size:12px;margin-bottom:4px">share of bases in ${st.chrom}:${d.start.toLocaleString()}-${d.end.toLocaleString()} (${fmtPos(span)}) by block type and UNKNOWN class; overlaps counted once; click a chip to highlight${note}</div>`;
  }

  // summary of block types and UNKNOWN classes in the loaded window; click a chip to highlight all of that kind
  function classSummary() {
    const d = st.data; const span = d.end - d.start;
    const groups = {};
    for (const b of d.blocks) { const a = Math.max(d.start, b.start), e = Math.min(d.end, b.end); if (e <= a) continue; const k = classOf(b); const g = groups[k] = groups[k] || {n: 0, bp: 0, iv: []}; g.n++; g.iv.push([a, e]); }
    for (const g of Object.values(groups)) { g.iv.sort((a, b) => a[0] - b[0]); let cur = null, bp = 0; for (const [a, b] of g.iv) { if (!cur || a > cur[1]) { if (cur) bp += cur[1] - cur[0]; cur = [a, b]; } else cur[1] = Math.max(cur[1], b); } if (cur) bp += cur[1] - cur[0]; g.bp = bp; }
    const keys = Object.keys(groups).sort((a, b) => groups[b].bp - groups[a].bp);
    const ucol = {interspersed_repeat: '#bf8700', interspersed_repeat_SINE: '#bf8700', tandem_repeat: '#e3b341', low_complexity: '#9a6700', long_orf: '#d1242f', satellite_array: '#7a5901', mixed_intergenic: '#6e7681', unique_intergenic: '#8b949e', promoter_like: '#1a7f37', centromere: '#0e8a8a', telomere: '#0e8a8a', gene_desert: '#3d444d', regulatory: '#2ea043'};
    const colorFor = k => k.startsWith('unknown:') ? (ucol[k.slice(8)] || '#57606a') : k.startsWith('ccre:') ? (CCRE_COLORS[k.slice(5)] || '#2ea043') : (COLORS[k] || '#999');
    const label = k => k.startsWith('unknown:') ? 'UNKNOWN · ' + k.slice(8) : k.startsWith('ccre:') ? 'ENCODE ' + k.slice(5) : k;
    $('#b-classes').innerHTML = keys.map(k => `<span class="chip${st.highlight === k ? ' on' : ''}" data-key="${k}" data-tip="Click to highlight every ${label(k)} block and mute the rest; click again to clear."><i style="background:${colorFor(k)}"></i>${label(k)} <b>${(groups[k].bp / span * 100).toFixed(groups[k].bp / span >= 0.1 ? 0 : 1)}%</b> <span class="muted">(${groups[k].n})</span></span>`).join('');
    $('#b-classes').querySelectorAll('.chip').forEach(c => c.onclick = () => { st.highlight = st.highlight === c.dataset.key ? null : c.dataset.key; classSummary(); draw(); });
  }

  function drawMini() {
    const W = mini.width, H = mini.height;
    mctx.clearRect(0, 0, W, H);
    mctx.fillStyle = getComputedStyle(document.documentElement).getPropertyValue('--line').trim(); mctx.fillRect(0, H / 2 - 2, W, 4);
    for (const b of st.data.blocks) if (b.type === 'gene') { mctx.fillStyle = COLORS.gene; mctx.fillRect(b.start / st.length * W, H / 2 - 5, Math.max(1, (b.end - b.start) / st.length * W), 10); }
    mctx.strokeStyle = '#fff'; mctx.lineWidth = 1.5; mctx.strokeRect(st.view[0] / st.length * W, 2, Math.max(2, (st.view[1] - st.view[0]) / st.length * W), H - 4);
  }
  const fmtPos = p => p >= 1e6 ? (p / 1e6).toFixed(2) + ' Mb' : p >= 1e3 ? (p / 1e3).toFixed(1) + ' kb' : p + ' bp';

  // ---- interaction
  function hit(x, y) { for (let i = st.hitboxes.length - 1; i >= 0; i--) { const h = st.hitboxes[i]; if (x >= h.x && x <= h.x + h.w && y >= h.y && y <= h.y + h.h) return h; } return null; }
  canvas.addEventListener('wheel', e => {
    e.preventDefault();
    const p = P(e.offsetX), f = e.deltaY > 0 ? 1.25 : 0.8;
    let a = p - (p - st.view[0]) * f, b = p + (st.view[1] - p) * f;
    a = Math.max(0, a); b = Math.min(st.length, b); if (b - a < 200) return;
    st.view = [a, b]; scheduleLoad();
  }, {passive: false});
  let pan = null;
  canvas.addEventListener('mousedown', e => {
    const h = hit(e.offsetX, e.offsetY);
    if (st.mode === 'edit' && h && h.b.type !== 'unknown') { st.drag = {block: h.b, y: h.y, x0: e.offsetX, newStart: h.b.start, ok: true, reasons: []}; return; }
    pan = {x: e.offsetX, view: [...st.view]};
  });
  canvas.addEventListener('mousemove', async e => {
    if (st.drag) {
      const span = st.view[1] - st.view[0];
      st.drag.newStart = Math.round(st.drag.block.start + (e.offsetX - st.drag.x0) / canvas.width * span);
      draw(); return;
    }
    if (pan) { const span = pan.view[1] - pan.view[0]; const d = (e.offsetX - pan.x) / canvas.width * span; let a = pan.view[0] - d, b = pan.view[1] - d; if (a < 0) { b -= a; a = 0; } if (b > st.length) { a -= b - st.length; b = st.length; } st.view = [a, b]; draw(); return; }
    const h = hit(e.offsetX, e.offsetY);
    canvas.style.cursor = h ? (st.mode === 'edit' ? 'grab' : 'pointer') : 'default';
    if (h !== st.hover) { st.hover = h; showInfo(h ? h.b : (st.sel ? st.byId[st.sel] : null), !!h && h.b.id !== st.sel); }
  });
  window.addEventListener('mouseup', async e => {
    if (st.drag) {
      const d = st.drag; st.drag = null;
      if (Math.abs(d.newStart - d.block.start) >= 1) {
        const r = await window.api('/api/blocks/check', {blocks: st.data.blocks, id: d.block.id, new_start: d.newStart, chrom_length: st.length});
        if (r.ok) { st.edits.push({id: d.block.id, name: d.block.name, type: d.block.type, from: d.block.start, to: d.newStart, delta: r.delta}); applyEdits(); layout(); renderEdits(); }
        else { $('#b-verdict').innerHTML = `<span class="err">not plausible: ${r.reasons.join('; ')}</span>`; }
      }
      draw(); return;
    }
    if (pan) { const moved = Math.abs(e.offsetX - pan.x) > 3; pan = null; if (moved) { scheduleLoad(); return; } }
    if (e.target !== canvas) return;
    const h = hit(e.offsetX, e.offsetY);
    st.sel = h ? h.b.id : null; showInfo(h ? h.b : null, false); draw();
  });
  mini.addEventListener('mousedown', e => { const span = st.view[1] - st.view[0]; const c = e.offsetX / mini.width * st.length; st.view = [Math.max(0, c - span / 2), Math.min(st.length, c + span / 2)]; scheduleLoad(); });
  let loadT = null;
  function scheduleLoad() { draw(); clearTimeout(loadT); loadT = setTimeout(() => { const d = st.data; const need = !d || d.coarse !== ((st.view[1] - st.view[0]) > 3e6) || st.view[0] < d.start || st.view[1] > d.end; if (need) load(true); else draw(); }, 150); }

  function applyEdits() {
    if (!st.data) return;
    const shift = (id, delta, byId, kids) => { const b = byId[id]; b.start += delta; b.end += delta; (kids[id] || []).forEach(k => shift(k.id, delta, byId, kids)); };
    const byId = Object.fromEntries(st.data.blocks.map(b => [b.id, b]));
    const kids = {}; for (const b of st.data.blocks) if (b.parent) (kids[b.parent] = kids[b.parent] || []).push(b);
    for (const ed of st.edits) if (byId[ed.id] && byId[ed.id].start === ed.from) shift(ed.id, ed.to - ed.from, byId, kids);
    for (const ed of st.edits) if (byId[ed.id]) { byId[ed.id].evidence = 'predicted'; byId[ed.id].confidence = 0.3; byId[ed.id].attrs.edited = true; }
  }
  function renderEdits() {
    $('#b-edits').innerHTML = st.edits.length ? '<table><tr><th>block</th><th>type</th><th class="num">from</th><th class="num">to</th><th class="num">shift</th></tr>' + st.edits.map(e => `<tr><td class="mono">${e.name}</td><td>${e.type}</td><td class="num">${fmtPos(e.from)}</td><td class="num">${fmtPos(e.to)}</td><td class="num">${e.delta >= 0 ? '+' : ''}${fmtPos(Math.abs(e.delta)).replace(/^/, e.delta < 0 ? '-' : '')}</td></tr>`).join('') + '</table>' : '<span class="muted">no edits yet</span>';
    $('#b-verdict').innerHTML = st.edits.length ? `<span class="pill predicted">predicted</span> <span class="muted">${st.edits.length} proposed move(s), confidence 0.3 — a design, not an observation</span>` : '';
  }
  function showInfo(b, hovering) {
    if (!b) { $('#b-info').innerHTML = '<span class="muted">Click a block. Hover shows a preview.</span>'; return; }
    const par = b.parent ? st.byId[b.parent] : null;
    const kids = (st.kids[b.id] || []);
    $('#b-info').innerHTML = `<div><b>${b.name}</b> <span class="pill ${b.evidence}">${b.evidence}</span> <span class="muted">conf ${b.confidence}</span>${hovering ? ' <span class="muted">(hover)</span>' : ''}</div>
      <div class="mono" style="font-size:12px">${b.type} · ${st.chrom}:${b.start.toLocaleString()}-${b.end.toLocaleString()} ${b.strand} · ${fmtPos(b.end - b.start)}</div>
      ${par ? `<div class="muted">inside ${par.type} <b>${par.name}</b></div>` : ''}
      ${kids.length ? `<div class="muted">contains ${Object.entries(kids.reduce((a, k) => (a[k.type] = (a[k.type] || 0) + 1, a), {})).map(([k, v]) => `${v} ${k}`).join(', ')}</div>` : ''}
      ${Object.keys(b.attrs).length ? `<div class="muted" style="font-size:12px">${Object.entries(b.attrs).map(([k, v]) => `${k}=${v}`).join(' · ')}</div>` : ''}
      ${b.type === 'unknown' && b.attrs.class ? `<div class="hint">Investigated: <b>${b.attrs.class}</b> (${b.evidence}, conf ${b.confidence}). ${b.attrs.patterns ? 'patterns: ' + b.attrs.patterns + '. ' : ''}${b.attrs.similar_to ? 'similar to ' + b.attrs.similar_to : ''}</div>` : ''}
      ${b.type === 'unknown' && !b.attrs.class ? '<div class="hint">UNKNOWN: no annotated gene here. In a human chromosome this space carries most of the regulation (enhancers, insulators) that we cannot yet read from sequence alone.</div>' : ''}
      <div class="row" style="margin-top:6px"><button class="ghost" id="b-zoomto"><i class="ic">🔍</i>Zoom to</button><button class="ghost" id="b-hl"><i class="ic">💡</i>${st.highlight === classOf(b) ? 'Clear highlight' : 'Highlight all ' + (b.type === 'unknown' && b.attrs.class ? b.attrs.class : b.type)}</button>${(b.type === 'gene' || par) ? `<button class="ghost" id="b-flow" data-tip="Open this gene in Flow: DNA → RNA → protein with its regulation"><i class="ic">🔁</i>Flow</button><button class="ghost" id="b-mol" data-tip="Open this gene in Molecules: isoforms, expression, protein definition, structure, graph"><i class="ic">🧫</i>Molecules</button>` : ''}</div>`;
    $('#b-zoomto').onclick = () => { const pad = (b.end - b.start) * 0.15; st.view = [Math.max(0, b.start - pad), Math.min(st.length, b.end + pad)]; scheduleLoad(); };
    const geneName = b.type === 'gene' ? b.name : (par && par.type === 'gene' ? par.name : (par && par.parent && st.byId[par.parent] ? st.byId[par.parent].name : null));
    const go = (tab, fill) => { fill(); const t = document.querySelector(`nav button[data-tab="${tab}"]`); if (t) t.click(); };
    if ($('#b-flow') && geneName) $('#b-flow').onclick = () => go('flow', () => { $('#f-gene').value = geneName; $('#f-chrom').value = st.chrom; setTimeout(() => $('#f-load').click(), 200); });
    if ($('#b-mol') && geneName) $('#b-mol').onclick = () => go('molecules', () => { $('#m-gene').value = geneName; $('#m-chrom').value = st.chrom; setTimeout(() => $('#m-load').click(), 200); });
    $('#b-hl').onclick = () => { st.highlight = st.highlight === classOf(b) ? null : classOf(b); classSummary(); draw(); showInfo(b, false); };
  }

  // ---- controls
  window.blocksInit = async function (files) {
    $('#b-file').innerHTML = files.genomes.map(g => `<option value="${g.path}">${g.path}</option>`).join('');
    const pref = files.genomes.find(g => g.path.includes('chr21')) || files.genomes.find(g => g.path.includes('chrM')) || files.genomes[0];
    if (pref) { $('#b-file').value = pref.path; $('#b-chrom').value = pref.path.includes('chr21') ? 'chr21' : pref.path.includes('chrM') ? 'chrM' : ''; }
    $('#b-load').onclick = () => { st.path = $('#b-file').value; st.chrom = $('#b-chrom').value; st.sel = null; st.edits = []; renderEdits(); load(false); };
    $('#b-search').oninput = e => { st.query = e.target.value.trim(); draw(); };
    $('#b-search').onkeydown = e => { if (e.key === 'Enter' && st.data) { const q = st.query.toLowerCase(); const g = st.data.blocks.find(b => b.type === 'gene' && b.name.toLowerCase() === q) || st.data.blocks.find(b => b.name.toLowerCase().includes(q)); if (g) { st.sel = g.id; const pad = (g.end - g.start) * 0.5; st.view = [Math.max(0, g.start - pad), Math.min(st.length, g.end + pad)]; showInfo(g, false); scheduleLoad(); } } };
    document.querySelectorAll('#b-filters input').forEach(cb => cb.onchange = () => { st.filter = new Set([...document.querySelectorAll('#b-filters input:checked')].map(x => x.value)); if (st.filter.size === document.querySelectorAll('#b-filters input').length) st.filter.clear(); draw(); });
    $('#b-mode').onchange = () => { st.mode = $('#b-mode').checked ? 'edit' : 'view'; $('#b-modelabel').textContent = st.mode === 'edit' ? 'edit mode: drag a block to propose a move' : 'view mode'; };
    $('#b-reset').onclick = () => { st.edits = []; renderEdits(); load(true); };
    $('#b-export').onclick = () => { const out = {chrom: st.chrom, edits: st.edits, evidence: 'predicted', confidence: 0.3, note: 'block layout proposal from the GenomeOS block editor'}; $('#b-export-out').textContent = JSON.stringify(out, null, 1); };
    $('#b-zoomout').onclick = () => { st.view = [0, st.length]; scheduleLoad(); };
    resize(); renderEdits();
    // deep link: #blocks?locus=chr21:25800000-26300000
    const m = ((window.initialHash || location.hash).split('?')[1] || '').match(/locus=([^:&]+):(\d+)-(\d+)/);
    if (m) { $('#b-chrom').value = m[1]; st.path = $('#b-file').value; st.chrom = m[1]; st.view = [+m[2], +m[3]]; st.length = 1e12; load(true).then(() => { st.view = [+m[2], +m[3]]; scheduleLoad(); }); }
    else if (pref) $('#b-load').click();
  };
})();
