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
  window.moleculesLoad = async function (gene, chrom) {
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
