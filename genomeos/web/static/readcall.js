// The reader's call on one gene in every biosample (/api/blocks/readcall), for the Blocks lane's Block card.
// "Read" is the promoter accessible or marked, and the card says nothing more about the gene than that.
// A call made by the histone marks (no DNase peak; H3K4me3 and H3K27ac over the promoter) is told apart
// from a call on the DNase file, and UNKNOWN is shown with its reason, never as closed.
// Loaded by blocks.js; window.readCallHtml is pure, so a test can run it in node.
(function () {
  const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c]));
  // glyph and colour match the dot blocks.js draws on each gene for the selected reader biosample
  const LOOK = {
    open: ['●', '#56d364', 'read · DNase'],
    read_by_marks: ['◐', '#39c5c5', 'read · by marks'],
    poised: ['●', '#d29922', 'poised · not read'],
    closed: ['○', '#8b949e', 'closed · not read'],
  };

  function row(c, d, selected) {
    // the biosample the reader selector shows on the map; index.html is not touched, so the mark is inline
    const sel = selected && c.result === `reader_${selected}_${d.chrom}`;
    const mark = sel ? ' style="background:rgba(127,127,127,.18)"' : '';
    const who = esc(c.biosample) + (c.in_reader_terms ? '' : ' <span class="muted" style="font-size:11px" data-tip="A reading the reader holds outside the thirteen biosamples compile and readergate gate on (READER_TERMS).">outside the 13 reader cells</span>');
    if (c.call === 'unknown') {
      return `<tr class="rc-row${sel ? ' rc-sel' : ''}"${mark} data-call="unknown"><td>${who}</td><td><span class="pill none">UNKNOWN</span></td><td class="muted">no call</td>` +
        `<td><span class="mono rc-reason">${esc(c.reason)}</span> <span class="muted">${esc(c.reason_text)}</span></td></tr>`;
    }
    const [glyph, colour, label] = LOOK[c.call] || ['?', '#8b949e', esc(c.call)];
    const detail = [
      c.call === 'open' && c.signal != null ? `DNase signal ${esc(c.signal)}` : '',
      c.call === 'open' && c.by_absence ? 'on none of the closed, poised or by-marks lists' : '',
      c.call === 'read_by_marks' ? 'no DNase peak; H3K4me3 and H3K27ac peaks' : '',
      c.call === 'poised' ? 'H3K27me3 without H3K27ac over an open promoter' : '',
      c.note ? esc(c.note) : '',
    ].filter(Boolean).join(' · ');
    return `<tr class="rc-row${sel ? ' rc-sel' : ''}"${mark} data-call="${esc(c.call)}"><td>${who}</td>` +
      `<td><span style="color:${colour}">${glyph}</span> ${label}</td><td>${esc(c.basis)}</td><td class="muted">${detail}</td></tr>`;
  }

  function readCallHtml(d, selected) {
    if (!d || !Array.isArray(d.calls)) return '';
    const n = d.counts || {};
    const head = d.calls.length
      ? `read in <b>${d.read}</b> of ${d.biosamples} biosamples (<b>${n.open || 0}</b> on DNase, <b>${n.read_by_marks || 0}</b> by marks) · ${n.poised || 0} poised · ${n.closed || 0} closed · ${n.unknown || 0} UNKNOWN`
      : 'the reader holds no biosample on this chromosome';
    return `<div class="rc" style="margin-top:6px">
      <div><b>Read by biosample</b> <span class="muted" style="font-size:12px">promoter (TSS ± 1 kb) accessible or marked</span></div>
      <div style="font-size:12px">${head}</div>
      <table class="rc-table" style="font-size:12px"><tr><th>biosample</th><th>call</th><th>rests on</th><th></th></tr>${d.calls.map(c => row(c, d, selected)).join('')}</table>
      <div class="muted" style="font-size:11px">${esc(d.evidence)}</div>
    </div>`;
  }

  window.readCallHtml = readCallHtml;
})();
