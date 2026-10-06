(function () {
  const $ = (id) => document.getElementById(id);
  const root = $('pdf-tools-root');
  if (!root) return;

  const pdfjsLib = window.pdfjsLib;
  const PDFLib = window.PDFLib;
  if (!pdfjsLib || !PDFLib) {
    root.insertAdjacentHTML('afterbegin', '<p class="digitizer-note">PDF libraries failed to load. Check your internet connection.</p>');
    return;
  }
  pdfjsLib.GlobalWorkerOptions.workerSrc = 'https://cdn.jsdelivr.net/npm/pdfjs-dist@3.11.174/build/pdf.worker.min.js';
  const { PDFDocument, StandardFonts, rgb } = PDFLib;

  const FONT_MAP = {
    Arial: { css: 'Helvetica, Arial, sans-serif', pdf: StandardFonts.Helvetica },
    'Times New Roman': { css: '"Times New Roman", Times, serif', pdf: StandardFonts.TimesRoman },
    'Courier New': { css: '"Courier New", Courier, monospace', pdf: StandardFonts.Courier },
  };
  const SIGNATURE_FONTS = {
    'Handwritten': '"Brush Script MT", "Segoe Script", "Lucida Handwriting", cursive',
    'Classic italic': 'italic "Times New Roman", serif',
    'Casual': 'italic "Comic Sans MS", "Comic Neue", cursive',
    'Modern italic': 'italic Arial, sans-serif',
    'Bold italic': 'italic bold Arial, sans-serif',
  };

  function setStatus(message, isError) {
    const box = $('pdf-status');
    box.innerHTML = '<strong>' + (isError ? 'Error' : 'Status') + '</strong><span></span>';
    box.lastChild.textContent = message;
  }

  function download(bytes, name, mime) {
    const url = URL.createObjectURL(new Blob([bytes], { type: mime }));
    const a = document.createElement('a');
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 5000);
  }

  const readFile = (file) => file.arrayBuffer().then((b) => new Uint8Array(b));
  const baseName = (file) => file.name.replace(/\.[^.]+$/, '');
  const hexToRgb = (hex) => rgb(
    parseInt(hex.slice(1, 3), 16) / 255,
    parseInt(hex.slice(3, 5), 16) / 255,
    parseInt(hex.slice(5, 7), 16) / 255
  );

  function parsePageGroups(text, count) {
    const groups = [];
    text.split(',').forEach((part) => {
      part = part.trim();
      if (!part) return;
      const m = /^(\d+)\s*(?:-\s*(\d+))?$/.exec(part);
      if (!m) throw new Error('Use page numbers such as 1, 3-5, 8.');
      const start = Number(m[1]);
      const end = Number(m[2] || m[1]);
      if (start < 1 || end < start || end > count) throw new Error('Page numbers must be between 1 and ' + count + '.');
      const group = [];
      for (let p = start; p <= end; p++) group.push(p);
      groups.push(group);
    });
    if (!groups.length) throw new Error('Enter at least one page number.');
    return groups;
  }

  function parsePages(text, count) {
    return Array.from(new Set(parsePageGroups(text, count).flat()));
  }

  // Sub-tabs
  const subTabs = root.querySelectorAll('[data-pdf-tab]');
  const subPanels = root.querySelectorAll('[data-pdf-panel]');
  subTabs.forEach((btn) => btn.addEventListener('click', () => {
    subTabs.forEach((b) => b.classList.toggle('is-active', b === btn));
    subPanels.forEach((p) => { p.hidden = p.dataset.pdfPanel !== btn.dataset.pdfTab; });
  }));

  // ---------- Edit ----------
  const edit = { bytes: null, doc: null, file: null, page: 1, pageCount: 0, scale: 1, pageW: 0, pageH: 0, items: [], selected: null, signature: null, stamp: null, seq: 0 };
  const stageEl = $('pdf-edit-stage');
  const pageCanvas = $('pdf-edit-canvas');
  const overlay = $('pdf-edit-overlay');

  async function renderPage() {
    const page = await edit.doc.getPage(edit.page);
    const base = page.getViewport({ scale: 1 });
    edit.pageW = base.width;
    edit.pageH = base.height;
    const maxW = Math.min(stageEl.clientWidth || 760, 900) - 4;
    edit.scale = maxW / base.width;
    const viewport = page.getViewport({ scale: edit.scale });
    const ratio = window.devicePixelRatio || 1;
    pageCanvas.width = viewport.width * ratio;
    pageCanvas.height = viewport.height * ratio;
    pageCanvas.style.width = viewport.width + 'px';
    pageCanvas.style.height = viewport.height + 'px';
    overlay.style.width = viewport.width + 'px';
    overlay.style.height = viewport.height + 'px';
    const ctx = pageCanvas.getContext('2d');
    await page.render({ canvasContext: ctx, viewport, transform: [ratio, 0, 0, ratio, 0, 0] }).promise;
    $('pdf-edit-page-label').textContent = 'Page ' + edit.page + ' of ' + edit.pageCount;
    refreshItems();
  }

  function refreshItems() {
    overlay.innerHTML = '';
    edit.items.filter((it) => it.page === edit.page).forEach((it) => overlay.appendChild(it.el));
  }

  const sizeInput = (kind) => $(kind === 'date' ? 'pdf-date-size' : 'pdf-font-size');

  function select(item) {
    edit.selected = item;
    edit.items.forEach((it) => it.el.classList.toggle('is-selected', it === item));
    if (item && item.text !== undefined) {
      sizeInput(item.kind).value = Math.round(item.fontSize * 10) / 10;
      $('pdf-font-family').value = item.family;
      $('pdf-text-color').value = item.color;
      if (item.kind === 'text') $('pdf-text').value = item.text;
      else if (/^\d{4}-\d{2}-\d{2}$/.test(item.text)) $('pdf-date').value = item.text;
    }
  }

  function setItemText(item, text) {
    item.text = text;
    item.el.firstChild.textContent = text;
  }

  function addItem(item) {
    item.id = ++edit.seq;
    item.page = edit.page;
    item.el.className = 'pdf-item pdf-item-' + item.kind;
    item.el.style.left = item.x + 'px';
    item.el.style.top = item.y + 'px';
    item.el.tabIndex = 0;
    item.el.addEventListener('pointerdown', (e) => {
      if (e.target.classList.contains('pdf-handle') || item.editing) return;
      select(item);
      const startX = e.clientX, startY = e.clientY, ox = item.x, oy = item.y;
      item.el.setPointerCapture(e.pointerId);
      const move = (ev) => {
        item.x = Math.max(0, ox + ev.clientX - startX);
        item.y = Math.max(0, oy + ev.clientY - startY);
        item.el.style.left = item.x + 'px';
        item.el.style.top = item.y + 'px';
      };
      const up = () => {
        item.el.removeEventListener('pointermove', move);
        item.el.removeEventListener('pointerup', up);
      };
      item.el.addEventListener('pointermove', move);
      item.el.addEventListener('pointerup', up);
      e.preventDefault();
    });
    item.el.addEventListener('keydown', (e) => {
      if (item.editing) return;
      if (e.key === 'Delete' || e.key === 'Backspace') { removeItem(item); e.preventDefault(); }
    });
    if (item.text !== undefined) {
      const label = item.el.firstChild;
      item.el.addEventListener('dblclick', () => {
        item.editing = true;
        label.contentEditable = 'true';
        label.focus();
      });
      label.addEventListener('blur', () => {
        item.editing = false;
        label.contentEditable = 'false';
        item.text = label.innerText.replace(/\n$/, '');
        if (edit.selected === item && item.kind === 'text') $('pdf-text').value = item.text;
      });
      const handle = document.createElement('span');
      handle.className = 'pdf-handle';
      item.el.appendChild(handle);
      handle.addEventListener('pointerdown', (e) => {
        e.stopPropagation();
        e.preventDefault();
        select(item);
        const startX = e.clientX, size0 = item.fontSize, w0 = item.el.offsetWidth;
        handle.setPointerCapture(e.pointerId);
        const move = (ev) => {
          item.fontSize = Math.min(120, Math.max(6, size0 * (w0 + ev.clientX - startX) / w0));
          styleTextItem(item);
          sizeInput(item.kind).value = Math.round(item.fontSize * 10) / 10;
        };
        const up = () => {
          handle.removeEventListener('pointermove', move);
          handle.removeEventListener('pointerup', up);
        };
        handle.addEventListener('pointermove', move);
        handle.addEventListener('pointerup', up);
      });
    }
    if (item.kind === 'signature' || item.kind === 'stamp') {
      const handle = document.createElement('span');
      handle.className = 'pdf-handle';
      item.el.appendChild(handle);
      handle.addEventListener('pointerdown', (e) => {
        e.stopPropagation();
        e.preventDefault();
        select(item);
        const startX = e.clientX, w0 = item.w, ratio = item.h / item.w;
        handle.setPointerCapture(e.pointerId);
        const move = (ev) => {
          item.w = Math.max(20, w0 + ev.clientX - startX);
          item.h = item.w * ratio;
          item.el.style.width = item.w + 'px';
          item.el.style.height = item.h + 'px';
        };
        const up = () => {
          handle.removeEventListener('pointermove', move);
          handle.removeEventListener('pointerup', up);
        };
        handle.addEventListener('pointermove', move);
        handle.addEventListener('pointerup', up);
      });
    }
    edit.items.push(item);
    refreshItems();
    select(item);
    item.el.focus();
  }

  function removeItem(item) {
    edit.items = edit.items.filter((it) => it !== item);
    if (edit.selected === item) edit.selected = null;
    refreshItems();
  }

  function styleTextItem(item) {
    const fam = FONT_MAP[item.family];
    item.el.style.fontFamily = fam.css;
    item.el.style.fontSize = item.fontSize * edit.scale + 'px';
    item.el.style.color = item.color;
    item.el.style.lineHeight = '1.2';
  }

  function addText(text, kind) {
    if (!edit.doc) return setStatus('Choose a PDF first.', true);
    if (!text.trim()) return setStatus('Enter some text first.', true);
    const item = {
      kind, text, x: edit.pageW * edit.scale * 0.08, y: edit.pageH * edit.scale * (kind === 'date' ? 0.9 : 0.12),
      fontSize: Number(sizeInput(kind).value) || 12, family: $('pdf-font-family').value, color: $('pdf-text-color').value,
      el: document.createElement('div'),
    };
    item.el.appendChild(document.createElement('span')).textContent = text;
    styleTextItem(item);
    addItem(item);
  }

  function addImage(kind, src, naturalW, naturalH, bytes) {
    if (!edit.doc) return setStatus('Choose a PDF first.', true);
    const w = Math.min(edit.pageW * edit.scale * (kind === 'signature' ? 0.3 : 0.22), naturalW);
    const h = w * naturalH / naturalW;
    const el = document.createElement('div');
    const img = document.createElement('img');
    img.src = src;
    img.draggable = false;
    el.appendChild(img);
    el.style.width = w + 'px';
    el.style.height = h + 'px';
    addItem({ kind, x: edit.pageW * edit.scale * 0.08, y: edit.pageH * edit.scale * (kind === 'signature' ? 0.76 : 0.5), w, h, bytes, el });
  }

  $('pdf-edit-file').addEventListener('change', async (e) => {
    const file = e.target.files[0];
    if (!file) return;
    try {
      edit.bytes = await readFile(file);
      edit.file = file;
      edit.doc = await pdfjsLib.getDocument({ data: edit.bytes.slice() }).promise;
      edit.pageCount = edit.doc.numPages;
      edit.page = 1;
      edit.items = [];
      edit.selected = null;
      await renderPage();
      setStatus('Loaded ' + file.name + ' (' + edit.pageCount + ' pages). Add items, drag them into place, then download.');
    } catch (err) {
      edit.doc = null;
      setStatus('Could not open this PDF: ' + err.message, true);
    }
  });

  $('pdf-edit-prev').addEventListener('click', () => { if (edit.doc && edit.page > 1) { edit.page--; renderPage(); } });
  $('pdf-edit-next').addEventListener('click', () => { if (edit.doc && edit.page < edit.pageCount) { edit.page++; renderPage(); } });
  $('pdf-add-text').addEventListener('click', () => addText($('pdf-text').value, 'text'));
  $('pdf-add-date').addEventListener('click', () => addText($('pdf-date').value, 'date'));
  $('pdf-delete-item').addEventListener('click', () => { if (edit.selected) removeItem(edit.selected); });

  ['text', 'date'].forEach((kind) => {
    sizeInput(kind).addEventListener('input', () => {
      const it = edit.selected;
      if (!it || it.kind !== kind) return;
      it.fontSize = Number(sizeInput(kind).value) || it.fontSize;
      styleTextItem(it);
    });
  });
  $('pdf-text-color').addEventListener('input', () => {
    const it = edit.selected;
    if (!it || it.text === undefined) return;
    it.color = $('pdf-text-color').value;
    styleTextItem(it);
  });
  $('pdf-font-family').addEventListener('change', () => {
    const it = edit.selected;
    if (!it || it.text === undefined) return;
    it.family = $('pdf-font-family').value;
    styleTextItem(it);
  });
  $('pdf-date').valueAsDate = new Date();
  $('pdf-text').addEventListener('input', () => {
    const it = edit.selected;
    if (it && it.kind === 'text') setItemText(it, $('pdf-text').value);
  });
  $('pdf-date').addEventListener('input', () => {
    const it = edit.selected;
    if (it && it.kind === 'date' && $('pdf-date').value) setItemText(it, $('pdf-date').value);
  });

  // Signature: draw or type
  const sigCanvas = $('pdf-sig-canvas');
  const sigCtx = sigCanvas.getContext('2d');
  let drawing = false;
  function sigPoint(e) {
    const r = sigCanvas.getBoundingClientRect();
    return [(e.clientX - r.left) * sigCanvas.width / r.width, (e.clientY - r.top) * sigCanvas.height / r.height];
  }
  function sigStyle() {
    sigCtx.lineWidth = 3;
    sigCtx.lineCap = 'round';
    sigCtx.lineJoin = 'round';
    sigCtx.strokeStyle = $('pdf-sig-color').value;
  }
  sigCanvas.addEventListener('pointerdown', (e) => {
    drawing = true;
    sigCanvas.setPointerCapture(e.pointerId);
    sigStyle();
    sigCtx.beginPath();
    const [x, y] = sigPoint(e);
    sigCtx.moveTo(x, y);
    sigCtx.lineTo(x + 0.01, y);
    sigCtx.stroke();
  });
  sigCanvas.addEventListener('pointermove', (e) => {
    if (!drawing) return;
    const [x, y] = sigPoint(e);
    sigCtx.lineTo(x, y);
    sigCtx.stroke();
  });
  ['pointerup', 'pointercancel'].forEach((n) => sigCanvas.addEventListener(n, () => { drawing = false; }));
  $('pdf-sig-clear').addEventListener('click', () => sigCtx.clearRect(0, 0, sigCanvas.width, sigCanvas.height));

  const sigMode = $('pdf-sig-mode');
  const syncSigMode = () => {
    const typed = sigMode.value === 'type';
    $('pdf-sig-draw-box').hidden = typed;
    $('pdf-sig-type-box').hidden = !typed;
  };
  sigMode.addEventListener('change', syncSigMode);
  syncSigMode();

  function cropToBlob(source) {
    const w = source.width, h = source.height;
    const data = source.getContext('2d').getImageData(0, 0, w, h).data;
    let minX = w, minY = h, maxX = -1, maxY = -1;
    for (let y = 0; y < h; y++) {
      for (let x = 0; x < w; x++) {
        if (data[(y * w + x) * 4 + 3] > 12) {
          if (x < minX) minX = x;
          if (x > maxX) maxX = x;
          if (y < minY) minY = y;
          if (y > maxY) maxY = y;
        }
      }
    }
    if (maxX < 0) return Promise.resolve(null);
    const pad = 4;
    minX = Math.max(0, minX - pad); minY = Math.max(0, minY - pad);
    maxX = Math.min(w - 1, maxX + pad); maxY = Math.min(h - 1, maxY + pad);
    const out = document.createElement('canvas');
    out.width = maxX - minX + 1;
    out.height = maxY - minY + 1;
    out.getContext('2d').drawImage(source, minX, minY, out.width, out.height, 0, 0, out.width, out.height);
    return new Promise((resolve) => out.toBlob((blob) => resolve({ blob, w: out.width, h: out.height }), 'image/png'));
  }

  async function addImageFromBlob(kind, result) {
    if (!result) return setStatus('Nothing to add yet.', true);
    const bytes = new Uint8Array(await result.blob.arrayBuffer());
    addImage(kind, URL.createObjectURL(result.blob), result.w, result.h, bytes);
  }

  $('pdf-sig-add').addEventListener('click', async () => {
    let source = sigCanvas;
    if (sigMode.value === 'type') {
      const name = $('pdf-sig-name').value.trim();
      if (!name) return setStatus('Type a name for the signature.', true);
      source = document.createElement('canvas');
      source.width = 900;
      source.height = 220;
      const c = source.getContext('2d');
      c.fillStyle = $('pdf-sig-color').value;
      c.textBaseline = 'middle';
      c.font = SIGNATURE_FONTS[$('pdf-sig-style').value].replace(/^((?:italic |bold )*)(.*)$/, '$1 88px $2');
      c.fillText(name, 20, 110);
    }
    await addImageFromBlob('signature', await cropToBlob(source));
  });

  $('pdf-stamp-file').addEventListener('change', async (e) => {
    const file = e.target.files[0];
    e.target.value = '';
    if (!file) return;
    try {
      let blob = file;
      if (file.type === 'application/pdf') {
        const doc = await pdfjsLib.getDocument({ data: await readFile(file) }).promise;
        const page = await doc.getPage(1);
        const vp = page.getViewport({ scale: 2 });
        const c = document.createElement('canvas');
        c.width = vp.width; c.height = vp.height;
        await page.render({ canvasContext: c.getContext('2d'), viewport: vp }).promise;
        const result = await cropToBlob(c);
        return addImageFromBlob('stamp', result);
      }
      const bitmap = await createImageBitmap(blob);
      const c = document.createElement('canvas');
      c.width = bitmap.width; c.height = bitmap.height;
      c.getContext('2d').drawImage(bitmap, 0, 0);
      await addImageFromBlob('stamp', await cropToBlob(c));
    } catch (err) {
      setStatus('Could not read stamp: ' + err.message, true);
    }
  });

  $('pdf-edit-download').addEventListener('click', async () => {
    if (!edit.doc) return setStatus('Choose a PDF first.', true);
    if (!edit.items.length) return setStatus('Add text, a date, a signature, or a stamp first.', true);
    try {
      const pdf = await PDFDocument.load(edit.bytes);
      const fonts = {};
      for (const key of Object.keys(FONT_MAP)) fonts[key] = await pdf.embedFont(FONT_MAP[key].pdf);
      const pages = pdf.getPages();
      for (const it of edit.items) {
        const page = pages[it.page - 1];
        const pageH = page.getHeight();
        const k = 1 / edit.scale;
        if (it.text !== undefined) {
          const font = fonts[it.family];
          const size = it.fontSize;
          const ascent = font.heightAtSize(size, { descender: false });
          it.text.split(/\r?\n/).forEach((line, i) => {
            page.drawText(line, {
              x: it.x * k,
              y: pageH - it.y * k - size * 0.1 - ascent - i * size * 1.2,
              size, font, color: hexToRgb(it.color),
            });
          });
        } else {
          const image = await pdf.embedPng(it.bytes);
          page.drawImage(image, { x: it.x * k, y: pageH - (it.y + it.h) * k, width: it.w * k, height: it.h * k });
        }
      }
      download(await pdf.save(), baseName(edit.file) + '_edited.pdf', 'application/pdf');
      setStatus('Edited PDF downloaded.');
    } catch (err) {
      setStatus('Could not apply changes: ' + err.message, true);
    }
  });

  window.addEventListener('resize', () => {
    if (!edit.doc) return;
    // Scale items to the new page width so placement stays aligned.
    const old = edit.scale;
    clearTimeout(window.__pdfResizeTimer);
    window.__pdfResizeTimer = setTimeout(async () => {
      await renderPage();
      const ratio = edit.scale / old;
      if (ratio === 1) return;
      edit.items.forEach((it) => {
        it.x *= ratio; it.y *= ratio;
        it.el.style.left = it.x + 'px'; it.el.style.top = it.y + 'px';
        if (it.text !== undefined) styleTextItem(it);
        else {
          it.w *= ratio; it.h *= ratio;
          it.el.style.width = it.w + 'px'; it.el.style.height = it.h + 'px';
        }
      });
    }, 250);
  });

  // ---------- Merge ----------
  let mergeFiles = [];
  function renderMergeList() {
    const list = $('pdf-merge-list');
    list.innerHTML = '';
    mergeFiles.forEach((file, i) => {
      const li = document.createElement('li');
      const name = document.createElement('span');
      name.textContent = file.name;
      li.appendChild(name);
      [['Up', -1], ['Down', 1], ['Remove', 0]].forEach(([label, dir]) => {
        const b = document.createElement('button');
        b.type = 'button';
        b.className = 'digitizer-button is-secondary';
        b.textContent = label;
        b.addEventListener('click', () => {
          if (dir === 0) mergeFiles.splice(i, 1);
          else if (i + dir >= 0 && i + dir < mergeFiles.length) [mergeFiles[i], mergeFiles[i + dir]] = [mergeFiles[i + dir], mergeFiles[i]];
          renderMergeList();
        });
        li.appendChild(b);
      });
      list.appendChild(li);
    });
  }
  $('pdf-merge-file').addEventListener('change', (e) => {
    mergeFiles = mergeFiles.concat(Array.from(e.target.files));
    e.target.value = '';
    renderMergeList();
  });
  $('pdf-merge-btn').addEventListener('click', async () => {
    if (!mergeFiles.length) return setStatus('Choose at least one PDF.', true);
    try {
      const out = await PDFDocument.create();
      for (const file of mergeFiles) {
        const src = await PDFDocument.load(await readFile(file));
        const copied = await out.copyPages(src, src.getPageIndices());
        copied.forEach((p) => out.addPage(p));
      }
      download(await out.save(), 'merged.pdf', 'application/pdf');
      setStatus('Merged ' + mergeFiles.length + ' PDFs.');
    } catch (err) {
      setStatus('Could not merge PDFs: ' + err.message, true);
    }
  });

  // ---------- Compress ----------
  $('pdf-compress-btn').addEventListener('click', async () => {
    const file = $('pdf-compress-file').files[0];
    if (!file) return setStatus('Choose a PDF first.', true);
    try {
      const [dpi, quality] = $('pdf-compress-preset').value.split(',').map(Number);
      const src = await pdfjsLib.getDocument({ data: await readFile(file) }).promise;
      const out = await PDFDocument.create();
      for (let n = 1; n <= src.numPages; n++) {
        setStatus('Compressing page ' + n + ' of ' + src.numPages + '...');
        const page = await src.getPage(n);
        const base = page.getViewport({ scale: 1 });
        const vp = page.getViewport({ scale: dpi / 72 });
        const c = document.createElement('canvas');
        c.width = Math.round(vp.width);
        c.height = Math.round(vp.height);
        const ctx = c.getContext('2d');
        ctx.fillStyle = '#fff';
        ctx.fillRect(0, 0, c.width, c.height);
        await page.render({ canvasContext: ctx, viewport: vp }).promise;
        const blob = await new Promise((r) => c.toBlob(r, 'image/jpeg', quality));
        const img = await out.embedJpg(new Uint8Array(await blob.arrayBuffer()));
        out.addPage([base.width, base.height]).drawImage(img, { x: 0, y: 0, width: base.width, height: base.height });
      }
      const bytes = await out.save();
      download(bytes, baseName(file) + '_compressed.pdf', 'application/pdf');
      setStatus('Compressed ' + (file.size / 1024).toFixed(0) + ' KB to ' + (bytes.length / 1024).toFixed(0) + ' KB.');
    } catch (err) {
      setStatus('Could not compress this PDF: ' + err.message, true);
    }
  });

  // ---------- Split ----------
  let splitFile = null;
  let splitCount = 0;
  $('pdf-split-file').addEventListener('change', async (e) => {
    splitFile = e.target.files[0] || null;
    if (!splitFile) return;
    try {
      const doc = await PDFDocument.load(await readFile(splitFile));
      splitCount = doc.getPageCount();
      $('pdf-split-pages').value = '1-' + splitCount;
      setStatus(splitFile.name + ' has ' + splitCount + ' pages.');
    } catch (err) {
      splitFile = null;
      setStatus('Could not open this PDF: ' + err.message, true);
    }
  });
  $('pdf-split-btn').addEventListener('click', async () => {
    if (!splitFile) return setStatus('Choose a PDF first.', true);
    try {
      const pages = parsePages($('pdf-split-pages').value, splitCount);
      const src = await PDFDocument.load(await readFile(splitFile));
      const extract = async (list) => {
        const out = await PDFDocument.create();
        (await out.copyPages(src, list.map((p) => p - 1))).forEach((p) => out.addPage(p));
        return out.save();
      };
      if ($('pdf-split-output').value === 'single') {
        download(await extract(pages), baseName(splitFile) + '_pages.pdf', 'application/pdf');
      } else {
        if (!window.JSZip) throw new Error('ZIP library failed to load.');
        const zip = new window.JSZip();
        for (const g of parsePageGroups($('pdf-split-pages').value, splitCount)) {
          const name = g.length === 1 ? 'page_' + g[0] : 'pages_' + g[0] + '-' + g[g.length - 1];
          zip.file(name + '.pdf', await extract(g));
        }
        download(await zip.generateAsync({ type: 'uint8array' }), baseName(splitFile) + '_split.zip', 'application/zip');
      }
      setStatus('Split complete.');
    } catch (err) {
      setStatus(err.message, true);
    }
  });
})();
