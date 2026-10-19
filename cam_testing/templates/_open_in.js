// The version page's "Open in" links (included inline by version.html; tests/test_open_in.py runs
// it under node).
//   a[data-open-file][data-open-url]: on click, fetches the file (relative to the page), appends
//     the standard base64 of its bytes to the app URL (which ends in "#") and opens that in a new
//     tab. If the file can't be read (a page opened from disk, in a browser that blocks reading
//     local files), the link's .open-msg says so and offers the file to download instead.
//   a[data-cuflynx-omex][data-cuflynx-open]: its href becomes the app URL followed by the
//     percent-encoded absolute URL of the archive, resolved against the page's own URL.
(function (root) {
  'use strict';
  var CHUNK = 0x8000;   // bytes per String.fromCharCode call: well inside any engine's argument limit
  // URLs much longer than this fail in some browsers (Chrome's limit is 2 MB)
  var MAX_URL = 2 * 1024 * 1024;

  function bytesToBase64(bytes) {
    var parts = [];
    for (var i = 0; i < bytes.length; i += CHUNK) {
      parts.push(String.fromCharCode.apply(null, bytes.subarray(i, i + CHUNK)));
    }
    return root.btoa(parts.join(''));
  }

  function cuflynxUrl(open, rel, base) {
    return open + encodeURIComponent(new URL(rel, base).href);
  }

  root.camOpenIn = { bytesToBase64: bytesToBase64, cuflynxUrl: cuflynxUrl };
  if (typeof document === 'undefined') { return; }

  function say(link, text, file) {
    var msg = link.parentNode.querySelector('.open-msg');
    if (!msg) { return; }
    msg.textContent = text;
    if (file) {
      msg.appendChild(document.createTextNode(' '));
      var a = document.createElement('a');
      a.href = file;
      a.setAttribute('download', '');
      a.textContent = 'Download ' + file.split('/').pop();
      msg.appendChild(a);
      msg.appendChild(document.createTextNode(' and open it in ' + link.dataset.app + ' instead.'));
    }
  }

  function open(event) {
    var link = event.currentTarget;
    var file = link.dataset.openFile;
    event.preventDefault();
    say(link, 'Reading ' + file.split('/').pop() + '…');
    // the tab is opened now, within the click, so a popup blocker allows it
    var tab = root.open('', '_blank');
    fetch(file).then(function (r) {
      if (!r.ok) { throw new Error('HTTP ' + r.status); }
      return r.arrayBuffer();
    }).then(function (buf) {
      var url = link.dataset.openUrl + bytesToBase64(new Uint8Array(buf));
      if (url.length > MAX_URL) {
        if (tab) { tab.close(); }
        say(link, 'The file is too large to open through a link.', file);
        return;
      }
      if (tab) {
        tab.opener = null;
        tab.location.href = url;
      } else {
        root.open(url, '_blank', 'noopener');
      }
      say(link, '');
    }).catch(function (e) {
      if (tab) { tab.close(); }
      say(link, 'Could not read the file from this page (' + e.message + '; a page opened from disk may not read local files).', file);
    });
  }

  function init() {
    document.querySelectorAll('a[data-open-file][data-open-url]').forEach(function (a) {
      a.addEventListener('click', open);
    });
    document.querySelectorAll('a[data-cuflynx-omex][data-cuflynx-open]').forEach(function (a) {
      a.href = cuflynxUrl(a.dataset.cuflynxOpen, a.dataset.cuflynxOmex, root.location.href);
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})(typeof window !== 'undefined' ? window : globalThis);
