"""Check document discovery, safe serving and the minimal PDF/Markdown/plain-text reader."""

from pathlib import Path
import subprocess

from dash import Dash, dcc, html, no_update

from dash_app.reports_tab import catalog, layout, tab
from dash_app.reports_tab.static import register_document_routes


def _document_root(tmp_path, monkeypatch):
    """Create ordinary documents alongside files the reader must ignore."""
    root = tmp_path / 'doc'
    root.mkdir()
    for relative, text in {
        'note.md': '# A note\n\n![Figure](figures/test.svg)\n\n<script>alert(1)</script>',
        'nested/note.md': '# Another note', 'equations.pdf': '%PDF-1.4\n',
        'z-last.pdf': '%PDF-1.4\n',
        'plain.TXT': '# Literal heading\n  <script>alert(1)</script>\n\n    indented text\n',
        'README.md': '# Root README', 'nested/README.md': '# Nested README',
        '.hidden/note.md': '# Hidden', 'analysis.ipynb': '{}',
        'old/report.html': '<h1>Old report</h1>',
        'figures/test.svg': '<svg xmlns="http://www.w3.org/2000/svg"/>',
    }.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    monkeypatch.setattr(catalog, 'DOC_ROOT', root)
    return root


def _callback(app, name):
    """Find the registered native callback so tests exercise its actual wiring."""
    return next(item['callback'].__wrapped__ for item in app.callback_map.values()
                if item.get('callback') and item['callback'].__name__ == name)


def test_discovery_groups_formats_and_includes_readmes(tmp_path, monkeypatch):
    """Both nested paths and simple filename searches use one document catalog."""
    root = _document_root(tmp_path, monkeypatch)
    (tmp_path / 'outside.md').write_text('# Outside')
    (root / 'escape.md').symlink_to(tmp_path / 'outside.md')
    documents = catalog.discover_documents()
    assert [item.path for item in documents] == [
        'equations.pdf', 'z-last.pdf', 'nested/note.md', 'nested/README.md',
        'note.md', 'README.md', 'plain.TXT',
    ]
    entries = layout.document_tabs(documents)
    assert [item.label for item in entries if getattr(item, 'disabled', False)] == ['PDFs', 'Markdown', 'Text']
    assert [item.label for item in entries if not getattr(item, 'disabled', False)] == [item.path for item in documents]
    assert layout.selected_value(documents) == 'equations.pdf'
    assert catalog.resolve_document('README.md') is not None
    assert catalog.resolve_document('nested/README.md') is not None
    assert [item.path for item in catalog.discover_documents(query='NESTED note')] == ['nested/note.md']
    assert [item.path for item in catalog.discover_documents(query='PLAIN txt')] == ['plain.TXT']
    assert catalog.discover_documents(query='not present') == ()
    for path in ['../outside.md', 'escape.md', '.hidden/note.md', 'old/report.html']:
        assert catalog.resolve_document(path) is None


def test_file_versions_detect_edits_and_catalog_additions(tmp_path, monkeypatch):
    """No manifest is needed to discover a new file or reload an edited note."""
    root = _document_root(tmp_path, monkeypatch)
    original = catalog.catalog_token(catalog.discover_documents())
    (root / 'note.md').write_text('# Updated note')
    assert catalog.catalog_token(catalog.discover_documents()) != original
    (root / 'fresh.pdf').write_bytes(b'%PDF-1.4\n')
    assert 'fresh.pdf' in {item.path for item in catalog.discover_documents()}


def test_pdf_and_markdown_use_native_renderers(tmp_path, monkeypatch):
    """PDFs stay in the browser; Markdown disables HTML and resolves local images."""
    _document_root(tmp_path, monkeypatch)
    docs = {item.path: item for item in catalog.discover_documents()}
    pdf = layout.document_view(docs['equations.pdf'])
    md = layout.document_view(docs['note.md'])
    assert isinstance(pdf, html.Iframe) and pdf.src.startswith('/documents/equations.pdf?v=')
    assert md.to_plotly_json()['props']['data-repository-path'] == 'doc/note.md'
    markdown = md.children
    assert isinstance(markdown, dcc.Markdown)
    assert not getattr(markdown, 'dangerously_allow_html', False)
    assert '/document-images/figures/test.svg' in markdown.children
    plain = layout.document_view(docs['plain.TXT'])
    assert isinstance(plain, html.Pre)
    assert plain.children == '# Literal heading\n  <script>alert(1)</script>\n\n    indented text\n'
    assert isinstance(layout.document_view(None), html.Div)


def test_static_routes_serve_documents_and_images_only(tmp_path, monkeypatch):
    """Traversal, missing files and executable formats return 404; READMEs are served."""
    _document_root(tmp_path, monkeypatch)
    app = Dash(__name__)
    app.layout = html.Div()
    register_document_routes(app)
    register_document_routes(app)
    client = app.server.test_client()
    pdf = client.get('/documents/equations.pdf')
    assert pdf.status_code == 200 and pdf.mimetype == 'application/pdf'
    assert client.get('/documents/note.md').status_code == 200
    plain = client.get('/documents/plain.TXT')
    assert plain.status_code == 200 and plain.mimetype == 'text/plain'
    assert plain.data.startswith(b'# Literal heading\n  <script>')
    assert client.get('/document-images/figures/test.svg').status_code == 200
    assert client.get('/documents/README.md').status_code == 200
    assert client.get('/documents/nested/README.md').status_code == 200
    for url in ['/documents/../outside.md',
                '/documents/old/report.html', '/documents/missing.pdf', '/document-images/note.md']:
        assert client.get(url).status_code == 404


def test_callbacks_search_refresh_retain_selection_and_handle_removal(tmp_path, monkeypatch):
    """Search changes the list; stable polls do not rebuild the current reader."""
    root = _document_root(tmp_path, monkeypatch)
    app = Dash(__name__, suppress_callback_exceptions=True)
    app.layout = tab.build_tab(app)
    app._setup_server()
    refresh = _callback(app, 'refresh_documents')
    read = _callback(app, 'display_document')
    handoff = _callback(app, 'select_agent_document')
    assert refresh(1, 'plots', '', 'old', None) == (no_update,) * 3
    entries, selection, token = refresh(1, 'reports', 'nested', 'old', 'equations.pdf')
    assert [entry.value for entry in entries if not getattr(entry, 'disabled', False)] == ['nested/note.md', 'nested/README.md']
    assert selection == 'nested/note.md'
    assert refresh(2, 'reports', 'nested', token, selection) == (no_update,) * 3
    assert isinstance(read(selection, token).children, dcc.Markdown)
    assert handoff({'tab': 'reports', 'operation': 'open_document', 'path': 'note.md'}) == ('note.md', '')
    assert handoff({'tab': 'reports', 'operation': 'open_document', 'path': 'plain.TXT'}) == ('plain.TXT', '')
    assert handoff({'tab': 'reports', 'operation': 'open_document', 'path': 'README.md'}) == ('README.md', '')
    (root / 'nested/note.md').unlink()
    (root / 'nested/README.md').unlink()
    entries, selection, _token = refresh(3, 'reports', 'nested', token, selection)
    assert entries == [] and selection is None
    assert isinstance(read('nested/note.md', _token), html.Div)


def test_existing_report_pdf_is_discoverable():
    """The migrated report is a PDF document without an index or HTML adapter."""
    path = 'sam-heatmap-color-signals.pdf'
    assert path in {item.path for item in catalog.discover_documents()}
    assert catalog.resolve_document(path).read_bytes().startswith(b'%PDF-')


def test_linked_sources_open_as_text_without_exposing_local_or_outside_files(tmp_path, monkeypatch):
    """Tracked Fortran/Python files are literal text; hidden, ignored and binary targets fail."""
    root = _document_root(tmp_path, monkeypatch)
    repo = root.parent
    source = repo / 'src'
    source.mkdir()
    for name, content in {
        'model.F90': b'program example\n  ! <script>literal source</script>\nend program\n',
        'helper.py': b'raise RuntimeError("must never execute")\n',
        'payload.html': b'<script>must stay text</script>',
        'binary.dat': b'abc\x00def',
    }.items():
        (source / name).write_bytes(content)
    subprocess.run(['git', 'init', '-q', str(repo)], check=True)
    subprocess.run(['git', 'add', 'src'], cwd=repo, check=True)
    (source / 'local.F90').write_text('untracked private file')
    (repo / '.secret').write_text('hidden file')
    outside = tmp_path.parent / (tmp_path.name + '-outside.F90')
    outside.write_text('outside repository')
    (source / 'escape.F90').symlink_to(outside)
    try:
        app = Dash(__name__)
        app.layout = html.Div()
        register_document_routes(app)
        client = app.server.test_client()
        for name in ['model.F90', 'helper.py', 'payload.html']:
            response = client.get('/repository-files/src/' + name)
            assert response.status_code == 200
            assert response.mimetype == 'text/plain'
            assert response.data == (source / name).read_bytes()
            assert response.headers['X-Content-Type-Options'] == 'nosniff'
        assert client.get('/repository-files/doc/note.md').status_code == 200
        assert client.get('/repository-files/doc/equations.pdf').mimetype == 'application/pdf'
        for path in ['src/local.F90', 'src/escape.F90', '.secret', '../outside.F90', 'src/missing.F90']:
            assert client.get('/repository-files/' + path).status_code == 404
        assert client.get('/repository-files/src/binary.dat').status_code == 415
    finally:
        outside.unlink()


def test_rendered_markdown_links_resolve_from_the_document_folder():
    """Native click, modifier/middle click and context-menu targets retain file-relative URLs."""
    asset = Path(__file__).parents[3] / 'assets' / '27_document_links.js'
    script = r"""
const fs = require('fs'), assert = require('assert');
const handlers = {};
global.window = {location: {origin: 'http://localhost:23404'}};
global.document = {addEventListener: (name, fn) => { handlers[name] = fn; }};
eval(fs.readFileSync(process.argv[1], 'utf8'));
function follow(href, path, event = 'click', inViewer = true) {
  const attrs = {href};
  const anchor = {
    getAttribute: name => attrs[name],
    setAttribute: (name, value) => { attrs[name] = value; },
    closest: () => inViewer ? {dataset: {repositoryPath: path}} : null
  };
  handlers[event]({target: {closest: () => anchor}});
  return attrs;
}
for (const event of ['click', 'pointerdown', 'focusin', 'auxclick', 'contextmenu']) {
  const link = follow('../src/CLUBB_core/mixing_length.F90#L10', 'doc/coherence_mixing_length.md', event);
  assert.equal(link.href, 'http://localhost:23404/repository-files/src/CLUBB_core/mixing_length.F90#L10');
  assert.equal(link.rel, 'noopener noreferrer');
}
assert.equal(follow('checked_functions.py', 'doc/triple_gaussian/README.md').href,
             'http://localhost:23404/repository-files/doc/triple_gaussian/checked_functions.py');
assert.equal(follow('../helper%20model.F90?raw=1', 'doc/nested/note.md').href,
             'http://localhost:23404/repository-files/doc/helper%20model.F90?raw=1');
assert.equal(follow('/src/clubb_driver.F90', 'doc/note.md').href,
             'http://localhost:23404/repository-files/src/clubb_driver.F90');
for (const href of ['https://example.org/paper', '//example.org/paper', 'mailto:author@example.org',
                    '/repository-files/src/model.F90', '/documents/note.md']) {
  assert.equal(follow(href, 'doc/note.md').href, href);
}
assert.equal(follow('helper.py', 'doc/note.md', 'click', false).href, 'helper.py');
"""
    subprocess.run(['node', '-e', script, str(asset)], check=True, capture_output=True, text=True)
