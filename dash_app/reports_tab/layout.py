"""Build the document list and PDF/Markdown/plain-text reading pane for doc/ files."""

from pathlib import Path
import posixpath
import re
from urllib.parse import quote, unquote, urlsplit

from dash import dcc, html

from .catalog import DOCUMENT_TYPES, catalog_token, discover_documents, resolve_document, resolve_image


def document_tabs(documents):
    """Group document paths under non-selectable PDF, Markdown and Text headings."""
    tabs = []
    for suffix, title in DOCUMENT_TYPES.items():
        group = [document for document in documents if Path(document.path).suffix.lower() == suffix]
        if not group:
            continue
        tabs.append(dcc.Tab(
            label=title, value=f'document-category-{suffix[1:]}', disabled=True,
            className='reports-category-heading', disabled_className='reports-category-heading',
        ))
        tabs.extend(dcc.Tab(label=document.label, value=document.path,
                            className='reports-directory-card',
                            selected_className='reports-directory-card-selected')
                    for document in group)
    return tabs


def selected_value(documents, current=None):
    """Retain a visible selection, otherwise choose the first match or no document."""
    if any(document.path == current for document in documents):
        return current
    return documents[0].path if documents else None


def _local_images(text, document_path):
    """Point ordinary inline Markdown images at their doc-relative static files."""
    def replace(match):
        """Leave external URLs and unresolved images as written in the document."""
        target = match.group(2)
        if urlsplit(target).scheme or target.startswith(('/', '#')):
            return match.group(0)
        relative = posixpath.normpath(posixpath.join(posixpath.dirname(document_path), unquote(target)))
        if resolve_image(relative) is None:
            return match.group(0)
        return match.group(1) + '/document-images/' + quote(relative, safe='/')
    return re.sub(r'(!\[[^\]]*\]\()([^\s)]+)', replace, text)


def document_view(document):
    """Display PDFs, render Markdown, and preserve plain-text spacing literally."""
    if document is None:
        return html.Div('No matching PDF, Markdown or text documents.', className='reports-empty-state')
    path = resolve_document(document.path)
    if path is None:
        return html.Div('This document is no longer available.', className='reports-empty-state')
    if path.suffix.lower() == '.pdf':
        return html.Iframe(src=document.url, title=document.path, className='reports-frame')
    try:
        text = path.read_text(encoding='utf-8')
    except (OSError, UnicodeError):
        return html.Div('This text document could not be read as UTF-8.', className='reports-empty-state')
    if path.suffix.lower() == '.txt':
        return html.Pre(text, className='reports-text')
    return html.Div(
        dcc.Markdown(_local_images(text, document.path), mathjax=True,
                     link_target='_blank', className='reports-markdown'),
        **{'data-repository-path': 'doc/' + document.path},
    )


def build_layout():
    """Build a searchable path list and reader without publication metadata."""
    documents = discover_documents()
    return html.Div([
        dcc.Interval(id='reports-catalog-poll', interval=3000),
        dcc.Store(id='reports-catalog-token', data=catalog_token(documents)),
        html.Aside([
            html.Div([
                html.H2('Documents'),
                html.P('PDF, Markdown and plain-text files in doc/.',
                       className='reports-directory-intro'),
                dcc.Input(id='reports-search', type='search', placeholder='Search filenames or paths',
                          value='', className='reports-search'),
            ], className='reports-directory-header'),
            dcc.Tabs(id='reports-pages', value=selected_value(documents), vertical=True,
                     persistence=True, persistence_type='local', className='reports-directory',
                     children=document_tabs(documents)),
        ], className='reports-directory-shell'),
        html.Section(id='reports-view', children=document_view(documents[0] if documents else None),
                     className='reports-view'),
    ], className='reports-tab-root')
