"""Serve documents and linked repository files without executing source code."""

from flask import Response, abort, send_file

from .catalog import resolve_document, resolve_image, resolve_linked_file


def register_document_routes(app):
    """Register bounded document, image and linked-file routes once."""
    if 'serve_document' in app.server.view_functions:
        return

    @app.server.get('/documents/<path:document_path>', endpoint='serve_document')
    def serve_document(document_path):
        """Serve a validated document with browser cache revalidation."""
        path = resolve_document(document_path)
        if path is None:
            abort(404)
        return send_file(path, max_age=0, conditional=True)

    @app.server.get('/document-images/<path:image_path>', endpoint='serve_document_image')
    def serve_document_image(image_path):
        """Serve only image formats needed by Markdown, confined to doc/."""
        path = resolve_image(image_path)
        if path is None:
            abort(404)
        response = send_file(path, max_age=0, conditional=True)
        response.headers['Content-Security-Policy'] = "sandbox"
        return response

    @app.server.get('/repository-files/<path:file_path>', endpoint='serve_linked_repository_file')
    def serve_linked_repository_file(file_path):
        """Open linked sources as literal text; PDFs/images use their native viewers."""
        path = resolve_linked_file(file_path)
        if path is None:
            abort(404)
        if path.suffix.lower() in {'.pdf', '.png', '.jpg', '.jpeg', '.gif', '.webp', '.svg'}:
            response = send_file(path, max_age=0, conditional=True)
            if path.suffix.lower() == '.svg':
                response.headers['Content-Security-Policy'] = 'sandbox'
        else:
            try:
                text = path.read_text(encoding='utf-8')
            except UnicodeError:
                abort(415, description='This file is not UTF-8 text.')
            except OSError:
                abort(404)
            if '\x00' in text:
                abort(415, description='This file is not a text document.')
            response = Response(text, mimetype='text/plain')
            response.headers['Cache-Control'] = 'no-cache'
        response.headers['X-Content-Type-Options'] = 'nosniff'
        return response
