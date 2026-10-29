"""Contracts for immutable static-report publication and Dash discovery."""

from __future__ import annotations

from dash import Dash, dcc
from PIL import Image

from dash_app.reports_tab.catalog import catalog_token, discover_reports, report_by_id, resolve_report_asset
from dash_app.reports_tab.layout import build_layout
from dash_app.reports_tab.publisher import ReportBuilder
from dash_app.reports_tab.static import register_static_report_routes
from dash_app.reports_tab.tab import _opened_report_url, build_tab


def _walk(component):
    yield component
    children = getattr(component, "children", None)
    if children is None:
        return
    if not isinstance(children, (list, tuple)):
        children = [children]
    for child in children:
        if hasattr(child, "to_plotly_json"):
            yield from _walk(child)


def test_checked_in_sam_report_is_discoverable_and_static():
    reports = discover_reports()
    report = report_by_id("sam-heatmap-coloring-examples-20h30")

    assert report is not None
    assert report.title == "SAM heatmap colors: ARM 20.5 h / 1900 m"
    assert report in reports
    assert report.url == "/static-reports/sam-heatmap-coloring-examples-20h30/report.html"
    assert resolve_report_asset("sam-heatmap-coloring-examples-20h30", "report.html") is not None
    assert resolve_report_asset(
        "sam-heatmap-coloring-examples-20h30", "figures/raw-sam-wrt-color-signals.png"
    ) is not None
    assert resolve_report_asset("sam-heatmap-coloring-examples-20h30", "../report.html") is None
    assert len(catalog_token(reports)) == 16


def test_report_builder_publishes_a_complete_bundle_atomically(tmp_path):
    builder = ReportBuilder(
        "test-report",
        "Test report",
        summary="A self-contained report used only for this test.",
        tags=("test",),
        root=tmp_path,
    )
    report = builder.heading("Result").paragraph("The report is ready.").code("print('ready')", language="python").publish()

    assert report.report_id == "test-report"
    assert (tmp_path / "test-report" / "report.html").is_file()
    assert (tmp_path / "test-report" / "manifest.json").is_file()
    assert not any((tmp_path / ".staging").iterdir())
    assert [entry.report_id for entry in discover_reports(tmp_path)] == ["test-report"]


def test_report_builder_replaces_a_published_bundle_in_place(tmp_path):
    ReportBuilder("test-report", "Original", summary="first", root=tmp_path).publish()
    replacement = ReportBuilder(
        "test-report", "Replacement", summary="second", root=tmp_path, replace=True
    ).publish()

    assert replacement.title == "Replacement"
    assert [entry.report_id for entry in discover_reports(tmp_path)] == ["test-report"]
    assert report_by_id("test-report", tmp_path).title == "Replacement"
    assert "Replacement" in (tmp_path / "test-report" / "report.html").read_text(encoding="utf-8")


def test_report_builder_caps_raster_figure_dimensions(tmp_path):
    source = tmp_path / "oversized.png"
    Image.new("RGBA", (2400, 600), "#38bdf8").save(source)

    ReportBuilder("test-report", "Test report", summary="raster cap", root=tmp_path).figure(
        source, caption="A deliberately oversized test image."
    ).publish()

    with Image.open(tmp_path / "test-report" / "figures" / "oversized.png") as image:
        assert image.size == (1200, 300)


def test_reports_layout_keeps_navigation_separate_from_the_iframe_viewer():
    layout = build_layout()
    components = list(_walk(layout))
    directory = next(component for component in components if getattr(component, "id", None) == "reports-pages")
    frame = next(component for component in components if getattr(component, "id", None) == "reports-frame")
    reports = discover_reports()

    assert isinstance(directory, dcc.Tabs)
    assert directory.vertical is True
    assert all(isinstance(tab.label, str) for tab in directory.children)
    assert frame.src == reports[0].url
    assert "sam-heatmap-coloring-examples-20h30" in {report.report_id for report in reports}
    assert "sam-heatmap-coloring-examples" not in {report.report_id for report in reports}


def test_agent_handoff_gets_a_unique_iframe_url_even_for_the_current_report():
    report = report_by_id("sam-heatmap-coloring-examples-20h30")

    assert report is not None
    assert _opened_report_url(report, 67) == "/static-reports/sam-heatmap-coloring-examples-20h30/report.html?open=67"


def test_reports_tab_and_static_route_build_with_dash():
    app = Dash(__name__, suppress_callback_exceptions=True, eager_loading=True)
    register_static_report_routes(app)
    tab = build_tab(app)
    app.layout = tab
    app._setup_server()

    response = app.server.test_client().get(
        "/static-reports/sam-heatmap-coloring-examples-20h30/report.html"
    )
    missing = app.server.test_client().get(
        "/static-reports/sam-heatmap-coloring-examples-20h30/../manifest.json"
    )

    assert tab.value == "reports"
    assert response.status_code == 200
    assert b"SAM heatmap color signals" in response.data
    assert missing.status_code == 404
