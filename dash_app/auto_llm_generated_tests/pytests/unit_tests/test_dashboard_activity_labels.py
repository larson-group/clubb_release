from dash_app.agent_integration.handoff import dashboard_tab_labels


def test_tab_labels_describe_work_in_each_broker_launching_tab():
    labels = dashboard_tab_labels(
        {
            "compile": {"state": "running", "job": {"kind": "rebuild"}},
            "run_summary": {"running": 1, "queued": 1, "stopping": 0},
            "profile": {"state": "running"},
            "pyplotgen": {"state": "running"},
            "tune": {"state": "running"},
            "loss_runs": {"loss-1": {"state": "running"}},
        }
    )

    assert labels == (
        "Compile · rebuilding",
        "Run · 1 running · 1 queued",
        "Profile · timing",
        "Tune · tuning · 1 result run",
        "Plots · generating",
    )
