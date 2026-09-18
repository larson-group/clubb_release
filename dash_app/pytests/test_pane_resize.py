"""Resize listeners must not retain old page elements after tab navigation."""

from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.mark.parametrize("asset,prefix", [
    ("30_run_ui_resize.js", "run"),
    ("31_plots_ui_resize.js", "plots"),
    ("36_tune_ui_resize.js", "tune"),
])
def test_resize_gesture_cleanup_and_saved_width(asset, prefix):
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is needed for the browser-event regression")
    script = r"""
    const assert = require('node:assert/strict');
    const fs = require('node:fs');
    const vm = require('node:vm');
    const prefix = process.argv[2];
    const listeners = new Map(), storage = new Map(), classes = new Set();
    let tick, layout, divider;
    const window = {
        innerWidth: 1200,
        localStorage: {getItem: key => storage.get(key), setItem: (key, val) => storage.set(key, val)},
        addEventListener(type, fn) {
            if (!listeners.has(type)) listeners.set(type, new Set());
            listeners.get(type).add(fn);
        },
        removeEventListener(type, fn) { listeners.get(type)?.delete(fn); },
        dispatchEvent(event) { for (const fn of [...(listeners.get(event.type) || [])]) fn(event); },
        requestAnimationFrame(fn) { fn(); },
        setTimeout(fn) { fn(); },
        setInterval(fn) { tick = fn; }
    };
    const document = {
        readyState: 'complete',
        getElementById: id => id === prefix + '-tab-layout' ? layout : divider,
        addEventListener() {},
        body: {classList: {add: name => classes.add(name), remove: name => classes.delete(name)}}
    };
    const idle = () => {
        for (const type of ['pointermove', 'pointerup', 'pointercancel', 'blur']) {
            assert.equal(listeners.get(type)?.size || 0, 0, type + ' leaked');
        }
        assert.equal(classes.size, 0);
    };
    const mount = () => {
        if (layout) layout.isConnected = false;
        layout = {dataset: {}, style: {}, isConnected: true,
                  getBoundingClientRect: () => ({width: 1200, right: 1200})};
        divider = {dataset: {}, addEventListener(type, fn) {this[type] = fn;}};
    };
    mount();
    vm.runInNewContext(fs.readFileSync(process.argv[1], 'utf8'), {window, document, Event});
    for (let i = 0; i < 100; i++) { mount(); tick(); idle(); }
    const start = () => divider.pointerdown({button: 0, preventDefault() {}});
    start();
    start(); // Repeated pointerdown must not attach multiple gesture listeners.
    assert.equal(listeners.get('pointermove').size, 1);
    window.dispatchEvent({type: 'pointermove', clientX: 700});
    assert.match(layout.style.gridTemplateColumns, /500px$/);
    window.dispatchEvent({type: 'pointerup', clientX: 750});
    assert.equal(storage.get(prefix + '-ui-right-width'), '450');
    idle();
    mount(); tick();
    assert.match(layout.style.gridTemplateColumns, /450px$/);
    for (const type of ['pointercancel', 'blur']) {
        start(); window.dispatchEvent({type}); idle();
        assert.equal(storage.get(prefix + '-ui-right-width'), '450');
    }
    start(); mount();
    window.dispatchEvent({type: 'pointermove', clientX: 500});
    idle();
    assert.equal(storage.get(prefix + '-ui-right-width'), '450');
    """
    path = Path(__file__).parents[1] / "assets" / asset
    subprocess.run([node, "-e", script, str(path), prefix], check=True, capture_output=True, text=True, timeout=10)
