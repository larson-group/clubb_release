from compile import update_install_aliases


def test_compile_updates_both_aliases_by_default(tmp_path):
    install = tmp_path / "install"
    install.mkdir()
    old = install / "old"
    new = install / "new"
    old.mkdir()
    new.mkdir()
    (install / "latest").symlink_to(old)
    (install / "selected").symlink_to(old)

    update_install_aliases(str(new), install_root=str(install))

    assert (install / "latest").resolve() == new
    assert (install / "selected").resolve() == new


def test_dash_compile_leaves_selected_alone(tmp_path):
    install = tmp_path / "install"
    install.mkdir()
    old = install / "old"
    new = install / "new"
    old.mkdir()
    new.mkdir()
    (install / "latest").symlink_to(old)
    (install / "selected").symlink_to(old)

    update_install_aliases(str(new), update_selected=False, install_root=str(install))

    assert (install / "latest").resolve() == new
    assert (install / "selected").resolve() == old
