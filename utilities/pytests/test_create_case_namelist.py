import pytest

from utilities.create_case_namelist import override_value


def test_bare_override_replaces_existing_value():
    result = override_value("dt_main=2.0", "&model_setting\n  dt_main = 1.0,\n/\n")

    assert "dt_main = 2.0" in result


def test_grouped_override_adds_missing_value():
    result = override_value(
        "microphysics_setting.lh_num_samples=8",
        "&microphysics_setting\n  microphys_scheme = 'none',\n/\n",
    )

    assert "lh_num_samples = 8," in result


def test_grouped_override_only_updates_its_namelist():
    result = override_value(
        "second_setting.value=3",
        "&first_setting\n  value = 1,\n/\n&second_setting\n  value = 2,\n/\n",
    )

    assert "&first_setting\n  value = 1," in result
    assert "&second_setting\n  value = 3," in result


def test_missing_bare_override_explains_how_to_add_a_key():
    with pytest.raises(ValueError, match="could not be matched for find-and-replace") as exc_info:
        override_value("lh_num_samples=8", "&microphysics_setting\n/\n")

    assert "[namelist.]key=value" in str(exc_info.value)
