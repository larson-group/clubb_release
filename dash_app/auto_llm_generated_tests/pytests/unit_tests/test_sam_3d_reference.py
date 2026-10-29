from netCDF4 import Dataset

from utilities.sam_3d_reference import inventory_run


def _write_snapshot(path, timestep, seconds_per_step):
    with Dataset(path / f"case_000000{timestep:04d}_micro.nc", "w") as dataset:
        dataset.createDimension("time", 1)
        dataset.createDimension("z", 2)
        time = dataset.createVariable("time", "f8", ("time",))
        time.units = "days"
        time[:] = timestep * seconds_per_step / 86400.0
        z = dataset.createVariable("z", "f8", ("z",))
        z[:] = [25.0, 75.0]


def test_inventory_derives_physical_time_from_snapshot_metadata(tmp_path):
    output = tmp_path / "OUT_3D"
    output.mkdir()
    for timestep in (1500, 1510, 3600):
        _write_snapshot(output, timestep, 6.0)

    inventory = inventory_run(tmp_path)

    assert inventory.steps_seconds == [9000, 9060, 21600]
    assert inventory.first_minute == 150.0
    assert inventory.last_minute == 360.0
