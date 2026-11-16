"""Exercise on-demand LES reuse, selective HTTP ZIP reads and loss preparation.

The existing clubb_python/-unit stage owns these tests. All network responses
come from an in-memory catalog/archive, with real ZIP CRC and NetCDF validation.
No model, external dataset or public-server access is required.
"""

from concurrent.futures import ThreadPoolExecutor
import errno
import io
import json
from pathlib import Path
import zipfile

from netCDF4 import Dataset
import numpy as np
import pytest

from input.scripts import confirm_les_data as acquisition


@pytest.fixture
def roots(tmp_path, monkeypatch):
    """Give every test its own input directory and shared data sources."""
    monkeypatch.setattr(acquisition, 'INPUT_ROOT', tmp_path / 'input' / acquisition.COLLECTION)
    monkeypatch.setattr(acquisition, 'SHARED_ROOTS', (tmp_path / 'shared', tmp_path / 'home-data'))
    monkeypatch.setattr(acquisition, 'BASE_URL', 'https://benchmarks.example/data/')
    return acquisition.INPUT_ROOT, *acquisition.SHARED_ROOTS


@pytest.fixture
def remote(tmp_path, monkeypatch, roots):
    """Publish one tiny SAM field in a ZIP containing a large unrelated member."""
    source = tmp_path / 'source.nc'
    with Dataset(source, 'w', format='NETCDF3_CLASSIC') as ds:
        ds.createDimension('time', 2)
        ds.createDimension('z', 3)
        t = ds.createVariable('time', 'f8', ('time',))
        t[:] = [0, 60]
        t.units = 'seconds'
        z = ds.createVariable('z', 'f8', ('z',))
        z[:] = [10, 100, 200]
        z.units = 'm'
        thl = ds.createVariable('THETAL', 'f8', ('time', 'z'))
        thl[:] = np.arange(6).reshape(2, 3) + 290
        thl.units = 'K'
        ds.createVariable('CLD', 'f8', ('time', 'z'))[:] = 0
    relative = 'sam_benchmark_runs/example/stats.nc'
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(acquisition.COLLECTION + '/' + relative, source.read_bytes())
        archive.writestr(acquisition.COLLECTION + '/unrelated.bin', np.random.default_rng(7).bytes(2 * 1024**2))
    state = {
        'relative': relative, 'source': source.read_bytes(), 'archive': buffer.getvalue(),
        'requests': [], 'transferred': 0, 'ignore_range': False, 'short_range': False,
        'changed_archive': False,
    }
    state['catalog'] = {'runs': [{
        'ready': True, 'url': 'downloads/run.zip', 'zip_bytes': len(state['archive']),
        'members': [{'path': relative, 'bytes': len(state['source'])}],
    }]}

    class Response(io.BytesIO):
        """Expose urllib's response metadata while serving bounded in-memory bytes."""
        def __init__(self, data, status=200, headers=None):
            super().__init__(data)
            self.status, self.headers = status, headers or {}

    def urlopen(request, timeout):
        """Honor Range and If-Match, recording bytes requested from the archive."""
        assert timeout > 0
        url = request if isinstance(request, str) else request.full_url
        method = 'GET' if isinstance(request, str) else request.get_method()
        state['requests'].append((url, method))
        if url.endswith('catalog.json'):
            return Response(json.dumps(state['catalog']).encode())
        assert url == acquisition.BASE_URL + 'downloads/run.zip'
        if method == 'HEAD':
            return Response(b'', headers={'Content-Length': str(len(state['archive'])), 'ETag': '"version1"'})
        assert request.get_header('If-match') == '"version1"'
        first, last = map(int, request.get_header('Range').removeprefix('bytes=').split('-'))
        data = state['archive'][first:last + 1]
        state['transferred'] += len(data)
        if state['short_range']:
            data = data[:-1]
        return Response(data, status=200 if state['ignore_range'] else 206, headers={
            'Content-Range': f'bytes {first}-{last}/{len(state["archive"])}',
            'ETag': '"version2"' if state['changed_archive'] else '"version1"',
        })

    monkeypatch.setattr(acquisition, 'urlopen', urlopen)
    return state


def test_existing_files_are_used_without_network_or_input_creation(roots, monkeypatch):
    """A file in any shared root is reused before checking destination writability."""
    requested = roots[2] / 'les_runs/existing.nc'
    requested.parent.mkdir(parents=True)
    requested.write_bytes(b'existing data is the converter\'s responsibility')
    monkeypatch.setattr(acquisition, 'urlopen', lambda *_a, **_k: pytest.fail('unexpected network'))
    assert acquisition.confirm_les_data(roots[0] / 'les_runs/existing.nc') == requested
    assert not roots[0].exists()
    assert not (roots[0].parent / f'.{acquisition.COLLECTION}.locks').exists()


def test_missing_file_is_downloaded_to_new_input_directory(roots, remote):
    """Only the chosen member is published, with no unrelated run data transferred."""
    actual = acquisition.confirm_les_data(remote['relative'])
    assert actual == roots[0] / remote['relative']
    assert actual.read_bytes() == remote['source']
    assert remote['transferred'] < len(remote['archive']) // 4
    assert not any(actual.parent.glob('*.part'))
    with Dataset(actual) as ds:
        np.testing.assert_array_equal(ds['THETAL'][:], np.arange(6).reshape(2, 3) + 290)


@pytest.mark.parametrize('input_symlink', [False, True])
def test_writable_shared_root_gets_just_the_missing_file(roots, remote, input_symlink):
    """Writable shared data accept the missing file without replacing a usable symlink."""
    roots[1].mkdir()
    if input_symlink:
        roots[0].parent.mkdir(parents=True)
        roots[0].symlink_to(roots[1], target_is_directory=True)
    actual = acquisition.confirm_les_data(remote['relative'])
    assert actual == roots[1] / remote['relative']
    assert actual.read_bytes() == remote['source']
    if input_symlink:
        assert roots[0].is_symlink() and roots[0].resolve() == roots[1]
    else:
        assert not roots[0].exists()


def test_unwritable_input_fails_without_using_another_cache(roots, remote, monkeypatch):
    """An unwritable local fallback fails before downloading instead of writing elsewhere."""
    def mkstemp(*args, **kwargs):
        """Model permission denial for the local destination."""
        raise PermissionError(errno.EACCES, 'input directory is not writable')

    monkeypatch.setattr(acquisition.tempfile, 'mkstemp', mkstemp)
    with pytest.raises(RuntimeError, match='must be writable'):
        acquisition.confirm_les_data(remote['relative'])
    assert not remote['requests']
    assert not (roots[0] / remote['relative']).exists()


@pytest.mark.parametrize('input_symlink', [False, True])
def test_unwritable_shared_root_uses_input_directory(roots, remote, monkeypatch, input_symlink):
    """Read-only shared data use a real input directory; existing shared files remain reusable."""
    roots[1].mkdir()
    existing = roots[1] / 'existing.nc'
    existing.write_bytes(b'untouched shared data')
    if input_symlink:
        roots[0].parent.mkdir(parents=True)
        roots[0].symlink_to(roots[1], target_is_directory=True)
    original = acquisition.tempfile.mkstemp

    def mkstemp(*args, **kwargs):
        """Model a shared read-only mount independently of the test user's privileges."""
        if Path(kwargs['dir']).resolve().is_relative_to(roots[1]):
            raise OSError(errno.EROFS, 'read-only shared data')
        return original(*args, **kwargs)

    monkeypatch.setattr(acquisition.tempfile, 'mkstemp', mkstemp)
    actual = acquisition.confirm_les_data(roots[0] / remote['relative'])
    assert actual == roots[0] / remote['relative']
    assert roots[0].is_dir() and not roots[0].is_symlink()
    assert actual.read_bytes() == remote['source']
    assert not (roots[1] / remote['relative']).exists()
    before = len(remote['requests'])
    assert acquisition.confirm_les_data(roots[0] / remote['relative']) == actual
    assert acquisition.confirm_les_data(roots[0] / 'existing.nc') == existing
    assert existing.read_bytes() == b'untouched shared data'
    assert len(remote['requests']) == before


def test_broken_input_symlink_is_replaced_with_local_directory(roots, remote):
    """A dangling input symlink becomes a real directory without creating its old target."""
    roots[0].parent.mkdir(parents=True)
    roots[0].symlink_to(roots[1], target_is_directory=True)
    actual = acquisition.confirm_les_data(roots[0] / remote['relative'])
    assert actual == roots[0] / remote['relative']
    assert roots[0].is_dir() and not roots[0].is_symlink()
    assert not roots[1].exists()


@pytest.mark.parametrize('input_symlink', [False, True])
def test_concurrent_workers_acquire_one_copy(roots, remote, input_symlink):
    """Two callers share one download even when acquisition replaces the input symlink."""
    if input_symlink:
        roots[0].parent.mkdir(parents=True)
        roots[0].symlink_to(roots[1], target_is_directory=True)
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(acquisition.confirm_les_data, [remote['relative']] * 2))
    assert results[0] == results[1]
    assert results[0].read_bytes() == remote['source']
    assert sum(url.endswith('catalog.json') for url, _method in remote['requests']) == 1


@pytest.mark.parametrize('failure', ['ignore_range', 'short_range', 'changed_archive', 'bad_catalog_size', 'invalid_netcdf'])
def test_failed_downloads_never_publish_partial_files(roots, remote, failure):
    """Range/version/format failures leave neither a target nor a temporary download."""
    if failure == 'bad_catalog_size':
        remote['catalog']['runs'][0]['zip_bytes'] += 1
    elif failure == 'invalid_netcdf':
        data = b'<html>This is not a NetCDF file</html>'
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr(acquisition.COLLECTION + '/' + remote['relative'], data)
        remote['archive'] = buffer.getvalue()
        remote['catalog']['runs'][0]['zip_bytes'] = len(remote['archive'])
        remote['catalog']['runs'][0]['members'][0]['bytes'] = len(data)
    else:
        remote[failure] = True
    with pytest.raises(RuntimeError, match='Could not acquire LES file'):
        acquisition.confirm_les_data(remote['relative'])
    assert not (roots[0] / remote['relative']).exists()
    assert not list(roots[0].rglob('*.part'))


def test_unpublished_path_fails_without_substituting_a_different_benchmark(roots, remote):
    """Missing catalog paths fail explicitly, even if the requested basename exists elsewhere."""
    with pytest.raises(RuntimeError, match='found 0'):
        acquisition.confirm_les_data('les_runs/stats.nc')
    assert not any(method == 'HEAD' for _url, method in remote['requests'])
    assert not list(roots[0].rglob('*.part'))


@pytest.mark.parametrize('path', ['../escape.nc', '/elsewhere/missing.nc', 'notes.txt'])
def test_unrelated_missing_paths_do_not_trigger_network(roots, monkeypatch, path):
    """Only collection-relative/configured NetCDF paths are eligible for acquisition."""
    monkeypatch.setattr(acquisition, 'urlopen', lambda *_a, **_k: pytest.fail('unexpected network'))
    with pytest.raises((ValueError, FileNotFoundError)):
        acquisition.confirm_les_data(path)
    assert not roots[0].exists()


def test_loss_normalization_uses_the_acquired_file(roots, remote, tmp_path):
    """The shared loss preparation converts downloaded SAM fields without changing case defaults."""
    from input.scripts.create_case_namelist import create_normalized_benchmark_file

    defaults = {'les_stats_file': str(roots[0] / remote['relative'])}
    original = dict(defaults)
    normalized = create_normalized_benchmark_file('synthetic', defaults, ['thlm'], tmp_path / 'job')
    assert defaults == original
    assert normalized.is_file()
    with Dataset(normalized) as ds:
        assert ds['thlm'].shape == (2, 3, 1, 1)
        np.testing.assert_array_equal(ds['thlm'][:, :, 0, 0], np.arange(6).reshape(2, 3) + 290)
