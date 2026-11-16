#!/usr/bin/env python3
"""Find or acquire LES benchmarks for loss evaluation and Dash plots.

find_les_data(file) returns an existing file or None without downloads or writes.
confirm_les_data(file) also downloads a missing NetCDF into a writable data
directory. file may be relative to the benchmark collection, or an
existing/configured absolute path. main() exposes acquisition from the command
line, for example:

    python3 input/scripts/confirm_les_data.py les_runs/wangara_rams.nc

Existing files are sought in input/les_and_clubb_benchmark_runs, /home/pub and
~/les_and_clubb_benchmark_runs. Missing files are added to an existing writable
data directory, falling back to input/les_and_clubb_benchmark_runs. When this
fallback is needed, an existing input symlink is replaced with a local directory;
its shared target is untouched. Unrelated explicit missing paths fail without
downloading.

The public catalog describes run ZIPs. HTTP byte ranges let ZipFile read the
archive index and the requested member, without transferring the other run
members. CRC/size and NetCDF checks precede atomic publication. A local per-file
lock beside the input collection serializes acquisition by workers in this
checkout, including while an input symlink is replaced. Existing files are
reused as-is; only new downloads are validated here. No source data is modified
by benchmark normalization, which remains a separate caller-owned step.
"""

from __future__ import annotations

import argparse
import errno
import fcntl
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import sys
import tempfile
from urllib.parse import urljoin, urlsplit
from urllib.request import Request, urlopen
import zipfile

if __name__ == '__main__':
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from utilities.setup_python_venv import ensure_python_venv

    ensure_python_venv()

COLLECTION = 'les_and_clubb_benchmark_runs'
REPO_ROOT = Path(__file__).resolve().parents[2]
INPUT_ROOT = REPO_ROOT / 'input' / COLLECTION
SHARED_ROOTS = (Path('/home/pub') / COLLECTION, Path.home() / COLLECTION)
BASE_URL = 'https://carson.math.uwm.edu/sites/benchmark-data/'
TIMEOUT_SECONDS = 30
CHUNK_BYTES = 1024 * 1024


class _RangeFile(io.RawIOBase):
    """Give ZipFile a seekable HTTP archive while transferring only requested ranges."""

    def __init__(self, url, expected_size):
        """Read archive size/version headers before any partial content requests."""
        super().__init__()
        self.url, self.position = url, 0
        with urlopen(Request(url, method='HEAD'), timeout=TIMEOUT_SECONDS) as response:
            self.size = int(response.headers['Content-Length'])
            self.etag = response.headers.get('ETag')
        if self.size != expected_size:
            raise ValueError('The catalog and archive sizes differ; the published data may be updating.')

    def readable(self):
        """Support ZipFile's binary reads."""
        return True

    def seekable(self):
        """Support ZIP index and member-offset lookups."""
        return True

    def tell(self):
        """Return the current archive byte offset."""
        return self.position

    def seek(self, offset, whence=os.SEEK_SET):
        """Move the logical cursor without downloading intervening archive content."""
        origins = {os.SEEK_SET: 0, os.SEEK_CUR: self.position, os.SEEK_END: self.size}
        if whence not in origins or origins[whence] + offset < 0:
            raise ValueError('Invalid archive seek')
        self.position = origins[whence] + offset
        return self.position

    def read(self, size=-1):
        """Require an exact partial response; never accept an ignored Range request."""
        count = max(0, self.size - self.position)
        if size >= 0:
            count = min(count, size)
        if not count:
            return b''
        end = self.position + count - 1
        headers = {'Range': f'bytes={self.position}-{end}'}
        if self.etag:
            headers['If-Match'] = self.etag
        with urlopen(Request(self.url, headers=headers), timeout=TIMEOUT_SECONDS) as response:
            expected_range = f'bytes {self.position}-{end}/{self.size}'
            if response.status != 206 or response.headers.get('Content-Range') != expected_range:
                raise ValueError('The download server must honor HTTP byte ranges; whole-run downloads are disabled.')
            if self.etag and response.headers.get('ETag') != self.etag:
                raise ValueError('The archive changed during acquisition; please try again.')
            data = response.read(count + 1)
        if len(data) != count:
            raise ValueError('Incomplete archive byte range')
        self.position += count
        return data


def _relative_file(file, roots):
    """Identify a collection member without interpreting an unrelated absolute path."""
    path = Path(file).expanduser()
    if path.is_absolute():
        for root in roots:
            try:
                relative = path.relative_to(root)
                break
            except ValueError:
                continue
        else:
            raise FileNotFoundError(f'LES input does not exist and is outside the benchmark data roots: {path}')
    else:
        relative = path
        for prefix in (Path('input') / COLLECTION, Path(COLLECTION)):
            if path.is_relative_to(prefix):
                relative = path.relative_to(prefix)
                break
    if '..' in relative.parts or relative.suffix.lower() != '.nc':
        raise ValueError('Request a NetCDF .nc path within les_and_clubb_benchmark_runs.')
    return PurePosixPath(*relative.parts).as_posix()


def _catalog_member(relative):
    """Match the exact requested path to one published run and its member size."""
    catalog_url = urljoin(BASE_URL, 'catalog.json')
    with urlopen(catalog_url, timeout=TIMEOUT_SECONDS) as response:
        catalog = json.load(response)
    matches = [(run, member) for run in catalog['runs'] for member in run['members']
               if member['path'] == relative]
    if len(matches) != 1:
        raise ValueError(f'Expected one published entry for {relative}, found {len(matches)} in {catalog_url}.')
    run, member = matches[0]
    if not run['ready']:
        raise ValueError(f'The published run for {relative} is not ready.')
    url = urljoin(BASE_URL, run['url'])
    if urlsplit(url).netloc != urlsplit(BASE_URL).netloc or not url.startswith(BASE_URL):
        raise ValueError('The catalog archive URL is outside the benchmark download site.')
    return url, int(run['zip_bytes']), int(member['bytes'])


def _destination(relative, roots):
    """Try existing data roots, then a local input directory without changing shared targets."""
    candidates = [root.resolve() for root in roots if root.is_dir()]
    candidates.append(INPUT_ROOT)
    seen = set()
    for root in candidates:
        if root == INPUT_ROOT and root.is_symlink():
            try:
                root.unlink()
            except (FileNotFoundError, IsADirectoryError):
                # Another acquisition may already have replaced this symlink.
                pass
        resolved = root.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        target = root / relative
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            descriptor, temporary = tempfile.mkstemp(prefix=f'.{target.name}.', suffix='.part', dir=target.parent)
            return target, Path(temporary), descriptor
        except OSError as exc:
            if exc.errno not in {errno.EACCES, errno.EPERM, errno.EROFS}:
                raise
    raise PermissionError(f'No writable LES destination for {relative}; {INPUT_ROOT} must be writable.')


def _download(relative, target, temporary, descriptor):
    """Stream one ZIP member, check its CRC/size and NetCDF header, then publish it."""
    try:
        with os.fdopen(descriptor, 'wb') as output:
            url, archive_bytes, file_bytes = _catalog_member(relative)
            print(f'[LES] Downloading {relative} -> {target}', flush=True)
            with _RangeFile(url, archive_bytes) as remote, zipfile.ZipFile(remote) as archive:
                member = archive.getinfo(f'{COLLECTION}/{relative}')
                if member.file_size != file_bytes:
                    raise ValueError('The catalog and ZIP member sizes differ.')
                written, next_progress = 0, 64 * CHUNK_BYTES
                with archive.open(member) as source:
                    while chunk := source.read(CHUNK_BYTES):
                        output.write(chunk)
                        written += len(chunk)
                        if written >= next_progress:
                            print(f'[LES] {relative}: {written / file_bytes:.0%}', flush=True)
                            next_progress += 64 * CHUNK_BYTES
                if written != file_bytes:
                    raise ValueError('The downloaded NetCDF size is incomplete.')
            output.flush()
            os.fsync(output.fileno())
        from netCDF4 import Dataset

        with Dataset(temporary, 'r'):
            pass
        temporary.chmod(0o644)
        os.replace(temporary, target)
        print(f'[LES] Ready: {target}', flush=True)
        return target.resolve()
    finally:
        temporary.unlink(missing_ok=True)


def find_les_data(file):
    """Return an existing requested file, or None, without creating files or downloading."""
    requested = Path(file).expanduser()
    if requested.is_file():
        return requested.resolve()
    roots = (INPUT_ROOT, *SHARED_ROOTS)
    relative = _relative_file(file, roots)
    for root in roots:
        if (root / relative).is_file():
            return (root / relative).resolve()
    return None


def confirm_les_data(file):
    """Return the requested LES file, acquiring it only when all configured sources lack it."""
    existing = find_les_data(file)
    if existing is not None:
        return existing
    roots = (INPUT_ROOT, *SHARED_ROOTS)
    relative = _relative_file(file, roots)

    lock_dir = INPUT_ROOT.parent / f'.{COLLECTION}.locks'
    lock_dir.mkdir(parents=True, exist_ok=True)
    lock_path = lock_dir / (hashlib.sha256(relative.encode()).hexdigest() + '.lock')
    with lock_path.open('a') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        # A sibling worker may have published this member while we waited.
        for root in roots:
            if (root / relative).is_file():
                return (root / relative).resolve()
        try:
            target, temporary, descriptor = _destination(relative, roots)
            return _download(relative, target, temporary, descriptor)
        except Exception as exc:
            raise RuntimeError(f'Could not acquire LES file {relative}: {exc}') from exc


def main(argv=None):
    """Ensure each requested NetCDF file exists and print its resolved location."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0], add_help=False, allow_abbrev=False)
    parser.add_argument('-h', '-help', action='help', help='Show this help and exit.')
    parser.add_argument('files', nargs='+', help='NetCDF paths relative to the collection, or configured input paths.')
    args = parser.parse_args(argv)
    try:
        for file in args.files:
            print(confirm_les_data(file))
    except (OSError, RuntimeError, ValueError) as exc:
        parser.exit(1, f'{exc}\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
