"""Exposure-specific, inspection-only handoffs from a nightly log to Qt."""
import json
from pathlib import Path
import shlex
from urllib.parse import parse_qs, urlencode, urlsplit


def issue_updates(run):
    """Read explicitly selected follow-up notes bound to this run and raw data."""
    from .inventory import sha256
    run=Path(run).resolve();qc=run.parent/'QC'
    path=qc/'issue-updates.json'
    if not path.exists():return []
    record=json.loads(path.read_text());parent_hash=sha256(run/'run.json')
    if record['parent_sha256']!=parent_hash:raise ValueError('issue updates belong to another reduction')
    for item in record['updates']:
        manifest=Path(item['manifest']).resolve()
        if not manifest.is_relative_to(qc) or sha256(manifest)!=item['manifest_sha256']:
            raise ValueError('issue update evidence changed or lies outside the revision')
        evidence=json.loads(manifest.read_text())
        if (evidence['parent_sha256']!=parent_hash or Path(evidence['parent_run']).resolve()!=run
                or evidence['raw_sha256']!=item['raw_sha256'] or evidence['target']!=item['object']
                or sha256(Path(item['raw']))!=item['raw_sha256']):
            raise ValueError('issue update input binding failed')
    return record['updates']


def review_url(run, raw):
    return 'efospec://review?' + urlencode(dict(run=str(Path(run).resolve()), raw=Path(raw).name))


def parse_review_url(value):
    url = urlsplit(value)
    query = parse_qs(url.query, strict_parsing=True)
    if (url.scheme != 'efospec' or url.netloc != 'review' or url.path not in ('', '/')
            or url.fragment or set(query) != {'run', 'raw'}
            or any(len(v) != 1 for v in query.values())):
        raise ValueError('Expected an efospec://review link with one run and raw exposure.')
    run, raw = query['run'][0], query['raw'][0]
    if not Path(run).is_absolute() or Path(raw).name != raw or raw in {'.', '..'}:
        raise ValueError('Review links require an absolute run and an exposure filename.')
    return run, raw


def review_command(run, raw):
    return shlex.join(['efospec', 'gui', str(Path(run).resolve()), '--raw', Path(raw).name])


def counts_delivery_status(counts, manifest):
    """Name the recorded blocking stage instead of implying a missing response."""
    if any(item.get('path') == counts['raw'] and item.get('code') == 'wavelength_validation_failed'
           for item in manifest.get('errors', [])):
        return 'Counts extracted; flux withheld — wavelength validation failed'
    return counts.get('delivery_status', 'Review counts; wavelength unvalidated')


def failures(manifest):
    items = manifest.get('errors', []) + manifest.get('skipped', [])
    recorded = {Path(item.get('path', '')).name for item in items}
    for counts in manifest.get('review_counts', []):
        if Path(counts['raw']).name not in recorded:
            items.append(dict(path=counts['raw'], object=counts['object'], setup=counts['setup'],
                code='calibration_pending', stage='counts_extracted', counts_review=counts['product'],
                reason='Counts extracted; flux calibration pending. '+' '.join(counts.get('pending_checks', []))))
    return items


def failure_reason(item):
    from .sky_quality import explained_sky_failure
    return explained_sky_failure(item) or item.get('reason') or item.get('code', 'Reduction failed')


def recovery_advice(item):
    code = item.get('code')
    if code == 'calibration_pending':
        return 'Review the saved counts and the listed calibration checks. Aperture edits do not validate a response or wavelength correction.'
    if code == 'background_support_required':
        return ('Inspect the final trace, aperture and rejected background pixels. Refit clean bands on both '
                'sides of the final trace, then rerun extraction and support checks in a new revision.')
    if code == 'wavelength_validation_failed':
        return ('Inspect the arc identification, sky-line fits and independent residuals. Background edits '
                'alone do not validate wavelengths. Review counts retain unshifted arc wavelengths.')
    if code == 'unsupported_setup':
        return ('This setup is outside the selected reduction route. Check matching arcs, flats and standards '
                'and the separate red/filtered routes; an aperture edit cannot add setup support.')
    return 'Inspect source identity, trace and calibration diagnostics before extracting a new revision.'


def failure_log(run, item):
    """Keep recorded reasons available even when an archived log has moved."""
    lines = [item.get('object', 'Unknown'), item.get('path', ''),
             item.get('code', 'Reduction failed'), failure_reason(item), recovery_advice(item)]
    if item.get('traceback'):
        path = Path(item['traceback'])
        try:
            path.resolve().relative_to(Path(run).resolve())
            lines += ['\nOriginal error log', path.read_text()]
        except (ValueError, OSError) as exc:
            lines += ['Original error log unavailable: ' + str(exc)]
    lines += ['\nSaved failure record', json.dumps(item, indent=2)]
    return '\n'.join(lines)
