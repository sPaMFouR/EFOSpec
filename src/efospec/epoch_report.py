"""Read-only nightly reports assembled from one epoch's saved evidence."""
import argparse
import base64
from collections import Counter
from datetime import date
from html import escape
import json
import math
from io import BytesIO
from pathlib import Path

from .inventory import sha256
from .naming import source_order
from .report import CSS
from .review_handoff import failure_reason, recovery_advice, review_url, review_command


def _text(value):
    return escape(str(value))


def _unvalidated(status):
    return any(word in str(status).lower() for word in
               ('fail', 'unvalidated', 'unsupported', 'unavailable', 'missing', 'unknown'))


def _table(headings, rows):
    return ('<div class="table-wrap"><table><thead><tr>'
            + ''.join('<th>'+_text(h)+'</th>' for h in headings)
            + '</tr></thead><tbody>'
            + ''.join('<tr>'+''.join('<td>'+str(c)+'</td>' for c in row)+'</tr>' for row in rows)
            + '</tbody></table></div>')


class EpochPage:
    def __init__(self, root):
        self.root = root.resolve()
        self.warnings = []
        self.images = 0

    def read(self, path):
        if not path.exists():
            return {}
        try:
            return json.loads(path.read_text())
        except (OSError, ValueError) as exc:
            self.warnings.append(f'{path.name}: {exc}')
            return {}

    def data(self, path, mime):
        # Embed only this date's saved artifacts, including relocated products.
        path.resolve().relative_to(self.root)
        return 'data:'+mime+';base64,'+base64.b64encode(path.read_bytes()).decode('ascii')

    def download(self, path, label, *, filename=None):
        if not path.is_file():
            return ''
        mime = {'.json':'application/json', '.csv':'text/csv', '.asci':'text/plain',
                '.log':'text/plain', '.pdf':'application/pdf'}.get(path.suffix, 'application/octet-stream')
        return f'<a download="{_text(filename or path.name)}" href="{self.data(path, mime)}">{_text(label)}</a>'

    def figure(self, path, title, *, full_width=False, aspect_ratio=None):
        if not path or not path.is_file():
            return '<p class="note">'+_text(title)+': figure unavailable.</p>'
        self.images += 1
        # Small web previews retain the existing figure layout and styling.
        # Publication PDF/JPG files are read-only and remain at full resolution.
        from PIL import Image
        path.resolve().relative_to(self.root)
        with Image.open(path) as source:
            preview = source.convert('RGB')
            preview.thumbnail((1600, 1600), Image.Resampling.LANCZOS)
            encoded = BytesIO()
            preview.save(encoded, format='JPEG', quality=85, optimize=True)
        url = 'data:image/jpeg;base64,'+base64.b64encode(encoded.getvalue()).decode('ascii')
        css_class = 'plot full-width' if full_width else 'plot'
        button_class = 'enlarge fixed-ratio' if aspect_ratio is not None else 'enlarge'
        sizing = f' style="aspect-ratio:{aspect_ratio:g}"' if aspect_ratio is not None else ''
        return (f'<figure class="{css_class}"><figcaption>{_text(title)}</figcaption>'
                f'<button class="{button_class}"{sizing} type="button" aria-label="Enlarge {_text(title)}">'
                f'<img loading="lazy" src="{url}" alt="{_text(title)}"></button></figure>')

    def failure_notice(self, run_dir, item):
        from .review_handoff import issue_updates
        update=next((u for u in issue_updates(run_dir) if u['raw']==item.get('path')),None)
        raw = Path(item.get('path', '')).name
        launch = (f'<button type="button" class="review-launch" data-review-url="{_text(review_url(run_dir, raw))}">Open in EFOSpec GUI</button>'
                  if raw else '')
        log = self.download(Path(item['traceback']), 'Error log') if item.get('traceback') else ''
        current=('<p class="qc-status-bad"><strong>Latest retry:</strong> '+_text(update['summary'])+'</p>'
                 if update else '')
        return (current+'<p class="qc-status-bad">'+('Original attempt — ' if update else '')+_text(item.get('code', 'Reduction failed'))
                +': '+_text(failure_reason(item))+'</p><p>'+_text(recovery_advice(item))+'</p>'
                +'<p>'+ ' · '.join(filter(None, [launch, log]))+'</p>'
                +('<details class="review-help"><summary>Desktop launch help</summary><p class="launch-status" role="status"></p><p>Requires EFOSpec on this computer '
                   'and access to the original run. If this viewer blocks desktop apps, open the report in Safari or run this command in Terminal:</p><pre>'
                   +_text(review_command(run_dir, raw))+'</pre></details>' if raw else ''))

    def exposure_time(self, result, frames):
        raw = Path(result.get('raw', result.get('path', ''))).name
        seconds = result.get('exptime_s', frames.get(raw, {}).get('metadata', {}).get('EXPTIME'))
        try:
            seconds = float(seconds)
            if not math.isfinite(seconds) or seconds < 0:return ''
        except (ValueError, TypeError):return ''
        return f' · {seconds:g} s'

    def spectra(self, run_dir, run, *, id_prefix='spectrum'):
        from .telluric_products import selected_pairs,overlay_record,status_text
        pairs=selected_pairs(run_dir)
        records = [self.read(p) for p in sorted((run_dir/'style-preview').glob('FitSpectrum_*.json'))]
        frames = {Path(f['path']).name: f for f in self.read(run_dir.parent/'observing-log.json').get('frames', [])}
        sections = []; section_order = []
        for i, result in sorted(enumerate(run.get('results', [])), key=lambda item: source_order(item[1])):
            product = run_dir/'products'/Path(result.get('product', '')).name
            digest = sha256(product) if product.is_file() else None
            matches = [r for r in records if digest and any(
                Path(p).name == product.name and h == digest for p, h in r.get('input_sha256', {}).items())]
            record = matches[-1] if matches else {}
            outputs = [run_dir/'style-preview'/Path(p).name for p in record.get('outputs', [])]
            pair=pairs.get(product.name);telluric_notice=''
            if pair:
                overlay=overlay_record(pair)
                if overlay:
                    record=overlay;outputs=[Path(p) for p in overlay['outputs']]
                else:
                    telluric_notice='<p class="qc-status-bad">Before/after telluric plot has not been generated.</p>'
                telluric_notice+='<p class="qc-status-bad">'+_text(status_text(pair))+'</p>'
                telluric_notice+='<p class="caption">Both reductions use the same flux scale. Shaded atmospheric regions are withheld from the corrected product.</p>'
            title = ' · '.join(filter(None, [result.get('object', 'Unknown'), result.get('grism'), result.get('slit')]))
            title += self.exposure_time(result, frames)
            cards = []
            for prefix, name in [('FitSpectrum_', 'Spectrum'), ('FitTrace_', '2D spectrum before and after sky subtraction'), ('FitSpatialProfile_', 'Spatial profile')]:
                images = [p for p in outputs if p.name.startswith(prefix) and p.suffix == '.jpg']
                cards.append(self.figure(images[0] if len(images) == 1 else None, title+' / '+name))
            downloads = ' · '.join(filter(None, [self.download(product, 'FITS'), self.download(product.with_suffix('.asci'), 'ASCII')]))
            if pair:
                downloads=' · '.join(self.download(Path(pair[key]),label,
                    filename=Path(pair[key]).stem+suffix+Path(pair[key]).suffix)
                    for key,label,suffix in [('uncorrected_product','Before telluric FITS','_before_telluric'),
                        ('uncorrected_ascii','Before telluric ASCII','_before_telluric'),
                        ('product','After telluric FITS','_after_telluric'),('ascii','After telluric ASCII','_after_telluric')])
            association = result.get('target', {}).get('association')
            identity = '' if association in {'historical_aperture', 'operator_selected'} else (
                ' · <span title="The extracted detector trace has not been independently matched to the named target; extraction itself succeeded.">Target identity unverified</span>')
            if result.get('target',{}).get('identity',{}).get('method') == 'operator_visual':
                identity = ' · C'+str(result['target']['identity']['candidate'])+' identified visually by operator'
            if result.get('status')=='PROVISIONAL':
                identity += ' · <strong class="qc-status-bad">Provisional: wavelength unvalidated</strong>'
                identity += ' · '+_text(result.get('interpretation',''))
            section_order.append(source_order(result))
            sections.append(f'<article class="exposure" id="{_text(id_prefix)}-{i}"><h3>{_text(title)}</h3>'
                f'<p class="caption">{_text(result.get("date_obs", ""))} UTC{identity} · {downloads}</p>'
                +telluric_notice+
                f'<div class="plots">{"".join(cards)}</div></article>')
        for i, result in enumerate(run.get('review_counts',[])):
            for path,digest in result['product_sha256'].items():
                if sha256(Path(path))!=digest:raise ValueError('review counts changed')
            plots=Path(result['product']).parent.parent/'plots'
            title=' · '.join(filter(None, [result['object'], result['setup']['grism'], result['setup'].get('filter'), result['setup'].get('slit'), 'Review counts']))
            title += self.exposure_time(result, frames)
            links=' · '.join([self.download(Path(result['product']),'Counts FITS'),self.download(Path(result['ascii']),'Counts ASCII')])
            section_order.append(source_order(result))
            pending=result.get('pending_checks', ['Wavelength validation failed; unshifted arc wavelengths, no flux calibration.'])
            from .review_handoff import failures, counts_delivery_status
            handoff=next((item for item in failures(run) if item.get('path')==result['raw']), None)
            sections.append(f'<article class="exposure" id="{_text(id_prefix)}-counts-{i}"><h3>'+_text(title)+'</h3>'
                '<p class="qc-status-bad">'+_text(counts_delivery_status(result,run))+'. '+_text(' '.join(pending))+'</p>'
                +(self.failure_notice(run_dir,handoff) if handoff else '')
                +'<p>'+links+'</p><div class="plots">'+''.join(self.figure(plots/(kind+'.png'),title+' / '+kind)
                    for kind in ('counts','traces','profiles'))+'</div></article>')
        for group in ('errors', 'skipped'):
            for i, item in sorted(enumerate(run.get(group, [])), key=lambda pair: source_order(pair[1])):
                setup = item.get('setup', {})
                title = ' · '.join(filter(None, [item.get('object', 'Unknown'), setup.get('grism'),
                                                setup.get('filter'), setup.get('slit')]))
                title += self.exposure_time(item, frames)
                section_order.append(source_order(item))
                from .failure_profiles import load_failure_profile
                saved_profile=load_failure_profile(run_dir,item.get('path',''))
                profile=''
                if saved_profile:
                    folder,preview=saved_profile
                    profile=self.figure(folder/'profiles.png',title+' / Spatial profile before subtraction')
                    if preview.get('caption'):
                        profile+='<p class="caption">'+_text(preview['caption'])+'</p>'
                sections.append(f'<article class="exposure" id="{_text(id_prefix)}-{group}-{i}">'
                    f'<h3>{_text(title)}</h3>'+self.failure_notice(run_dir, item)
                    +profile
                    +'<p class="caption">'+_text(Path(item.get('path', '')).name)+'</p></article>')
        return ''.join(html for _, html in sorted(zip(section_order, sections), key=lambda pair: pair[0])) or '<p>No spectra recorded in this revision.</p>'

    def observations(self, revision, run, provisional=None, redshifts=None, manuals=()):
        targets = self.read(revision/'night.json').get('targets')
        outcomes = {}
        for group, label in [('results', 'Spectrum'), ('errors', 'Failed'), ('skipped', 'Unsupported')]:
            for i, result in enumerate(run.get(group, [])):
                raw = Path(result.get('raw', result.get('path', ''))).name
                reason = {'sky_flexure_unreliable':'Sky-line check', 'source_selection_required':'Choose source',
                          'background_support_required':'Review background',
                          'native_calibration_not_implemented':'Calibration route unavailable',
                          'unsupported_setup':'Unsupported setup'}.get(result.get('code'), label)
                outcomes[raw] = f'<a href="#spectrum-{i}">Spectrum ↓</a>' if group == 'results' else f'<a class="qc-status-bad" href="#spectrum-{group}-{i}">{_text(reason)} ↓</a>'
        for standard in run.get('standards', []):
            raw = Path(standard['frame']['path']).name
            outcomes.setdefault(raw, 'Standard counts; response unvalidated')
        for i,result in enumerate(run.get('review_counts',[])):
            raw = Path(result['raw']).name
            from .review_handoff import counts_delivery_status
            label=counts_delivery_status(result,run)
            outcomes[raw] = f'<a class="qc-status-bad" href="#spectrum-counts-{i}">{_text(label)} ↓</a>'
        for i,result in enumerate((provisional or {}).get('results',[])):
            outcomes[Path(result['raw']).name]=f'<a class="qc-status-bad" href="#rescue-spectrum-{i}">Provisional flux; wavelength unvalidated ↓</a>'
        for j,(_,manual) in enumerate(manuals):
            for i,result in enumerate(manual['results']):
                raw=Path(result['raw']).name
                outcomes[raw]=f'<a href="#manual-{j}-spectrum-{i}">Manual flux recovery (review) ↓</a><br>Original attempt: '+outcomes.get(raw,'')
        rows = []
        for frame in sorted(self.read(revision/'observing-log.json').get('frames', []), key=source_order):
            if frame.get('role') not in {'science', 'standard_candidate'} or 'duplicate_of' in frame:
                continue
            h, setup = frame['metadata'], frame['setup']
            raw = Path(frame['path']).name
            start = str(h.get('DATE-OBS', '')).replace('T', ' ')
            try:
                seconds = round(float(h.get('EXPTIME') or 0))
            except (ValueError, TypeError, OverflowError):
                seconds = h.get('EXPTIME', 'Unknown')
            from .redshift_report import redshift_cells
            cells=redshift_cells((redshifts or {}).get(raw),standard=frame['role']=='standard_candidate',
                                 previous_revision=raw in getattr(self,'previous_redshift_raws',set()))
            if frame['role']!='standard_candidate':
                from .typing_report import record_raw, ned_cell
                annotation=next((r for r in getattr(self,'typing_annotations',{}).values() if record_raw(r)==raw),None)
                if annotation:
                    catalogue_value=ned_cell(annotation)
                    if catalogue_value is not None:cells[0]=catalogue_value
            target = _text(h.get('OBJECT', ''))
            rows.append([target, _text(start),
                'Standard' if frame['role'] == 'standard_candidate' else 'Science',
                _text(setup.get('grism', ''))+' / '+_text(setup.get('filter', '')),
                _text(setup.get('slit', '')), _text(seconds),
                outcomes.get(raw, 'Outside selected targets' if targets and h.get('OBJECT') not in targets else 'Not reduced'),*cells])
        note = {
            'not_run':'Host-emission measurements are a separate post-reduction step; no saved run is available for this observing night.',
            'previous_revision':'Saved host-emission measurements belong to an earlier reduction revision and have not been reused. Current NED values may still be available from the saved catalogue associations.'
        }.get(getattr(self,'redshift_status',None), '')
        return (_table(['Target', 'UTC start', 'Role', 'Setup', 'Slit', 'Seconds', 'Result', 'z (NED)', 'z (Host Emission)'], rows)
                + ('<p class="caption">'+_text(note)+'</p>' if note else ''))

    def problems(self, run, run_dir=None):
        from .review_handoff import failures
        groups = {}
        for item in failures(run):
            code = item.get('code', 'unknown')
            groups.setdefault(code, []).append(item)
        rows = []
        for code, items in groups.items():
            names = Counter(x.get('object', 'Unknown') for x in sorted(items, key=source_order))
            targets = ', '.join(k+(f' ×{v}' if v > 1 else '') for k, v in names.items())
            detail = {'sky_flexure_unreliable':'Sky lines did not pass the flexure checks.',
                'source_selection_required':'Multiple traces; target selection needs review.',
                'background_support_required':'Insufficient clean two-sided background; inspect the saved trace geometry.',
                'wavelength_validation_failed':'Science sky-line validation failed; any review counts retain unshifted arc wavelengths.',
                'unsupported_setup':'The selected reduction route does not support the recorded setup.'}.get(code, items[0].get('reason', items[0].get('error', code)))
            rows.append([_text(code), len(items), _text(targets), _text(detail)])
        content=_table(['Issue', 'Exposures', 'Targets', 'Meaning'], rows) if rows else '<p>No recorded reduction failures.</p>'
        if run_dir is not None:
            from .review_handoff import issue_updates
            updates=issue_updates(run_dir)
            if updates:
                latest=''.join('<p class="qc-status-bad"><strong>'+_text(u['object'])+
                    ' — latest retry:</strong> '+_text(u['summary'])+'</p>' for u in updates)
                content=latest+'<h3>Original nightly attempts</h3>'+content
        if run_dir is not None:
            if run.get('flux_calibration') == 'same-night':
                from .sky_quality import explained_sky_failure
                for item in run.get('errors', []):
                    content+='<p class="qc-status-bad">'+_text(item.get('object','Exposure')+': '+explained_sky_failure(item))+'</p>'
                    if item.get('traceback'):
                        content += '<p>'+self.download(Path(item['traceback']),
                            str(item.get('object','Exposure'))+' failure log')+'</p>'
                    if item.get('diagnostics',{}).get('review_file'):
                        content+='<p>'+self.download(Path(item['diagnostics']['review_file']),'Sky-line measurements')+'</p>'
                for item in run.get('review_counts_errors',[]):
                    content+='<p class="qc-status-bad">'+_text(item['object']+' review extraction: '+item['error'])+'</p>'
                content += '<p>'+self.download(run_dir/'run.json','Daily response and extraction record')+'</p>'
            content+=self.wavelength_checks(run_dir)
            readiness=run_dir/'setup-checks.json'
            if readiness.exists():
                content+='<p>'+self.download(readiness,'Raw calibration readiness by exposure')+'</p>'
        return content

    def wavelength_checks(self, run_dir):
        content=''
        records=[(path,self.read(path)) for path in (run_dir/'checks').glob('*/check.json')]
        for path,check in sorted(records,key=lambda item:source_order({
                'object':item[1].get('metadata',{}).get('object',''), 'path':str(item[0])})):
            title=check.get('metadata',{}).get('object',path.parent.name)+' / '+path.parent.name
            status=check.get('status','unknown')
            line=check.get('line_check',{})
            reason=', '.join(line.get('failed_gates',[])) or (
                'Anchor checks passed' if status=='passed' else 'Wavelength correction is not validated')
            if status=='passed' and 'held_out_6363_pass' in line and line['held_out_6363_pass'] is not True:
                status='anchors passed; 6363 check unvalidated'
            comparison=check.get('correlation',{})
            description=('PypeIt correlation: '+str(comparison.get('status','not recorded'))
                +'; comparison only, not an applied correction.')
            css=' class="qc-unvalidated"' if _unvalidated(status) else ''
            content+='<details'+css+'><summary>'+_text(title)+' — '+_text(status)+'</summary>'
            content+='<p>'+_text(reason)+'. '+_text(description)+'</p>'
            content+=self.figure(path.parent/'wavelength.png',title+' / wavelength checks')
            content+='<p>'+self.download(path,'Fit measurements and provenance')+'</p></details>'
        return content

    def qc(self, revision, run):
        packets = {}
        # Current QC wins; retained original packets supply reviewed experiments.
        for folder in [self.root/'QC', revision/'QC']:
            for path in sorted(folder.glob('*/review.json')):
                review = self.read(path)
                if review.get('raw_name'):
                    packets[review['raw_name']] = (path.parent, review)
        current = {Path(e.get('path', '')).name for e in run.get('errors', [])}
        current.update(Path(r['raw']).name for r in run.get('results', []) if r.get('target',{}).get('reviewed_source'))
        sections, backgrounds = [], []
        for raw, (folder, review) in sorted(packets.items(), key=lambda item: source_order(item[1][1])):
            if raw not in current:
                continue
            title = review.get('object', raw)
            historical = Path(review.get('parent_run', '')).resolve() != (revision/'reduction').resolve()
            provenance = ('Retained historical candidate geometry; IDs do not refer to the current reduction. '
                          if historical else 'Candidate geometry from this reduction. ')
            original = self.root/'QC'/Path(raw).stem
            cards = []
            for suffix, label in [('science', '2D Spectrum'), ('traces', 'Sky-subtracted 2D spectrum'),
                                  ('profiles', 'Spatial Profile and Trace Identification'),
                                  ('acquisition-1', 'Acquisition 1'), ('acquisition-2', 'Acquisition 2')]:
                matches = sorted(folder.glob(f'FitQC_*_{suffix}.jpg'))
                if matches:
                    ratio = 5/3 if suffix in ('science','traces') else None
                    cards.append(self.figure(matches[0], title+' / '+label,
                                             full_width=suffix=='profiles', aspect_ratio=ratio))
            candidates = ', '.join(f'C{c.get("number", i+1)}: x={c.get("anchor_raw_x", c.get("raw_spatial_pixel", 0)):.2f}'
                                   for i, c in enumerate(review.get('candidates', [])))
            profile_info = self.read(folder/'profile-display.json')
            summed = profile_info.get('profile_statistic', '').startswith('sum ')
            profile_description = (
                'The spatial profiles sum background-subtracted electrons over the stated y rows at each fixed x. '
                'Background subtraction uses a saved fit; the zero level is not imposed on the curves. '
                if summed else
                'A fitted background is unavailable for these profiles; they show median electrons including sky across the stated y rows at each fixed x. ')
            if profile_info.get('background_source'):
                profile_description += 'The profiles use a retained two-sided PypeIt fit; the 2D panels above retain the parent reduction state. '
            if profile_info.get('profile_background_x_limits'):
                profile_description += 'Curves are shown only within the supported background interval; candidate markers can lie outside it. '
            current_layout=profile_info.get('layout')=='three wavelength panels; before subtraction only'
            if current_layout:
                profile_description=('Three equal air-wavelength intervals show summed data before subtraction, '
                    'with the saved sky model overlaid where available. No new background is fitted. '
                    'Candidate labels retain the saved review IDs and do not establish target identity. ')
            sections.append((source_order(review), f'<article class="exposure"><h3>{_text(title)}</h3>'
                f'<p class="caption">{_text(raw)} · Identity {_text(review.get("identity", "unverified").lower())}</p>'
                f'<p class="caption">{_text(provenance)}</p>'
                f'<p>{_text(candidates)} <span class="caption">(raw detector pixels)</span></p>'
                '<p class="caption">Raw detector: full frame and extraction window. Candidate traces: saved source paths. '
                'Sky residual: background-subtracted window, shown without and with those paths. '
                'Where no sky fit exists, the plots are labeled as counts including sky. '
                + profile_description
                +('' if current_layout else 'Q2 and Q3 are the second and third quarters of the extraction window along y; '
                'the legend gives inclusive raw-pixel ranges. Line-style markers locate candidates at the middle row.')+'</p>'
                f'<div class="qc-plots">{"".join(cards)}</div>'
                + ('<p class="caption">Candidate line styles: solid, dashed, dash-dot, dotted. '
                   'Candidate IDs remain fixed. '
                   + self.download(folder/'profile-display.json', 'Profile display details')+'</p>'
                   if (folder/'profile-display.json').is_file() else '')+'</article>'))
            # An explicit current review supersedes old candidate drawings.
            # Keep historical evidence accessible without presenting its rejected
            # markers as the current source geometry.
            error=next((e for e in run.get('errors', []) if Path(e.get('path','')).name==raw),None)
            if error:
                from .failure_profiles import load_failure_profile
                selected=load_failure_profile(revision/'reduction',error['path'])
                if selected and selected[1].get('selection'):
                    selected_folder,display=selected
                    current_profile=(f'<article class="exposure"><h3>{_text(title)}</h3>'
                        '<p class="caption">Current reviewed aperture geometry.</p>'
                        '<div class="qc-plots">'+self.figure(selected_folder/'profiles.png',
                            title+' / Reviewed spatial profile',full_width=True)+'</div>'
                        '<p class="caption">'+_text(display.get('caption',''))+'</p></article>')
                    order,older=sections[-1]
                    sections[-1]=(order,current_profile+'<details><summary>Earlier candidate diagnostics (superseded)</summary>'
                                  +older+'</details>')
            # Discover saved polynomial experiments; never invent a global C3 recipe.
            rows, background_cards, experiments = [], [], {}
            for path in sorted(set(original.glob('*/background.json')) | set(folder.glob('*/background.json'))):
                result, recipe = self.read(path), self.read(path.parent/'recipe.json')
                if not recipe or recipe.get('input_sha256', {}).get('raw') != review.get('input_sha256', {}).get('raw'):
                    continue
                key = json.dumps([recipe.get('objects'), recipe.get('sky')], sort_keys=True)
                # Named later revisions supersede an identical recipe; retain other choices.
                experiments[key] = (path, result, recipe)
            for path, result, recipe in experiments.values():
                obj = recipe['objects'][0];sky = recipe['sky']
                bands = ', '.join(f'[{a:g}, {b:g})' for a, b in sky.get('include', []))
                rows.append([_text(obj['label']), f'{obj["spat"]:.2f}', _text(sky.get('degree')),
                    _text(bands), f'{100*result.get("valid_spectrum_fraction", 0):.1f}%'])
                if sky.get('degree') == 2:
                    for image in path.parent.glob('FitQC_*_polynomial-background.jpg'):
                        background_cards.append(self.figure(image, title+' / Quadratic background'))
            if rows:
                for image in original.glob('FitQC_*background-order.jpg'):
                    background_cards.insert(0, self.figure(image, title+' / Polynomial order comparison'))
                backgrounds.append(f'<article class="exposure"><h3>{_text(title)}</h3>'
                    '<p class="note">Provisional source selection for this exposure only. Background-subtracted counts; no flux calibration or verified host separation.</p>'
                    + _table(['Candidate', 'Raw x', 'Degree', 'Background x bands', 'Valid spectrum'], rows)
                    + '<div class="qc-plots">'+''.join(background_cards)+'</div></article>')
        from .catalogue_report import catalogue_qc
        sections.extend(catalogue_qc(self, run))
        from .report_evidence import background_evidence
        return (self.trace_geometry_review(revision)+''.join(html for _, html in sorted(sections, key=lambda item:item[0])) or '<p>No source-review packets for this revision.</p>',
                background_evidence(run)+''.join(backgrounds)+self.local_background_checks(revision,run))

    def trace_geometry_review(self, revision):
        """Render an explicitly selected, hash-bound diagnostic revision."""
        pointer=self.read(self.root/'QC/trace-geometry-current.json')
        if not pointer or not pointer.get('enabled', False):return ''
        manifest=revision/'reduction/run.json'
        if pointer.get('parent_manifest_sha256')!=sha256(manifest):return ''
        directory=(self.root/pointer['directory']).resolve()
        directory.relative_to(self.root.resolve())
        audit=directory/'audit.json'
        if not audit.is_file() or sha256(audit)!=pointer.get('audit_sha256'):
            self.warnings.append('Trace geometry review changed; not displayed')
            return ''
        sections=[]
        for target in self.read(audit):
            name=target['object'];folder=directory/name;rows=[];cards=[]
            for candidate in target['candidates']:
                check=candidate['calibrated_pixel_check']
                rows.append([f'C{candidate["candidate"]}',f'{candidate["raw_x"]:.2f}',
                    f'{check["positive_bins"]} / {len(check["bins"])}',
                    'Persistent local peak; identity unverified' if check['accepted'] else 'Not supported as a persistent trace'])
            for filename,digest in target.get('figures',{}).items():
                image=folder/filename
                if image.is_file() and sha256(image)==digest:
                    cards.append(self.figure(image,name+' / '+filename.removesuffix('.jpg').split('_')[-1]))
                else:self.warnings.append(name+': trace review figure changed')
            sections.append('<article class="exposure"><h3>'+_text(name)+' / Trace geometry review</h3>'
                '<p class="note">Historical C-labels are retained here. Their paths were checked against '
                'the current bias/flat-corrected pixels in six independent spectral sections. '
                'At least three sections need local contrast S/N ≥3. A passing check is not target identification. '
                'The historical trial FWHM of 4 pixels is used for this comparison.</p>'
                +_table(['Historical candidate','Raw x','Positive sections','Evidence'],rows)
                +'<p class="caption">The saved historical curves are shown unchanged: solid C1, dashed C2, '
                'dash-dot C3, with thicker lines. These diagnostics do not establish on-slit centring.</p>'
                +'<div class="qc-plots">'+''.join(cards)+'</div></article>')
            error=next((e for e in self.read(manifest).get('errors',[]) if e.get('object')==name),None)
            if error:
                from .failure_profiles import load_failure_profile
                selected=load_failure_profile(revision/'reduction',error['path'])
                if selected and selected[1].get('selection'):
                    sections[-1]='<details><summary>'+_text(name)+' / Earlier trace geometry checks</summary>'+sections[-1]+'</details>'
        return ''.join(sections)

    def local_background_checks(self, revision, run):
        """Keep geometry-only checks separate from extraction/identity claims."""
        manifest=revision/'reduction/run.json'
        if not manifest.is_file():return ''
        digest=sha256(manifest)
        raw_hashes={Path(r.get('raw',r.get('path',''))).name:
                    r.get('raw_sha256',r.get('diagnostics',{}).get('raw_sha256'))
                    for r in run.get('results',[])+run.get('errors',[])}
        sections=[];seen=set()
        for path in sorted((self.root/'QC').glob('local-background-*/local-background.json'),reverse=True):
            record=self.read(path)
            if record.get('night')!=self.root.name or record.get('parent_manifest_sha256')!=digest:continue
            for target in record.get('targets',[]):
                expected_raw=raw_hashes.get(target.get('raw_name'))
                if not expected_raw or target.get('input_sha256',{}).get('raw')!=expected_raw:continue
                if target['raw_name'] in seen:continue
                seen.add(target['raw_name'])
                name=target['object'];directory=path.parent/name
                directory.resolve().relative_to(self.root)
                rows=[]
                for case in target.get('cases',[]):
                    bands=', '.join(f'[{a}, {b})' for a,b in case.get('midrow_bands',[])) or 'Unavailable'
                    status=(('Fits expanded window' if case.get('expanded_window') else 'Fits saved window') if case['status']=='geometry_pass' else 'Insufficient two-sided space')
                    rows.append([f'C{case["candidate"]}',f'{case["raw_x"]:.2f}',case.get('selection',{}).get('actual_width_spatial_pixels',case['width_pixels']),
                                 _text(bands),status])
                figures=[];extractions=[];downloads=[]
                for case in target.get('cases',[]):
                    for trial in case.get('extractions',[]):
                        if trial['status']=='extracted_unverified':
                            result=trial['results'][0];fit=result['processing']['background_fit']
                            extractions.append([f'C{case["candidate"]}',trial['degree'],'Extracted · identity unverified',
                                f'{100*result["valid_flux_fraction"]:.1f}%',str(fit['minimum_side_samples'])])
                            recipe_path=path.parent/trial['recipe'];trial_dir=recipe_path.with_suffix('')
                            sidecar=self.read(trial_dir/'qa/local-background-figure.json')
                            saved=list((trial_dir/'native').glob('*_images.npz'))
                            if len(saved)==1 and sidecar.get('images_sha256')==sha256(saved[0]):
                                for filename,digest_image in sidecar.get('figures',{}).items():
                                    image=trial_dir/'qa'/Path(filename).name
                                    if image.is_file() and sha256(image)==digest_image:
                                        figures.append(self.figure(image,f'{name} / C{case["candidate"]}, degree {trial["degree"]}'))
                            product=trial_dir/'products'/Path(result['product']).name
                            downloads.append(self.download(product,f'C{case["candidate"]} degree {trial["degree"]} FITS'))
                        else:
                            extractions.append([f'C{case["candidate"]}',trial['degree'],_text(trial.get('error','Failed')),'—','—'])
                for filename,expected in target.get('figures',{}).items():
                    image=directory/Path(filename).name
                    if not image.is_file() or sha256(image)!=expected:
                        self.warnings.append(f'{name}: local-background figure is missing or changed')
                        continue
                    figures.append(self.figure(image,name+' / Local background geometry'))
                note=('Geometry checks only; test widths are spatial pixels. No new extraction or default change.'
                      if record.get('mode')=='geometry_preflight' else
                      'PypeIt extraction trials; source identities remain unverified. See trial record for failures.')
                sections.append(f'<article class="exposure"><h3>{_text(name)} · Bounded background</h3>'
                    f'<p class="note">{note} Neighbour exclusions use the saved candidates and requested source guard.</p>'
                    +_table(['Candidate','Raw x','Pixels per side','Bands at middle row','Geometry'],rows)
                    +(_table(['Candidate','Degree','Extraction','Valid flux bins','Minimum accepted / side'],extractions) if extractions else '')
                    +'<div class="qc-plots">'+''.join(figures)+'</div>'
                    +'<p>'+' · '.join(filter(None,downloads+[self.download(path,'Trial record')]))+'</p></article>')
        return ''.join(sections)

    def calibrations(self, revision=None):
        rows, seen = [], set()
        paths=sorted((self.root/'calibrations').glob('*/run.json'))
        selected=self.read(revision/'night.json').get('native_calibrations') if revision else None
        if selected:
            selected_path=Path(selected).resolve()/'run.json'
            paths=[selected_path]+[p for p in paths if p.resolve()!=selected_path]
        for path in paths:
            for group in self.read(path).get('groups', []):
                if not group.get('id') or not group.get('setup'):
                    continue
                if group['id'] in seen:
                    continue
                seen.add(group['id']);s = group.get('setup', {})
                wave = group.get('wave_qa', [])
                rms = ', '.join(f'{w["held_out_rms_pixel"]:.3f}' for w in wave if 'held_out_rms_pixel' in w)
                status=group.get('status','unknown').replace('_',' ')
                if any(w.get('independent_science_validation') is False for w in wave):
                    status += '; science unvalidated'
                label=_text(status)
                if _unvalidated(status):label='<span class="qc-status-bad">'+label+'</span>'
                rows.append([_text(s.get('grism', '')), _text(s.get('filter', '')), _text(s.get('slit', '')),
                             label, rms or '—'])
        validations = sorted((self.root/'archive/native-validation').glob('*/validation.json'))
        note = ''
        if validations:
            v = self.read(validations[-1]);counts = Counter(x.get('status', 'unknown') for x in v.get('results', []))
            note = '<p class="caption">Latest retained native diagnostic: '+_text(', '.join(f'{n} {k.replace("_", " ")}' for k, n in counts.items()))
            note += f'; {len(v.get("errors", []))} failed checks. These retained diagnostics are separate from this revision’s calibration provenance.</p>'
        historical = (_table(['Grism', 'Filter', 'Slit', 'Recorded status', 'Held-out RMS (px)'], rows)+note
                      if rows else '<p>No retained native calibration groups.</p>')
        if not revision:
            return historical
        run = self.read(revision/'reduction/run.json')
        if run.get('flux_calibration') != 'same-night':
            return historical
        from .report_evidence import calibration_overview, response_status
        responses = []
        for entry in run.get('standards', []):
            frame, fit = entry['frame'], entry['fit']
            setup = frame['setup']; response = fit.get('response')
            download = 'No fitted response'
            if response:
                path = Path(response['path'])
                if sha256(path) != response['sha256']:
                    raise ValueError('daily response changed since reduction')
                download = self.download(path, 'Response FITS')
            responses.append([_text(frame['metadata']['OBJECT']), _text(frame['metadata']['DATE-OBS']),
                _text(setup['grism']+' / '+setup['filter']+' / '+setup['slit']),
                response_status(fit), download])
        associations = [[_text(r['object']), _text(r['sensitivity'].get('standard_name','')),
                         _text(Path(r['sensitivity']['standard_raw']).name)] for r in run.get('results', [])]
        current = calibration_overview(run)
        from .telluric_products import selected_pairs,status_text
        pairs=selected_pairs(revision/'reduction')
        if pairs:
            current+='<h3>Paired telluric reductions</h3><p>The response scope below describes the original reduction. '
            current+='The Spectra section also retains separate atmospheric corrections with both downloads and overlaid curves.</p>'
            current+=_table(['Transient','Corrected bands','Transfer status'],[
                [_text(r['object']),_text(', '.join(pairs[Path(r['product']).name]['bands'])),
                 '<span class="qc-status-bad">'+_text(status_text(pairs[Path(r['product']).name]))+'</span>']
                for r in sorted(run['results'],key=source_order) if Path(r['product']).name in pairs])
            current+='<p>'+self.download(Path(next(iter(pairs.values()))['manifest']),'Telluric fits and transfer checks')+'</p>'
        current += (_table(['Standard','UTC start','Setup','Status','Download'],responses)
                    if responses else '<p>No contemporary response was built; see failed exposures.</p>')
        current += _table(['Spectrum','Reference standard','Observed standard exposure'], associations)
        current += ('<p>'+str(len(run.get('held_out',[])))+' held-out exposure transfers recorded; '
                    +str(len(run.get('held_out_errors',[])))+' transfer checks failed. '
                    'These test internal repeatability, not an independent absolute scale.</p>')
        transfers=[]
        for check in run.get('held_out',[]):
            setup=check['setup'];ratio=check.get('median_ratio')
            bands='; '.join(f'{b["band_air_angstrom"][0]}–{b["band_air_angstrom"][1]} Å: {b["median_ratio"]:.3f}'
                           for b in check.get('continuum_bands',[]) if b.get('median_ratio') is not None)
            transfers.append([_text(setup['grism']+' / '+setup['filter']+' / '+setup['slit']),
                _text(Path(check['training_raw']).name),_text(Path(check['raw']).name),
                'Unavailable' if ratio is None else f'{ratio:.3f}',_text(bands)])
        if transfers:
            current+='<h3>Independent-exposure response checks</h3><p>Recovered flux divided by the physical standard reference; unity is closure. Slit losses and transparency can change these ratios.</p>'
            current+=_table(['Setup','Response exposure','Test exposure','Median ratio','Continuum bands'],transfers)
        return current+historical

    def archive(self, revision):
        paths = list((self.root/'archive').glob('*/run.json'))
        paths += list((self.root/'archive/reductions').glob('*/run.json'))
        paths += [p for p in (self.root/'nightly').glob('*/reduction/run.json') if p.parent.parent != revision]
        rows = []
        for path in sorted(paths):
            run = self.read(path)
            if not any(k in run for k in ['results', 'errors', 'skipped']):
                continue
            rows.append([_text(path.parent.relative_to(self.root)), len(run.get('results', [])),
                         len(run.get('errors', [])), len(run.get('skipped', []))])
        return _table(['Retained reduction', 'Spectra', 'Failed', 'Unsupported'], rows) if rows else '<p>No earlier reductions.</p>'

    def archive_folders_and_logs(self, revision):
        current = revision.relative_to(self.root)
        folders = [
            (current/'reduction/products', 'Current reduced spectra (FITS and ASCII)'),
            (current/'QC', 'Current source and background inspection'),
            (current/'logs', 'Current workflow, reduction and report logs'),
            (Path('nightly'), 'Current and previous processing revisions'),
            (Path('QC'), 'Retained source and background experiments'),
            (Path('calibrations'), 'Native masters and group selections'),
            (Path('archive'), 'Earlier reductions, diagnostics and logs'),
        ]
        if (self.root/'filtered').is_dir():
            folders.insert(-1, (Path('filtered'), 'Separate filtered-standard reductions and reports'))
        content = '<h3>Folders</h3><p>Night folder: <code>'+_text(str(self.root)+'/')+'</code></p>'
        content += _table(['Folder', 'Contents'], [
            ['<code>'+_text(str(path)+'/')+'</code>', label] for path, label in folders])
        downloads = [self.download(p, p.name) for p in sorted((revision/'logs').glob('*.log'))]
        content += '<details><summary>Download logs</summary>'
        content += ('<ul>'+''.join('<li>'+d+'</li>' for d in downloads if d)+'</ul>'
                    if any(downloads) else '<p>No saved logs.</p>')
        return content+'</details>'


SCRIPT = """
document.querySelectorAll('.review-launch').forEach(button => button.addEventListener('click', () => {
  const help = button.closest('.exposure').querySelector('.review-help');
  help.open = true;
  help.querySelector('.launch-status').textContent = 'Desktop launch requested. If no window appears, use the command below.';
  // Keep the log in place even when the embedding viewer rejects custom protocols.
  let bridge = document.getElementById('desktop-launch-bridge');
  if (!bridge) {
    bridge = document.createElement('iframe'); bridge.id = 'desktop-launch-bridge';
    bridge.hidden = true; bridge.title = 'Desktop application launch'; document.body.appendChild(bridge);
  }
  bridge.src = button.dataset.reviewUrl;
}));
const dialog = document.querySelector('dialog');
document.querySelectorAll('.enlarge').forEach(button => button.addEventListener('click', () => {
  const source = button.querySelector('img');
  dialog.querySelector('img').src = source.src;
  dialog.querySelector('img').alt = source.alt;
  dialog.showModal();
}));
dialog.querySelector('button').addEventListener('click', () => dialog.close());
dialog.addEventListener('click', event => { if(event.target === dialog) dialog.close(); });
"""


EPOCH_CSS = """
    main{max-width:1480px}.tabs{display:flex;flex-wrap:wrap;gap:8px;margin:20px 0}.tabs a{padding:7px 12px;border:1px solid var(--line);border-radius:5px;text-decoration:none;background:var(--panel)}
    .target>h2{margin-bottom:18px}.exposure{margin:24px 0 32px}.exposure~.exposure{border-top:1px solid #bcc5b8;padding-top:28px}
    .review-launch{font:inherit;color:var(--forest);background:var(--panel);border:1px solid var(--line);border-radius:5px;padding:7px 12px;cursor:pointer}.review-help pre{white-space:pre-wrap;overflow-wrap:anywhere}.qc-plots{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}
    #archive code{overflow-wrap:anywhere}
    .qc-plots>.full-width{grid-column:1 / -1}
    .enlarge{display:block;padding:0;border:0;background:transparent;width:100%;cursor:zoom-in}.enlarge:focus-visible{outline:3px solid var(--forest)}
    .enlarge.fixed-ratio{background:white}.enlarge.fixed-ratio img{width:100%;height:100%;object-fit:contain}
    .plots .plot-image{aspect-ratio:auto}td{white-space:normal}dialog{max-width:95vw;max-height:95vh;padding:10px;border:1px solid var(--line);border-radius:8px}dialog img{display:block;max-width:100%;height:auto}dialog button{position:sticky;top:0;padding:10px 18px;cursor:pointer}dialog::backdrop{background:#172c43bb}
    @media(max-width:800px){.plots{grid-template-columns:1fr;overflow:visible}.qc-plots{grid-template-columns:1fr}.plot{min-width:0}h1{font-size:30px}.intro{display:block}}
    """


def build_epoch_page(night, project_root='.', *, revision_dir=None):
    """Rebuild presentation only; saved science products and recipes stay unchanged."""
    if date.fromisoformat(night).isoformat() != night:
        raise ValueError('observing date must be YYYY-MM-DD')
    root = Path(project_root).resolve()/'runs'/night
    page = EpochPage(root)
    if revision_dir is None:
        pointer = page.read(root/'nightly/current.json')
        revision = (root/'nightly'/pointer.get('revision', '')).resolve()
        if pointer.get('night') != night or revision.parent != root/'nightly':
            raise ValueError('invalid dated revision pointer')
        destination = root/'index.html'
    else:
        revision=Path(revision_dir).resolve()
        if not any(revision.is_relative_to(root/sub) for sub in ('nightly','QC')):
            raise ValueError('review revision must belong to this observing night')
        destination=revision/'index.html'
    state = page.read(revision/'night.json')
    if state.get('night') != night:
        raise ValueError('workflow belongs to a different observing date')
    run = page.read(revision/'reduction/run.json')
    rescue=None;rescue_path=None
    pointer_path=root/'rescue/current.json'
    if pointer_path.exists():
        pointer=page.read(pointer_path);rescue_path=(root/pointer['path']).resolve()
        if not rescue_path.is_relative_to(root/'rescue') or pointer.get('night')!=night:
            raise ValueError('invalid provisional rescue pointer')
        if sha256(rescue_path/'run.json')!=pointer['run_sha256']:
            raise ValueError('provisional rescue manifest changed')
        candidate=page.read(rescue_path/'run.json')
        parent_manifest=revision/'reduction/run.json'
        if parent_manifest.is_file() and candidate.get('parent_sha256')==sha256(parent_manifest):
            if candidate.get('wavelength_policy')!='provisional_arc_only' or candidate.get('night')!=night:
                raise ValueError('invalid provisional rescue policy')
            failures={e['path'] for e in run['errors'] if e.get('code')=='wavelength_validation_failed'}
            for result in candidate['results']:
                if result.get('status')!='PROVISIONAL' or result['raw'] not in failures or result.get('wavelength_validated') is not False:
                    raise ValueError('rescue result cannot replace an accepted product')
                for path,digest in result['product_sha256'].items():
                    if sha256(Path(path))!=digest:raise ValueError('provisional product changed')
            rescue=candidate
    from .typing_report import records_for_run, HOST_CSS
    page.typing_run = revision/'reduction'
    page.typing_annotations = records_for_run(page.typing_run) if (page.typing_run/'typing-metadata.json').exists() else {}
    from .redshift_report import load_redshifts,redshift_section
    from .manual_recovery_report import selected_manual_recoveries
    manuals=selected_manual_recoveries(revision/'reduction')
    redshifts=load_redshifts(page,revision/'reduction',rescue_path if rescue else None,manuals=[p for p,_ in manuals])
    page.redshift_records=redshifts
    qc, background = page.qc(revision, run)
    if rescue:
        audit=page.read(rescue_path/'verification.json')
        if audit and (audit.get('run_sha256')!=sha256(rescue_path/'run.json') or audit.get('status')!='passed'):
            raise ValueError('current rescue needs a matching product audit')
        geometry=[]
        for item in sorted(audit.get('spectra',[]),key=source_order):
            sides=item.get('background_gap_from_aperture_edge_pixels')
            if not sides:continue
            gaps=[f"{s['median_gap']:.1f} (range {s['minimum_gap']:.1f}–{s['maximum_gap']:.1f})" for s in sides]
            geometry.append([_text(item['object']),*gaps,item['background_minimum_side_samples']])
        if geometry:
            retained_background=background
            background='<h3>Current provisional backgrounds</h3><p>Nearby finite bands follow the source trace, with guards for source wings and trace movement. '
            background+='Distances below run from the final fitted aperture edge to the nearest sampled pixel centre, in detector pixels. '
            background+='The local background can include extended galaxy light; this is not validated host deblending. All four fits retain two-sided support without aperture overlap.</p>'
            background+=_table(['Target','Left gap (pixels)','Right gap (pixels)','Minimum retained pixels per side'],geometry)
            background+='<p>'+page.download(rescue_path/'verification.json','Provisional product and background audit')+'</p>'
            background+='<h3>Other retained background diagnostics</h3>'+retained_background
    # Redshift evidence is its own tab: it has its own inputs, figures and
    # provenance, and reads as a separate measurement from source QC.
    redshift_content=redshift_section(page,redshifts)
    titles = [('observations', 'Observing log'), ('spectra', 'Spectra'), ('issues', 'Issues'),
              ('qc', 'Source QC'), ('background', 'Background'), ('calibrations', 'Calibrations'),
              ('archive', 'Reduction Archive')]
    if redshift_content:
        titles.insert(4, ('redshifts', 'Redshift evidence'))
    calibration_content=page.calibrations(revision)
    native=revision/'native-validation'
    if (native/'validation.json').exists():
        calibration_content+='<h3>Contemporary raw-calibration checks</h3><p>These checks use the selected native masters. They do not change the archive-calibrated spectra above.</p>'
        calibration_content+=page.download(native/'validation.json','Native validation record')
        calibration_content+=page.wavelength_checks(native)
    spectra=page.spectra(revision/'reduction',run)
    for j,(manual_path,manual_run) in enumerate(manuals):
        spectra+='<h3>Checked manual recovery — retained for review</h3><p>Native PypeIt extraction with explicit two-sided background and passing wavelength checks. Original attempts remain listed separately; absolute flux and host subtraction remain unvalidated.</p>'
        if manual_run.get('validation_policy')=='two-anchor':
            spectra+='<p>Wavelength acceptance: 5577 and 6300 Å quality and agreement checks. The 6363 Å line is diagnostic only; its measured discrepancy is retained in the sky-check download.</p>'
        spectra+=page.spectra(manual_path,manual_run,id_prefix=f'manual-{j}-spectrum')
    if rescue:
        spectra+='<h3>Provisional rescue spectra</h3><p class="qc-status-bad">Same-night flux calibration; science wavelengths remain unvalidated. No flexure correction applied. These do not replace accepted products.</p>'
        spectra+=page.spectra(rescue_path,rescue,id_prefix='rescue-spectrum')
    issues=page.problems(run,revision/'reduction')
    if manuals:
        issues='<p>Failures from the original nightly attempts are retained here. Later checked manual recoveries are listed separately in Spectra and the observing log.</p>'+issues
    content = [page.observations(revision, run, rescue, redshifts, manuals), spectra, issues, qc,
               background, calibration_content, page.archive(revision)+page.archive_folders_and_logs(revision)]
    if redshift_content:
        content.insert(4, redshift_content)
    from .iraf_report import comparison_section, write_comparison_page
    iraf_page=EpochPage(root if revision_dir is None else revision)
    iraf_content=comparison_section(iraf_page,revision)
    page.warnings.extend(iraf_page.warnings)
    if iraf_content:
        write_comparison_page(iraf_page,revision,iraf_content)
        iraf_content='<p><a href="iraf-comparison.html">Open separate IRAF comparison report</a>. The nightly spectra above show PypeIt only.</p>'
        titles.insert(2,('iraf-comparison','IRAF comparison'))
        content.insert(2,iraf_content)
    sections = ''.join(f'<section class="target" id="{anchor}"><h2>{title}</h2>{body}</section>'
                       for (anchor, title), body in zip(titles, content))
    stats = ' · '.join(f'{len(run.get(k, []))} {label}' for k, label in [('results','spectra'), ('errors','failed'), ('skipped','unsupported')])
    if run.get('review_counts'):stats+=' · '+str(len(run['review_counts']))+' review counts'
    if manuals:stats+=' · '+str(sum(len(r['results']) for _,r in manuals))+' checked manual recovery spectra'
    if rescue:stats+=' · '+str(len(rescue['results']))+' provisional flux spectra'
    warning = '<p class="note">'+_text('; '.join(page.warnings))+'</p>' if page.warnings else ''
    scope = ('<p class="caption">Selected targets: '+_text(', '.join(state['targets']))
             +'. Other exposures are listed for context; earlier results remain in Reduction Archive.</p>') if state.get('targets') else ''
    if revision_dir is not None:
        scope += '<p class="caption">Review of processing revision '+_text(revision.name)+'. The selected nightly reduction is unchanged.</p>'
    html = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Nightly LOG - NTT/EFOSC2 (ePESSTO+) · {night} · EFOSpec</title><style>{CSS}{EPOCH_CSS}{HOST_CSS}</style></head><body><main>
<header class="masthead"><span class="brand">Nightly LOG - NTT/EFOSC2 (ePESSTO+)</span><span class="eyebrow">EFOSC2 · {night}</span></header>
<div class="intro"><div><span class="eyebrow">Observing night</span><h1>{night}</h1></div><p>{stats}</p></div>
{scope}
<nav class="tabs" aria-label="Report sections">{''.join(f'<a href="#{a}">{t}</a>' for a,t in titles)}</nav>
{warning}{sections}<footer>EFOSpec · {night} · Saved figures and measurements</footer></main>
<dialog aria-label="Enlarged figure"><button type="button">Close</button><img alt=""></dialog><script>{SCRIPT}</script></body></html>'''
    from .report_assets import write_report_bundle
    bundle=write_report_bundle(destination,html)
    return dict(bundle,revision=revision.name,embedded_figures=page.images,warnings=page.warnings)


def build_filtered_page(run_dir):
    """A portable supplement using the nightly report's style and QC renderer."""
    root = Path(run_dir).resolve(); page = EpochPage(root)
    run = page.read(root/'run.json')
    if run.get('method') not in {'native_filtered_standard_extraction','native_red_counts_diagnostic',
                               'native_optical_counts_diagnostic','native_nightly_counts_diagnostic'} or run.get('status') == 'running':
        raise ValueError('counts report requires a completed separate native run')
    red=run['method']=='native_red_counts_diagnostic'
    purpose=('Nightly calibration review' if run['method']=='native_nightly_counts_diagnostic' else 'Red spectra' if red else 'Wavelength review counts' if
             run['method']=='native_optical_counts_diagnostic' else 'Filtered standards')
    roles={f['path']:f['role'] for f in page.read(root/'inputs.json').get('frames',[])}
    night = date.fromisoformat(run['night']).isoformat()
    rows, spectra, qc, backgrounds, calibrations = [], [], [], [], []
    for result in sorted(run['results'], key=lambda row: (row['object'].casefold(), row['date_obs'])):
        directory = Path(result['product']).parent.parent
        directory.resolve().relative_to(root)
        for filename, digest in result['product_sha256'].items():
            path = Path(filename); path.resolve().relative_to(root)
            if sha256(path) != digest: raise ValueError('filtered report product changed')
        plots = directory/'plots'
        drawing = page.read(plots/'plots.json')
        if drawing.get('product_sha256') != result['product_sha256']:
            raise ValueError('filtered figures do not match the extracted products')
        for name, digest in drawing['figures_sha256'].items():
            path = plots/name; path.resolve().relative_to(plots.resolve())
            if sha256(path) != digest: raise ValueError('filtered figure changed')
        title = _text(result['object']+' / '+result['date_obs'])
        flux = result['flux_response']
        if 'response' in flux and sha256(Path(flux['response']['path'])) != flux['response']['sha256']:
            raise ValueError('filtered response changed before reporting')
        rows.append([_text(result['object']), _text(result['date_obs']),
                     'Science' if roles.get(result['raw'])=='science' else 'Standard',
                     _text(result['setup']['grism']+' / '+result['setup']['filter']), _text(result['setup']['slit']),
                     _text(result.get('delivery_status','Counts extracted; flux unvalidated'))])
        links = [page.download(Path(result[key]), label) for key, label in
                 [('product','Counts FITS'), ('ascii','Counts ASCII'), ('native_spec1d','PypeIt spec1d')]]
        links += [page.download(p, p.suffix[1:].upper()+' figure') for p in sorted(plots.glob('FitQC_*_counts.*'))
                  if p.suffix in {'.pdf', '.jpg'}]
        spectra.append('<article class="exposure"><h3>'+title+'</h3><p>'+' · '.join(links)+'</p>'
            +page.figure(plots/'counts.png', result['object']+' / Extracted counts')+'</article>')
        pending=result.get('pending_checks',['Target identity unverified; exposure-specific flexure unmeasured.'])
        qc.append('<article class="exposure"><h3>'+title+'</h3>'
            '<p class="qc-status-bad">'+_text(' '.join(pending))+'</p>'
            '<div class="qc-plots">'+page.figure(plots/'science.png',result['object']+' / 2D Spectrum')
            +page.figure(plots/'traces.png',result['object']+' / Sky-subtracted 2D spectrum')
            +page.figure(plots/'profiles.png',result['object']+' / Spatial Profile and Trace Identification',full_width=True)
            +'</div></article>')
        bg = result['background']
        backgrounds.append([title, str(bg['spatial_degree']), str(bg['minimum_side_samples']),
                            f'{bg["reduced_chi2"]:.3f}', page.download(directory/'extraction.json','Fit and provenance')])
        wave = result['wavelength']
        detail = (page.download(Path(flux['response']['path']), 'Experimental response FITS')
                  if 'response' in flux else _text(flux.get('error', 'No physical flux reference')))
        calibrations.append('<article class="exposure"><h3>'+title+'</h3><p>'
            +f'Arc RMS {wave["rms_pixel"]:.3f} pixel; held-out RMS {wave["held_out_rms_pixel"]:.3f} pixel. '
            +'Wavelength support: '+'–'.join(f'{x:.1f}' for x in wave['valid_wave_vacuum_angstrom'])
            +' vacuum Å.</p><p class="qc-status-bad">'+_text(flux['status'].replace('_',' '))
            +'; no independent flux validation.</p><p>'+detail+'</p></article>')
    issues = '<ul>'+''.join('<li>'+_text(s)+'</li>' for s in run['limitations'])+'</ul>'
    issues += ''.join('<p class="qc-status-bad">'+_text(e['object']+': '+e['error'])+'</p>' for e in run['errors'])
    titles = [('observations','Observing log'),('spectra','Spectra'),('issues','Issues'),('qc','Source QC'),
              ('background','Background'),('calibrations','Calibrations'),('archive','Archive'),('records','Run record')]
    bodies = [_table(['Source','UTC','Role','Setup','Slit','Outcome'],rows), ''.join(spectra), issues,
              ''.join(qc), _table(['Exposure','Spatial degree','Minimum samples per side','Reduced chi2','Details'],backgrounds),
              ''.join(calibrations), '<p>Earlier reductions and the main nightly report remain separate and unchanged.</p>',
              '<p>Native counts; no flux calibration. Each exposure records its wavelength checks.</p><p>'
              +page.download(root/'run.json','Run manifest')+' · '+page.download(root/'setup-checks.json','Setup checks')+'</p>']
    sections = ''.join(f'<section class="target" id="{a}"><h2>{t}</h2>{b}</section>'
                       for (a,t),b in zip(titles,bodies))
    html = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Nightly LOG - NTT/EFOSC2 (ePESSTO+) · {night} · {purpose}</title><style>{CSS}{EPOCH_CSS}</style></head><body><main>
<header class="masthead"><span class="brand">Nightly LOG - NTT/EFOSC2 (ePESSTO+)</span><span class="eyebrow">EFOSC2 · {night}</span></header>
<div class="intro"><div><span class="eyebrow">{purpose}</span><h1>{night}</h1></div><p>{len(run['results'])} extracted · {len(run['errors'])} failed</p></div>
<p class="qc-status-bad">Experimental counts. Review the remaining calibration checks for each exposure.</p>
<nav class="tabs" aria-label="Report sections">{''.join(f'<a href="#{a}">{t}</a>' for a,t in titles)}</nav>
{sections}<footer>EFOSpec · {night} · {purpose}</footer></main>
<dialog aria-label="Enlarged figure"><button type="button">Close</button><img alt=""></dialog><script>{SCRIPT}</script></body></html>'''
    destination = root/'index.html'; destination.write_text(html)
    return destination


def main(argv=None):
    parser = argparse.ArgumentParser(description='Rebuild one portable epoch report from saved results')
    parser.add_argument('night')
    parser.add_argument('--project-root', default=str(Path(__file__).resolve().parents[2]))
    parser.add_argument('--revision', help='Render this saved dated revision beside its own files without changing the selected nightly log')
    args = parser.parse_args(argv)
    print(json.dumps(build_epoch_page(args.night, args.project_root, revision_dir=args.revision), indent=2))


if __name__ == '__main__':
    main()
