"""Shared PypeIt extraction for automatic reductions and saved manual recipes."""
from contextlib import contextmanager

import numpy as np
from astropy import units as u


@contextmanager
def native_qa_context():
    """Allow literal slit-number hashes in PypeIt 3 QA while retaining TeX.

    Upstream renders the title before saving, even with show=False. Only the
    engine call receives this TeX catcode; restore the caller's rcParams.
    """
    import matplotlib
    settings={}
    if matplotlib.rcParams['text.usetex']:
        preamble=matplotlib.rcParams['text.latex.preamble']
        settings['text.latex.preamble']=preamble+'\n'+r'\AtBeginDocument{\catcode35=12}'
    with matplotlib.rc_context(settings):
        yield



class SourceSelectionError(ValueError):
    """Keep detected positions available when target identity is unresolved."""
    code = 'source_selection_required'
    action = 'Verify the target against acquisition data and supply an explicit extraction recipe.'

    def __init__(self, message, found, target_position):
        super().__init__(message)
        self.diagnostics = {'expected_spatial_pixel': target_position,
            'coordinate_system': 'zero_based_extraction_window',
            'candidates': [{'name': str(s.NAME), 'spatial_pixel': float(s.SPAT_PIXPOS)} for s in found],
            'review_recommendation': {'workflow':'aperture_and_background_review',
                'reason':'Multiple or unassociated traces may include field sources or host structure.',
                'suggested_background':'two_sided_polynomial_preview',
                'host_contamination_confirmed':False,'automatic_subtraction':False}}


def select_sources(found, target_position=None, manual=None):
    """Map stable recipe IDs to found objects without silently swapping sources."""
    if len(found) == 0:
        raise SourceSelectionError('no source detected in target window', found, target_position)
    if manual is None:
        if target_position is None or not np.isfinite(target_position):
            raise ValueError('automatic selection requires a finite target position')
        # Observer policy: strongest positive source peak within the existing
        # 15-pixel central search window. This selects geometry, not identity.
        candidates = [s for s in found if abs(s.SPAT_PIXPOS-target_position) <= 15
                      and getattr(s, 'sign', 1) != -1]
        ranked = []
        for source in candidates:
            peak = getattr(source, 'smash_peakflux', None)
            snr = getattr(source, 'smash_snr', None)
            if peak is not None and (not np.isfinite(peak) or peak <= 0):
                continue
            if snr is not None and (not np.isfinite(snr) or snr < 5):
                continue
            ranked.append((source, peak))
        if len(ranked) == 1:
            return {'primary': str(ranked[0][0].NAME)}
        if ranked and all(peak is not None for _, peak in ranked):
            ranked.sort(key=lambda item: item[1], reverse=True)
            if not np.isclose(ranked[0][1], ranked[1][1], rtol=1e-8, atol=0):
                return {'primary': str(ranked[0][0].NAME)}
        raise SourceSelectionError(f'no distinct positive source peak near expected position {target_position}',
                                   found, target_position)
    selected = {}
    for obj in manual:
        matches = [s for s in found if s.hand_extract_flag
                   and abs(s.hand_extract_spat-obj['spat']) < 1e-6
                   and abs(s.hand_extract_spec-obj['spec']) < 1e-6]
        if len(matches) != 1:
            raise ValueError(f'manual source {obj["id"]} was not found uniquely')
        selected[obj['id']] = str(matches[0].NAME)
    if len(set(selected.values())) != len(manual):
        raise ValueError('manual sources do not have distinct extraction identities')
    return selected


def validated_candidate_packet(selection, found, evidence):
    """Publish supported/manual candidates; retain all proposals only as audit data.

    Never recompute a saved sky mask using the filtered catalogue.
    """
    import json
    keep = np.array([bool(s.hand_extract_flag) or evidence.get(str(s.NAME), {}).get('accepted', False)
                     for s in found], dtype=bool)
    rows = selection['sky_input_mask'].shape[0]
    traces = np.array([s.TRACE_SPAT for s in found], dtype=float).reshape(len(found), rows)
    names = np.array([str(s.NAME) for s in found], dtype='U80')
    return dict(selection, found_names=names[keep], found_traces=traces[keep],
                rejected_names=names[~keep], rejected_traces=traces[~keep],
                candidate_evidence_json=np.array(json.dumps(evidence))), keep


@native_qa_context()
def extract_spectra(data, waveair, detector, *, exposure_time, read_noise,
                    output, basename, target_position=None, manual=None, sky_input_mask=None, processed=None,
                    sky_model=None, sky_trace=None, bounded_background=False):
    """Jointly model neighbouring objects and return requested spectra and images.

    Data and read noise are electrons. Coordinates, manual FWHM and boxcar
    radii are pixels in this image, after binning. Air wavelengths are Angstroms.
    A supplied sky mask controls the global fit, held during optimal extraction.
    PypeIt diagnostics are written below output; callers own product delivery.
    """
    import pypeit
    from pypeit.images.pypeitimage import PypeItImage
    from pypeit.slittrace import SlitTraceSet
    from pypeit.find_objects import FindObjects
    from pypeit.extraction import Extract
    from pypeit.manual_extract import ManualExtractionObj
    from pypeit.core.wave import airtovac
    from .pypeit_adapter import EFOSpecEFOSC2Spectrograph
    from .sky_regions import capture_sky_selection

    if pypeit.__version__ != '3.0.0':
        raise ValueError('extraction is validated only with PypeIt 3.0.0')
    data, waveair = np.asarray(data, dtype=float), np.asarray(waveair, dtype=float)
    if (data.ndim != 2 or data.shape != waveair.shape or min(data.shape) < 3
            or not np.isfinite(waveair).all() or np.any(waveair <= 0)
            or not np.all(np.diff(waveair, axis=0) > 0)):
        raise ValueError('science image and increasing wavelength map must have matching 2D geometry')
    if not np.isfinite([exposure_time, read_noise]).all() or exposure_time <= 0 or read_noise <= 0:
        raise ValueError('exposure time and read noise must be positive and finite')
    if sky_input_mask is not None:
        sky_input_mask = np.asarray(sky_input_mask, dtype=bool)
        if sky_input_mask.shape != data.shape or np.any(sky_input_mask.sum(axis=1) < 8):
            raise ValueError('sky mask must match the image and retain eight pixels per row')
    if bounded_background and (sky_model is not None or sky_input_mask is not None or
                               (manual is not None and (len(manual)!=1 or manual[0]['id']!='primary'))):
        raise ValueError('automatic bounded background cannot ignore an explicit recipe')
    if sky_model is not None and (sky_input_mask is None or manual is None or len(manual) != 1):
        raise ValueError('bounded PypeIt background needs one explicit source and two selected bands')
    inst = EFOSpecEFOSC2Spectrograph()
    par = inst.profile_parameters('rapid')
    par['flexure']['spec_method'] = 'skip'
    from .trace_geometry import TRACE_DEGREE, trace_evidence, measured_aperture
    par['reduce']['findobj']['trace_npoly'] = TRACE_DEGREE
    par['reduce']['findobj']['maxnumber_sci'] = 3
    par['reduce']['findobj']['snr_thresh'] = 5
    par['reduce']['findobj']['skip_second_find'] = True
    par['reduce']['skysub']['bspline_spacing'] = 1.0
    par['rdx']['redux_path'] = str(output)
    if sky_input_mask is not None or bounded_background:
        par['reduce']['skysub']['no_local_sky'] = True
    ny, nx = data.shape
    tilts = (waveair-waveair.min())/(waveair.max()-waveair.min())
    rn = np.full(data.shape, float(read_noise)**2)
    if processed is None:
        sci = PypeItImage(data, ivar=1/(np.maximum(data, 0)+rn), rn2img=rn, base_var=rn,
                         detector=detector, exptime=exposure_time, noise_floor=0., shot_noise=True,
                         bpm=~np.isfinite(data))
    else:
        required = {'ivar', 'rn2img', 'base_var', 'img_scale', 'bad_pixels'}
        if set(processed) != required or any(np.shape(processed[k]) != data.shape for k in required):
            raise ValueError('processed variance and mask planes must match the image')
        if any(np.any(~np.isfinite(processed[k])) or np.any(processed[k] < 0)
               for k in required - {'bad_pixels'}):
            raise ValueError('processed variance planes must be finite and nonnegative')
        sci = PypeItImage(data, ivar=processed['ivar'], rn2img=processed['rn2img'],
                         base_var=processed['base_var'], img_scale=processed['img_scale'],
                         bpm=processed['bad_pixels'], detector=detector, exptime=exposure_time,
                         noise_floor=0., shot_noise=True)
    gain = float(detector['gain'][0])
    sci.build_mask(saturation=detector['saturation']*gain, mincounts=-1e10)
    slits = SlitTraceSet(np.zeros(ny), np.full(ny, nx-1.), 'MultiSlit',
                        detname='DET01', nspat=nx, PYP_SPEC=inst.name, binspec=2, binspat=2,
                        specmin=np.array([-.5]), specmax=np.array([ny-.5]))
    instruction = None
    if manual is not None:
        if not manual:
            raise ValueError('manual extraction needs at least one source')
        instruction = ManualExtractionObj(
            frame=basename, spec=np.array([o['spec'] for o in manual], dtype=float),
            spat=np.array([o['spat'] for o in manual], dtype=float),
            fwhm=np.array([o['fwhm'] for o in manual], dtype=float),
            boxcar_rad=np.array([o['radius'] for o in manual], dtype=float),
            detname=np.array(['DET01']*len(manual)), neg=np.zeros(len(manual), dtype=bool))
    finder = FindObjects.get_instance(sci, slits, inst, par, 'science', tilts=tilts,
                                     show=False, basename=basename, manual=instruction,
                                     initial_skymask=sky_input_mask)

    def save_detection_failure(exc, found, sky, selection=None):
        from pathlib import Path
        from tempfile import NamedTemporaryFile
        from .pypeit_diagnostics import write_spec2d
        if selection is None:
            selection=capture_sky_selection(finder,found)
            selection.update(sky_input_mask=np.zeros_like(data,dtype=bool),
                sky_mask_stage=np.array('source_finding_only_no_accepted_background_fit'))
        selection, keep = validated_candidate_packet(selection, found, evidence)
        exc.diagnostics.update(SourceSelectionError(str(exc),found[keep],target_position).diagnostics)
        exc.diagnostics['rejected_candidates'] = selection['rejected_names'].tolist()
        directory=Path(output)/'native';directory.mkdir(parents=True,exist_ok=True)
        with NamedTemporaryFile(dir=directory,prefix=basename+'_selection_failed_',
                                suffix='.npz',delete=False) as stream:
            np.savez_compressed(stream,data=data,waveair=waveair,sky=sky,**selection)
            diagnostic=Path(stream.name)
        exc.diagnostics['candidate_evidence'] = dict(evidence)
        if aperture_info is not None:
            exc.diagnostics['aperture'] = aperture_info
        exc.diagnostics['diagnostic_file']=str(diagnostic.resolve())
        exc.diagnostics['native_spec2d']=write_spec2d(
            directory/('spec2d_'+diagnostic.stem+'.fits'),sci,slits,
            airtovac(waveair*u.AA).value,tilts,sky,
            stage='source_selection_failed' if isinstance(exc,SourceSelectionError) else 'background_support_failed')

    fitted_background = None
    # An automatic selection already fitted its winning configuration; reuse it
    # rather than repeating the same PypeIt solve.
    selected_background = None
    background_radius=manual[0]['radius'] if sky_model is not None else None
    auto_info=None
    unresolved=False
    detection_ivar=sci.ivar
    detection_support=None
    evidence = {}
    aperture_info = None

    def validate_candidates(found, image):
        for source in found:
            check = trace_evidence(image, sci.ivar, source.TRACE_SPAT,
                                   float(np.nanmedian(source.FWHM)),
                                   support=sci.select_flag(invert=True) & detection_support
                                   if detection_support is not None else sci.select_flag(invert=True))
            evidence[str(source.NAME)] = check
        # Explicit recipes remain explicit; automatic neighbours need evidence.
        return found[np.array([bool(s.hand_extract_flag) or evidence[str(s.NAME)]['accepted']
                               for s in found], dtype=bool)]

    def find_in_supported_background(image):
        # PypeIt's collapsed spatial profile uses the image mask as well as
        # inverse variance. Zero ivar alone does not exclude extrapolated sky.
        saved_mask=sci.fullmask.mask.copy()
        try:
            if detection_support is not None:
                sci.update_mask('BPM',indx=~detection_support)
            return finder.find_objects(image,detection_ivar)
        finally:
            sci.fullmask.mask[:]=saved_mask
    if bounded_background:
        from .local_sky import nearby_mask, BackgroundSupportError
        # Find traces without making a full-slit sky model. Ambiguity is still
        # handled by the diagnostic-writing selection gate below.
        found,finder.nobj=finder.find_objects(data,sci.ivar)
        try:
            first_names=select_sources(validate_candidates(found,data),target_position,manual)
        except SourceSelectionError:
            unresolved=True;sky=np.zeros_like(data)
        if not unresolved:
            primary=next(s for s in found if s.NAME==first_names['primary'])
            # PypeIt object finding expects a sky-subtracted frame. A raw
            # collapsed background pedestal can masquerade as a very broad
            # neighbour and veto every sky band. Bootstrap with finite bands,
            # then find sources again before choosing the final sampling bands.
            # This provisional model never becomes the delivered background.
            from .local_sky import preliminary_background
            initial_found = [{'window_pixel':float(s.SPAT_PIXPOS),
                              'fwhm_pixels':float(np.nanmedian(s.FWHM))} for s in found]
            # Test the initial detections against the full calibrated frame.
            # A validated outside source keeps its guard even if the later
            # restricted search cannot see it. An unmeasurable raw pedestal is
            # retained as evidence, not assigned an unbounded exclusion width.
            initial_neighbours = []
            for source, record in zip(found, initial_found):
                if source.NAME == primary.NAME:
                    continue
                check = trace_evidence(data, sci.ivar, source.TRACE_SPAT,
                    record['fwhm_pixels'], support=sci.select_flag(invert=True))
                record['local_contrast'] = check
                if check['accepted']:
                    initial_neighbours.append(source)
                    record['guard_retained'] = True
                else:
                    record['guard_retained'] = False
                    record['guard_reason'] = ('raw peak lacks measurable flanks; retained for review'
                        if not all(b['supported'] for b in check['bins']) else
                        'raw proposal fails positive-peak validation')
            seed_trace=primary.TRACE_SPAT
            seed_radius=max(15.,1.5*float(np.nanmedian(primary.FWHM)))
            try:
                seed_mask,seed=preliminary_background(data,
                    np.where(sci.select_flag(invert=True),sci.ivar,0.),tilts,seed_trace,
                    radius=seed_radius,fwhm=float(np.nanmedian(primary.FWHM)))
                # A polynomial sky model is only constrained between its
                # sampling bands. Never detect neighbours in its extrapolation.
                x=np.arange(nx)[None,:]
                left=np.min(np.where(seed_mask,x,nx),axis=1)
                right=np.max(np.where(seed_mask,x,-1),axis=1)
                detection_support=(x>=left[:,None])&(x<=right[:,None])
                detection_ivar=np.where(detection_support,sci.ivar,0.)
                found,finder.nobj=find_in_supported_background(data-seed['sky'])
                accepted = validate_candidates(found, data-seed['sky'])
                first_names=select_sources(accepted,target_position,manual)
                primary=next(s for s in found if s.NAME==first_names['primary'])
            except (SourceSelectionError,BackgroundSupportError) as exc:
                save_detection_failure(exc,found,np.zeros_like(data))
                raise
            sky_trace=primary.TRACE_SPAT
            fwhm=float(np.nanmedian(primary.FWHM))
            background_radius=(float(manual[0]['radius']) if manual is not None else
                float(par['reduce']['extraction']['boxcar_radius'])/(float(detector['platescale'])*2))
            neighbours=[(s.TRACE_SPAT,float(np.nanmedian(s.FWHM)),background_radius)
                        for s in accepted if s.NAME!=primary.NAME]
            for source in initial_neighbours:
                fw = float(np.nanmedian(source.FWHM))
                # Keep the broader original guard even if a wing is detected
                # as a narrower source after the provisional subtraction.
                neighbours.append((source.TRACE_SPAT, fw, background_radius))
            if manual is None:
                aperture_info = dict(policy='configured_aperture_with_diagnostic_wings',
                    adopted_radius_pixels=background_radius, total_flux_guaranteed=False)
                try:
                    measured_radius, wing_info = measured_aperture(
                        data-seed['sky'], np.where(sci.select_flag(invert=True),sci.ivar,0.),
                        sky_trace, fwhm, neighbours=neighbours)
                    aperture_info.update(wing_measurement_status='measured',
                        measured_radius_pixels=measured_radius, wing_measurement=wing_info)
                except ValueError as error:
                    aperture_info.update(wing_measurement_status='unavailable',
                        wing_measurement_error=str(error))
                # Wing growth is a diagnostic, not a total-flux requirement.
                # Keep the configured radius for extraction and validate its
                # actual two-sided background support below.
            # Let the sampled residuals choose the band width and spatial
            # degree. A steeply varying background under the slit is not
            # described by the default wide band and quadratic term.
            from .local_sky import select_bounded_background
            try:
                sky_input_mask,selected_background,auto_selection=select_bounded_background(
                    data,np.where(sci.select_flag(invert=True),sci.ivar,0.),tilts,sky_trace,
                    radius=background_radius,fwhm=fwhm,neighbours=neighbours,
                    allowed_mask=sci.select_flag(invert=True),max_gap_pixels=60,sigma_clip=4.)
            except BackgroundSupportError as exc:
                exc.diagnostics['preliminary_detection_background']=seed['summary']
                save_detection_failure(exc,found,np.zeros_like(data))
                raise
            finder.initial_skymask=sky_input_mask
            # The final finite bands define the supported detection interval.
            left=np.min(np.where(sky_input_mask,x,nx),axis=1)
            right=np.max(np.where(sky_input_mask,x,-1),axis=1)
            detection_support=(x>=left[:,None])&(x<=right[:,None])
            detection_ivar=np.where(detection_support,sci.ivar,0.)
            sky_model={'degree':auto_selection['selected_spatial_degree'],'sigma_clip':4.}
            auto_info={'requested_width_spatial_pixels':30,'minimum_width_spatial_pixels':15,
                       'actual_width_spatial_pixels':auto_selection['selected_width_spatial_pixels'],
                       'selection':'nearest_clean_two_sided_bands',
                       'automatic_background_selection':auto_selection,
                       'host_identity_inferred':False,
                       'preliminary_detection_background':seed['summary'],
                       'preliminary_source_candidates':initial_found,
                       'detection_sequence':'raw seed trace, bounded preliminary sky, source finding, bounded final sky'}
    if sky_model is None and not unresolved:
        sky, found = finder.run()
    elif sky_model is not None:
        from .local_sky import fit_background, BackgroundSupportError
        if set(sky_model) != {'degree', 'sigma_clip'}:
            raise ValueError('bounded PypeIt background requires degree and sigma_clip')
        trace = (np.full(ny, manual[0]['spat']) if sky_trace is None else np.asarray(sky_trace, float))
        if selected_background is not None:
            fitted_background = selected_background
        else:
            try:
                fitted_background = fit_background(data, np.where(sci.select_flag(invert=True),sci.ivar,0.), tilts,
                    sky_input_mask, trace,
                    radius=background_radius, **sky_model)
            except BackgroundSupportError as exc:
                if bounded_background:save_detection_failure(exc,found,np.zeros_like(data))
                raise
        if auto_info is not None:fitted_background['summary']['selection']=auto_info
        sky = fitted_background['sky']
        sci.update_mask('BPM', indx=np.broadcast_to(~fitted_background['valid_rows'][:, None], data.shape))
        # FindObjects.run would fit the sky again. Use its public object-finding
        # stage with our audited PypeIt model, then hold it during extraction.
        found, finder.nobj = find_in_supported_background(data-sky)
    accepted = validate_candidates(found, data-sky)
    if aperture_info is not None:
        for source in accepted:
            if abs(source.SPAT_PIXPOS-primary.SPAT_PIXPOS) <= fwhm:
                source.BOX_R_PIX = background_radius
    selection = capture_sky_selection(finder, found)
    selection, _ = validated_candidate_packet(selection, found, evidence)
    if unresolved:
        selection.update(sky_input_mask=np.zeros_like(data,dtype=bool),
                         sky_mask_stage=np.array('source_finding_only_no_background_fit'))
    if fitted_background is not None:
        selection.update(sky_input_mask=fitted_background['input_mask'],
                         sky_accepted_mask=fitted_background['accepted_mask'],
                         sky_valid_rows=fitted_background['valid_rows'],
                         sky_side_counts=fitted_background['side_counts'],
                         sky_mask_stage=np.array('bounded_input_and_post_rejection_masks'),
                         sky_spatial_degree=np.array(sky_model['degree']))
    try:
        names = select_sources(accepted, target_position, manual)
    except SourceSelectionError as exc:
        save_detection_failure(exc,found,sky,selection)
        raise
    extractor = Extract.get_instance(sci, slits, accepted, inst, par, 'science', global_sky=sky,
                                     tilts=tilts, waveimg=airtovac(waveair*u.AA).value,
                                     show=False, basename=basename)
    native_arrays = extractor.run()
    # PypeIt 3 adds the separate A-frame sky model at index 1. Preserve
    # EFOSpec's existing [sky, object, ivar, mask, spectra, ...] contract.
    arrays = (native_arrays[0], *native_arrays[2:])
    if fitted_background is not None:
        from .local_sky import BackgroundSupportError, side_counts
        for obj in arrays[4]:
            if obj.NAME not in names.values():
                continue
            supported = np.all(side_counts(fitted_background['accepted_mask'],
                                          obj.TRACE_SPAT, background_radius) >= 5, axis=0)
            # A drifting fitted trace must not silently cross a sampling band.
            x = np.arange(nx)[None, :]
            overlap = np.any(fitted_background['input_mask'] &
                             (abs(x-obj.TRACE_SPAT[:, None]) <= background_radius), axis=1)
            if np.any((~supported | overlap) & fitted_background['valid_rows']):
                exc=BackgroundSupportError('Fitted source trace moved outside the two-sided background support; review bands.',
                    initial_trace_median=float(np.median(trace)),fitted_trace_median=float(np.median(obj.TRACE_SPAT)),
                    maximum_trace_change_pixels=float(np.max(abs(obj.TRACE_SPAT-trace))),
                    unsupported_rows=np.flatnonzero((~supported|overlap)&fitted_background['valid_rows']).tolist())
                save_detection_failure(exc,found,sky,selection)
                raise exc
    # The profile fitter can adjust a low-order seed with splines. Reject
    # large non-polynomial corrections rather than cosmetically smoothing a
    # trace after extracting flux along a different path.
    for obj in arrays[4]:
        if obj.NAME not in names.values():
            continue
        y = np.linspace(-1., 1., ny)
        low_order = np.polynomial.legendre.legval(y, np.polynomial.legendre.legfit(y, obj.TRACE_SPAT, 4))
        if np.max(abs(obj.TRACE_SPAT-low_order)) > .5:
            exc = SourceSelectionError('Extracted trace departs from degree <=4 by more than half a pixel; review required.', found, target_position)
            save_detection_failure(exc, found, sky, selection)
            raise exc
    extracted = arrays[4]
    objects = {}
    for identity, name in names.items():
        selected = np.array([i for i, obj in enumerate(extracted) if obj.NAME == name])
        if len(selected) != 1:
            raise ValueError(f'target {identity} lost during extraction')
        objects[identity] = extracted[selected]
    result = {'objects': objects, 'found': accepted, 'arrays': arrays, 'selection': selection,
              'diagnostic_context':{'science':sci,'slits':slits,'tilts':tilts},
              'sky_policy': 'selected_global_sky' if sky_input_mask is not None else 'automatic_global_and_local'}
    if manual is None:
        result['source_selection'] = {
            'method':'brightest_positive_central_peak', 'half_width_pixels':15,
            'expected_window_pixel':float(target_position), 'identity_confirmed':False,
            'candidates':[{'name':str(s.NAME),'window_pixel':float(s.SPAT_PIXPOS),
                           'peak':float(s.smash_peakflux) if s.smash_peakflux is not None else None,
                           'snr':float(s.smash_snr) if s.smash_snr is not None else None,
                           'trace_evidence':evidence.get(str(s.NAME))} for s in accepted]}
    result['trace_geometry'] = dict(trace_polynomial_degree=TRACE_DEGREE, candidate_evidence=evidence,
                                    aperture=aperture_info, maximum_nonpolynomial_departure_pixels=.5)
    if fitted_background is not None:
        result.update(sky_policy='bounded_two_sided_pypeit', background_fit=fitted_background['summary'])
    return result
