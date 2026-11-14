"""Two-sided spatial sampling and the actual PypeIt background kernel."""
import numpy as np
import pytest

from efospec.local_sky import BackgroundSupportError, nearby_mask, fit_background, side_counts


def test_nearest_short_clean_run_wins_over_distant_full_width():
    from efospec.local_sky import closest_mask
    trace=np.full(20,100.)
    # The neighbour leaves an 18-pixel clean run next to the source guard.
    mask,bounds,info=closest_mask((20,240),trace,radius=7.,fwhm=4.,
        neighbours=[(trace+36,6.,7.)])
    assert np.all(bounds[:,1,0]==108)
    assert np.all(bounds[:,1,1]<=127)
    assert 15<=info['minimum_actual_width_pixels']<30
    assert np.all(side_counts(mask,trace,7.)>=15)


def test_closest_bands_follow_curvature_and_fail_if_one_side_crowded():
    from efospec.local_sky import closest_mask
    trace=95+3*np.sin(np.arange(100)/25.)
    mask,bounds,_=closest_mask((100,240),trace,radius=7.,fwhm=4.)
    assert np.ptp(bounds[:,0,0])>2
    assert not np.any(mask&(abs(np.arange(240)-trace[:,None])<=7.5))
    with pytest.raises(BackgroundSupportError,match='two-sided'):
        closest_mask((100,115),trace,radius=7.,fwhm=4.)


def test_bands_follow_curved_trace_avoid_neighbour_and_keep_fixed_width():
    trace = 95+3*np.sin(np.arange(100)/25.)
    neighbour = trace+22
    mask, bounds = nearby_mask((100, 240), trace, radius=6., fwhm=4., width_pixels=40,
                              neighbours=[(neighbour, 4., 6.)])
    assert np.all(side_counts(mask, trace, 6.) == 40)
    assert np.all(bounds[:, :, 1]-bounds[:, :, 0] == 40)
    assert not mask[:, :3].any() and not mask[:, -3:].any()
    x = np.arange(240)
    assert not np.any(mask & (abs(x-neighbour[:, None]) <= 6.5))
    assert np.all(bounds[:, 1, 0] > neighbour+6.5)
    assert np.ptp(bounds[:, 0, 0]) > 2


def test_crowding_or_window_edge_never_falls_back_to_one_side_or_whole_slit():
    trace = np.full(40, 70.)
    with pytest.raises(BackgroundSupportError, match='both sides'):
        nearby_mask((40, 141), trace, radius=7., fwhm=5., width_pixels=40,
                    neighbours=[(trace+23., 6., 7.), (trace+50., 5., 7.)])


def test_authorized_adaptive_width_stops_at_15_pixels():
    from types import SimpleNamespace
    from efospec.local_sky import exposure_sky_inputs
    from efospec.recipes import DetectorWindow
    trace=np.full(40,35.)
    exposure=SimpleNamespace(window=DetectorWindow((40,85),(0,0)),
        images={'found_traces':np.array([trace,trace+40])})
    recipe={'objects':[dict(spec=20.,spat=35.,radius=7.,fwhm=4.)],
            'sky':dict(mode='pypeit_local',include=[],exclude=[],width_pixels=30,minimum_width_pixels=15)}
    mask,selected,info=exposure_sky_inputs(exposure,recipe)
    assert 15<=info['actual_width_spatial_pixels']<30
    assert np.all(side_counts(mask,selected,7.)==info['actual_width_spatial_pixels'])


def curved_background():
    y, x = np.indices((160, 180))
    trace = 89.+.01*(np.arange(160)-80.)
    true_sky = 250.+1.2*(x-89)+.028*(x-89)**2+60.*np.exp(-.5*((y-70.)/3.)**2)
    source = 3000.*np.exp(-.5*((x-trace[:, None])/1.7)**2)
    rng = np.random.default_rng(316)
    variance = true_sky+source+7.8**2
    data = true_sky+source+rng.normal(size=x.shape)*np.sqrt(variance)
    mask, _ = nearby_mask(data.shape, trace, radius=7., fwhm=4., width_pixels=35)
    data[62, 55] += 10000.
    tilts = (y+.004*x)/(159+.004*179)
    return data, variance, tilts, mask, trace, true_sky


def test_detector_defects_leave_holes_without_relocating_background_bands():
    from efospec.local_sky import select_bounded_background
    data, variance, tilts, _, trace, _ = curved_background()
    allowed=np.ones(data.shape,bool)
    allowed[30:33,60:70]=False
    allowed[90:93,105:115]=False
    mask,fit,selection=select_bounded_background(data,1/variance,tilts,trace,
        radius=7.,fwhm=4.,allowed_mask=allowed,widths=(25,),degrees=(2,))
    expected,_=nearby_mask(data.shape,trace,
        radius=selection['selected_inner_edge_spatial_pixels'],fwhm=4.,width_pixels=25)
    np.testing.assert_array_equal(mask,expected)
    assert not np.any(fit['input_mask'] & ~allowed)
    assert not np.any(fit['accepted_mask'] & ~allowed)
    assert np.all(fit['side_counts'][:,fit['valid_rows']]>=5)


def test_quadratic_pypeit_model_recovers_curved_background_and_rejects_cosmic():
    data, variance, tilts, mask, trace, truth = curved_background()
    fit = fit_background(data, 1/variance, tilts, mask, trace, radius=7., degree=2)
    center = np.round(trace).astype(int)
    rows = np.arange(15, 145)
    error = fit['sky'][rows, center[rows]]-truth[rows, center[rows]]
    assert abs(np.median(error)) < 2.
    assert np.sqrt(np.mean(error**2)) < 5.
    assert not fit['accepted_mask'][62, 55]
    assert not fit['accepted_mask'][~mask].any()
    assert np.all(fit['side_counts'][:, fit['valid_rows']] >= 5)
    assert fit['summary']['legendre_terms'] == 3
    assert fit['summary']['pypeit_exit_status'] == 0
    assert fit['summary']['constant_fallback'] is False


def test_one_sided_mask_rejected_before_fit():
    data, variance, tilts, mask, trace, _ = curved_background()
    mask[:, 90:] = False
    with pytest.raises(BackgroundSupportError, match='each side'):
        fit_background(data, 1/variance, tilts, mask, trace, radius=7.)


def test_bad_detector_rows_are_masked_on_both_sides_without_one_sided_fits():
    data,variance,tilts,mask,trace,truth=curved_background()
    ivar=1/variance;ivar[80,:90]=0
    result=fit_background(data,ivar,tilts,mask,trace,radius=7.,degree=2)
    assert not result['valid_rows'][80]
    assert not result['accepted_mask'][80].any()
    assert not result['input_mask'][80].any()
    assert result['summary']['pre_fit_invalid_rows']==[80]
    assert result['valid_rows'].mean()>.98
    ivar[50:80,:90]=0
    with pytest.raises(BackgroundSupportError,match='10%'):
        fit_background(data,ivar,tilts,mask,trace,radius=7.,degree=2)


def test_nonconvergence_never_changes_polynomial_order(monkeypatch):
    from pypeit.core import fitting
    calls = []
    def failed(*args, **kwargs):
        calls.append(kwargs['basis'].shape)
        return None, np.ones_like(args[0], bool), np.zeros_like(args[0]), 1., 1
    monkeypatch.setattr(fitting, 'iterative_bspline_fit', failed)
    data, variance, tilts, mask, trace, _ = curved_background()
    with pytest.raises(BackgroundSupportError, match='no lower-order fallback'):
        fit_background(data, 1/variance, tilts, mask, trace, radius=7., degree=2)
    assert len(calls) == 1


def test_rejection_removing_one_side_stops_extraction(monkeypatch):
    from pypeit.core import fitting
    def one_sided(*args, **kwargs):
        accepted = kwargs['gpm'].copy()
        # The second Legendre term encodes the spatial coordinate.
        accepted &= kwargs['basis'][:, 1] < 0
        return None, accepted, np.ones_like(args[0]), 1., 0
    monkeypatch.setattr(fitting, 'iterative_bspline_fit', one_sided)
    data, variance, tilts, mask, trace, _ = curved_background()
    with pytest.raises(BackgroundSupportError, match='after rejection'):
        fit_background(data, 1/variance, tilts, mask, trace, radius=7., degree=2)


def test_pypeit_optimal_extraction_holds_bounded_quadratic_model(tmp_path):
    from astropy.io import fits
    from test_pypeit_adapter import write_raw
    from efospec.pypeit_adapter import EFOSpecEFOSC2Spectrograph
    from efospec.extraction import extract_spectra
    data, variance, tilts, mask, trace, truth = curved_background()
    raw = write_raw(tmp_path,'raw.fits')
    with fits.open(raw) as hdus:
        detector = EFOSpecEFOSC2Spectrograph().get_detector_par(1,hdus)
    (tmp_path/'QA/PNGs').mkdir(parents=True)
    result = extract_spectra(data, 5000.+800.*tilts, detector,
        exposure_time=600., read_noise=7.8, output=tmp_path, basename='local-quadratic',
        manual=[dict(id='primary',spec=80.,spat=89.,fwhm=4.,radius=7.)],
        sky_input_mask=mask, sky_trace=trace, sky_model={'degree':2,'sigma_clip':4.})
    obj = result['objects']['primary'][0]
    rows = obj.OPT_MASK & (np.arange(len(data))>25) & (np.arange(len(data))<135)
    assert rows.sum() > 90
    true_flux = 3000.*np.sqrt(2*np.pi)*1.7
    assert abs(np.median(obj.OPT_COUNTS[rows])/true_flux-1) < .02
    assert result['sky_policy'] == 'bounded_two_sided_pypeit'
    assert result['background_fit']['spatial_degree'] == 2
    assert not result['selection']['sky_accepted_mask'][~mask].any()
    # A later full-slit/constant refit would erase the strong curvature.
    y,x = np.indices(data.shape)
    aperture = (abs(x-trace[:,None])<7) & rows[:,None]
    assert abs(np.median(result['arrays'][0][aperture]-truth[aperture])) < 3.


@pytest.mark.parametrize('background_curvature',[0.,.12])
@pytest.mark.parametrize('wing_measurement',['measured','unavailable','raw_pedestal'])
def test_automatic_bounded_default_recovers_flux_without_full_slit_fit(tmp_path,monkeypatch,background_curvature,wing_measurement):
    from astropy.io import fits
    from test_pypeit_adapter import write_raw
    from test_extraction import synthetic_image
    from efospec.pypeit_adapter import EFOSpecEFOSC2Spectrograph
    from efospec.extraction import extract_spectra
    from pypeit.core import skysub
    from efospec import trace_geometry
    if wing_measurement == 'raw_pedestal':
        from copy import deepcopy
        from pypeit.find_objects import MultiSlitFindObjects
        original=MultiSlitFindObjects.find_objects
        calls=[]
        def with_raw_pedestal(self,*args,**kwargs):
            objects,count=original(self,*args,**kwargs)
            if not calls:
                pedestal=deepcopy(objects[0])
                pedestal.NAME='UNVALIDATED_RAW_PEDESTAL'
                pedestal.TRACE_SPAT=np.full_like(pedestal.TRACE_SPAT,110.)
                pedestal.SPAT_PIXPOS=110.
                pedestal.FWHM=90.
                objects.add_sobj(pedestal);count+=1
            calls.append(True)
            return objects,count
        monkeypatch.setattr(MultiSlitFindObjects,'find_objects',with_raw_pedestal)
    if wing_measurement == 'unavailable':
        def unavailable(*args, **kwargs):
            raise ValueError('no uncontaminated room to measure both aperture wings')
        monkeypatch.setattr(trace_geometry, 'measured_aperture', unavailable)
    def forbidden(*args,**kwargs):raise AssertionError('full-slit global sky called')
    monkeypatch.setattr(skysub,'global_skysub',forbidden)
    data,wave,truth=synthetic_image()
    pedestal=background_curvature*(np.arange(data.shape[1])-70.)**2
    data=data+pedestal[None,:]+np.random.default_rng(42).normal(0,np.sqrt(pedestal),data.shape)
    raw=write_raw(tmp_path,'raw.fits')
    with fits.open(raw) as hdus:detector=EFOSpecEFOSC2Spectrograph().get_detector_par(1,hdus)
    (tmp_path/'QA/PNGs').mkdir(parents=True)
    result=extract_spectra(data,wave,detector,exposure_time=600.,read_noise=7.8,
        output=tmp_path,basename='bounded-auto',target_position=48.,bounded_background=True)
    obj=result['objects']['primary'][0];rows=obj.OPT_MASK&(np.arange(160)>25)&(np.arange(160)<135)
    assert abs(np.median(obj.OPT_COUNTS[rows]/truth['primary'][rows])-1)<.02
    assert 15<=result['background_fit']['selection']['actual_width_spatial_pixels']<=30
    aperture=result['trace_geometry']['aperture']
    configured=EFOSpecEFOSC2Spectrograph().default_pypeit_par()['reduce']['extraction']['boxcar_radius']/(float(detector['platescale'])*2)
    assert obj.BOX_R_PIX == pytest.approx(configured)
    assert aperture['adopted_radius_pixels'] == pytest.approx(configured)
    assert aperture['wing_measurement_status'] == ('measured' if wing_measurement=='raw_pedestal' else wing_measurement)
    assert aperture['total_flux_guaranteed'] is False
    assert result['background_fit']['minimum_side_samples']>=5
    if wing_measurement=='raw_pedestal':
        proposals=result['background_fit']['selection']['preliminary_source_candidates']
        rejected=next(p for p in proposals if p['fwhm_pixels']==90.)
        assert not rejected['guard_retained']
        assert not rejected['local_contrast']['accepted']
        assert any(p['guard_retained'] for p in proposals if p.get('local_contrast',{}).get('accepted'))


def test_reviewed_source_can_recompute_automatic_bounded_bands(tmp_path):
    from astropy.io import fits
    from test_pypeit_adapter import write_raw
    from test_extraction import synthetic_image
    from efospec.pypeit_adapter import EFOSpecEFOSC2Spectrograph
    from efospec.extraction import extract_spectra
    data,wave,truth=synthetic_image()
    raw=write_raw(tmp_path,'raw.fits')
    with fits.open(raw) as hdus:detector=EFOSpecEFOSC2Spectrograph().get_detector_par(1,hdus)
    (tmp_path/'QA/PNGs').mkdir(parents=True)
    result=extract_spectra(data,wave,detector,exposure_time=600.,read_noise=7.8,
        output=tmp_path,basename='reviewed',target_position=98.,bounded_background=True,
        manual=[dict(id='primary',spec=80.,spat=48.,fwhm=4.,radius=7.)])
    obj=result['objects']['primary'][0]
    rows=obj.OPT_MASK&(np.arange(160)>25)&(np.arange(160)<135)
    assert np.median(obj.TRACE_SPAT)==pytest.approx(48.,abs=.2)
    assert abs(np.median(obj.OPT_COUNTS[rows]/truth['primary'][rows])-1)<.02
    assert result['background_fit']['minimum_side_samples']>=5


def steep_background(curvature, scale=44.):
    """A smooth slit profile whose curvature a wide-band quadratic must arch over."""
    y, x = np.indices((160, 180))
    trace = 89.+.01*(np.arange(160)-80.)
    # A broad, flat-topped host hump, as under AT2026agvs: a quadratic through
    # two distant bands bulges above it across the unsampled aperture gap, while
    # a band narrow enough to sit beside the aperture is locally almost linear.
    true_sky = 700.+curvature*np.exp(-.5*((x-89.)/scale)**4)
    source = 2500.*np.exp(-.5*((x-trace[:, None])/1.7)**2)
    rng = np.random.default_rng(2026)
    variance = true_sky+source+7.8**2
    data = true_sky+source+rng.normal(size=x.shape)*np.sqrt(variance)
    tilts = (y+.004*x)/(159+.004*179)
    return data, 1/variance, tilts, trace


def haloed_source():
    """A source whose power-law wings never die inside the slit, as a standard star.

    The saved window is narrow enough that the profile is still falling at its
    edge, so no ring of columns inside it carries background without starlight.
    """
    y, x = np.indices((160, 90))
    trace = 44.+np.zeros(160)
    separation = (x-trace[:, None])/6.
    source = 90000./(1.+separation**2)**1.1
    true_sky = 40.+np.zeros_like(x, dtype=float)
    rng = np.random.default_rng(7)
    variance = true_sky+source+7.8**2
    data = true_sky+source+rng.normal(size=x.shape)*np.sqrt(variance)
    tilts = (y+.004*x)/(159+.004*89)
    return data, 1/variance, tilts, trace


def test_a_flat_background_keeps_the_validated_default_geometry():
    from efospec.local_sky import select_bounded_background
    data, ivar, tilts, trace = steep_background(0.)
    mask, fit, selection = select_bounded_background(
        data, ivar, tilts, trace, radius=12., fwhm=4.)
    assert selection['selected_width_spatial_pixels'] == 30
    assert selection['selected_spatial_degree'] == 2
    assert selection['escalated'] is False
    assert selection['geometry_unresolved'] is False
    # Nothing beat the default by a worthwhile margin, and the measured error
    # beside the source is a small fraction of a single pixel's noise.
    assert abs(selection['holdout_residual_pixel_sigma']) < .1
    assert 'no offered width or spatial degree' in selection['reason']
    assert fit['summary']['spatial_degree'] == 2
    assert np.all(fit['side_counts'][:, fit['valid_rows']] >= 5)


def test_a_curved_background_is_measured_and_the_geometry_corrected():
    from efospec.local_sky import select_bounded_background, fit_background, nearby_mask
    data, ivar, tilts, trace = steep_background(700.)
    mask, fit, selection = select_bounded_background(
        data, ivar, tilts, trace, radius=12., fwhm=4.)
    assert selection['escalated'] is True
    default = next(t for t in selection['trials']
                   if t['width_spatial_pixels'] == 30 and t['spatial_degree'] == 2
                   and t['holdout_residual_counts'] is not None)
    # The choice is justified by a measurement, not by a preference for a shape.
    assert (abs(selection['holdout_residual_counts'])
            < abs(default['holdout_residual_counts']))
    assert (abs(default['holdout_residual_counts'])-abs(selection['holdout_residual_counts'])
            > selection['required_improvement_pixel_sigma']*selection['pixel_noise_counts'])
    # Every trial stayed bounded and two-sided.
    assert all(5 <= t['width_spatial_pixels'] <= 50 for t in selection['trials'])
    assert all(t['spatial_degree'] in (1, 2, 3) for t in selection['trials']
               if t['spatial_degree'] is not None)
    assert np.all(fit['side_counts'][:, fit['valid_rows']] >= 5)
    # And it really does subtract the background beside the source better.
    wide, _ = nearby_mask(data.shape, trace, radius=12., fwhm=4., width_pixels=30)
    old = fit_background(data, ivar, tilts, wide, trace, radius=12., degree=2)
    x = np.arange(data.shape[1])
    distance = np.abs(x[None, :]-trace[:, None])
    shoulder = (distance > 8.) & (distance <= 12.5) & fit['valid_rows'][:, None]
    truth = np.broadcast_to(700.+700.*np.exp(-.5*((x-89.)/44.)**4), data.shape)
    assert (abs(np.median((fit['sky']-truth)[shoulder]))
            < abs(np.median((old['sky']-truth)[shoulder])))


def test_the_withheld_ring_is_scored_unfitted_then_returned_to_the_final_fit():
    """Selection must not see the ring; the delivered model must still use it."""
    from efospec.local_sky import select_bounded_background, nearby_mask
    data, ivar, tilts, trace = steep_background(700.)
    mask, fit, selection = select_bounded_background(
        data, ivar, tilts, trace, radius=12., fwhm=4.)
    low, high = selection['holdout_ring_pixels']
    assert high == low+3 and low == selection['clean_background_radius_pixels']
    distance = np.abs(np.arange(data.shape[1])[None, :]-trace[:, None])
    ring = (distance >= low) & (distance <= high)
    scored = [t for t in selection['trials'] if t.get('role') == 'trial']
    assert scored and all(t['holdout_residual_counts'] is not None for t in scored
                          if t['status'] == 'fitted')
    # Selection never sampled the ring, and the delivered fit is not penalised
    # for it: the returned mask is the whole band that geometry is entitled to.
    full, _ = nearby_mask(data.shape, trace, radius=float(
        selection['selected_inner_edge_spatial_pixels']), fwhm=4.,
        width_pixels=selection['selected_width_spatial_pixels'])
    assert np.array_equal(mask, full)
    assert np.all(fit['side_counts'][:, fit['valid_rows']] >= 5)


def test_a_source_with_no_clean_ring_keeps_the_default_and_says_so():
    """A halo that never dies must not be chased by the geometry selector."""
    from efospec.local_sky import select_bounded_background, clean_background_radius
    data, ivar, tilts, trace = haloed_source()
    radius, _ = clean_background_radius(data, ivar, trace)
    assert radius is None
    # The narrow window only has room for the smaller permitted bands.
    _, fit, selection = select_bounded_background(
        data, ivar, tilts, trace, radius=12., fwhm=4., widths=(20, 15))
    assert selection['escalated'] is False
    assert selection['geometry_unresolved'] is True
    assert selection['selected_width_spatial_pixels'] == 20
    assert selection['selected_spatial_degree'] == 2
    assert selection['clean_background_radius_pixels'] is None
    assert selection['holdout_residual_counts'] is None
    assert 'no ring of columns' in selection['reason']
    assert fit['summary']['spatial_degree'] == 2


def test_selection_rejects_impossible_geometry_without_a_full_slit_fallback():
    from efospec.local_sky import select_bounded_background
    data, ivar, tilts, trace = steep_background(0.)
    with pytest.raises(BackgroundSupportError, match='No bounded two-sided background'):
        select_bounded_background(data, ivar, tilts, trace, radius=12., fwhm=4.,
                                  max_gap_pixels=0.)
    for kwargs in (dict(widths=()), dict(degrees=()), dict(holdout_width=1),
                   dict(holdout_width=4.), dict(improvement_sigma=0.),
                   dict(improvement_sigma=float('nan'))):
        with pytest.raises(ValueError):
            select_bounded_background(data, ivar, tilts, trace, radius=12., fwhm=4., **kwargs)


def test_the_clean_radius_grows_with_the_brightness_of_a_compact_source():
    """A brighter core stays above the noise further out and must be respected."""
    from efospec.local_sky import clean_background_radius
    radii = []
    for amplitude in (300., 30000.):
        y, x = np.indices((160, 180))
        trace = 89.+np.zeros(160)
        source = amplitude*np.exp(-.5*((x-trace[:, None])/2.6)**2)
        true_sky = 700.+np.zeros_like(x, dtype=float)
        rng = np.random.default_rng(11)
        variance = true_sky+source+7.8**2
        data = true_sky+source+rng.normal(size=x.shape)*np.sqrt(variance)
        radius, noise = clean_background_radius(data, 1/variance, trace)
        assert radius is not None and noise > 0
        radii.append(radius)
    assert radii[1] > radii[0]


def test_fractional_aperture_survives_final_background_refit(monkeypatch):
    from efospec import local_sky
    data, ivar, tilts, trace = steep_background(0.)
    monkeypatch.setattr(local_sky, 'clean_background_radius', lambda *args: (None, 1.))
    mask, fit, selection = local_sky.select_bounded_background(
        data, ivar, tilts, trace, radius=12.75, fwhm=4.)
    assert selection['selected_inner_edge_spatial_pixels'] == 12.75
    assert not np.any(mask & (abs(np.arange(data.shape[1])-trace[:, None]) <= 13.25))
    np.testing.assert_array_equal(mask, fit['input_mask'])


def test_background_holdout_never_scores_inside_measured_aperture(monkeypatch):
    from efospec import local_sky
    data, ivar, tilts, trace = steep_background(0.)
    monkeypatch.setattr(local_sky, 'clean_background_radius', lambda *args: (13, 1.))
    _, _, selection = local_sky.select_bounded_background(
        data, ivar, tilts, trace, radius=23.5, fwhm=4.)
    ring = selection['holdout_ring_pixels']
    assert ring is None or ring[0] >= 24
    if ring is None: assert selection['holdout_residual_counts'] is None


def test_geometry_failure_reports_that_solver_was_not_called(monkeypatch):
    from efospec import local_sky
    data, ivar, tilts, trace = steep_background(0.)
    def unexpected(*args, **kwargs): pytest.fail('impossible bands must fail before the solver')
    monkeypatch.setattr(local_sky, 'fit_background', unexpected)
    with pytest.raises(BackgroundSupportError, match='sky solver was not called'):
        local_sky.select_bounded_background(data, ivar, tilts, trace, radius=80., fwhm=40.)


def test_broad_neighbour_outside_bootstrap_search_remains_guarded(tmp_path, monkeypatch):
    from astropy.io import fits
    from test_pypeit_adapter import write_raw
    from efospec.pypeit_adapter import EFOSpecEFOSC2Spectrograph
    from efospec.extraction import extract_spectra
    from efospec import local_sky
    y,x=np.indices((160,201));trace=102.6+.002*(np.arange(160)-80)
    neighbour=trace+58.
    expected=150.+8000.*np.exp(-.5*((x-trace[:,None])/4.5)**2)
    expected+=2000.*np.exp(-.5*((x-neighbour[:,None])/7.85)**2)
    data=expected+np.random.default_rng(728).normal(size=x.shape)*np.sqrt(expected+7.8**2)
    wave=4500.+25.*y+np.zeros_like(x)
    raw=write_raw(tmp_path,'raw.fits')
    with fits.open(raw) as hdus:detector=EFOSpecEFOSC2Spectrograph().get_detector_par(1,hdus)
    (tmp_path/'QA/PNGs').mkdir(parents=True)
    original=local_sky.select_bounded_background;guarded=[]
    def select(*args,**kwargs):
        guarded.extend(kwargs['neighbours'])
        return original(*args,**kwargs)
    monkeypatch.setattr(local_sky,'select_bounded_background',select)
    try:
        result=extract_spectra(data,wave,detector,exposure_time=600.,read_noise=7.8,
            output=tmp_path,basename='broad-neighbour',target_position=102.6,bounded_background=True)
    except BackgroundSupportError:
        pass  # Refusing a contaminated extraction is a valid outcome.
    else:
        mask=result['selection']['sky_input_mask']
        assert not np.any(mask & (abs(x-neighbour[:,None])<=1.5*18.48))
    assert any(abs(np.median(t-neighbour))<1. and fw>16. for t,fw,radius in guarded)
