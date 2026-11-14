"""Bounded, two-sided background fitting with the pinned PypeIt B-spline solver.

The polynomial is spatial; the B-spline coordinate follows the wavelength map.
All widths here are *spatial detector pixels*, never Angstroms. No full-slit or
constant-background fallback is allowed. The selected and rejected pixels are
returned separately, including rows that lose support on either side.
"""
import numpy as np


class BackgroundSupportError(ValueError):
    code = 'background_support_required'
    action = 'Review both local background bands and neighbouring traces; no full-slit fallback was used.'

    def __init__(self, message, **diagnostics):
        super().__init__(message)
        self.diagnostics = diagnostics


def side_counts(mask, trace, radius):
    x = np.arange(mask.shape[1])[None, :]
    return np.array([(mask & (x < trace[:, None]-radius)).sum(axis=1),
                     (mask & (x > trace[:, None]+radius)).sum(axis=1)])


def validate_bounded_mask(mask, trace, radius, *, maximum_width=50):
    """Validate two finite local bands, including holes for bad pixels."""
    x = np.arange(mask.shape[1])
    for row, center in enumerate(trace):
        if np.any(mask[row] & (abs(x-center) <= radius)):
            raise BackgroundSupportError('Background selection includes the source aperture.', row=row)
        for side in (x < center-radius, x > center+radius):
            selected = np.flatnonzero(mask[row] & side)
            if len(selected) < 5:
                raise BackgroundSupportError('Background needs at least five samples on each side.', row=row)
            if selected[-1]-selected[0]+1 > maximum_width:
                raise BackgroundSupportError(f'Each local background band must span at most {maximum_width} spatial pixels.', row=row)
            if np.min(abs(selected-center))-radius > 2*maximum_width:
                raise BackgroundSupportError('Background band is too far from the source aperture.', row=row)


def closest_mask(shape, trace, *, radius, fwhm, width_pixels=30, minimum_width=15,
                 neighbours=(), max_gap_pixels=5., aperture_padding=.5):
    """Prioritize proximity, then width, for two trace-following clean bands.

    The protected core follows the existing 1.5-FWHM/aperture rule. A short
    nearby clean run wins over a wider distant run. No one-sided fallback.
    """
    trace=np.asarray(trace,float);ny,nx=shape
    if trace.shape!=(ny,) or not np.isfinite(trace).all():raise ValueError('invalid background trace')
    if not 15<=minimum_width<=width_pixels<=50:raise ValueError('background widths must be 15–50 pixels')
    if not np.isfinite([radius,fwhm,max_gap_pixels]).all() or min(radius,fwhm)<=0 or max_gap_pixels<0:
        raise ValueError('invalid nearest-background geometry')
    if not np.isfinite(aperture_padding) or aperture_padding<.5:raise ValueError('aperture padding must be at least half a pixel')
    x=np.arange(nx);guard=max(radius+aperture_padding,1.5*fwhm)
    clean=np.broadcast_to((x>=3)&(x<nx-3),shape).copy()
    clean &= abs(x[None,:]-trace[:,None])>guard
    for other,fw,rad in neighbours:
        other=np.asarray(other,float)
        if other.shape!=(ny,) or not np.isfinite(other).all() or min(fw,rad)<=0:
            raise ValueError('invalid neighbouring source')
        clean &= abs(x[None,:]-other[:,None])>max(rad+.5,1.5*fw)
    mask=np.zeros(shape,bool);bounds=np.full((ny,2,2),-1,int);gaps=np.full((ny,2),np.nan)
    for row,center in enumerate(trace):
        for side in (0,1):
            allowed=clean[row]&(x<center-guard if side==0 else x>center+guard)
            indices=np.flatnonzero(allowed)
            runs=np.split(indices,np.flatnonzero(np.diff(indices)>1)+1)
            choices=[r for r in runs if len(r)>=minimum_width]
            if not choices:continue
            chosen=choices[-1] if side==0 else choices[0]
            chosen=chosen[-width_pixels:] if side==0 else chosen[:width_pixels]
            gap=center-guard-chosen[-1] if side==0 else chosen[0]-center-guard
            if gap>max_gap_pixels:continue
            mask[row,chosen]=True;bounds[row,side]=chosen[0],chosen[-1]+1;gaps[row,side]=gap
    if np.any(bounds<0):
        raise BackgroundSupportError('Cannot place nearby two-sided bands without crossing source guards.',
            unsupported_rows=np.flatnonzero(np.any(bounds[:,:,0]<0,axis=1)).tolist())
    validate_bounded_mask(mask,trace,radius)
    return mask,bounds,dict(policy='nearest clean run first, then widest within that run',
        source_guard_pixels=float(guard),core_fwhm_pixels=float(fwhm),
        aperture_padding_pixels=float(aperture_padding),
        maximum_gap_beyond_guard_pixels=float(np.max(gaps)),
        minimum_actual_width_pixels=int(np.min(bounds[:,:,1]-bounds[:,:,0])),
        maximum_actual_width_pixels=int(np.max(bounds[:,:,1]-bounds[:,:,0])))


def nearby_mask(shape, trace, *, radius, fwhm, width_pixels, neighbours=(),
                max_gap_pixels=None, allowed_mask=None):
    """Find the nearest clean contiguous band of the requested width per side.

    Bands follow the trace, avoiding 1.5 FWHM (at least the boxcar aperture)
    around each source. A neighbouring source can displace a band, but never
    expands it. No narrower, one-sided or detector-edge fallback is made.
    Neighbours are (trace array, FWHM, boxcar radius) tuples in local pixels.
    """
    trace = np.asarray(trace, float)
    if (len(shape) != 2 or trace.shape != (shape[0],) or not np.isfinite(trace).all()
            or type(width_pixels) is not int or not 5 <= width_pixels <= 50
            or not np.isfinite([radius, fwhm]).all() or radius <= 0 or fwhm <= 0):
        raise ValueError('local background requires a finite trace and an explicit width of 5–50 spatial pixels')
    max_gap = 2*width_pixels if max_gap_pixels is None else float(max_gap_pixels)
    if not np.isfinite(max_gap) or max_gap < 0:
        raise ValueError('maximum background gap must be finite and nonnegative')
    x = np.arange(shape[1])[None, :]
    guard = max(radius+.5, 1.5*fwhm)
    clean = np.broadcast_to((x >= 3) & (x < shape[1]-3), shape).copy()
    if allowed_mask is not None:
        if np.shape(allowed_mask) != shape:
            raise ValueError('allowed background mask has the wrong shape')
        clean &= np.asarray(allowed_mask, bool)
    clean &= abs(x-trace[:, None]) > guard
    for other, other_fwhm, other_radius in neighbours:
        other = np.asarray(other, float)
        if (other.shape != trace.shape or not np.isfinite(other).all()
                or not np.isfinite([other_fwhm, other_radius]).all()
                or other_fwhm <= 0 or other_radius <= 0):
            raise ValueError('invalid neighbouring trace geometry')
        clean &= abs(x-other[:, None]) > max(other_radius+.5, 1.5*other_fwhm)
    mask = np.zeros(shape, bool)
    bounds = np.full((shape[0], 2, 2), -1, int)
    for row, center in enumerate(trace):
        for side in (0, 1):
            allowed = clean[row] & (x[0] < center-guard if side == 0 else x[0] > center+guard)
            starts = np.flatnonzero(np.convolve(allowed.astype(int), np.ones(width_pixels, int), 'valid') == width_pixels)
            if not len(starts):
                continue
            start = int(starts[-1] if side == 0 else starts[0])
            gap = center-guard-(start+width_pixels) if side == 0 else start-center-guard
            if gap > max_gap:
                continue
            mask[row, start:start+width_pixels] = True
            bounds[row, side] = start, start+width_pixels
    invalid = np.flatnonzero(np.any(bounds[:, :, 0] < 0, axis=1))
    if len(invalid):
        raise BackgroundSupportError(
            f'Cannot place {width_pixels}-pixel local bands on both sides in {len(invalid)}/{shape[0]} rows; '
            'review neighbours or prepare a wider calibrated window.',
            width_spatial_pixels=width_pixels, unsupported_rows=invalid.tolist(),
            window_shape=list(shape), max_gap_pixels=max_gap,
            source_guard_pixels=float(guard), source_fwhm_pixels=float(fwhm),
            source_radius_pixels=float(radius), source_trace_median=float(np.median(trace)),
            neighbour_guards=[dict(trace_median=float(np.median(t)),
                fwhm_pixels=float(fw), radius_pixels=float(rad),
                guard_pixels=float(max(rad+.5,1.5*fw))) for t,fw,rad in neighbours],
            unsupported_left_rows=int(np.sum(bounds[:,0,0]<0)),
            unsupported_right_rows=int(np.sum(bounds[:,1,0]<0)))
    return mask, bounds


def preliminary_background(data, ivar, tilts, trace, *, radius, fwhm):
    """Bootstrap source finding with bounded bands, narrowing before giving up.

    This provisional quadratic fit is never the delivered background. The
    final selector still measures width/degree alternatives after source finding.
    """
    trials=[]
    for width in (30,25,20,15):
        try:
            mask,_=nearby_mask(np.shape(data),trace,radius=radius,fwhm=fwhm,
                               width_pixels=width,max_gap_pixels=60)
        except BackgroundSupportError as exc:
            trials.append(dict(width_spatial_pixels=width,status='geometry_unavailable',reason=str(exc)))
            continue
        try:
            fit=fit_background(data,ivar,tilts,mask,trace,radius=radius,degree=2,sigma_clip=4.)
        except BackgroundSupportError as exc:
            trials.append(dict(width_spatial_pixels=width,status='fit_failed',reason=str(exc)))
            continue
        trials.append(dict(width_spatial_pixels=width,status='fitted'))
        fit['summary'].update(selected_width_spatial_pixels=width,width_trials=trials)
        return mask,fit
    raise BackgroundSupportError('Preliminary background failed with 30, 25, 20 and 15-pixel two-sided bands.',
        phase='preliminary_source_detection',trials=trials)


def exposure_sky_inputs(exposure, recipe):
    """Translate a hash-bound recipe into a trace-following local selection.

    Saved candidate traces guide the geometry, never the astrophysical identity.
    Legacy records do not store every neighbour's FWHM: use the requested source
    guard for them and record that limitation instead of inferring a host class.
    """
    from .recipes import sky_mask
    obj, sky = recipe['objects'][0], recipe['sky']
    y0, x0 = exposure.window.origin
    traces = np.asarray(exposure.images.get('found_traces', []), float)
    if not len(traces) and 'trace' in exposure.images:
        traces = np.asarray([exposure.images['trace']], float)
    if traces.ndim != 2 or traces.shape[1] != exposure.window.shape[0] or not len(traces):
        raise BackgroundSupportError('Saved exposure has no trace geometry for local background selection.')
    row = obj['spec']-y0
    anchors = np.array([np.interp(row, np.arange(len(t)), t) for t in traces])
    index = int(np.argmin(abs(anchors-(obj['spat']-x0))))
    if abs(anchors[index]-(obj['spat']-x0)) > obj['radius']:
        raise BackgroundSupportError('Requested source is too far from every saved trace; retrace before selecting background.')
    trace = traces[index]+obj['spat']-x0-anchors[index]
    neighbours = [(t,obj['fwhm'],obj['radius']) for i,t in enumerate(traces) if i != index]
    allowed = np.ones(exposure.window.shape, bool)
    x = np.arange(exposure.window.shape[1])+x0
    for low, high in sky['exclude']:
        allowed[:, (x >= low) & (x < high)] = False
    if not sky['include']:
        for actual_width in range(sky['width_pixels'],sky.get('minimum_width_pixels',sky['width_pixels'])-1,-1):
            try:
                mask, bounds = nearby_mask(exposure.window.shape, trace, radius=obj['radius'],
                    fwhm=obj['fwhm'], width_pixels=actual_width, neighbours=neighbours,
                    allowed_mask=allowed, max_gap_pixels=2*sky['width_pixels'])
                break
            except BackgroundSupportError:
                if actual_width==sky.get('minimum_width_pixels',sky['width_pixels']):raise
    else:
        mask = sky_mask(recipe, exposure.window) & allowed
        local_x = np.arange(exposure.window.shape[1])[None, :]
        for other, fwhm, radius in [(trace,obj['fwhm'],obj['radius'])]+neighbours:
            mask &= abs(local_x-other[:, None]) > max(radius+.5,1.5*fwhm)
        validate_bounded_mask(mask,trace,obj['radius'])
        bounds = None
    info = dict(coordinates='zero_based_extraction_window', raw_spatial_origin=x0,
                requested_raw_bands=sky['include'],
                selection='nearest_clean_two_sided_bands' if bounds is not None else 'explicit_two_sided_bands',
                width_spatial_pixels=sky.get('width_pixels'),
                actual_width_spatial_pixels=actual_width if bounds is not None else None,
                neighbour_guard='requested source FWHM/radius applied to saved candidates; inspect QC',
                source_guard_pixels=max(obj['radius']+.5,1.5*obj['fwhm']),
                saved_candidate_index=index+1)
    return mask, trace, info


# Automatic selection ladder. Width 30 with spatial degree 2 reproduces the
# previously validated default exactly and is always tried first; narrower
# bands sit closer to the aperture because each band is placed on the nearest
# clean run on its side. Degree 1 is offered because a freer spatial term is
# unconstrained across the unsampled aperture gap.
BACKGROUND_WIDTHS = (30, 25, 20, 15, 35, 40, 45, 50)
BACKGROUND_DEGREES = (2, 1, 3)
# Columns withheld from every trial fit, immediately outside the measured
# source extent. Four is enough to average down the noise without reaching
# back into the wings or out past the inner band edge.
HOLDOUT_WIDTH = 4
# How much a rival geometry must reduce the measured error beside the source,
# in units of the per-pixel noise, before the validated default is replaced.
# Measured over 2026-10-28: every well-behaved exposure gains at most 0.09,
# while the two genuinely misfitted ones gain 0.59 and 0.72.
REQUIRED_IMPROVEMENT_SIGMA = 0.3
# A correctly described background leaves normalized residuals of about one
# robust sigma over its sampled bands.
ADEQUATE_BAND_RESIDUAL_SIGMA = 1.15


def band_residual_sigma(data, ivar, sky, mask, valid_rows):
    """Robust normalized residual of a fitted model over its sampled bands.

    Uses the input selection rather than the post-rejection mask: an underfit
    gradient shows up exactly where sigma clipping would discard it. The median
    absolute deviate is scaled so that Gaussian noise returns about 1.0, which
    keeps cosmic rays from deciding the spatial degree.
    """
    data, ivar, sky = (np.asarray(v, float) for v in (data, ivar, sky))
    use = (np.asarray(mask, bool) & np.asarray(valid_rows, bool)[:, None]
           & np.isfinite(data) & np.isfinite(sky) & np.isfinite(ivar) & (ivar > 0))
    if use.sum() < 20:
        raise ValueError('too few sampled background pixels to assess a fit')
    normalized = (data[use]-sky[use])*np.sqrt(ivar[use])
    return float(1.4826*np.median(np.abs(normalized)))


def radial_background_profile(data, ivar, trace):
    """Row-median spatial profile about the trace, with its own noise scale.

    Returns the profile, the matching signed offset from the trace and the
    expected noise on that median. Measuring against the row median rather
    than a single row keeps a cosmic ray or a sky line out of the geometry
    decision.
    """
    data, ivar, trace = (np.asarray(v, float) for v in (data, ivar, trace))
    good = np.isfinite(data) & np.isfinite(ivar) & (ivar > 0)
    if good.sum() < 100:
        raise ValueError('too few usable pixels to measure a spatial profile')
    columns = np.arange(data.shape[1])
    offset = np.median(columns[None, :]-trace[:, None], axis=0)
    profile = np.array([np.median(data[good[:, c], c]) if good[:, c].any() else np.nan
                        for c in columns])
    rows = max(int(np.median(np.maximum(good.sum(0), 1))), 1)
    noise = float(np.median(1/np.sqrt(ivar[good]))/np.sqrt(rows))
    return profile, offset, noise


def _core_radius(profile, offset, noise):
    """Radius at which a Gaussian core of the measured width falls under the noise."""
    outer = (np.abs(offset) >= 16) & (np.abs(offset) <= 24) & np.isfinite(profile)
    if outer.sum() < 4:
        return 0
    base = np.interp(offset, [-20., 20.],
                     [float(np.mean(profile[outer & (offset < 0)])) if (outer & (offset < 0)).any() else 0.,
                      float(np.mean(profile[outer & (offset > 0)])) if (outer & (offset > 0)).any() else 0.])
    excess = profile-base
    near = np.abs(offset) <= 8
    if not near.any() or not np.isfinite(excess[near]).any():
        return 0
    peak = int(np.nanargmax(np.where(near, excess, -np.inf)))
    amplitude = excess[peak]
    if not np.isfinite(amplitude) or amplitude <= 0:
        return 0
    half = amplitude/2
    left = np.arange(max(peak-12, 0), peak+1)
    right = np.arange(peak, min(peak+13, len(offset)))
    if len(left) < 2 or len(right) < 2:
        return 0
    low = np.interp(half, excess[left], offset[left])
    high = np.interp(-half, -excess[right], offset[right])
    width = high-low
    if not np.isfinite(width) or width <= 0:
        return 0
    sigma = width/2.3548
    ratio = max(amplitude/max(0.5*noise, 1e-9), 1.0001)
    return int(np.ceil(sigma*np.sqrt(2*np.log(ratio))))


def _wing_radius(profile, offset, noise, *, window=16, tolerance=3.):
    """Radius beyond which the profile stops climbing above a LOCAL straight line.

    Sky, and a host galaxy, vary smoothly, so over a short span either is close
    to linear; a source wing is convex and rises above that line. Nothing global
    is assumed, so a galaxy gradient cannot be mistaken for source light nor the
    reverse. Returns None when no radius qualifies -- a bright halo that never
    dies inside the saved window.
    """
    found = 0
    for side in (-1, 1):
        same = np.sign(offset) == side
        limit = int(np.nanmax(np.abs(offset[same])))-window-2 if same.any() else 0
        radius = None
        for test in range(3, max(limit, 4)):
            fit = same & (np.abs(offset) >= test+2) & (np.abs(offset) <= test+2+window) & np.isfinite(profile)
            near = same & (np.abs(offset) >= test) & (np.abs(offset) < test+2) & np.isfinite(profile)
            if fit.sum() < 8 or not near.any():
                continue
            line = np.polyfit(np.abs(offset[fit]), profile[fit], 1)
            scatter = max(float(np.std(profile[fit]-np.polyval(line, np.abs(offset[fit])))), noise)
            if float(np.mean(profile[near]-np.polyval(line, np.abs(offset[near])))) < tolerance*scatter:
                radius = test
                break
        if radius is None:
            return None
        found = max(found, radius)
    return found


def clean_background_radius(data, ivar, trace):
    """Smallest offset beyond which source light no longer biases the background.

    Two independent estimates are combined and the larger taken. The first
    extrapolates a Gaussian core of the measured FWHM down to half the profile
    noise, which catches a bright compact source. The second asks where the
    profile stops rising above a local straight line, which catches the
    power-law halo of a standard star that a Gaussian badly underestimates.
    Returns None when no such radius exists inside the saved window.
    """
    profile, offset, noise = radial_background_profile(data, ivar, trace)
    wing = _wing_radius(profile, offset, noise)
    if wing is None:
        return None, noise
    return max(int(wing), int(_core_radius(profile, offset, noise))), noise


def holdout_residual(data, ivar, sky, holdout, valid_rows):
    """Mean model error over source-free columns that the fit never sampled.

    This is the only direct measurement available of whether a background model
    extrapolates correctly towards the source. It is signed: negative means the
    model sits above the data and the background is being over-subtracted.
    """
    data, ivar, sky = (np.asarray(v, float) for v in (data, ivar, sky))
    use = (np.asarray(holdout, bool) & np.asarray(valid_rows, bool)[:, None]
           & np.isfinite(data) & np.isfinite(sky) & np.isfinite(ivar) & (ivar > 0))
    if use.sum() < 50:
        return None, None
    deviation = data[use]-sky[use]
    return float(np.mean(deviation)), float(np.std(deviation)/np.sqrt(deviation.size))


def select_bounded_background(data, ivar, tilts, trace, *, radius, fwhm, neighbours=(),
                              allowed_mask=None, max_gap_pixels=60, sigma_clip=4.,
                              widths=BACKGROUND_WIDTHS, degrees=BACKGROUND_DEGREES,
                              holdout_width=HOLDOUT_WIDTH,
                              improvement_sigma=REQUIRED_IMPROVEMENT_SIGMA, maxiter=500):
    """Choose the band width and spatial degree by measuring, not by assumption.

    Two bands constrain the background either side of the source and nothing
    measures it underneath, so the band geometry decides how the model behaves
    exactly where it is used and cannot be checked. It can be checked here,
    because between the source wings and the inner band edge there is a ring of
    columns carrying background but no source light. Those columns are withheld
    from every trial fit and used to score it: the mean of data minus model over
    them is a direct measurement of the error the model is making beside the
    source, with a negative value meaning over-subtraction.

    Where that ring sits is measured per exposure, not assumed. A faint
    supernova on a bright host is clean a few pixels out; the power-law halo of
    a standard star is not clean anywhere inside the slit, and that exposure
    keeps its default geometry rather than letting the selector chase the halo.

    The previously validated default is kept unless some other geometry actually
    reduces the measured error by a worthwhile margin, so a well-behaved
    exposure is never churned. When nothing achieves that, the default is kept
    and the residual error is recorded rather than hidden: an asymmetric host,
    or a neighbour inside the bands, is a reduction for a human to look at, not
    something a width and a degree can repair.

    Every trial stays bounded and two-sided; no one-sided, full-slit or constant
    fallback is reachable. This selects a background description, not a host
    model, and a small score is evidence about the ring that was measured, not
    proof that the background under the aperture is known.
    """
    if not widths or not degrees:
        raise ValueError('background selection needs at least one width and degree')
    if type(holdout_width) is not int or holdout_width < 2:
        raise ValueError('the withheld ring must be at least two pixels wide')
    if not np.isfinite(improvement_sigma) or improvement_sigma <= 0:
        raise ValueError('the required improvement must be positive')
    # Detector defects remove fitting samples, not geometric band locations.
    # Keep source/neighbour guards in nearby_mask; retain holes in the actual
    # fit mask and let fit_background enforce two-sided support per row.
    if allowed_mask is not None:
        if np.shape(allowed_mask) != np.shape(data):
            raise ValueError('allowed background mask has the wrong shape')
        ivar = np.where(np.asarray(allowed_mask, bool), ivar, 0.)
    finite = np.isfinite(ivar) & (ivar > 0)
    if not finite.any():
        raise BackgroundSupportError('No pixel carries a usable inverse variance.')
    pixel_noise = float(np.median(1/np.sqrt(np.asarray(ivar, float)[finite])))
    clean, profile_noise = clean_background_radius(data, ivar, trace)
    columns = np.arange(np.shape(data)[1])
    offset = columns[None, :]-np.asarray(trace, float)[:, None]
    trials, fits, masks = [], {}, {}

    def place(width, inner):
        key = (width, inner)
        if key not in masks:
            try:
                masks[key] = nearby_mask(data.shape, trace, radius=float(inner), fwhm=fwhm,
                    width_pixels=int(width), neighbours=neighbours,
                    max_gap_pixels=max_gap_pixels)[0]
                # Distance/source-overlap failures are geometry failures too.
                # Classify before fitting so native window expansion can try
                # more calibrated room without relaxing any support gate.
                validate_bounded_mask(masks[key], trace, radius)
            except BackgroundSupportError as exc:
                masks[key] = None
                trials.append(dict(width_spatial_pixels=int(width), spatial_degree=None,
                                   status='geometry_unavailable', reason=str(exc)))
        return masks[key]

    def attempt(width, degree, inner, holdout, label):
        key = (width, degree, inner, label)
        if key in fits:
            return fits[key]
        fits[key] = None
        mask = place(width, inner)
        if mask is None:
            return None
        sampled = mask & ~holdout if holdout is not None else mask
        try:
            fit = fit_background(data, ivar, tilts, sampled, trace, radius=radius,
                                 degree=int(degree), sigma_clip=sigma_clip, maxiter=maxiter)
            band = band_residual_sigma(data, ivar, fit['sky'], sampled, fit['valid_rows'])
        except (BackgroundSupportError, ValueError) as exc:
            trials.append(dict(width_spatial_pixels=int(width), spatial_degree=int(degree),
                               role=label, status='fit_failed', reason=str(exc)))
            return None
        score, error = (holdout_residual(data, ivar, fit['sky'], holdout, fit['valid_rows'])
                        if holdout is not None else (None, None))
        record = dict(width_spatial_pixels=int(width), spatial_degree=int(degree), role=label,
                      status='fitted', band_residual_sigma=band,
                      inner_edge_spatial_pixels=float(inner),
                      holdout_residual_counts=score,
                      holdout_residual_error_counts=error,
                      holdout_residual_pixel_sigma=(None if score is None
                                                    else float(score/pixel_noise)),
                      reduced_chi2=fit['summary']['reduced_chi2'],
                      minimum_side_samples=fit['summary']['minimum_side_samples'])
        trials.append(record)
        fits[key] = (mask, fit, record)
        return fits[key]

    # A 15--29 pixel clean band beside the source beats a distant 30-pixel
    # band. Rank offered widths by their worst inner-edge distance first.
    proximity = []
    for width in widths:
        candidate = place(width, radius)
        if candidate is not None:
            distances = np.where(candidate, abs(offset), np.inf)
            left = np.min(np.where(offset < 0, distances, np.inf), axis=1)
            right = np.min(np.where(offset > 0, distances, np.inf), axis=1)
            proximity.append((float(np.max([left, right])), abs(width-30), width))
    if proximity:
        preferred = min(proximity)[2]
        widths = (preferred,) + tuple(w for w in widths if w != preferred)
    default_width, default_degree = widths[0], degrees[0]
    # All scored models use a common, source-free holdout outside the measured
    # aperture. Their bands begin beyond it; the default width/degree competes
    # on the same support as the alternatives.
    aperture = int(np.ceil(radius))
    measurable = clean is not None
    inner = None
    holdout = None
    if measurable:
        clean = max(int(clean), int(np.ceil(radius+.5)))
        # A band may never start inside the aperture, nor inside the ring.
        inner = max(aperture, int(clean)+int(holdout_width))
        holdout = np.broadcast_to((np.abs(offset) >= clean)
                                  & (np.abs(offset) < clean+holdout_width), np.shape(data)).copy()
        for other, fw, other_radius in neighbours:
            guard = max(1.5*fw, other_radius+.5)
            holdout &= abs(columns[None, :]-np.asarray(other)[:, None]) > guard
        if allowed_mask is not None:
            holdout &= allowed_mask
        holdout &= finite
        if any(np.any(np.sum(holdout & (offset*side > 0), axis=1) < 2) for side in (-1, 1)):
            measurable = False
        # Leave room for the narrowest permitted band outside the withheld ring.
        if inner+min(int(w) for w in widths) > int(np.max(np.abs(offset))):
            measurable = False

    default_inner = None
    if measurable:
        placed = place(default_width, inner)
        if placed is None:
            measurable = False
        else:
            reach = np.abs(offset)[placed]
            default_inner = int(np.floor(reach.min())) if reach.size else 0
            if int(clean)+int(holdout_width) > default_inner:
                measurable = False

    scored, chosen, reason, unresolved = [], None, None, False
    if measurable:
        attempt(default_width, default_degree, inner, holdout, 'default')
        for width in widths:
            for degree in degrees:
                attempt(width, degree, inner, holdout, 'trial')
        scored = [v for v in fits.values()
                  if v is not None and v[2]['holdout_residual_counts'] is not None]
    if not scored:
        # Nothing could be measured, so keep the previously validated behaviour:
        # the widest offered band that can actually be placed and fitted, at the
        # default spatial degree. A crowded slit has always narrowed this way.
        selected = next((v for v in (attempt(width, default_degree, radius, None, 'default')
                                     for width in widths) if v is not None), None)
        if selected is None:
            fitted_any = any(t['status'] != 'geometry_unavailable' for t in trials)
            raise BackgroundSupportError(
                ('No bounded two-sided background converged at the default spatial degree.'
                 if fitted_any else 'No bounded two-sided background could be placed: '
                 'source/neighbour guards and calibrated support leave insufficient room; '
                 'the sky solver was not called.'), trials=trials)
        chosen, unresolved = selected, True
        reason = ('default retained: no ring of columns inside this slit carries background '
                  'without source light, so no geometry could be measured'
                  if clean is None else
                  'default retained: no two-sided holdout outside the measured aperture '
                  'and neighbour guards fits beside the finite bands; width and degree '
                  'could not be compared on independent background pixels')
    else:
        best = min(scored, key=lambda v: abs(v[2]['holdout_residual_counts']))
        default = next((v for v in scored if v[2]['role'] == 'default'), None)
        if default is None:
            chosen = best
            reason = ('the default geometry could not be fitted; the smallest measured error '
                      'beside the source was taken')
        else:
            gain = abs(default[2]['holdout_residual_counts'])-abs(best[2]['holdout_residual_counts'])
            if gain <= improvement_sigma*pixel_noise:
                chosen = default
                reason = ('default retained: no offered width or spatial degree reduced the '
                          'measured error beside the source by a worthwhile margin')
                unresolved = abs(default[2]['holdout_residual_counts']) > improvement_sigma*pixel_noise
            else:
                # Among geometries that measure as well as the best, keep the one
                # closest to the validated default rather than the extreme of the grid.
                tolerance = max(3*(best[2]['holdout_residual_error_counts'] or 0.), .3)
                near = [v for v in scored
                        if abs(v[2]['holdout_residual_counts'])
                        <= abs(best[2]['holdout_residual_counts'])+tolerance]
                chosen = min(near, key=lambda v: (
                    abs(v[2]['width_spatial_pixels']-default_width),
                    abs(v[2]['spatial_degree']-default_degree),
                    v[2]['inner_edge_spatial_pixels'],
                    abs(v[2]['holdout_residual_counts'])))
                reason = ('measured: this width and spatial degree reduced the model error over '
                          'the withheld source-free ring beside the source')
    mask, fit, record = chosen
    # The withheld ring is ordinary background; only the choice needed it held
    # out. Refit the selected geometry using every column it is entitled to.
    final = attempt(record['width_spatial_pixels'], record['spatial_degree'],
                    record['inner_edge_spatial_pixels'], None, 'final')
    if final is not None:
        mask, fit = final[0], final[1]
    selection = dict(policy='bounded two-sided background chosen by withheld-ring measurement',
        selected_width_spatial_pixels=record['width_spatial_pixels'],
        selected_spatial_degree=record['spatial_degree'],
        selected_inner_edge_spatial_pixels=record['inner_edge_spatial_pixels'],
        clean_background_radius_pixels=(None if clean is None else int(clean)),
        default_band_inner_edge_spatial_pixels=default_inner,
        holdout_ring_pixels=(None if not measurable else [int(clean), int(clean)+holdout_width-1]),
        holdout_residual_counts=record['holdout_residual_counts'],
        holdout_residual_error_counts=record['holdout_residual_error_counts'],
        holdout_residual_pixel_sigma=record['holdout_residual_pixel_sigma'],
        band_residual_sigma=record['band_residual_sigma'],
        pixel_noise_counts=pixel_noise, profile_noise_counts=profile_noise,
        required_improvement_pixel_sigma=float(improvement_sigma),
        geometry_unresolved=bool(unresolved),
        reason=reason,
        offered_widths=[int(w) for w in widths], offered_degrees=[int(d) for d in degrees],
        residual_statistic='1.4826 x median |(data-sky) sqrt(ivar)| over the input bands and supported rows',
        holdout_statistic='mean of data minus fitted model over withheld source-free columns '
                          'between the source wings and the inner band edge; negative means the '
                          'background is being over-subtracted beside the source',
        holdout_scope='measures the model error over the ring that was withheld, which is the '
                      'closest measurable proxy for the unsampled aperture gap; it is not proof '
                      'that the background under the aperture is known',
        clean_radius_statistic='larger of a Gaussian core extrapolated to half the profile noise '
                               'and the radius where the profile stops rising above a local line',
        escalated=(record['width_spatial_pixels'] != default_width
                   or record['spatial_degree'] != default_degree),
        identity_inferred=False, trials=trials)
    return mask, fit, selection


def fit_background(data, ivar, tilts, mask, trace, *, radius, degree=1,
                   sigma_clip=4., bsp=1., maxiter=500):
    """Use PypeIt's B-spline × Legendre fit, retaining its convergence evidence.

    ``degree=2`` supplies three Legendre terms. Unlike global_skysub in 1.14,
    failure never retries a spatial constant. We deliberately omit its positive
    log-sky prefit: that spatially constant prefit can clip a real gradient.
    Returned valid rows have >=5 accepted pixels *on each side*. Callers must
    mask invalid rows before extraction. Model covariance is not propagated by
    PypeIt's fixed-sky optimal extraction and must not be claimed as included.
    """
    import pypeit
    from pypeit.core import basis, fitting
    if pypeit.__version__ != '3.0.0':
        raise ValueError('local background adapter requires PypeIt 3.0.0')
    data, ivar, tilts = (np.asarray(v, float) for v in (data, ivar, tilts))
    mask, trace = np.asarray(mask, bool), np.asarray(trace, float)
    if (data.ndim != 2 or data.shape != ivar.shape or data.shape != tilts.shape
            or mask.shape != data.shape or trace.shape != (len(data),)
            or not np.isfinite(tilts).all() or not np.isfinite(trace).all()):
        raise ValueError('background arrays and trace geometry disagree')
    if type(degree) is not int or not 1 <= degree <= 3:
        raise ValueError('PypeIt local spatial degree must be 1, 2 or 3')
    if not np.isfinite([radius, sigma_clip, bsp]).all() or radius <= 0 or sigma_clip < 3 or bsp <= 0:
        raise ValueError('invalid background aperture, clipping or B-spline spacing')
    if type(maxiter) is not int or maxiter < 1:
        raise ValueError('background iteration limit must be a positive integer')
    good = mask & np.isfinite(data) & np.isfinite(ivar) & (ivar > 0)
    validate_bounded_mask(mask, trace, radius)
    counts = side_counts(good, trace, radius)
    initial_valid=np.all(counts>=5,axis=0)
    if initial_valid.mean()<.9:
        raise BackgroundSupportError('More than 10% of rows lack five usable background pixels on each side before fitting.',
                                     unsupported_rows=np.flatnonzero(~initial_valid).tolist())
    # Losing a few rows to detector masks must not reject a whole exposure.
    # Remove both sides of each unsupported row and never deliver it as valid.
    good &= initial_valid[:,None]
    pixels = (tilts*(len(data)-1)).ravel()
    order = np.argsort(pixels)
    spatial = np.broadcast_to(np.linspace(-1., 1., data.shape[1]), data.shape).ravel()
    # Fit all coordinates with zero weight outside the bounded sampling mask.
    # This also returns the interpolated model under the source aperture.
    _, accepted, model, chi2, status = fitting.iterative_bspline_fit(
        pixels[order], np.where(np.isfinite(data), data, 0).ravel()[order],
        np.where(good, ivar, 0).ravel()[order],
        basis=basis.flegendre(spatial[order], degree+1), gpm=good.ravel()[order],
        nord=4, upper=sigma_clip, lower=sigma_clip, maxiter=maxiter,
        kwargs_knots={'spacing':bsp}, kwargs_reject={'groupbadpix':False, 'maxrej':10})
    if status != 0 or not np.isfinite(model).all():
        raise BackgroundSupportError(f'PypeIt spatial degree {degree} did not converge (exit status {status}); no lower-order fallback.',
                                     spatial_degree=degree, pypeit_exit_status=int(status))
    sky, used = np.empty_like(pixels), np.empty_like(pixels, dtype=bool)
    sky[order], used[order] = model, accepted
    sky, used = sky.reshape(data.shape), used.reshape(data.shape) & good
    final_counts = side_counts(used, trace, radius)
    valid = initial_valid & np.all(final_counts >= 5, axis=0)
    if valid.mean() < .9:
        raise BackgroundSupportError('More than 10% of rows lost two-sided background support after rejection.',
                                     unsupported_rows=np.flatnonzero(~valid).tolist(), spatial_degree=degree)
    summary = dict(engine='pypeit.core.fitting.iterative_bspline_fit', pypeit_version=pypeit.__version__,
                   spatial_degree=degree, legendre_terms=degree+1, pypeit_exit_status=int(status),
                   sigma_clip=float(sigma_clip), bspline_spacing=float(bsp),
                   maximum_iterations=maxiter,
                   pre_fit_invalid_rows=np.flatnonzero(~initial_valid).tolist(),
                   reduced_chi2=float(chi2), invalid_rows=np.flatnonzero(~valid).tolist(),
                   accepted_fraction=float(used.sum()/good.sum()),
                   minimum_side_samples=int(final_counts[:, valid].min()),
                   no_local_sky=True, constant_fallback=False, positive_log_prefit=False,
                   uncertainty_scope='PypeIt fixed-sky extraction; background-model covariance not propagated')
    return dict(sky=sky, input_mask=good, accepted_mask=used, valid_rows=valid,
                side_counts=final_counts, summary=summary)
