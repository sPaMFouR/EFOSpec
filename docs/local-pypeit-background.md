# Bounded PypeIt backgrounds

New nightly reductions and the desktop editor use bounded, two-sided PypeIt
backgrounds. The confirmed unit is spatial detector pixels: prefer 30 pixels
per side, narrowing to 15 when crowded. GUI widths can be set from 15 to 50.
The default spatial model is quadratic; degrees 1–3 are selectable. Existing
`automatic` recipes retain the old algorithm through explicit legacy replay.
New automatic recipes use `bounded`; anchored recipes use `pypeit_local`.

- Background is spatial: detector pixels or arcseconds. Å describes the
  wavelength direction. Both controls must have separate units in the GUI.
- One finite band on each side, outside the aperture and source wings. The
  nearest clean band moves past saved neighbours when necessary; its width
  never expands to the detector edge. Insufficient room raises a recorded
  `background_support_required` error.
- Degrees 1, 2 and 3 select a slope, quadratic or cubic spatial model. A smooth,
  uneven host background can motivate curvature; nearby peaks alone do not
  establish that it is a galaxy. Compare residuals before accepting an order.
- Nightly reductions choose the width and spatial degree automatically; see
  [automatic selection](#automatic-width-and-degree-selection).
- The actual fit uses PypeIt 1.14.0 `bspline_profile`: spectral B-splines times
  `degree+1` spatial Legendre terms, followed by PypeIt optimal extraction with
  `no_local_sky=True`. This mode is distinct from the existing NumPy polynomial
  boxcar experiment.
- No silent constant-order fallback or positive log-sky prefit. The latter's
  constant spatial model can reject a real gradient. Non-convergence stops the
  trial and records the requested degree and PypeIt exit code.
- At least five accepted pixels on each side after rejection. Unsupported rows
  are masked; losing more than 10% stops extraction. A fitted source trace
  crossing the selected background also stops extraction.
- Saved NPZ files retain the input and accepted masks, side counts, valid rows
  and spatial degree. Product QA records convergence, settings and limitations.
- Fixed-sky PypeIt uncertainties do not propagate fitted background covariance.
  This is not validated SN/host deblending. Saved candidate guards are
  conservative: old records lack neighbour FWHM measurements.

## Automatic width and degree selection

`efospec.local_sky.select_bounded_background` chooses the band width and the
spatial degree from the exposure instead of fixing 30 pixels and a quadratic.

Two bands constrain the background on either side of the source; nothing
measures it underneath. A quadratic fitted to two wide, distant bands can arch
across that gap and over-subtract beside the source while still fitting its own
sampled pixels perfectly. A statistic computed over the fitted bands cannot see
this, because the error lives exactly where no pixel was used. The selection
therefore measures the model where it is not fitted:

- **Clean background radius**: the larger of a Gaussian core extrapolated out to
  half the profile noise and the radius at which the row-median spatial profile
  stops rising above a local straight line. The second test catches power-law
  halos that a Gaussian core underestimates, as a bright standard star has.
- **Withheld ring**: `HOLDOUT_WIDTH` columns immediately outside that radius, on
  both sides, excluded from every trial fit. Each trial is scored by the mean of
  `data - sky` there. Negative means the background is being over-subtracted
  beside the source. The chosen geometry is then refitted with the ring included.
- **Band residual**: `1.4826 x median |(data-sky) sqrt(ivar)|` over the input
  bands and supported rows, as a secondary check that the model describes what it
  did sample. The input selection is used, not the post-rejection mask, so an
  underfit gradient shows up where sigma clipping would have discarded it.

The default 30-pixel band with degree 2 is fitted and scored first. It is
replaced only when some other offered geometry reduces the measured ring error
by more than `REQUIRED_IMPROVEMENT_SIGMA` (0.3) times the per-pixel noise;
among geometries within noise of the best, the one closest to the default is
preferred. Widths 30/25/20/15/35/40/45/50 and degrees 2/1/3 are offered. Every
trial stays bounded and two-sided; no one-sided, full-slit or constant fallback
is reachable.

`background_fit.selection` records the chosen width, degree and inner edge, the
clean radius, the ring, every statistic and threshold, the reason and all
trials.

### When no ring can be measured

The ring must lie in the gap the default band leaves: the clean radius plus the
ring width must fall inside the default band's inner edge. When it does not, no
geometry can be compared against anything, `geometry_unresolved` is recorded and
the default is retained. If the default width cannot be placed at all, the width
ladder is walked at the default degree until one fits.

**This is currently the case for every exposure in the nightly path.** In the
28 October reduction the bands are placed much closer to the source than the
clean radius: for L745a Gr13 slit1.0 the clean radius is 14 pixels while the
default band's inner edge is 9 and the fitted band begins at 6. The source is
still bright where the background is being measured, so no clean ring exists and
the selector returns the default unmeasured. The band placement — driven by
`aperture_padding_pixels` in the recipe — has to be reconciled with the measured
clean radius before this selection does anything in production.

Measured off-line on the 2026-10-28 exposures with bands placed outside the
clean radius, the ring is measurable for all 14 and 12 keep width 30 degree 2
with ring errors of -0.40 to +0.31 per-pixel sigma. AT2026agvs Gr13 slit1.5
moves to width 15 degree 1 (-0.02 sigma, from a default that over-subtracts),
and vma2 Gr13 slit1.5 to degree 3.

A small ring error is not evidence that the background under the aperture is
known, and this selection chooses a background description, not a host model.

## One date, explicit targets

Inspect the saved geometry and create house-style QC and candidate recipes:

```sh
python scripts/test_local_background.py 2026-10-28 --targets 2026ayt 2026agwm
python scripts/build_epoch_fieldnotes.py 2026-10-28
```

This checks 30, 40 and 50 **spatial pixels** as alternatives without selecting
a default or extracting new spectra. Outputs live under
`runs/2026-10-28/QC/local-background-<revision>/`.

An explicit `--width-pixels N` runs linear
and quadratic PypeIt candidate trials in fresh revisions; `--degrees 1 2 3`
adds a cubic comparison. Source identities remain unverified, including every
SN2026ayt candidate. No C3 choice from SN2026aebl transfers to these exposures.

The completed 2026-10-28 trial invocation was:

```sh
python scripts/test_local_background.py 2026-10-28 --targets 2026ayt 2026agwm --width-pixels 30 --minimum-width-pixels 15 --padding-pixels 50
```

Padding is an explicit archive-bridge operation, recorded in the recipe and
hashed window geometry. It must reproduce every original data and wavelength
pixel exactly. It cannot expand already processed native-calibration data.
The new pixels still need QC review for contaminants outside the saved traces.

## Evidence and limits

On 2026-10-28, the saved 141-pixel windows allow 30-pixel bands for all three
SN2026ayt candidate hypotheses. The 40/50-pixel alternatives fail. For
AT2026agwm, all 30/40/50-pixel alternatives fail when the other saved candidates
are excluded. Expanding this archived-calibration window by 50 pixels per side
allowed two 30-pixel bands and successful linear/quadratic extractions.
SN2026ayt C2 also extracted with both degrees. C1 and C3 failed the fitted-trace
background-overlap gate. C2's numerical success does not verify its identity.
All products remain experimental, in new dated QC revisions; the accepted
epoch spectra and current-revision pointer were not replaced.

Synthetic tests cover curved traces, neighbouring sources, cosmic-ray
rejection, quadratic-background recovery, PypeIt optimal flux recovery,
non-convergence, loss of one-sided support and hash-bound revision replay.

References: [PypeIt 1.14 sky subtraction](https://pypeit.readthedocs.io/en/1.14.0/skysub.html),
and the installed 1.14.0 `pypeit/core/skysub.py` and `core/fitting.py` source.
The adapter calls the underlying solver directly to retain its rejection mask
and convergence status, which `global_skysub` does not return.

## Detector masks and smooth band geometry

Automatic finite bands follow the source trace and avoid neighbouring sources.
Detector validity is applied through zero inverse variance after choosing those
bands, rather than searching for a replacement contiguous band around each bad
pixel. Saved fitting masks therefore contain holes without displaced band edges.
Existing two-sided support and rejection checks remain active; unsupported rows
cannot become valid by moving their bands elsewhere. Explicit legacy recipe
exclusions retain their existing geometry semantics.
