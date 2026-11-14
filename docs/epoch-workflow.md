# One observing epoch at a time

Always name the observing date. One invocation processes one date only. For the next epoch, change the date; there is no automatic multi-night loop.

The user-approved **2026-10-29** page is the final daily-log template for both October nights and future epochs. Keep its approved section structure and use the masthead **ePESSTO+ / EFOSpec**. Calibration checks appear in AT, SN, then other-source order; failed or unvalidated checks are red. The optional `$webpage-spectralatlas` coding-agent skill is installed in the user's agent configuration, outside this repository. The pipeline does not require a skill installation to generate its report.

Activate the existing EFOSpec Python environment and run from the EFOSpec checkout:

```sh
efospec night 2026-10-28
# Equivalent checkout script:
python scripts/reduce_epoch.py 2026-10-28
```

To inspect each stage separately:

```sh
python scripts/01_log_epoch.py 2026-10-28
python scripts/02_reduce_epoch.py 2026-10-28
python scripts/03_report_epoch.py 2026-10-28
```

- CLI equivalents: `efospec night DATE --stage log`, `reduce`, or `report`.
- Numbered scripts use the checkout as project root; the CLI uses the current directory.
- `--project-root PATH` selects another project containing `data/` and `runs/`.
- The engine and reference files always come from the code checkout.
- Downloads remain explicit: `efospec download DATE --with-calibrations`.

## Stages

| Stage | Work |
| --- | --- |
| Observing log | Inventory `data/DATE/`; record FITS roles, dates and hashes in CSV/JSON |
| Reduction | Verify unchanged inputs; build native raw calibrations and same-night standard responses |
| Final report | Build the nightly log, plots and candidate QC; retain partial results |

- Observing date spans noon UTC to noon the next day.
- Science/standard exposures from another night stop the workflow.
- Calibration/acquisition frames retain their actual dates and roles.
- Complete worker output goes into dated logs.

## Dated outputs

Open **`runs/2026-10-28/index.html`**.

Use this main nightly HTML for status and issues. Do not generate separate
diagnostic HTML pages unless explicitly requested. Show diagnostic plots directly
in the conversation for aperture/background review; retain their data and
provenance in dated QC folders. Later retry notes must be hash-bound to the
selected parent reduction and raw exposure before appearing in the main report.

The opening page embeds lightweight WebP figure previews, so figures also work
in HTML viewers without access to neighbouring files. Full-resolution figures
and scientific products remain unchanged. Downloads use the sibling
`report-assets/` folder: open the report in Chrome or Safari with that folder
alongside to download files. For a single-file offline copy including downloads,
use `index-portable.html`, also linked from the footer. It can exceed 100 MB;
open it in a browser rather than a chat file preview.
Both versions are written by the same renderer without re-extracting spectra.

For links opened inside T3, use the returned `viewer_path` (`index-view-<hash>.html`).
This identical copy stays beside `report-assets/` and is never overwritten.
Rebuilding `index.html` can invalidate T3's cached file-identity URL; a new report
gets a new viewer filename, while existing viewer links keep working.

Spatial-profile figures show only the data before subtraction and the saved
PypeIt sky model in three vertically stacked wavelength bins. Successful and
failed-review profiles share their renderer, colours, fonts and layout. Failed
reviews show only validated candidates or an explicitly selected aperture.
Visible apertures are numbered consecutively from left to right: a single
retained aperture is always C1. The original candidate identifiers, selection
and evidence remain unchanged in provenance, with an explicit display-label map.
An absent sky fit is explicitly labelled and
never replaced with a zero-valued model curve. The curves use identical native
pixels, with linear electron units. The spatial degree comes from the saved
fit; no new polynomial is fitted to the sums. Model curves stop at the common
background-band envelope rather than displaying unsupported extrapolation. The
figures carry no explanatory paragraphs or per-panel scale annotations: the
summation rule, display limits and excluded columns stay in the saved
diagnostics and the report caption.
The y-axis follows each panel's local profile range, with 15% headroom above
the retained peak. Positive backgrounds keep a positive lower limit; a
background median within three robust scatter estimates of zero gets negative
space. Real negative values remain visible. Legends use `Observed Data` and
`Background Sampling`; the panel heading is `Before Subtraction: Observed data
and fitted background` when a saved model is available.
Use `scripts/preview_spectrum_style.py --run RUN --profiles-only` to refresh
profiles without recomputing spectra. For a separate comparison plot directory,
add `--plot-directory DIRECTORY` and `--target NAME` as needed. Previous profile
figures and their record are archived; refreshed figure hashes remain verifiable.

The per-source 2D card shows both sky-subtraction stages on one wavelength axis,
above and below, so the fitted sky can be judged against its residual.

Recovery-QC spatial profiles exclude the outermost three unprotected columns on
each side from the limit estimate. Slit-boundary columns have no fitted
background and run an order of magnitude above the sources, and the spike rule
deliberately never masks a run touching the array end. Those samples remain
plotted; only the viewport changes. Protected candidate and aperture columns are
never dropped, and `profile-display.json` records the margin.

Acquisition image displays use `RdBu`; spectral 2D and raw spectral detector
images use grayscale (`gray_r`). Redshift-evidence cut-outs use `RdBu_r` about
the clipped median of each panel, so line emission above the local background
reads red and absorption or over-subtraction reads blue. All retain the arcsinh
stretch. The 2D spectral panels take their noise width from the sky-sampling
bands as before, but their upper limit from the 99.5th data percentile: a
sky-subtracted trace sits hundreds of background sigma above zero, and the
earlier five-sigma limit saturated it so only noise texture stayed visible. Three-sigma
clipping (median/MAD, up to five iterations) estimates the
colour scale only: all finite image samples remain in the display and numerical
products are untouched. Residual scales are symmetric about zero; raw/acquisition
scales about the clipped median, with limits at five robust sigma and a one-sigma
arcsinh width. Acquisition full/zoom views share one normalization. A yellow
outline on each full image marks the exact displayed bounds of the right-hand zoom.
Implementation references: [Astropy sigma clipping](https://docs.astropy.org/en/stable/api/astropy.stats.sigma_clip.html)
and [Matplotlib invertible colour normalization](https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.FuncNorm.html).
Per-date `report-assets/sky-model-arcsinh-verification.json` records unchanged
inputs; prior plots are retained under each run's `archive/before-sky-model-arcsinh/`.

```text
runs/2026-10-28/nightly/<processing-UTC>-<unique-id>/
  observing-log.csv
  observing-log.json
  night.json
  index.html
  logs/2026-10-28-workflow.log
  logs/2026-10-28-reduction.log
  logs/2026-10-28-report.log
  reduction/                 # run.json, spectra, plots, nightly report
  QC/<raw-exposure>/         # candidate images and recipes
```

- Every log/all-stages invocation creates a new processing revision.
- `current.json` selects the revision used by separate reduction/report steps.
- Repeating extraction in an existing revision is rejected; start with `log` again.
- The report can be regenerated without re-extraction.
- Raw files and accepted products remain unchanged.

## Independent redshift candidates

After reducing one night, run the automated 2D emission-line search and recorded
NED cross-match into a fresh QC revision:

```sh
python scripts/measure_host_redshifts.py 2026-10-28 \
  --output runs/2026-10-28/QC/host-redshifts-NEW --update-report
```

The current accepted and rescue products are selected and hash-checked. The
observing table keeps NED and spectral redshifts separate, with catalogue methods,
candidate associations and provisional wavelength flags. No NED prior enters
the spectral fit. Missing or ambiguous line systems remain without an estimate;
this emission-line workflow does not cover absorption-only hosts. Measurements,
plots and input hashes remain in the revision; `QC/redshift-current.json` selects
the current results. An incompatible parent revision prevents activation.
See [the delivered method and validation record](rescue-2026-11-09.md).

## Portable nightly log

- The date's `index.html` contains observations, spectra, issues, source QC,
  redshift evidence, background QC, calibrations and archive summaries.
- **Redshift evidence** is its own tab when measurements are current: it has its
  own inputs, figures and provenance and reads as a separate measurement from
  source QC. Spectral redshifts are quoted as `z ± error`, where the error is
  the larger of the inverse-variance fit error and the standard error of the
  group-to-group scatter over every consistent narrow emission-line group. It
  excludes wavelength, background and resampling systematics.
- The final section is **Reduction Archive**: retained reductions, folder locations,
  then **Download logs** (log files only). There is no separate Run record section
  or archive-calibration/status/revision caption beneath the date. Manifests and
  provenance remain saved on disk.
- `index.html` embeds its figure previews; downloads load from `report-assets/` offline.
  `index-portable.html` embeds them all and can be copied as one self-contained file.
- Figure enlargement and navigation stay within the page; existing PDF/JPG exports are preserved.
- `nightly/index.html` remains an entry to revision-specific records.
- Rebuild from saved results without re-extraction or new plots:

```sh
python -m efospec.epoch_report 2026-10-28
```

For explicitly selected targets, [bounded PypeIt background trials](local-pypeit-background.md)
save their geometry checks, candidate recipes and new extraction revisions inside
the date's `QC/` folder. These trials do not replace `nightly/current.json`.

[Acquisition catalogue QC](acquisition-catalogues.md) optionally solves acquisition
images locally and adds NED/SDSS/AllWISE overlays to Source QC. Its dated command
is separate from reduction; rebuilding the report uses only saved results.

## Retained results

| Folder under `runs/DATE/` | Contents |
| --- | --- |
| `nightly/` | Current workflow and processing revisions |
| `QC/` | Retained manual source/background experiments |
| `calibrations/` | Retained masters and group selections |
| `archive/` | Earlier reductions, diagnostics, inventories and logs |
| `report-assets/` | Byte-identical figure/download payloads for the small opening page |

- Use the date's `index.html`/`README.md` when available, or `nightly/index.html`.
- Write new results inside dated revisions.
- Older compatibility links are read-only destinations for existing recipes/records.

## QC cleanup after the PypeIt 3 reruns

On 11 November 2026, the authorized cleanup removed 45 superseded generated QC
folders (4,994,157,498 bytes) across the two October nights. This includes early
daily-reduction trials, obsolete sky/background probes, redundant report backups
and the old single-engine overview plots. It does not remove pipeline validation
code or weaken rejection thresholds.

| Observing date | Removed folders | Removed files | Bytes removed | Retained files hash-verified |
| --- | ---: | ---: | ---: | ---: |
| 2026-10-28 | 28 | 1,517 | 984,363,405 | 3,652 |
| 2026-10-29 | 17 | 1,634 | 4,009,794,093 | 2,202 |

Each date has `archive/qc-cleanup-20261111/manifest.json`, containing the exact
removed-file inventory, sizes, historical documentation references, retained
dependencies and protected SHA-256 values. The scripts and original comparison
baseline hashes accompany it. Historical debugging notes retain their original
paths as provenance; a removed diagnostic is no longer an available output.

All nightly and rescue revisions, selected-source recipes, active catalogue and
redshift dependencies, current filtered/red count diagnostics, the PypeIt migration
benchmark and final `QC/pypeit3-overlay-20261110-v1/` comparisons remain. In
particular, catalogue parent records and backgrounds referenced by retained rescue
products are dependencies, even when their folder names look old. The current
reports, spectra, comparison baselines and retained QC files passed before/after
hash checks. Raw data and the preserved legacy checkout were outside deletion scope.

## Status and limits

- [Integrated reduction checks](reduction-checks.md) now save sky-line fits,
  PypeIt correlation comparisons, calibration readiness and native 2D diagnostics.
  `--native-calibrations RUN` adds contemporary raw-master checks to the same
  dated invocation, with their own status and report section.

| Exit | Meaning |
| --- | --- |
| `0` | Requested stage complete; all three stages complete for a successful full run |
| `1` | Partial reduction or report/QC failures |
| `2` | Operational/input error |

- Reporting continues after partial reduction. For separate steps, run `report` even after exit `1`.
- `night.json` records stage outcomes; missing/incomplete workers cannot become successes.
- A per-date lock prevents concurrent operations. Inspect its PID before removing a stale `.running` lock.
- New revisions default to `--flux-calibration same-night`: native raw bias/flat/arc processing and contemporary responses. `--flux-calibration archive` explicitly replays the archived bridge. Existing revisions without a recorded policy retain archive replay. Policy changes require a new revision.
- Candidate QC prepares inspection materials; it does not identify or extract a target automatically.
- Recovery recipes remain bound to exposure/input hashes. A previous source's bands are never global defaults.

## GUI revisions and nightly publication

- **Save recipe** and **Extract new revision** create separate inspection products;
  neither replaces `nightly/current.json` nor automatically updates the daily log.
- `nightly/current.json` selects a workflow revision, not an individual GUI spectrum.
  Do not point it at a GUI extraction folder or copy a product over a nightly export.
- To use reviewed geometry in a new nightly reduction, prepare a `--reviewed-sources`
  JSON record using the contract in `src/efospec/reviewed_sources.py`: raw basename,
  raw SHA-256, exact OBJECT, absolute parent-run/recipe paths, recipe SHA-256 and
  explicit operator identification/background-review evidence. The native route
  currently accepts one manually selected source with a `pypeit_local` recipe in
  unchanged raw detector coordinates; arbitrary GUI recipes are not promotable.
- Run `efospec night DATE --reviewed-sources /absolute/path/reviewed-sources.json`.
  This creates and selects a new workflow revision, re-extracts on native calibrated
  pixels and regenerates the report. Review its failures and scientific status;
  selection is not scientific acceptance. Earlier revisions remain available.
- For a report refresh of the **already selected** revision, run
  `python -m efospec.epoch_report DATE`. This reads saved products and does not
  re-extract, change selection, or import standalone GUI revisions.

## Optional IRAF comparison

- `efospec night DATE --iraf-comparison` enables a read-only comparison after
  reduction. Add `--stage report` to compare an existing selected reduction.
- Search is recursive within `data/DATE/`. Only header-identified PESSTO/IRAF
  flux-calibrated FITS multispectra (`FILETYPE=22208`) are considered. Counts,
  2D images, ASCII files and ZIP archives are not treated as flux spectra.
- Match a unique raw exposure using provenance plus OBJECT, setup, exposure time
  and observation time. Without raw identifiers, all metadata must still identify
  one exposure. Byte-identical files are deduplicated with aliases retained;
  distinct reductions are compared separately. Combined or ambiguous exposures
  require review and are not silently paired.
- Preserve delivered linear wavelength grids, physical flux units and native
  validity masks. No shift, normalization, interpolation or science acceptance
  decision is made. IRAF wavelength medium/frame omissions remain explicit.
- Results, input hashes and PDF/JPG overlays live in new revision-local
  `iraf-comparison/<processing-id>/` folders. Its `current.json` selects the
  comparison for that workflow revision; the report rejects stale parent hashes.
- The existing renderer adds an **IRAF comparison** tab and a separate linked
  `iraf-comparison.html` (with `iraf-comparison-portable.html` for single-file export).
  The main nightly spectrum panels show PypeIt only. Comparison y-limits include
  both spectra; explicitly requested legacy spectrum panels also include the
  native uncertainty band. Missing native products
  and unsupported comparisons stay visible. A failed comparison makes the report
  stage partial without changing reduction products or its scientific checks.

## Reviewed targets and targeted reruns

- `--target SN2026ayt` selects that exact FITS `OBJECT`; the same-night route also includes matching standard exposures needed for calibration.
  The nightly log states this scope; other exposures remain visible as outside
  the selection. Earlier full-night reductions remain in the archive.
- `--reviewed-sources FILE.json` records explicit visual identifications and
  exposure-specific bounded PypeIt recipes. The selection file is copied and
  hashed inside the revision. Each entry identifies the raw basename/hash,
  object, parent run, recipe path/hash, candidate number and operator evidence.
- The recipe's original input hashes are verified. Fresh raw electrons, detector
  window and wavelength map must match the reviewed parent exactly before the
  source geometry can be reused. The sky fit and extraction are recomputed.
  A mismatch fails explicitly; identification is never transferred by name alone.
- Successful reviewed targets retain their original candidate IDs in Source QC,
  with new residual images and summed background-subtracted profiles. Visual
  identification remains distinct from automatic astrometric identification.
- Confirmed choices take precedence: SN2026ayt C2 and SN2026aadx C3. Otherwise
  choose the strongest positive PypeIt peak with detection S/N >=5 inside the
  existing +/-15-pixel central window; ties or missing detections remain explicit
  failures. This geometry rule alone does not confirm target identity.
- An `operator_background` review permits an exposure-specific bounded recipe
  without promoting its nominal source identity to visually confirmed status.
  Local fits allow up to 100 rejection iterations while retaining the same
  clipping, convergence and two-sided support gates.
- AT/SN targets appear first in observation tables and spectrum/QC sections,
  grouped by target and time. Sorting never changes header-based frame roles.

## Recorded verification

- Night **2026-10-28**: observing log/report complete; reduction explicitly partial.
- Candidate QC generated for that date; the next epoch was not started.
- All **193** raw hashes unchanged; report links and 300-dpi PDF/JPG pairs checked.
- Original workflow validation: **174 tests passed, 90.69% coverage**; wheel build passed.
- Evidence: the dated `nightly/current.json`, its revision's `verification.json`, and `runs/epoch-workflow-full-tests.log`.
- [Current desktop validation](validation-desktop-2026-11.md) · [Desktop guide](desktop.md).

## Contemporary daily standards

- Implementation and limits: [daily response calibration](daily-standard-calibration.md).
- Same observing night and identical grism/filter/slit/binning/readout are required.
  Among successful fits, choose the nearest observation in MJD; ties fail explicitly.
- Standards and science use the same native raw processing and bounded extraction.
  By default, science requires both sky anchors and the independent 6363 Å check.
  An explicitly requested native manual replay can use `background-extract
  --native-flux --sky-validation two-anchor`: 5577/6300 must retain their quality,
  agreement (3 Å), and shift (20 Å) gates; 6363 remains a measured diagnostic.
  The policy is recorded in new manifests/FITS/ASCII and respected by report,
  telluric and redshift validation. Earlier products are never relabelled. Short standards
  remain labelled arc-only; their wavelength status accompanies the applied response.
- Missing/failed standards never trigger an archived-response fallback. Counts from
  a successfully extracted standard survive a response-fit failure.
- New nightly reductions include Gr11/Gr13 Free/GG495 and Gr16/OG530 native
  routes. Counts survive response failures; Gr16 retains counts and response
  diagnostics pending independent red-wavelength/flexure validation. Remaining
  checks appear per exposure in the CLI, Qt Source review and nightly log.
  No filtered response substitutes for Free. See [integration details](daily-standard-calibration.md#nightly-integration-12-november-2026).
- Archive-grid reviewed recipes fail explicitly in the native route: they require
  new native geometry review, including previously identified C2/C3 sources.
- The Calibration section shows actual standard exposures, response downloads and
  the number of held-out transfers. Absolute flux accuracy remains unvalidated.
