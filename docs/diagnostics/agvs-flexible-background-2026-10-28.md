# AT2026agvs: flexible background trials

- Removed the fixed 30-pixel preliminary sky requirement. The source-finding bootstrap tries 30, 25, 20 and 15 pixels per side, retaining each failure or selected width. This provisional model never becomes the delivered background.
- The final selector retains 15–50-pixel alternatives and independently supported width/degree comparisons. Source/neighbour guards, two-sided support, detector masks and distance limits remain unchanged.
- Distance violations are classified as geometry failures before fitting, allowing the native route to try wider calibrated windows.
- Fresh 2026-10-28 replay: 201-pixel window cannot place even the 15-pixel preliminary bands. In 301- and 401-pixel windows the preliminary fit completes, but final source/neighbour guards push every available permitted band too far from the aperture. No counts or flux spectrum was delivered.
- Separate wavelength failure persists: the 5577/6300 offsets disagree by approximately 3.15 Å; the two-anchor policy still rejects it.
- Evidence and review HTML: `runs/2026-10-28/nightly/20261113T141208Z-6fee4b/QC/agvs-flexible-background-20261114T155037Z/`. `final-replay/run.json` retains all attempts, preliminary-fit summaries, diagnostic FITS and exact source hashes. Earlier raw/science files and nightly selections remain unchanged.
