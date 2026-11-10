# Neto reference sources

Source: https://github.com/Pedro-Miguel-Neto/PBS-MINI-0.1L
Paper: Neto et al. (2022), https://doi.org/10.3390/bioengineering9080386

The upstream repository does not declare a software license. Its MATLAB source files and data table are therefore not relicensed or bundled here. `provenance.json` records the original URLs, commit and checksums. `python publication/fetch_reference.py` retrieves and checksum-checks the one original diffusivity table required for the oxygen calculation. The publication reproduction command performs this step automatically. This file is an external input and is outside the release's Apache/CC BY license declarations.

Our Julia oxygen calculation implements the published mathematical correlation and water-property equations with attribution; it does not require running the upstream MATLAB code.
