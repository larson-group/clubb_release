# Experimental PDFs

This is an archive of two CLUBB PDF prototypes:

- PDF-9 explores parcel-based, vertically coherent cloudy transport.
- PDF-10 explores a moment-driven trivariate two-Gaussian closure.

The directory keeps their core Fortran code, tests, and short research notes.
They are archived for reference and are not part of the active build. Bergmann
PDF-8 is a separate, retained closure.

Start with `doc/trivariate_transport_method_report.md` for the main overview.

The `python/` directory keeps the related retired Dash labs, diagnostic
utilities, and PDF-10 report scripts. The obsolete Dash pytest copies were
removed after their imports no longer matched the active app; their history
remains in Git. Generated figures, HTML, and
data are intentionally omitted; these scripts may depend on the old Dash/SAM
environment.
