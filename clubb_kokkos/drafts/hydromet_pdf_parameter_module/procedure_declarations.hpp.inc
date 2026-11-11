// Merge inside clubb_kokkos after the generated types.
void init_hydromet_pdf_params(
    hydromet_pdf_parameter& hydromet_pdf_params // InOut
);
void init_precip_fracs_api(
    int nzt, int ngrdcol, // In
    precipitation_fractions& precip_fracs // Out
);
void zero_precip_fracs_api(
    precipitation_fractions& precip_fracs // InOut
);
// Native repeated per-cell initialization through the existing SoA owner.
inline void init_hydromet_pdf_params(
    hydromet_pdf_parameter_collection& hydromet_pdf_params // InOut
) {
  // Reset the existing case-owned collection; its schema owns all 18 fields.
  hydromet_pdf_params.zero();
}
