// Merge into the existing clubb_api_module.hpp owner after its type includes.
inline void init_pdf_params_api(int nz, int ngrdcol,       // In
                                pdf_parameter& pdf_params) // Out
{
  // Description: Allocate and initialize PDF parameters through their Native owner.
  // nz: Number of vertical grid levels [-].
  // ngrdcol: Number of grid columns [-].
  // pdf_params: Caller-owned PDF parameters [units vary].
  init_pdf_params(nz, ngrdcol, // In
                  pdf_params); // Out
}
