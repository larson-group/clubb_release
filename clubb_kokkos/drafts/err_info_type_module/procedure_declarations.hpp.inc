// Merge inside clubb_kokkos after the generated types.
void init_default_err_info_api(
    int ngrdcol, // In
    err_info_type& err_info // InOut
);
void init_err_info_api(
    int ngrdcol, int chunk_idx_in, int mpi_rank_in, Real1D lat_in, Real1D lon_in, // In
    err_info_type& err_info // InOut
);
void set_err_info_values_api(
    int ngrdcol, // In
    err_info_type& err_info, // InOut
    std::optional<int> chunk_idx_in = std::nullopt, // Optional In
    std::optional<int> mpi_rank_in = std::nullopt, // Optional In
    std::optional<Real1D> lat_in = std::nullopt, // Optional In
    std::optional<Real1D> lon_in = std::nullopt // Optional In
);
void cleanup_err_info_api(
    err_info_type& err_info // InOut
);
template <class Lat, class Lon>
inline void init_err_info_api(
    int ngrdcol, int chunk_idx_in, int mpi_rank_in, Lat lat_in, Lon lon_in, // In
    err_info_type& err_info // InOut
) {
  // Stage host/strided inputs once through the existing array owner.
  init_err_info_api(
      ngrdcol, chunk_idx_in, mpi_rank_in, to_execution_space(lat_in), // In
      to_execution_space(lon_in), // In
      err_info // InOut
  );
}
template <class Lat = Real1D, class Lon = Real1D>
inline void set_err_info_values_api(
    int ngrdcol, // In
    err_info_type& err_info, // InOut
    std::optional<int> chunk_idx_in = std::nullopt, // Optional In
    std::optional<int> mpi_rank_in = std::nullopt, // Optional In
    std::optional<Lat> lat_in = std::nullopt, // Optional In
    std::optional<Lon> lon_in = std::nullopt // Optional In
) {
  // Optional HostSpace/strided Views retain independent presence semantics.
  std::optional<Real1D> lat_exec;
  std::optional<Real1D> lon_exec;

  if (lat_in.has_value())
    lat_exec = to_execution_space(*lat_in);

  if (lon_in.has_value())
    lon_exec = to_execution_space(*lon_in);
  set_err_info_values_api(
      ngrdcol, // In
      err_info, // InOut
      chunk_idx_in, mpi_rank_in, lat_exec, lon_exec // Optional In
  );
}
template <class Lat, class Lon,
          std::enable_if_t<Kokkos::is_view_v<Lat> && Kokkos::is_view_v<Lon>, int> = 0>
inline void set_err_info_values_api(
    int ngrdcol, // In
    err_info_type& err_info, // InOut
    std::optional<int> chunk_idx_in, std::optional<int> mpi_rank_in, Lat lat_in, // Optional In
    Lon lon_in // Optional In
) {
  // Existing driver supplies both coordinate arrays as ordinary Views.
  set_err_info_values_api(
      ngrdcol, // In
      err_info, // InOut
      chunk_idx_in, mpi_rank_in, std::optional<Lat>(lat_in), // Optional In
      std::optional<Lon>(lon_in) // Optional In
  );
}
