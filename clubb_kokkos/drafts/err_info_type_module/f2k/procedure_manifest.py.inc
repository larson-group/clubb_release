# Merge into PROCEDURES.update in the existing transient procedure generator.
# The distinct C symbol selects the PUBLIC Native setter with all optionals absent.
"format_err_info_headers": ("set_err_info_values_api", set()),

# Merge into existing PUBLIC_NAMES to retain the Native public C++ overload name.
"format_err_info_headers": "set_err_info_values_api",

# Private C++ update_headers calls f2k::set_err_info_values_api(ngrdcol, snapshot).
# After permanent C++ caller replacement, remove the whole-init and lat/lon-only
# setter entries. Keep existing private converter accessors transient.
