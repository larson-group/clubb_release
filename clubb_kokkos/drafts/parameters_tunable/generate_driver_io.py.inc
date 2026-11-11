    # Insert before the existing individual-control section in generate_driver_io.main.
    # Read initialized parameter arrays using the exact Native namelist declaration.
    params_source = root / "src/CLUBB_core/parameters_tunable.F90"
    params_text = params_source.read_text()
    params_body = re.search(
        r"(?ims)^  subroutine init_clubb_params_api\b.*?^  end subroutine init_clubb_params_api\b",
        params_text,
    )[0]
    params_code = fp.statements(params_body)
    params_declarations = {}
    for statement in params_code:
        params_declarations.update(fp.declaration(statement))
    params_namelist = next(
        statement for statement in params_code
        if re.match(r"namelist\s*/clubb_params_nl/", statement, re.I)
    )
    params_names = split_top(params_namelist.split("/", 2)[2])
    params_args = [
        fp.Argument("runfile", "character", (), "in", False, False, "*"),
        dataclasses.replace(params_declarations["ngrdcol"], intent="in"),
        dataclasses.replace(params_declarations["iunit"], intent="in"),
    ] + [
        dataclasses.replace(params_declarations[name.lower()], intent="inout")
        for name in params_names
    ] + [fp.Argument("read_status", "integer", (), "out", False, False, "")]
    params_proc = fp.Procedure(
        "read_clubb_params_nl", "f2k_driver_io",
        "src/CLUBB_core/parameters_tunable.F90", params_args, None,
    )
    cpp.append('extern "C" void f2k_proc_read_clubb_params_nl_c(const WireField*, WireField*);')
    cpp += [
        line.replace(
            "Native contract: src/CLUBB_core/parameters_tunable.F90:read_clubb_params_nl.",
            "Read &clubb_params_nl declared in parameters_tunable::init_clubb_params_api.",
        ) for line in gp.cpp_wrapper("read_clubb_params_nl", params_proc)
    ]
    params_wrapper = gp.f_wrapper("read_clubb_params_nl", params_proc)
    params_wrapper = [
        line.replace(
            "Convert one call to src/CLUBB_core/parameters_tunable.F90:read_clubb_params_nl.",
            "Read &clubb_params_nl using initialized caller arrays and a call-local unit.",
        ) for line in params_wrapper if "native_call =>" not in line
    ]
    params_wrapper = [re.sub(r"\bv_(\w+)", r"\1", line) for line in params_wrapper]
    params_implicit = next(i for i, line in enumerate(params_wrapper) if "integer :: n," in line)
    params_wrapper[params_implicit + 1:params_implicit + 1] = [
        "    integer :: unit",
        "    ! Required iunit is retained, but NEWUNIT follows the existing call-local I/O owner.",
        "    namelist /clubb_params_nl/ " + ", &\n      ".join(params_names),
    ]
    params_call = next(i for i, line in enumerate(params_wrapper) if "call native_call(" in line)
    params_wrapper[params_call:params_call + 1] = [
        '    open(newunit=unit,file=trim(runfile),status="old",action="read",iostat=read_status)',
        "    if (read_status == 0) then",
        "      read(unit,nml=clubb_params_nl,iostat=read_status)",
        "      close(unit)",
        "    end if",
    ]
    fort += params_wrapper
