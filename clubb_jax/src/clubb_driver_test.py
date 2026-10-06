"""Minimalist CLUBB frontend adapted from src/clubb_driver_test.F90.

Read the same generated namelist as the standalone. Initialize, clean up and
initialize again, then advance without statistics, reset initial conditions
and advance with statistics. Compare the resulting NetCDF output with a
separate normal JAX standalone run, as the native driver Jenkins job does.

Run through run_scripts/run_scm.py -jax -driver_test CASE, or the managed
launcher with -module=clubb_jax.src.clubb_driver_test NAMELIST. The default
namelist is clubb.in, matching the Fortran program.

JAX adaptations: a state dictionary replaces driver globals, exceptions retain
cleanup through try/finally, and successful execution returns the native code
6. The labeled output-batch loop extends the native program to write every
requested parameter column. Additional window/batch assertions are kept
separately in clubb_jax/tests/run_driver_extensions_test.py.
"""

import sys

from clubb_jax.src.clubb_case_initalization import (
    init_clubb_case,
    set_case_initial_conditions,
    clean_up_clubb,
)
from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end


def main():
    # Constant parameters.
    success_code = 6
    l_stdout = True

    # Set default namelist values.
    namelist_filename = "clubb.in"
    if len(sys.argv) >= 2 and sys.argv[1]:
        namelist_filename = sys.argv[1].rstrip(' ')

    # ------------------------------ Test Sections ------------------------------

    print("======================== REINITIALIZATION TEST ========================",
          file=sys.stderr, flush=True)
    print("This section ensures that everything allocated in init_clubb_case\n"
          "will be deallocated in clean_up_clubb. This may cause a runtime\n"
          "error if there is a mismatch between\n"
          "allocate/deallocate statements, but could",
          file=sys.stderr, flush=True)

    # Read in model parameter values by initializing, cleaning up, and initializing again.
    state = init_clubb_case(namelist_filename)
    clean_up_clubb(state)
    state = init_clubb_case(namelist_filename)

    try:
        print("======================== DOUBLE TIMESTEP RUN ========================",
              file=sys.stderr, flush=True)
        print("Calling advance_clubb_to_end, then set_case_initial_conditions,\n"
              "then advance_clubb_to_end again. This could result in a runtime\n"
              "error if an allocation statement appears in these routines, since\n"
              "we don't call clean_up_clubb in between the advance_ calls.",
              file=sys.stderr, flush=True)

        # Run once with stats off. Turning stats on writes output to disk,
        # and we only want that for the second standalone comparison run.
        advance_clubb_to_end(state, l_stdout, l_suppress_stats=True)

        # Reset model back to initial conditions.
        set_case_initial_conditions(state)

        # JAX extension: the native test writes its initialized active columns.
        # Visit the remaining batches as run_clubb does so the NetCDF comparison
        # covers every requested parameter column. Batch 1 was just reset above.
        num_batches = state['total_param_sets'] // state['ngrdcol']
        for batch_idx in range(1, num_batches + 1):
            if batch_idx > 1:
                set_case_initial_conditions(state, batch_num=batch_idx)

            # Turn stats on and run. Differences from run_clubb's output mean
            # set_case_initial_conditions is not correctly resetting the state.
            advance_clubb_to_end(state, l_stdout)
            if state['err_info'].is_fatal():
                raise RuntimeError("Fatal error in clubb_driver_test")

        print("======================== RUN OVER ========================",
              file=sys.stderr, flush=True)
        print("WARNING: The double timestep test is not complete until the\n"
              "output is compared with the output from calling clubb_standalone.\n"
              "This driver should produce BFB output with clubb_standalone.\n"
              "Save the output, rerun with clubb_standalone, then compare.\n"
              "If there are differences, set_case_initial_conditions is likely\n"
              "not resetting everything it needs to.",
              file=sys.stderr, flush=True)
    finally:
        clean_up_clubb(state)

    # Clean up err_info. JAX owns this value in state rather than an allocation.
    if state['err_info'].is_fatal():
        raise RuntimeError("Fatal error in clubb_driver_test")
    print("Program exited normally", file=sys.stderr, flush=True)
    return success_code


if __name__ == '__main__':
    raise SystemExit(main())
