<p align="center">
# Standards and Coding Practices for CLUBB

## Larson Group UW-Milwaukee
</p>

The following rules apply only to the parts of CLUBB that we have
written. We try not to modify code from external contributors such as
BUGSrad, Numerical Recipes, and so forth, because such external code is
less familiar to us, and we might introduce bugs. Furthermore, our
changes would be lost if we were to update to a new version of the
external code.

1.  When adding functions and subroutines create a header with a short
    description of their purpose and reference if there is one. Include
    in the header a description of each variable and its respective
    units. Use the following standardized layout:

        !-----------------------------------------------------------------------
        function xyzzy( etherx, ethery, etherz )

        ! Description:
        !   Evaluates the xyzzy function. 

        ! References:
        !   Eqn. 33 on p. 4551 of  
        !   Tweedle-dee and Tweedle-dum (1842), J. Imaginary Atmospheric 
        !   Sci., Vol. 23, pp. 4550--4556
        !-----------------------------------------------------------------------

          ! Included Modules
          use grid_dimensions, only: ndim ! Variable
         
          implicit none

          ! External Calls
          complex, external :: &
            xyxy_fnc ! (Description of function goes here)

          external :: &
            zyzy_sub ! (Description of subroutine goes here)

          ! It is a good idea to declare intrinsics, in case another 
          ! function or variable is accidentally given the same name.
          intrinsic :: &
            matmul, exp, log, dble

          ! If we had any derived types to declare, then 
          ! we would put them here before the constants.

          ! Local Constants
          real, parameter :: &
            local_tol = 1.0e-6 ! Local tolerance [no units]

          ! Input Variables
          complex, dimension(ndim), intent(in) :: &
            etherx,  & ! Ether in the x dimension  [kg/kg]
            ethery,  & ! Ether in the y dimension  [kg/kg]
            etherz     ! Ether in the z dimension  [kg/kg]

          ! Output Variables
          complex, dimension(ndim,ndim,ndim), intent(out) :: &
            xyzzy    ! Ether in the xyyzy dimension  [kg/kg]

          ! Local Variables
          integer :: i, j, k ! Loop iterator

        !-----------------------------------------------------------------------

          !----- Begin Code -----

            ...

          return
        end function xyzzy
        !-----------------------------------------------------------------------

2.  In functions and subroutines put a line with
    `!----- Begin Code -----` after the declarative portion of the code.
    This is to make it clear where the algorithmic portion of the
    procedure starts.

3.  Be generous with code comments. Make sure your code is
    understandable to an outsider looking at it for the first time. When
    writing a code comment, leave your name at the end of the comment.
    This will allow other users to contact you with any questions about
    the section of code. Use the same name identifier, e.g. dschanen or
    Vince Larson, everywhere. Then everyone can use `grep` to search for
    all code comments you’ve made.

4.  Minimize the number of variables and the amount of data that can be
    modified by subroutines and functions. This includes variables
    passed as an argument and those coming from a `use` statement.
    Always specify the `intent` attributes of arguments, and use
    `intent(in)` and `intent(out)` rather than `intent(inout)` whenever
    possible. This makes the code more clear to new users by allowing
    them to trace where a given variable is changed in the code and spot
    at a glance what are the inputs and outputs of each subroutine.
    Furthermore, this makes it easier for the compiler to optimize,
    because modern computers can vectorize floating point instructions,
    thereby providing some parallelization when your calculation does
    not rely on previous results. Note however, that an `intent(out)`
    variable that is not set in the scope of the subroutine is *not
    defined*! It often preserves the previous value, but it may not.
    When using modules with public variables, limit the number of
    variables with an `only:` statement. This makes code easier to read,
    because it is then obvious which variables from the module are being
    used. Even if all variables from the module are used, listing all
    variables after the `only:` statement is still worthwhile because it
    provides documentation about the inputs. E.g.

        !-----------------------------------------------------------------------
        subroutine mulmul( x, y, z )

        ! Description:
        !   Bogus example

        ! References:
        !   None
        !-----------------------------------------------------------------------

          ! Here we limit our use statement to just the needed variables
          use constants, only: &
            pi,         & ! Ratio of radii to their circumferance [no units]
            euler_const   ! Made up constant [no units]

          use dimensions, only: &
            nzmax ! Number of data points of our model [count]

          implicit none

          ! First list input variables
          real, dimension(nzmax), intent(in) :: &
            x, &! First factor   [no units]
            y   ! Second factor  [no units]

          ! Then list inout variables, and then outputs 
          real, dimension(nzmax), intent(out) :: &
            z  ! output         [no units]
        !-----------------------------------------------------------------------

          !----- Begin Code -----

          ! Note that every element of z could be calculated 
          ! independently of other elements of the subroutine
          ! allowing for parallelism if possible.

          z(1:nzmax) = ( x(1:nzmax) / pi ) * ( y(1:nzmax) / euler_const )

          return
        end subroutine mulmul
        !-----------------------------------------------------------------------

5.  Do not use implicitly saved variables in the scope of a subroutine.
    E.g.

          subroutine set_fluxes( i )
            implicit none

            logical :: l_use_fixed_fluxes = .false.

            ! ---- Begin Code ----
            if ( i < 9 ) then
              l_use_fixed_fluxes = .true.
            end if

            return
          end subroutine set_fluxes

    In this example `l_use_fixed_fluxes` will be saved until the next
    call to `set_fluxes`, but it’s not obvious this is the case unless a
    user looks at the declaration and knows that the variable will be
    implicitly saved. If you wish to save a variable either declare it
    in the common block of the module, bring it in via a use statement,
    or make it an argument for the subroutine.

6.  The general convention in CLUBB is to have the dummy arguments of a
    subroutine ordered with input first, input/output second, and output
    arguments third. The exception to this is optional variables, which
    must occur last. In that case, standard variables should be ordered
    as above; below this, the optional variables should be listed with
    input first, then input/output, and then output variables. E.g.


          ! The wrong way
          subroutine sub( input_variable_1, output_variable_1, &
                          input_variable_2, input_output_variable_1 )

          ! The right way
          subroutine sub( input_variable_1, input_variable_2, &
                          input_output_variable_1, output_variable_1 )

7.  Keep on the same lines the arguments of both the originating call
    and its corresponding subroutine declaration. This makes it easier
    to check the consistency of arguments of subroutines and functions
    with long argument lists. E.g.


          ! The wrong way
          call sub( roast_beef, ham, tomato, cheese, &
                    lettuce, onion, mayo )

          subroutine sub( roast_beef, ham, &
                          tomato, cheese,  &
                          lettuce, onion, mayo )

           ! The right way
           call sub( roast_beef, ham, tomato, cheese, &
                     lettuce, onion, mayo )

           subroutine sub( roast_beef, ham, tomato, cheese, &
                           lettuce, onion, mayo )

8.  When breaking up an arithmetic expression over multiple lines, do so
    according to operator precedence. This makes things more readable at
    a glance. E.g. Do this:

          x = t + y*x  &
            + g*v + t*h

    Not this:

          x = t + y*x + g &
            *v + t*h

9.  Do not use assumed size or assumed shape arrays in CLUBB. Generally
    speaking, the grid class defines the length of model arrays, and the
    use of assumed size arrays introduces another level of complexity
    that can contribute to introducing bugs. An exception to this rule
    is that should you need to interface with code developed by external
    contributors that contains assumed size or shape arrays you may have
    to declare your arrays the same way.

10. When a logical expression or a do loop extends more than ten lines,
    put a small comment after the end statement with the name of the
    loop or logical variable. E.g.

          if ( l_code_enabled ) then

            ...many lines of code...

          end if ! l_code_enabled

    or

          if ( .not. l_code_enabled ) then

            ...many lines of code...

          end if ! ~l_code_enabled

    A long do statement

          do i = 1, 100, 1

          ...many lines of code...

          end do ! i = 1..100

11. When making a logical expression, leave spaces between the logical
    operators and use the Fortran 90 style logical operators. E.g.

          ! Like this
          if ( ( x == y ) .and. ( y >= z ) .or. ( z /= 0.0 ) ) then

            ...

          end if

          ! Not like this
          if ( x.eq.y.and.y.ge.z.or.z.ne.0.0 ) then

            ...

          end if

    The F90 style operators are:

        == instead of .eq.
        /= instead of .ne.
        >  instead of .gt.
        >= instead of .ge.
        <  instead of .lt.
        <= instead of .le.

12. When creating a model flag variable of type `logical`, precede the
    name with a `l_` to make it easy to see that it is a logical.

13. When writing flow control statements (` if, else if, select case `),
    keep the logic as straightforward and clear as possible. If a set of
    `if...then` expressions grows to more than three cases, then rewrite
    the logic as a `select...case` expression instead. E.g.

          ! Rewrite this:
          if ( trim( run ) == "BOMEX" ) then
            execute_stuff( )

          else if ( trim( run ) == "FIRE" ) then
            execute_other_stuff( )

          else if ( trim( run ) == "ARM" ) then
            execute_other_other_stuff( ) 

          else if ( trim( run ) == "DAVE'S_CAMPING_TRIP" ) then
            execute_other_other_other_stuff( ) 

          else
            write(unit=0,fmt=*) "Cannot determine what to do for "//trim( run )

          end if ! run

          ! Like this:
          select case( trim( run ) )
          case( "BOMEX" )
            execute_stuff( )

          case( "FIRE" )
            execute_other_stuff( )

          case( "ARM" )
            execute_other_other_stuff( )
         
          case ( "DAVE'S_CAMPING_TRIP" ) 
            execute_other_other_other_stuff( ) 

          case default
            write(unit=0,fmt=*) "Cannot determine what to do for "//trim( run )

          end select ! run

14. When mixing scalar variables with array variables in a single
    mathematical statement, it is often more clear to use indices, even
    though Fortran 90 allows you to omit them. E.g. Do this:

         xyzzy(1:gr%nz) = variable_arr1(1:gr%nz) * variable_scl1 &
                            * variable_arr2(1:gr%nz)**2

    Not this:

         xyzzy = variable_arr1 * variable_scl1  &
                 * variable_arr2**2

15. Use a consistent character case for variables names. E.g. if you use
    `thlm` in function a, then do not use `Thlm` in function b. This
    makes the code easier to read and search using your text editor and
    the `grep` command.

16. Always use `implicit none` in your code. It will make debugging much
    easier, because it prevents the accidental use of variables of the
    wrong type.

17. When adding a variable that isn’t a loop iterator and occurs in more
    than one context, use a variable name that is in some way unique.
    This will make searching for the variable in the code much easier.
    There several ways to do this:

    1.  Use mixed case. E.g. `Rd` for the Dry air gas constant rather
        than `R` or `rd`

    2.  Use of a number in the name when the mathematical equations
        contain them. E.g. `p0` for $P_0$.

    3.  Use 3 or more characters in the variable name. E.g. CLUBB uses
        `grav` rather than `g` or `gr` for the gravitational constant.

    Avoid using a variable name that overlaps with the text of a
    reserved word or Fortran directive as well (e.g. `inte` overlaps
    with `intent` and `integer` and would not be a good variable name).
    Try using a command like `grep` if you’re uncertain if your variable
    is a semi-unique string.

18. Do not use magic numbers, i.e. undefined numbers rather than
    variables. Use of magic numbers makes it harder to read the code and
    to change the value parameters, e.g. grav = 9.8, consistently
    throughout the code.

    This is hard for a non-meteorologist to understand:

          rtm = rtm * 1000. 

    Instead, do this:

          g_per_kg = 1000.
          rtm = rtm * g_per_kg

    If you absolutely need to use a magic number, perhaps because a case
    specification uses one, include a comment on that line saying:

          ! known magic number

    This will prevent the check\_for\_errors.py script from reporting
    that there is a magic number on this line.

19. Do not use magic flags, i.e. undefined input arguments to
    subroutines or functions. Instead, define the argument just above
    the call to the subroutine.

    To understand this, you’d need to go to the subroutine:

          call complex_function( wp2, wp3, .true. ) 

    This is better documented:

          l_enable_damping = .true.
          call complex_function( wp2, wp3, l_enable_damping )

    If you absolutely need to use a magic flag, include a comment on
    that line saying:

          ! known magic flag

    This will prevent the check\_for\_errors.py script from reporting
    that there is a magic flag on this line.

20. Do not use the prognostic or diagnostic variable modules anywhere
    they are not used already. Currently the module
    `prognostic_variables` is used in `numerical_check` and
    `clubb_driver`, while `diagnostic_variables` is used in
    `clubb_driver`, `clubb_core`, and `stats_subs`. Global variables
    like these make debugging the code more difficult, and complicate
    implementing CLUBB in a host model. The proper way to use these
    variables in a new section of code is to pass them as an argument.

21. Put all functions and subroutines within a module. This allows for
    compile time checking of input arguments.

22. When putting several subroutines or functions in a module, limit the
    number of public interfaces to those that are called from outside
    the module.\
    E.g.

          ! module gcss

            ...

          public ::
            cloud_rad, atex_tndcy, atex_sfclyr, fire_tndcy, &
            wangara_tndcy, wangara_sfclyr, &
            astex_tndcy, astex_sfclyr, &
            arm_tndcy, arm_sfclyr, &
            bomex_tndcy, bomex_sfclyr, &
            dycoms2_rf01_tndcy, dycoms2_rf01_sfclyr, &
            dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr, &
            nov11_altocu_tndcy 

          ! Note that all these are called from 
          ! subroutines within gcss;  This declaration
          ! makes them unavailable to outside routines

          private ::
            diag_ustar, arm_sfcflx, Diff_denom, &
            altocu_icedf

              ...

23. In order to keep CLUBB threadsafe, all variables in a module prior
    to the `contains` statement (i.e. common block variables) should be
    declared as `$omp threadprivate` in a commented line just after the
    variable declaration. The only exception to this rule is a variable
    with the `parameter` attribute, since these hold a constant value.
    Below is a code example:\

        module declination
          real, public :: angle_in_deg = 30. ! Degree of the angle

        !$omp threadprivate(angle_in_deg)
          contains
        ...
        end module declination

24. Better yet, don’t use module variables at all unless they have the
    parameter attribute. That will avoid the need for `threadprivate`
    directives and will make the flow of information easier to
    understand in a global model that calls CLUBB.

25. Good use of whitespace can make code easier to read and debug. One
    simple convention used in CLUBB is to put an extra space between
    parenthesis for functions and subroutines to differentiate them from
    array variables. E.g.

          ! Choosing less obscure names for arrays and functions 
          ! will also help with this problem.

          f( x ) ! This is the function f evaluated with argument x.
          f(x)   ! This is the one dimensional array f indexed at x.

    This applies to intrinsic functions as well. E.g.

          mycalc = exp( x ) * log( y )  &
                 + dble( g ) * mod( f, 5 )

26. For clarity, functions, subroutines, modules, and programs should
    end with their proper Fortran 90 style end statement. E.g.

          ! Like this
          subroutine foo

             ...

          end subroutine foo

          ! Not like this
          subroutine foo

            ...

          end

27. Separate end statements with a single space. E.g.

          subroutine sub( )

            ...

          do k = 1, gr%nz, 1

                   ...

          end do ! not 'enddo'

                   ...

          where ( rtm < rttol )
            ...
          end where ! not 'endwhere'

          subroutine sub

                   ...

            return
          end subroutine sub ! not 'endsubroutine'

28. For consistency use lowercase letters for Fortran keywords. E.g.

        program HAL ! lowercase program

            ... 

          if ( l_TurnOff ) then ! lowercase if, lowercase then

            write(unit=*,fmt=*) "I know that you and Frank were " &
              //"trying to turn me off Dave." ! lowercase write

            stop "Daisy, Daisy..." ! lowercase stop

          end if ! lowercase end, lowercase if

            ...

        end program HAL ! lowercase end, lowercase program 

29. For the sake of consistency and code readability, always use a 2
    column indent after all directives. Also, indent continuation lines
    by 2 columns. E.g.\
    `module dynamics`\
    `  ` <span>2 column indent within the scope of `module`</span>\
    ``\
    `  implicit`` none`\
    ``\
    `  contains`\
    ``\
    `  subroutine`` advect_2D( ... )`\
    `    use`` grid, only: & ` <span>2+2 columns within the scope of
    `dynamics/advect_2D`</span>\
    `      nzmax,`` & ` <span>2+2+2 columns since we’re continuing the
    `use` statement</span>\
    `      nx`\
    ``\
    `    if`` ( l_do_something ) then`\
    `      ` <span>2+2+2 columns within the scope of
    `dynamics/advect_2D/if..then`</span>\
    ``\
    `      ``call something( ... )`\
    ``\
    `    end`` if`\
    ``\
    `    return`\
    `  end`` subroutine advect_2D`\

    `end module dynamics`

30. Wrap all code lines at 100 columns. Fortran free format source code
    allows for 132, but the number 100 was chosen to so that we can view
    the code browser without word-wrap on one monitor.

31. Do not use deprecated or obsolescent features of Fortran, and do not
    use extensions to the Fortran 95 standard. This includes but is not
    limited to `pause`, `equivalence`, Cray style pointers, arithmetic
    `if`, alternate `return`, and `float()` for real data type
    conversions.

32. All files added to the CLUBB code should use a `.F90` file
    extension. This is needed for a kluge used with `mkmf` that only
    enables compiler warnings when the file is pre-processed Fortran 90
    source.

33. Any new files added to the repository should have matching filenames
    and module names. This is because some host models in which CLUBB is
    implemented, e.g. SAM and CAM, use make/compile scripts that assume
    all filenames and module names are the same. Also, if a new module
    name requires an extension, use `<module name>_module` rather than
    `<module name>_mod`, for clarity and consistency with the rest of
    the code.

34. CLUBB uses the convention that logical variables should have a
    prefix of `l_`. For example, use `l_cloud_sed` rather than
    `cloud_sed`.

35. All non-fatal print or write statements should recast as
    `clubb_debug( 1 or 2, "message")` or should be enclosed in
    `clubb_at_least_debug_level( 1 or 2 )` conditionals. All fatal
    errors should be write statements to stderr
    (`clubb_constants fstderr`) with no `clubb_at_least_debug_level`
    conditional.

36. When using allocatable pointers or derived types they should be
    deallocated in the same function or subroutine in which they were
    allocated. This is because unlike arrays, the memory used by
    pointers is not automatically freed without a deallocate statement,
    and we want to prevent memory leaks. This includes derived types
    that contain allocatable pointer components. E.g.\

          function func

              ...

            allocate( file%z )

              ...

            deallocate( file%z )

            return
          end function func

37. Variable locations in CLUBB:

    <span>**A. Thermodynamic-level Variables**</span>

    CLUBB generally defines mean fields (thlm, rtm, rcm, um, vm, etc.)
    and third-order moments (wp3, wp2thlp, wprtp2, etc.) to be
    thermodynamic-level variables. Furthermore, pressure variables
    (p\_in\_Pa and exner) are thermodynamic-level variables. Their
    \`\`home" is on thermodynamic levels.

    Furthermore, most variables that are diagnosed based purely on
    thermodynamic-level variables are thermodynamic-level variables
    themselves. A good example of this is Lscale. Lscale is computed in
    subroutine compute\_length based on the following
    thermodynamic-level inputs: thvm, thlm, rtm, rcm, p\_in\_Pa, and
    exner. There is only one momentum-level variable input into
    compute\_length, and that is em (as a side note, em is a rare
    exception to the \`\`mean field" rule because it is diagnosed based
    purely on three momentum-level variables, wp2, up2, and vp2). Thus,
    em is interpolated to thermodynamic-levels within compute\_length.
    Another good example of this rule is thvm, which is diagnosed based
    on thermodynamic-level variables thlm, rtm, exner, and rcm.

    When thermodynamic-level variables are at their \`\`home“
    thermodynamic levels, they simply use their \`\`base” name (thlm,
    wp2rcp, wp3, Lscale, etc.). However, when thermodynamic-level
    variables are interpolated to momentum levels (or computed on
    momentum levels in some instances), the general naming convention is
    to tack a \`\`\_zm“ onto the end of their \`\`base” name. For
    example, thlm is thlm at thermodynamic levels and thlm\_zm when
    interpolated to momentum levels.

    <span>**B. Momentum-level Variables**</span>

    CLUBB generally defines second-order moments (wpthlp, rtp2, wp2,
    sclrpthlp, etc.) and fourth-order moments (wp4) to be momentum-level
    variables. Their \`\`home" is on momentum levels.

    Furthermore, most variables that are diagnosed based purely on
    momentum-level variables are momentum-level variables themselves. A
    good example of this is sigma\_sqd\_w, which is diagnosed based
    purely on momentum-level variables wpthlp, wp2, thlp2, wprtp, and
    rtp2. Another example of this, as referenced above in the section on
    \`\`Thermodynamic-level variables", is em, which is diagnosed based
    purely on momentum-level variables wp2, up2, and vp2.

    When momentum-level variables are at their \`\`home“ momentum
    levels, they simply use their \`\`base” name (wprtp, up2,
    sigma\_sqd\_w, rtpthlp, etc.). However, when momentum-level
    variables are interpolated to thermodynamic levels (or computed on
    thermodynamic levels in some instances), the general naming
    convention is to tack a \`\`\_zt“ onto the end of their \`\`base”
    name. For example, rtp2 is rtp2 at momentum levels and rtp2\_zt when
    interpolated to thermodynamic levels. Likewise, sigma\_sqd\_w is
    sigma\_sqd\_w at momentum levels and sigma\_sqd\_w\_zt when
    interpolated to thermodynamic levels.

    <span>**C. Homeless Variables**</span>

    There are some variables in the CLUBB code that don’t have a
    definitive \`\`home". Rather, they just have numerous cardboard
    boxes scattered amongst the various different vertical levels.

    Generally, these are variables that are diagnosed based on nearly
    equal contributions from momentum-level variables and
    thermodynamic-level variables. A great example of this is Skewness
    of w. The skewness of w, Skw, equals `wp3 / wp2^1.5`. It is based on
    one thermodynamic-level variable, wp3, and one momentum-level
    variable, wp2. Thus, when calculating Skw on thermodynamic levels,
    `Skw = wp3 / wp2_zt^1.5`, and when calculating Skw on momentum
    levels, `Skw = wp3_zm / wp2^1.5`. Two more great examples of this
    are eddy diffusivity, Kh, and time-scale, tau. Both are based on one
    thermodynamic-level variable, Lscale, and one momentum-level
    variable, em.

    When a homeless variable is at its cardboard box on momentum levels,
    a \`\`\_zm“ extension is tacked onto it, and when it is at its
    cardboard box on thermodynamic levels, a \`\`\_zt” extension is
    tacked onto it. Thus, skewness of w is called Skw\_zt at
    thermodynamic levels and Skw\_zm at momentum levels. Likewise, eddy
    diffusivity and time-scale are Kh\_zt and tau\_zt, respectively, at
    thermodynamic levels, and Kh\_zm and tau\_zm, respectively, at
    momentum levels.

38. When using find/replace to make changes to the code, such as
    renaming a function or variable, mistakes can be made and
    overlooked. In particular, unwanted changes to comments frequently
    occur. To alleviate this, whenever using find/replace to make
    changes to the code, you must check all changes made to comments.

    A regular expression can easily be used with the git diff and egrep
    commands to find all changes made to comments in a local revision.
    Simply run the command from the base directory of your CLUBB
    checkout. The command is as follows:

        git diff | egrep '(^(\+|-).*!|^---|^\+\+\+|^@@|^===)' -

    Each pipe in this egrep command represents an OR statement. The
    first section includes all lines from the diff that include an
    exclamation mark, which represents a comment. The rest of the lines
    include the lines that show information about the filenames and line
    numbers.

39. CLUBB uses a generalized precision that can be specified with a
    compiler flag. This compiler flag is CLUBB\_REAL\_TYPE. This maps to
    the variable in clubb\_precision.F90 called core\_rknd. When using
    any real variables, please use a real of kind core\_rknd instead of
    a standard real, unless it is absolutely necessary to use single or
    double precision all the time.

    It is particularly important to specify core\_rknd when passing
    variables between CLUBB and an external model, such as WRF, which
    uses single precision by default, or MG microphysics, which uses
    double precision by default. If core\_rknd is not explicitly
    specified in an argument list from the external model, e.g.

        real( exciting_WRF_variable(k), kind = core_rknd ) ,

    then CLUBB will fail to compile when CLUBB’s precision is changed,
    i.e. CLUBB\_REAL\_TYPE is changed.

    Also, when assigning values to these variables, you must specify the
    precision of the value using \_core\_rknd (E.g. 1.23\_core\_rknd).
    Here is an example of declaring a variable as core\_rknd and
    assigning a value to it:

        use clubb_api_module, only &
          core_rknd ! generalized precision

        real( kind = core_rknd ) :: &
          variable

        variable = 0.0_core_rknd

40. We should avoid adding new use of command line applications, e.g.
    mktemp, to our compile or run scripts unless we’re certain they’re
    portable between Linux/BSD/MacOS/Unix. See
    <http://carson.math.uwm.edu/trac/clubb/ticket/531#comment:24>

41. There are certain nuances to nested WHERE statements that can cause
    runtime crashes on certain compilers. The general cause is that the
    logical evaluation of an inner WHERE statement would cause a runtime
    error if the code between between the outer and inner WHERE
    statements is not run. For example, the following code has the
    potential to cause a runtime error:

            where ( a > 0 )    <--- Outer where statement 

                b = 1 / a      <--- Code between outer and inner where statements
                                    This is the first time b is defined 

                    where ( b > 5 )        <--- Inner where statement
                                                b should be defined at this point
                        b = 5

                    end where

            end where

    This seems perfectly safe, but compilers may compile this code in a
    way similiar to the following:

            where ( a > 0 .and. b > 5 )     <--- b has NOT been defined yet
                                            This could cause a floating point exception
                b = 5

            elsewhere ( a > 0 ) 

                b = 1 / a      

            end where

    The statements are logically equivalent, but because b is defined
    within the WHERE statement, the evaluation of the logical value may
    cause a floating point exception. If floating point trapping is not
    enabled, the output behavior would be exactly as intended, but with
    floating point trapping the run could stop. It’s also possible for
    the inner WHERE statement to contain a statement that could cause a
    runtime error regardless of debugging option, which could be even
    more dangerous.

    It’s unclear whether or not this behavior is a bug with certain
    compilers or if it is the intended behavior as defined by Fortran
    standards. To be safe, we should avoid this potential error, either
    by refactoring the WHERE statement in some way, or by replacing the
    WHERE statement altogether with a DO loop containing IF statements.
    This issue is explained in more depth in CLUBB issue 859 on github,
    https://github.com/larson-group/clubb/issues/859.

42. When including the `clubb_config_flags` or the `silhs_config_flags`
    derived types or any of their components in argument lists only pass
    the complete `clubb_config_flags` to the subroutine
    `advance_clubb_core` and only pass the complete `silhs_config_flags`
    to the subroutine `generate_silhs_sample`. For all other subroutines
    and functions where any component of one of the derived types is
    needed just pass the specific components along in the argument list.
    This increases code readability and makes debugging easier because
    one can easily see which components are used where.


