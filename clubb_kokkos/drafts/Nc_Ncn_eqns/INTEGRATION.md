This draft ports all six Native elemental functions in `Nc_Ncn_eqns.F90`.
The four public scalar interfaces keep their original names and arguments;
the two private PDF integrals remain device-callable in their own namespace.
All Native branches, derivations, tolerances and ERFC/SQRT expressions are retained.

The four additional View overloads implement the actual elemental array boundary
without another physics routine. They accept conformable homogeneous rank-1 or
rank-2 Views under the same Native names and arguments, and return an owned
`LayoutLeft` View in the first input's memory space. Their explicit return types
support CUDA extended lambdas. Execution follows that input's execution space;
other inputs must be accessible from it. Host and strided inputs therefore need
no forced CUDA dispatch or per-element host bridge. All-scalar calls select the
Native scalar overload. Mixed scalar/View broadcasting and other ranks have no
current repository caller and are not introduced as a new C++ API.

For the compatible setup PDF batch, include `Nc_Ncn_eqns.hpp` and replace its
`f2k::Nc_in_cloud_to_Ncnm(nzt,ngrdcol,...)` with
`Nc_in_cloud_to_Ncnm(...)`, passing the original ten inputs. Remove the transient
shape adapter/declaration when that replacement is integrated. The existing
shared array owner supplies storage types; no scientific lower bridge, new
constant, duplicate type or mutable global state is required.

The full 953-line Native and 484-line JAX files were read before the outline.
JAX is scope context: its eager branch masks and safe denominators are not
imported. Native branches evaluate ERFC/division only where called. No FMA,
association barrier, identity special case or GPU math emulation was added.

`SOURCE_AUDIT.json` records the owner source review of all 57 arguments,
21 ordered assignments, ten predicates and four array kernels.
`CLEANUP_AUDIT.json` records the separate argument-group/comment/whitespace pass,
including the six nested scientific function calls. Independent full-file
review, compilation and numerical evidence remain pending. No BFB claim is made;
CPU zero-tolerance validation belongs to integration of the compatible group.
