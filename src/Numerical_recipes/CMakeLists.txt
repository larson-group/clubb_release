SET(NUMERICAL_RECIPE_SOURCES
                    amebsa.f90
                    amoeba.f90
                    nr.f90
                    nrtype.f90
                    nrutil.f90
                    ran1_s.f90
                    ran1_v.f90
                    ran_state.f90
)

add_library(numerical_recipes_lib ${NUMERICAL_RECIPE_SOURCES})
