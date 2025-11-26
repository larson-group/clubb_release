SET(NUMERICAL_RECIPE_SOURCES
                    amebsa.f90
                    amoeba.f90
                    nr.f90
                    nrtype.f90
                    nrutil.f90
                    ran1_s.f90
                    ran1_v.f90
                    ran_state.f90
                    enhanced_simann.F90
                    quicksort.F90
)

add_library(numerical_recipes_lib ${NUMERICAL_RECIPE_SOURCES})

# Only needed for mt95 in old version of algorithm in enhanced_simann.F90
# which would be nice to remove
target_link_libraries(numerical_recipes_lib PRIVATE silhs_lib )
