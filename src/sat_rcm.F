        ! This function was created by Brian Griffin.
        ! If you have information about your initial profile in thlm and rtm,
        ! and your profile is saturated at some point, you need to use an
        ! iterative method in order to accurately find the value of rcm.

        FUNCTION sat_rcm( thlm, rtm, p, exner )

        USE constants

        real, external:: rsat

        REAL, INTENT(IN):: thlm
        REAL, INTENT(IN):: rtm
        REAL, INTENT(IN):: p
        REAL, INTENT(IN):: exner

        REAL:: sat_rcm

        REAL:: theta
        REAL:: answer, too_low, too_high

        INTEGER:: iteration

        REAL, PARAMETER:: tolerance = 0.001


        theta = thlm
        iteration = 0

        DO

           iteration = iteration + 1

           answer = theta - (Lv/(Cp*exner))
     .                     *(MAX( rtm - rsat(p,theta*exner), 0.0 ))

           IF ( ABS(answer - thlm) <= tolerance ) THEN
              EXIT
           ELSEIF ( answer - thlm > tolerance ) THEN
              too_high = theta
           ELSEIF ( thlm - answer > tolerance ) THEN
              too_low = theta
           ENDIF

           ! For the first timestep, be sure to set a "too_high"
           ! that is "way too high."
           IF ( iteration == 1 ) THEN
              too_high = theta + 20.0
           ENDIF

           theta = (too_low + too_high)/2.0

        ENDDO

        sat_rcm = MAX( rtm - rsat(p,theta*exner), 0.0 )

        PRINT *, theta, thlm, rtm, sat_rcm

        END FUNCTION sat_rcm
