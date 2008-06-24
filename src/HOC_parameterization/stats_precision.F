!-----------------------------------------------------------------------
! $Id: stats_precision.F,v 1.1 2008-06-24 18:50:54 faschinj Exp $
        module stats_precision

        implicit none

        public :: stat_nknd, stat_rknd, time_prec

        private ! Default scope

        ! The precisions below are arbitrary, and could be adjusted as
        ! needed for long simulations or time averaging.  Note that on
        ! most machines 12 digits of precision will use a data type
        ! which is 8 bytes long.
        integer, parameter :: 
     .    stat_nknd = selected_int_kind( 8 ),
     .    stat_rknd = selected_real_kind( p=12 ),
     .    time_prec = selected_real_kind( p=12 )

        end module stats_precision
!-----------------------------------------------------------------------
