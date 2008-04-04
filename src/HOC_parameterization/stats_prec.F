!-----------------------------------------------------------------------
! $Id: stats_prec.F,v 1.1 2008-04-04 22:02:49 dschanen Exp $
        module stats_prec

        implicit none

        public

        ! The precisions below are arbitrary, and could be adjusted as
        ! needed for long simulations or time averaging.  Note that on
        ! most machines 12 digits of precision will use a data type
        ! which is 8 bytes long.
        integer, parameter :: 
     .    stat_nknd = selected_int_kind( 12 ),
     .    stat_rknd = selected_real_kind( p=12 ),
     .    time_prec = selected_real_kind( p=12 )

        end module stats_prec
!-----------------------------------------------------------------------
