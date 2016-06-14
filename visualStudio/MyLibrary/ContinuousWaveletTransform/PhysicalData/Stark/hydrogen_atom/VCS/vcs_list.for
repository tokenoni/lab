        program vcs_list
C ***
C ******************************************************************************
C ***
C ***	Fortran illustration to read VCS tables.
C ***
C ***	From A&AS paper `Extended VCS Stark broadening tables for hydrogen'
C ***	by Michael Lemke (michael@astro.as.utexas.edu,
C ***                     ai26@a400.sternwarte.uni-erlangen.de)
C ***
C ***	VCS tables are read from unit 5,
C ***	    result is written to unit 6,
C ***	            errors go to unit 0.
C ***
C ***	 2-AUG-1996 18:00:10.47		ML / Bamberg.  
C ***	11-NOV-1996 21:44:38.12		ML / Bamberg.  Increased pmp.
C ***
C ******************************************************************************
C ***
        implicit none

        integer pmne, pmp, pmt, mline
        parameter( pmne = 17,                ! max # of n_e
     >             pmt = 7,                  ! max # of T
     >             pmp = 65,                 ! max # of profile points
     >             mline = 21                ! max # of H lines
     >           )

        real svcs( pmt, pmne, 0:pmp, mline ) ! VCS profiles

        real log_alpha0(mline),              ! log Delta_alpha_min
     >       log_ne0(mline),                 ! log n_e_min
     >       log_t0(mline)                   ! log T_min

        real log_alpha_inc(mline),           ! Delta log Delta_alpha
     >       log_t_inc(mline),               ! Delta log n_e
     >       log_ne_inc(mline)               ! Delta log T

        integer nl(mline), nu(mline), mp(mline), mne(mline), mt(mline)
        character*1 null(pmt)

        integer i, j, k, line, nline, nnl, nnu
        real ne, f0


        real exp10, x
        exp10( x ) = exp( 2.30258509299405E0*x )
c
C       READ IN VCS ARRAYS
c
        read( 5, * ) nline
        if( nline .gt. mline ) then
           write( 6, * ) 'Table too big.  Not more than ', mline,
     >             ' lines.'
           call exit
        end if
        do i = 1, nline
           read( 5, * ) nl(i), nu(i),
     >                  log_alpha0(i), log_ne0(i), log_t0(i),
     >                  log_alpha_inc(i), log_ne_inc(i), log_t_inc(i),
     >                  mp(i), mne(i), mt(i)
           if( mp(i)  .gt. pmp .or.
     >         mne(i) .gt. pmne .or.
     >         mt(i)  .gt. pmt ) then
              write( 6, * ) 'Table too big in one of these:'
              write( 6, * ) 'mp:', mp(i), ' >', pmp
              write( 6, * ) 'mne:', mne(i), ' >', pmne
              write( 6, * ) 'mt:', mt(i), ' >', pmt
              call exit(1)
           end if
        end do
        do line = 1, nline
           read( 5, '(3x,i2,4x,i2)' ) nnl, nnu
           if( nnl .ne. nl(line) .or. nnu .ne. nu(line) ) then
              write( 0, * ) 'Inconsistency in table for', nl(line),
     >                      ' ->', nu(line)
              call exit(1)
           end if
           read( 5, * ) (( (svcs(i,j,k,line), k = 0, mp(line)),
     >                        i = 1, mt(line)),
     >                          j = 1, mne(line))
        end do

        do line = 1, nline
           write( 6, * ) 'nl = ', nl(line), ';  nu = ', nu(line)
           do j = 1, mne(line)
              ne = exp10( log_ne0(line) + log_ne_inc(line)*(j - 1) )
              write( 6, * ) ne, ' cm^-3'
              write( 6, '(a5,1x,a8,1x,7(:'' ('',i1,'')'',i6))' )
     >                'alpha', 'lambda',
     >                (abs(int(svcs(i+1,j,0,line))),
     >                 nint(exp10(log_t0(line) + log_t_inc(line)*i)),
     >                 i = 0, mt(line)-1 )
              f0 = 1.25e-9*ne**(2./3.)
              do k = 1, mp(line)
                x = log_alpha0(line) + log_alpha_inc(line)*(k - 1)
                do i = 1, mt(line)
                   null(i) = ' '
*                   if( svcs(i,j,1,line) - svcs(i,j,k,line) .le. 5. )
*     >                    null(i) = '*'
                end do
                write( 6, '(f5.1,1x,1p,e8.3e1,1x,0p,7(f9.5,a))' ) x,
     >               exp10(x)*f0,
     >               (real(svcs(i,j,k,line)), null(i), i = 1, mt(line) )
              end do
              write( 6, * )
           end do
        end do

        end
  
