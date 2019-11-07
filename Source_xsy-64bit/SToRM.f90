PROGRAM SToRM
  USE parameters
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!                _/_/_/  _/_/_/_/_/    _/_/    _/_/_/    _/      _/           !
!             _/            _/      _/    _/  _/    _/  _/_/  _/_/            !
!              _/_/        _/      _/    _/  _/_/_/    _/  _/  _/             !
!                 _/      _/      _/    _/  _/    _/  _/      _/              !
!          _/_/_/        _/        _/_/    _/    _/  _/      _/               !
!                                                                             !
!          Program SToRM  -  System for Transport and River Modeling          !
!                                                                             !
!                            Francisco J.M. Simões                            !
!                          National Research Program                          !
!                            U.S. Geological Survey                           !
!             USGS Geomorphology and Sediment Transport Laboratory            !
!                       4620 Technology Drive, Suite 400                      !
!                         Golden, CO 80403  -  U.S.A.                         !
!             303-278-7956, 303-279-4165 (fax), frsimoes@usgs.gov             !
!                                                                             !
!                                                                             !
!  This program and its source code were developed by the U.S. Geological     !
!  Survey (USGS).  The USGS does not guarantee the performance of this        !
!  program.  The USGS assumes no responsibility for the correct use of this   !
!  software or any of its derivatives, and makes no warranties concerning     !
!  the accuracy, completeness, reliability, usability, or suitability for     !
!  any particular purpose of this software.  The USGS shall not be liable     !
!  for any special, collateral, incidental, or consequential damages in       !
!  connection with the use of this software.                                  !
!                                                                             !
!  SToRM is an open source project, which means that the set of ASCII files   !
!  containing the complete program source codes are available to the public.  !
!  However, unauthorized changes to any part of the contents of those files,  !
!  no matter how small those changes may be, will invalidate and nullify any  !
!  claims that the resulting programs and source codes have any affinity,     !
!  association, relationship, or kinship of any type, direct or implicit,     !
!  with the original USGS program named SToRM.                                !
!                                                                             !
!  For a complete list of user's rights see the file Notice.txt.              !
!                                                                             !
!-----------------------------------------------------------------------------!
! Date of last change: 19 June 2018, F. Simoes

! Local variables.
  INTEGER :: elapsedt,i
  REAL(KIND=mp) :: r, ttime !rmcd mod 'ttime'
  CHARACTER (LEN=250) :: filename
  CHARACTER (LEN=160) :: title

! Print splash screen.
  PRINT *,''
  PRINT *,' _______ _______  _____   ______ _______'
  PRINT *,' |______    |    |     | |_____/ |  |  |'
  PRINT *,' ______|    |    |_____| |    \_ |  |  |'
  PRINT *,''
  PRINT *,' System for Transport and River Modeling'
  PRINT *,' U.S. Geological Survey - GSTL'
  PRINT *,' Golden, CO 80403 - USA'
  PRINT *,' Version 0.4  - 2018'
  PRINT *,''

  CALL read_data(filename,title)

  CALL solver(i,ttime,r,filename,title)

  CALL output_data(i,ttime,r,filename,title)
  !CALL output_dataV(i,r,title)

! Print elapsed time.  Non-standard Fortran.
  elapsedt = MCLOCK()
  WRITE (*,'(A,ES9.2,A)')"Elapsed time:", &
    REAL(elapsedt,mp)/1000.0_mp/60.0_mp," minutes."

  CALL byebye('')

END PROGRAM SToRM
