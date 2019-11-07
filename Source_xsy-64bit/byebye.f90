SUBROUTINE byebye(error_msg)
  USE geometry
  USE dep_vars
  USE vbc_arrays
  USE io
  USE RKparams
  USE options
  USE write_cgns_cblock
  IMPLICIT NONE

!-----------------------------------------------------------------------------!
!                                                                             !
!  This subroutine is used to stop SToRM.  The argument "error_msg" is        !
!  supposed to contain an error message that is printed in the standard       !
!  console, but it may left blank.  All the main SToRM arrays are             !
!  DEALLOCATEd upon exit.                                                     !
!                                                                             !
!  F. Simoes, November 2005                                                   !
!  Last updated (mm-dd-yyyy): 06-02-2017 by F. Simoes                         !
!                                                                             !
!-----------------------------------------------------------------------------!

! Dummy argument.
  CHARACTER(LEN=*), INTENT(IN) :: error_msg

! Local variables.
  INTEGER :: ierror

  IF (LEN_TRIM(error_msg) > 0) WRITE (*,'(A)') TRIM(error_msg)

  ! Temporary arrays from 'tecplot'.
  IF (ALLOCATED(cd_tcp)) DEALLOCATE(cd_tcp)
  IF (ALLOCATED(h_tcp)) DEALLOCATE(h_tcp)
  IF (ALLOCATED(u_tcp)) DEALLOCATE(u_tcp)
  IF (ALLOCATED(v_tcp)) DEALLOCATE(v_tcp)
  IF (ALLOCATED(z_tcp)) DEALLOCATE(z_tcp)

  ! From 'find_edges'.
  IF (ALLOCATED(edges)) DEALLOCATE(edges)

  ! From 'ibc_arrays'.
  IF (ALLOCATED(pnormals)) DEALLOCATE(pnormals)
  IF (ALLOCATED(l_inflow)) DEALLOCATE(l_inflow)
  IF (ALLOCATED(h_inflow)) DEALLOCATE(h_inflow)
  IF (ALLOCATED(a_inflow)) DEALLOCATE(a_inflow)
  IF (ALLOCATED(u_inflow)) DEALLOCATE(u_inflow)
  IF (ALLOCATED(v_inflow)) DEALLOCATE(v_inflow)
  IF (ALLOCATED(hnormals)) DEALLOCATE(hnormals)
  IF (ALLOCATED(wtangs)) DEALLOCATE(wtangs)
  IF (ALLOCATED(wtangs1)) DEALLOCATE(wtangs1)
  IF (ALLOCATED(hbar)) DEALLOCATE(hbar)
  IF (ALLOCATED(qbar)) DEALLOCATE(qbar)
  IF (ALLOCATED(bcedges)) DEALLOCATE(bcedges)
  IF (ALLOCATED(bcedgesinflow)) DEALLOCATE(bcedgesinflow)
  IF (ALLOCATED(bcedgesoutflow)) DEALLOCATE(bcedgesoutflow)

  ! From 'ioflow'.
  IF (ALLOCATED(n_qin)) DEALLOCATE(n_qin)
  IF (ALLOCATED(qin_nodes)) DEALLOCATE(qin_nodes)
  IF (ALLOCATED(n_hbc)) DEALLOCATE(n_hbc)
  IF (ALLOCATED(hbc_nodes)) DEALLOCATE(hbc_nodes)
  IF (ALLOCATED(hbc)) DEALLOCATE(hbc)
  IF (ALLOCATED(htype)) DEALLOCATE(htype)
  IF (ALLOCATED(timeloch)) DEALLOCATE(timeloch)
  IF (ALLOCATED(timelocq)) DEALLOCATE(timelocq)
  IF (ALLOCATED(timeserq)) DEALLOCATE(timeserq)
  IF (ALLOCATED(timeserq2)) DEALLOCATE(timeserq2)
  IF (ALLOCATED(timeserh)) DEALLOCATE(timeserh)
  IF (ALLOCATED(cqin)) DEALLOCATE(cqin)
  IF (ALLOCATED(cqin2)) DEALLOCATE(cqin2)
  IF (ALLOCATED(cqin3)) DEALLOCATE(cqin3)
  IF (ALLOCATED(cqout)) DEALLOCATE(cqout)
  IF (ALLOCATED(cqout2)) DEALLOCATE(cqout2)
  IF (ALLOCATED(cqout3)) DEALLOCATE(cqout3)

  IF (ALLOCATED(n_channel)) DEALLOCATE(n_channel)
  IF (ALLOCATED(channel)) DEALLOCATE(channel)
  IF (ALLOCATED(chanedges)) DEALLOCATE(chanedges)
  IF (ALLOCATED(chantrigs)) DEALLOCATE(chantrigs)
  IF (ALLOCATED(chant)) DEALLOCATE(chant)

  ! Culvert variables allocated in 'ioculvert'.
  IF (ALLOCATED(cvptoin)) DEALLOCATE(cvptoin)
  IF (ALLOCATED(cvptoout)) DEALLOCATE(cvptoout)
  IF (ALLOCATED(ncvtable)) DEALLOCATE(ncvtable)
  IF (ALLOCATED(cvhead)) DEALLOCATE(cvhead)
  IF (ALLOCATED(cvq)) DEALLOCATE(cvq)
  IF (ALLOCATED(cvtrigin)) DEALLOCATE(cvtrigin)
  IF (ALLOCATED(cvtrigout)) DEALLOCATE(cvtrigout)
  IF (ALLOCATED(cvrin)) DEALLOCATE(cvrin)
  IF (ALLOCATED(cvsrc)) DEALLOCATE(cvsrc)

  ! From 'node_db'.
  IF (ALLOCATED(cv_area)) DEALLOCATE(cv_area)
  IF (ALLOCATED(cv_perim)) DEALLOCATE(cv_perim)
  IF (ALLOCATED(n2t)) DEALLOCATE(n2t)
  IF (ALLOCATED(n2t2)) DEALLOCATE(n2t2)
  IF (ALLOCATED(n2n)) DEALLOCATE(n2n)

  ! From 'elemnt_db'.
  IF (ALLOCATED(t2t)) DEALLOCATE(t2t)
  IF (ALLOCATED(t2tHD)) DEALLOCATE(t2tHD)
  IF (ALLOCATED(t2t3)) DEALLOCATE(t2t3)

  ! From 'solver'.
  IF (ALLOCATED(residual)) DEALLOCATE(residual)
  IF (ALLOCATED(phi)) DEALLOCATE(phi)
  IF (ALLOCATED(phi0)) DEALLOCATE(phi0)
  IF (ALLOCATED(flux)) DEALLOCATE(flux)
  IF (ALLOCATED(source)) DEALLOCATE(source)
  IF (ALLOCATED(e_est)) DEALLOCATE(e_est)
  IF (ALLOCATED(u_avg)) DEALLOCATE(u_avg)
  IF (ALLOCATED(u_mag)) DEALLOCATE(u_mag)
  IF (ALLOCATED(wavec)) DEALLOCATE(wavec)
  IF (ALLOCATED(delh)) DEALLOCATE(delh)
  IF (ALLOCATED(delu)) DEALLOCATE(delu)
  IF (ALLOCATED(delv)) DEALLOCATE(delv)
  IF (ALLOCATED(delz)) DEALLOCATE(delz)
  IF (ALLOCATED(delz2)) DEALLOCATE(delz2)
  IF (ALLOCATED(z2)) DEALLOCATE(z2)
  IF (ALLOCATED(rkalpha)) DEALLOCATE(rkalpha)
  IF (ALLOCATED(rkbeta)) DEALLOCATE(rkbeta)
  IF (ALLOCATED(rkalpha)) DEALLOCATE(rkalpha)
  IF (ALLOCATED(rkh)) DEALLOCATE(rkh)
  IF (ALLOCATED(rku)) DEALLOCATE(rku)
  IF (ALLOCATED(rkv)) DEALLOCATE(rkv)
  IF (ALLOCATED(bdrycell)) DEALLOCATE(bdrycell)
  IF (ALLOCATED(drycell)) DEALLOCATE(drycell)
  IF (ALLOCATED(partdry)) DEALLOCATE(partdry)
  IF (ALLOCATED(csortedz)) DEALLOCATE(csortedz)
  IF (ALLOCATED(areap)) DEALLOCATE(areap)

  ! From 'tecplot'.
  IF (ALLOCATED(nodes)) DEALLOCATE(nodes)
  IF (ALLOCATED(grid)) DEALLOCATE(grid)
  IF (ALLOCATED(h)) DEALLOCATE(h)
  IF (ALLOCATED(hvtx)) DEALLOCATE(hvtx)
  IF (ALLOCATED(zeta)) DEALLOCATE(zeta)
  IF (ALLOCATED(zetavtx)) DEALLOCATE(zetavtx)
  IF (ALLOCATED(u)) DEALLOCATE(u)
  IF (ALLOCATED(uvtx)) DEALLOCATE(uvtx)
  IF (ALLOCATED(v)) DEALLOCATE(v)
  IF (ALLOCATED(vvtx)) DEALLOCATE(vvtx)
  IF (ALLOCATED(z)) DEALLOCATE(z)
  IF (ALLOCATED(zvtx)) DEALLOCATE(zvtx)
  IF (ALLOCATED(cd)) DEALLOCATE(cd)
  IF (ALLOCATED(rough)) DEALLOCATE(rough)
  IF (ALLOCATED(rcoef1)) DEALLOCATE(rcoef1)
  IF (ALLOCATED(rcoef2)) DEALLOCATE(rcoef2)
  IF (ALLOCATED(zbx)) DEALLOCATE(zbx)
  IF (ALLOCATED(zby)) DEALLOCATE(zby)
  IF (ALLOCATED(wcoef1)) DEALLOCATE(wcoef1)
  IF (ALLOCATED(wcoef2)) DEALLOCATE(wcoef2)

  ! From 'walls'.
  IF (ALLOCATED(walledg)) DEALLOCATE(walledg)
  IF (ALLOCATED(flowedg)) DEALLOCATE(flowedg)
  IF (ALLOCATED(walledg)) DEALLOCATE(walledg1)
  IF (ALLOCATED(flowedg)) DEALLOCATE(flowedg1)
  IF (ALLOCATED(bpolygon)) DEALLOCATE(bpolygon)
  IF (ALLOCATED(wall_pts)) DEALLOCATE(wall_pts)

  ! From 'read_options'.
  IF (ALLOCATED(iout)) DEALLOCATE(iout)
  IF (ALLOCATED(probepts)) DEALLOCATE(probepts)

  ! From 'lsq_setup' and 'hdlsq_setup'.
  IF (ALLOCATED(ai)) DEALLOCATE(ai)
  IF (ALLOCATED(bi)) DEALLOCATE(bi)
  IF (ALLOCATED(ci)) DEALLOCATE(ci)
  IF (ALLOCATED(det)) DEALLOCATE(det)
  IF (ALLOCATED(dx)) DEALLOCATE(dx)
  IF (ALLOCATED(dy)) DEALLOCATE(dy)
  IF (ALLOCATED(lsqweight)) DEALLOCATE(lsqweight)

  ! From 'c2v_setup'.
  IF (ALLOCATED(c2v_weights)) DEALLOCATE(c2v_weights)
  IF (ALLOCATED(c2v_weights_id)) DEALLOCATE(c2v_weights_id)

  IF (ALLOCATED(RCNumVals)) DEALLOCATE(RCNumVals)
  IF (ALLOCATED(TSNumVals)) DEALLOCATE(TSNumVals)
  IF (ALLOCATED(TSIndexUsed)) DEALLOCATE(TSIndexUsed)
  IF (ALLOCATED(RCVals)) DEALLOCATE(RCVals)
  IF (ALLOCATED(TSVals)) DEALLOCATE(TSVals)
  IF (ALLOCATED(TSType)) DEALLOCATE(TSType)
  IF (ALLOCATED(RCType)) DEALLOCATE(RCType)

  ! CGNS.
  IF(ALLOCATED(SolNames)) DEALLOCATE(SolNames)
  IF(ALLOCATED(CellSolNames)) DEALLOCATE(CellSolNames)
  IF(ALLOCATED(TimeIncrements)) DEALLOCATE(TimeIncrements)
  IF(ALLOCATED(CellTimeIncrements)) DEALLOCATE(CellTimeIncrements)
  IF (cgns) CALL cg_close_f(findex,ierror)  ! close cgns file.

  ! From 'write_cgns'.
  IF (ALLOCATED(ibc)) DEALLOCATE(ibc)
  IF (ALLOCATED(ibcell)) DEALLOCATE(ibcell)
  IF (ALLOCATED(wetv)) DEALLOCATE(wetv)
  IF (ALLOCATED(deltaubx)) DEALLOCATE(deltaubx)
  IF (ALLOCATED(deltauby)) DEALLOCATE(deltauby)
  IF (ALLOCATED(depth)) DEALLOCATE(depth)
  IF (ALLOCATED(dShields)) DEALLOCATE(dShields)
  IF (ALLOCATED(dtau)) DEALLOCATE(dtau)
  IF (ALLOCATED(dtauv)) DEALLOCATE(dtauv)
  IF (ALLOCATED(tau)) DEALLOCATE(tau)
  IF (ALLOCATED(taucell)) DEALLOCATE(taucell)
  IF (ALLOCATED(taubx)) DEALLOCATE(taubx)
  IF (ALLOCATED(tauby)) DEALLOCATE(tauby)
  IF (ALLOCATED(tauv)) DEALLOCATE(tauv)
  IF (ALLOCATED(wselev)) DEALLOCATE(wselev)

  !PAUSE  ! For cgns output console.
  STOP

END SUBROUTINE byebye
