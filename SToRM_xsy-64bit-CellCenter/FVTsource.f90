SUBROUTINE add_array_ptoH(i,n,iarray)
IMPLICIT NONE
INTEGER, INTENT(IN) :: i
INTEGER, INTENT(INOUT) :: n
INTEGER, DIMENSION(n+1), INTENT(INOUT) :: iarray
INTEGER :: j
DO j = n+1,2,-1
iarray(j) = iarray(j-1)
END DO
iarray(1) = i
n = n + 1
END SUBROUTINE add_array_ptoH
SUBROUTINE add_array_ptoT(i,n,iarray)
IMPLICIT NONE
INTEGER, INTENT(IN) :: i
INTEGER, INTENT(INOUT) :: n
INTEGER, DIMENSION(n+1), INTENT(INOUT) :: iarray
n = n + 1
iarray(n) = i
END SUBROUTINE add_array_ptoT
SUBROUTINE adjust_bed_elevation(zchange)
USE parameters
USE dep_vars
USE geometry
USE options
USE constants
IMPLICIT NONE
LOGICAL, INTENT(OUT) :: zchange
INTEGER :: i,j
zchange = .FALSE.
DO i = 1,n_elems
IF (h(i) < h_dry) CYCLE
IF (zeta(i) < zvtx(csortedz(i)%p(3))) THEN
zchange = .TRUE.
DO j = 1,3
zvtx(grid(i)%vertex(j)) = MIN(zeta(i),zvtx(grid(i)%vertex(j)))
END DO
z(i) = (zvtx(grid(i)%vertex(1)) + zvtx(grid(i)%vertex(2)) + &
zvtx(grid(i)%vertex(3)))*one_third
h(i) = zeta(i) - z(i)
DO j = 1,3
hvtx(grid(i)%vertex(j)) = zetavtx(grid(i)%vertex(j)) - &
zvtx(grid(i)%vertex(j))
END DO
END IF
END DO
END SUBROUTINE adjust_bed_elevation
SUBROUTINE adjust_water_depth
USE parameters
USE dep_vars
USE geometry
USE constants
IMPLICIT NONE
REAL (KIND=mp) :: delta
INTEGER :: i
DO i = 1,n_elems
delta = zstore(i) - z(i)
h(i) = h(i) - delta
END DO
WHERE (h < zero) h = zero
z = zstore
zvtx = zvtxstore
END SUBROUTINE adjust_water_depth
SUBROUTINE adjustZb(flag,e,gradzb,gradzb2)
USE parameters
USE constants
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER, INTENT(IN) :: e
LOGICAL, INTENT(IN) :: flag
TYPE(vector), INTENT(OUT) :: gradzb,gradzb2
INTEGER :: i,j
REAL (KIND=mp) :: a,b(3),c(3),hp(3),x(3),y(3),zp(3)
LOGICAL :: adjust
adjust = .FALSE.
DO i = 1,3
j = grid(e)%vertex(i)
hp(i) = h(e) + delh(e)%x*(nodes(j)%x - grid(e)%xc) + &
delh(e)%y*(nodes(j)%y - grid(e)%yc)
IF (hp(i) < zero) adjust = .TRUE.
END DO
IF (adjust) THEN
DO i = 1,3
j = grid(e)%vertex(i)
hp(i) = MIN(zero,hp(i))
zp(i) = zvtx(j) + hp(i)
x(i) = nodes(j)%x
y(i) = nodes(j)%y
END DO
b(1) = y(2) - y(3) ; c(1) = x(3) - x(2)
b(2) = y(3) - y(1) ; c(2) = x(1) - x(3)
b(3) = y(1) - y(2) ; c(3) = x(2) - x(1)
a = half/grid(e)%area
gradzb%x = a*(b(1)*zp(1) + b(2)*zp(2) + b(3)*zp(3))
gradzb%y = a*(c(1)*zp(1) + c(2)*zp(2) + c(3)*zp(3))
ELSE
gradzb%x = delz(e)%x
gradzb%y = delz(e)%y
END IF
IF (flag) THEN
gradzb2%x = 2.0_mp*z(e)*gradzb%x
gradzb2%y = 2.0_mp*z(e)*gradzb%y
END IF
END SUBROUTINE adjustZb
FUNCTION azimuth2angle(x)
USE parameters
USE constants
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: x
REAL (KIND=mp) :: azimuth2angle
REAL (KIND=mp) :: a
IF (x < zero .OR. x > 360_mp) THEN
PRINT *,"Wrong argument in wind direction: value < 0 or > 360 degrees."
CALL byebye("Program SToRM stopped.")
END IF
a = (360_mp - x) + 90_mp
IF (a < zero) THEN
a = a + 360_mp
ELSE IF (a >= 360_mp) THEN
a = a - 360_mp
END IF
azimuth2angle = a*pi/180.0_mp
END FUNCTION azimuth2angle
SUBROUTINE bed_slope_DFB
USE parameters
USE constants
USE geometry
USE dep_vars
USE options
IMPLICIT NONE
REAL (KIND=mp) :: eta,hr,ksign,s0(3)
REAL (KIND=mp), PARAMETER :: three = 3.0_mp
INTEGER :: i,j,k,p1,p2
SELECT CASE (opt_solver)
CASE (2)
DO i = 1,n_elems
IF (drycell(i)) CYCLE
eta = z(i) + h(i)
s0 = zero
DO j = 1,3
k = grid(i)%edge(j)
ksign = SIGN(one,REAL(k,mp))
k = ABS(k)
p1 = edges(k)%p(1)
p2 = edges(k)%p(2)
hr = eta - half*(zvtx(p1) + zvtx(p2))
IF (hr < h_dry) CYCLE
s0(2) = s0(2) + ksign*hr*hr*edges(k)%normal(1)*edges(k)%length
s0(3) = s0(3) + ksign*hr*hr*edges(k)%normal(2)*edges(k)%length
END DO
phi(i,2) = phi(i,2) + half*g*s0(2)
phi(i,3) = phi(i,3) + half*g*s0(3)
END DO
CASE (3)
DO i = 1,n_elems
IF (drycell(i)) CYCLE
s0 = zero
DO j = 1,3
k = grid(i)%edge(j)
ksign = SIGN(one,REAL(k,mp))
k = ABS(k)
p1 = edges(k)%p(1)
p2 = edges(k)%p(2)
hr = zeta(i) - half*(zvtx(p1) + zvtx(p2))
IF (hr < h_dry) CYCLE
s0(2) = s0(2) + ksign*hr*hr*edges(k)%normal(1)*edges(k)%length
s0(3) = s0(3) + ksign*hr*hr*edges(k)%normal(2)*edges(k)%length
END DO
phi(i,2) = phi(i,2) + half*g*s0(2)
phi(i,3) = phi(i,3) + half*g*s0(3)
END DO
CASE DEFAULT
PRINT *,"ERROR: invalid option in bed_slope_DFB."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE bed_slope_DFB
SUBROUTINE bed_slope_ETA(solver)
USE parameters
USE constants
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER, INTENT(IN) :: solver
REAL (KIND=mp) :: eta
TYPE(vector) :: gradzb,gradzb2
INTEGER :: i
SELECT CASE (solver)
CASE (2)
DO i = 1,n_elems
IF (drycell(i)) CYCLE
eta = z(i) + h(i)
phi(i,2) = phi(i,2) - g*(eta*delz(i)%x - half*delz2(i)%x)*grid(i)%area
phi(i,3) = phi(i,3) - g*(eta*delz(i)%y - half*delz2(i)%y)*grid(i)%area
END DO
CASE (3)
DO i = 1,n_elems
IF (drycell(i)) CYCLE
phi(i,2) = phi(i,2) - g*(zeta(i)*delz(i)%x - half*delz2(i)%x)*grid(i)%area
phi(i,3) = phi(i,3) - g*(zeta(i)*delz(i)%y - half*delz2(i)%y)*grid(i)%area
END DO
CASE DEFAULT
PRINT *,"ERROR: invalid option in bed_slope_ETA."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE bed_slope_ETA
SUBROUTINE bed_slope_STD
USE parameters
USE constants
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER :: i,j,k
TYPE(vector) :: dummy,gradzb
LOGICAL :: set_grad
DO i = 1,n_elems
IF (drycell(i)) CYCLE
CALL adjustZb(.FALSE.,i,gradzb,dummy)
phi(i,2) = phi(i,2) - g*h(i)*gradzb%x*grid(i)%area
phi(i,3) = phi(i,3) - g*h(i)*gradzb%y*grid(i)%area
END DO
END SUBROUTINE bed_slope_STD
SUBROUTINE c2v_setup(type)
USE parameters
USE geometry
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: type
INTEGER :: i,ierror,j,k,l
REAL (KIND=mp) :: d,ixx,ixy,iyy,lx,ly,rx,ry,tdist,xk,yk
LOGICAL :: clipping
k = SIZE(n2t2,2)
ALLOCATE(c2v_weights(n_pts,k-1),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(8*n_pts*(k - 1))
SELECT CASE (type)
CASE (1)
CASE (2)
DO i = 1,n_pts
rx = zero
ry = zero
ixx = zero
ixy = zero
iyy = zero
k = n2t2(i,1)
DO j = 2,k+1
l = n2t2(i,j)
xk = grid(l)%xc - nodes(i)%x
yk = grid(l)%yc - nodes(i)%y
rx = rx + xk
ry = ry + yk
ixx = ixx + xk*xk
ixy = ixy + xk*yk
iyy = iyy + yk*yk
END DO
d = ixx*iyy - ixy*ixy
lx = (ixy*ry - iyy*rx)/d
ly = (ixy*rx - ixx*ry)/d
DO j = 2,k+1
l = n2t2(i,j)
xk = grid(l)%xc - nodes(i)%x
yk = grid(l)%yc - nodes(i)%y
c2v_weights(i,j-1) = (one + lx*xk + ly*yk)/(REAL(k,mp) + lx*rx + ly*ry)
END DO
END DO
CASE (3)
ALLOCATE(c2v_weights_id(n_pts,k-1),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(8*n_pts*(k - 1))
DO i = 1,n_pts
tdist = zero
DO j = 2,n2t(i,1)+1
k = n2t(i,j)
xk = grid(k)%xc - nodes(i)%x
yk = grid(k)%yc - nodes(i)%y
c2v_weights(i,j-1) = one/SQRT(xk*xk + yk*yk)
c2v_weights_id(i,j-1) = c2v_weights(i,j-1)
tdist = tdist + c2v_weights(i,j-1)
END DO
DO j = 2,n2t(i,1)+1
c2v_weights(i,j-1) = c2v_weights(i,j-1)/tdist
END DO
END DO
CASE DEFAULT
PRINT *,"ERROR: invalid option in c2v_setup."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE c2v_setup
SUBROUTINE cde_t2t(t2tw,i,dim1,dim2)
IMPLICIT NONE
INTEGER, INTENT(IN) :: i,dim1,dim2
INTEGER, INTENT(INOUT), DIMENSION(dim1,dim2) :: t2tw
INTEGER :: j,k,l,m,n,nn
n = t2tw(i,1)
l = 1
DO
j = t2tw(i,l+1)
DO m = l+1,n
k = t2tw(i,m+1)
IF (j == k) THEN
n = n - 1
DO nn = m,n
t2tw(i,nn+1) = t2tw(i,nn+2)
END DO
t2tw(i,1) = n
END IF
END DO
l = l + 1
IF (l > n) EXIT
END DO
END SUBROUTINE cde_t2t
SUBROUTINE chan_arrays
USE parameters
USE constants
USE geometry
USE dep_vars
USE memory
IMPLICIT NONE
INTEGER :: d1,d2,e,i,i0,ierror,j,t1,t2
INTEGER, EXTERNAL :: is_edge
REAL (KIND=mp) :: length,rlength
d1 = SIZE(channel,1)
d2 = SIZE(channel,2)
IF (d2 /= n_chanbdr) THEN
PRINT *,"ERROR in 1D channel arrays:,"
PRINT *,"Dimension does not match count."
CALL byebye('Program SToRM stopped.')
END IF
ALLOCATE(chanedges(0:d1,d2),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(4*(d1 + 1)*d2)
length = 0
DO i0 = 1,n_chanbdr
chanedges(0,i0) = n_channel(i0) - 1
DO i = 1,n_channel(i0)-1
e = is_edge(channel(i,i0),channel(i+1,i0))
IF (e == 0) THEN
PRINT *,"ERROR in 1D channel node assignement,"
PRINT *,"nodes",channel(i,i0),channel(i+1,i0)
PRINT *,"These nodes are not connected by an edge."
PRINT *,"Please revise your channel link data."
CALL byebye('Program SToRM stopped.')
ELSE IF (e < 0) THEN
PRINT *,"ERROR in 1D channel node assignement,"
PRINT *,"nodes",channel(i,i0)," or ",channel(i+1,i0),"(or both)."
PRINT *,"Nodes out of range.  Please revise your channel link data."
CALL byebye('Program SToRM stopped.')
END IF
chanedges(i,i0) = e
length = length + edges(e)%length
END DO
END DO
ALLOCATE(chantrigs(0:d1,d2),chant(d1,d2),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(4*(d1 + 1)*d2 + 8*d1*d2)
DO j = 1,n_chanbdr
chantrigs(0,j) = chanedges(0,j)
rlength = zero
DO i = 1,chanedges(0,j)
t1 = edges(chanedges(i,j))%e(1)
t2 = edges(chanedges(i,j))%e(2)
IF (t1 < 0) THEN
chantrigs(i,j) = t2
ELSE IF (t2 < 0) THEN
chantrigs(i,j) = t1
ELSE
IF (z(t1) < z(t2)) THEN
chantrigs(i,j) = t1
ELSE
chantrigs(i,j) = t2
END IF
END IF
rlength = rlength + half*(edges(chanedges(i,j))%length)
IF (i > 1) rlength = rlength + half*(edges(chanedges(i-1,j))%length)
chant(i,j) = rlength/length
END DO
END DO
END SUBROUTINE chan_arrays
SUBROUTINE correct_vtx(depth,wselev,wetv,n)
USE parameters
USE constants
USE geometry
USE dep_vars
USE options
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
REAL(8), INTENT(OUT) :: depth(n),wselev(n)
LOGICAL, INTENT(OUT) :: wetv(n)
INTEGER :: counter,i,j,k
REAL (KIND=mp) :: a,id,vmag,zetaref
DO k = 1,n_pts
depth(k) = hvtx(k)
wselev(k) = zetavtx(k)
END DO
wetv = .FALSE.
DO k = 1,n_pts
vmag = sqrt((uvtx(k)*uvtx(k))+(vvtx(k)*vvtx(k)))
IF (hvtx(k) > h_wet) THEN
wetv(k) = .TRUE.
ELSE IF (hvtx(k) > h_dry .AND. vmag > 1.0e-5) THEN
wetv(k) = .TRUE.
END IF
END DO
DO i = 1,n_pts
id = zero
a = zero
counter = zero
DO j = 2,n2t(i,1)+1
k = n2t(i,j)
IF (drycell(k)) CYCLE
IF (zeta(k) < zvtx(csortedz(k)%p(3))) CYCLE
counter = counter + 1
id = id + c2v_weights_id(i,j-1)
a = a + c2v_weights_id(i,j-1)*zeta(k)
END DO
IF (counter > 0) THEN
wselev(i) = a/id
ELSE
wselev(i) = zvtx(i)
END IF
END DO
depth = wselev - zvtx
DO i = 1,n_pts
zetaref = zvtx(i)
DO j = 2,n2t(i,1) + 1
k = n2t(i,j)
IF (drycell(k)) CYCLE
zetaref = MAX(zeta(k),zetaref)
END DO
wselev(i) = MIN(wselev(i),zetaref)
depth(i) = wselev(i) - zvtx(i)
IF (depth(i) > h_wet) THEN
wetv(i) = .TRUE.
ELSE
wetv(i) = .FALSE.
END IF
END DO
END SUBROUTINE correct_vtx
SUBROUTINE ctr2vtx(type,carray,varray)
USE parameters
USE constants
USE geometry
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: carray(n_elems)
REAL (KIND=mp), INTENT(OUT) :: varray(n_pts)
INTEGER, INTENT(IN) :: type
INTEGER :: i,j,k
REAL (KIND=mp) :: a,tarea,weight
SELECT CASE (type)
CASE (1)
DO i = 1,n_pts
tarea = zero
a = zero
DO j = 2,n2t(i,1)+1
k = n2t(i,j)
weight = grid(k)%area
tarea = tarea + weight
a = a + weight*carray(k)
END DO
varray(i) = a/tarea
END DO
CASE (2)
DO i = 1,n_pts
a = zero
DO j = 2,n2t2(i,1)+1
k = n2t2(i,j)
a = a + c2v_weights(i,j-1)*carray(k)
END DO
varray(i) = a
END DO
CASE (3)
DO i = 1,n_pts
a = zero
DO j = 2,n2t(i,1)+1
k = n2t(i,j)
a = a + c2v_weights(i,j-1)*carray(k)
END DO
varray(i) = a
END DO
CASE DEFAULT
PRINT *,"ERROR: invalid option in ctr2vtx."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE ctr2vtx
SUBROUTINE ctr2vtxWSE(type,carray,varray)
USE parameters
USE constants
USE geometry
USE dep_vars
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: carray(n_elems)
REAL (KIND=mp), INTENT(OUT) :: varray(n_pts)
INTEGER, INTENT(IN) :: type
INTEGER :: counter,i,j,k
REAL (KIND=mp) :: a,id,tarea,weight
SELECT CASE (type)
CASE (1)
DO i = 1,n_pts
tarea = zero
a = zero
counter = 0
DO j = 2,n2t(i,1)+1
k = n2t(i,j)
IF (drycell(k)) CYCLE
counter = counter + 1
weight = grid(k)%area
tarea = tarea + weight
a = a + weight*carray(k)
END DO
IF (counter > 0) THEN
varray(i) = a/tarea
ELSE
varray(i) = zvtx(i)
END IF
END DO
CASE (2)
DO i = 1,n_pts
a = zero
DO j = 2,n2t2(i,1)+1
k = n2t2(i,j)
a = a + c2v_weights(i,j-1)*carray(k)
END DO
varray(i) = a
END DO
CASE (3)
DO i = 1,n_pts
id = zero
a = zero
counter = zero
DO j = 2,n2t(i,1)+1
k = n2t(i,j)
IF (drycell(k)) CYCLE
counter = counter + 1
id = id + c2v_weights_id(i,j-1)
a = a + c2v_weights_id(i,j-1)*carray(k)
END DO
IF (counter > 0) THEN
varray(i) = a /id
ELSE
varray(i) = zvtx(i)
END IF
END DO
CASE DEFAULT
PRINT *,"ERROR: invalid option in ctr2vtx."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE ctr2vtxWSE
SUBROUTINE cvdisch
USE dep_vars
USE constants
USE options
USE geometry
IMPLICIT NONE
INTEGER :: i,j
REAL (KIND=mp) :: head,t
cvsrc = zero
DO i = 1,nculvert
IF (drycell(cvtrigin(i))) CYCLE
head = zeta(cvtrigin(i)) + (delh(cvtrigin(i))%x + &
delz(cvtrigin(i))%x)*cvrin(1,i) + (delh(cvtrigin(i))%y + &
delz(cvtrigin(i))%y)*cvrin(2,i)
head = head - cvptoin(3,i)
IF (head < h_dry) CYCLE
IF (head < cvhead(1,i)) THEN
cvsrc(i) = cvq(1,i)
CYCLE
ELSE IF (head >= cvhead(ncvtable(i),i)) THEN
cvsrc(i) = cvq(ncvtable(i),i)
CYCLE
ELSE
DO j = 1,ncvtable(i)-1
IF (head < cvhead(j+1,i)) THEN
t = (head - cvhead(j,i))/(cvhead(j+1,i) - cvhead(j,i))
cvsrc(i) = cvq(j+1,i)*t + cvq(j,i)*(one - t)
EXIT
END IF
END DO
END IF
END DO
END SUBROUTINE cvdisch
SUBROUTINE cvprep
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER :: i,j
LOGICAL, EXTERNAL :: pto_in_triangle
DO i = 1,nculvert
cvtrigin(i) = 0
DO j = 1,n_elems
IF (pto_in_triangle(grid(j),cvptoin(1,i),cvptoin(2,i))) THEN
cvtrigin(i) = j
EXIT
END IF
END DO
IF (cvtrigin(i) == 0) THEN
WRITE(*,'(" ERROR: inlet of culvert ",I5, &
" not in computational domain.")') i
CALL byebye(' Program SToRM stopped.')
END IF
cvrin(1,i) = cvptoin(1,i) - grid(j)%xc
cvrin(2,i) = cvptoin(2,i) - grid(j)%yc
cvtrigout(i) = 0
DO j = 1,n_elems
IF (pto_in_triangle(grid(j),cvptoout(1,i),cvptoout(2,i))) THEN
cvtrigout(i) = j
EXIT
END IF
END DO
IF (cvtrigout(i) == 0) THEN
WRITE(*,'(" ERROR: outlet of culvert ",I5, &
" not in computational domain.")') i
CALL byebye(' Program SToRM stopped.')
END IF
END DO
END SUBROUTINE cvprep
SUBROUTINE dirichlet(i,rkstep)
USE geometry
USE constants
USE dep_vars
USE RKparams
USE options
IMPLICIT NONE
INTEGER, INTENT(IN) :: i,rkstep
INTEGER :: k,l,total
REAL (KIND=mp) :: utotal,vtotal,ztotal
total = 0
utotal = zero
vtotal = zero
ztotal = zero
DO k = 1,3
l = t2t3(k,i)
IF (l < 1) CYCLE
IF (.NOT. drycell(l)) THEN
total = total + 1
utotal = utotal + u(l)
vtotal = vtotal + v(l)
ztotal = ztotal + zeta(l)
END IF
END DO
IF (total > 0) THEN
zeta(i) = ztotal/REAL(total,mp)
h(i) = zeta(i) - z(i)
IF (h(i) > h_dry) THEN
u(i) = utotal/REAL(total,mp)
v(i) = vtotal/REAL(total,mp)
rku(i,rkstep) = h(i)*u(i)
rkv(i,rkstep) = h(i)*v(i)
rkh(i,rkstep) = h(i)
ELSE
u(i) = zero
v(i) = zero
h(i) = zero
zeta(i) = z(i)
rku(i,rkstep) = zero
rkv(i,rkstep) = zero
rkh(i,rkstep) = zero
drycell(i) = .TRUE.
bdrycell(i) = .FALSE.
END IF
ELSE
u(i) = zero
v(i) = zero
h(i) = zero
zeta(i) = z(i)
rku(i,rkstep) = zero
rkv(i,rkstep) = zero
rkh(i,rkstep) = zero
drycell(i) = .TRUE.
bdrycell(i) = .FALSE.
END IF
END SUBROUTINE dirichlet
SUBROUTINE eddy_viscosity
USE parameters
USE geometry
USE dep_vars
USE constants
USE options
IMPLICIT NONE
INTEGER :: i,j,k
REAL (KIND=mp) :: a,weight,tarea,ustar
rough = -rough
DO i = 1,n_pts
IF (hvtx(i) < h_dry) CYCLE
tarea = zero
a = zero
DO j = 2,n2t(i,1)+1
k = n2t(i,j)
IF (h(k) > h_dry) THEN
weight = grid(k)%area
tarea = tarea + weight
a = a + weight*ABS(rough(i))*u_mag(i)
END IF
END DO
ustar = SQRT(a/tarea)
ev(i) = omega*hvtx(i)*ustar
END DO
END SUBROUTINE eddy_viscosity
SUBROUTINE enforce_bcs(h0,u0,v0,z0,ndim)
USE constants
USE geometry
USE dep_vars
USE options
USE vbc_arrays
IMPLICIT NONE
INTEGER, INTENT(IN) :: ndim
REAL (KIND=mp), INTENT(INOUT), DIMENSION(ndim) :: h0,u0,v0,z0
INTEGER :: i,i0,j
REAL (KIND=mp) :: dotprod,hnot,incr,q
REAL (KIND=mp), EXTERNAL :: vbc_by_h
DO i0 = 1,n_outflowbdr
SELECT CASE (htype(i0))
CASE (1)
DO i = 1,n_hbc(i0)
j = hbc_nodes(i,i0)
hnot = MAX(zero,hbc(i0) - z0(j))
incr = alpha*(hnot - h0(j))
IF (ABS(incr) < kappa) incr = hnot - h0(j)
h0(j) = h0(j) - incr
j = hbc_nodes(i,i0)
h0(j) = MAX(zero,hbc(i0) - z0(j))
END DO
IF (hvalve) THEN
DO i = 1,n_hbc(i0)
dotprod = u0(hbc_nodes(i,i0))*hnormals(i,i0)%x + &
v0(hbc_nodes(i,i0))*hnormals(i,i0)%y
IF (dotprod < zero) THEN
u0(hbc_nodes(i,i0)) = zero
v0(hbc_nodes(i,i0)) = zero
END IF
END DO
END IF
CASE (0)
DO i = 1,n_hbc(i0)
dotprod = u0(hbc_nodes(i,i0))*hnormals(i,i0)%x + &
v0(hbc_nodes(i,i0))*hnormals(i,i0)%y
IF (dotprod < zero) THEN
u0(hbc_nodes(i,i0)) = zero
v0(hbc_nodes(i,i0)) = zero
END IF
END DO
CASE (2)
CASE DEFAULT
PRINT *,''
PRINT *,'ERROR: invalid value in htype.'
CALL byebye('Program SToRM stopped.')
END SELECT
END DO
IF (btype == 0) THEN
DO i = 1,n_wall
u0(wall_pts(i)) = zero
v0(wall_pts(i)) = zero
END DO
ELSE
DO i = 1,n_wall
dotprod = u0(wall_pts(i))*wtangs(i)%x + v0(wall_pts(i))*wtangs(i)%y
u0(wall_pts(i)) = dotprod*wtangs(i)%x
v0(wall_pts(i)) = dotprod*wtangs(i)%y
END DO
END IF
DO i = 1,n_pts
IF (h0(i) < h_dry) THEN
u0(i) = zero
v0(i) = zero
END IF
END DO
DO i0 = 1,n_inflowbdr
SELECT CASE (vtype(i0))
CASE (0)
DO i = 1,n_qin(i0)
IF (h0(qin_nodes(i,i0)) > h_dry) THEN
u0(qin_nodes(i,i0)) = qin(i0)*pnormals(i,i0)%x
v0(qin_nodes(i,i0)) = qin(i0)*pnormals(i,i0)%y
END IF
END DO
CASE (1,5)
q = vbc_by_h(i0,h0,ndim)
DO i = 1,n_qin(i0)
IF (h0(qin_nodes(i,i0)) > h_dry) THEN
u0(qin_nodes(i,i0)) = v_inflow(i)*pnormals(i,i0)%x
v0(qin_nodes(i,i0)) = v_inflow(i)*pnormals(i,i0)%y
END IF
END DO
CASE (2)
PRINT *,"VTYPE = 2: option not implemented yet."
CALL byebye('Program SToRM stopped.')
CASE (3)
DO i = 1,n_qin(i0)
IF (h0(qin_nodes(i,i0)) > h_dry) THEN
u0(qin_nodes(i,i0)) = qin(i0)
v0(qin_nodes(i,i0)) = qin2(i0)
END IF
END DO
CASE (4)
CASE DEFAULT
PRINT *,''
PRINT *,'ERROR: invalid value in vtype.'
CALL byebye('Program SToRM stopped.')
END SELECT
END DO
END SUBROUTINE enforce_bcs
FUNCTION equals(a,b)
USE parameters
USE constants
IMPLICIT NONE
LOGICAL :: equals
REAL(KIND=mp), INTENT(IN) :: a,b
REAL(KIND=mp) :: eps
eps = ABS(a)*fmachp
IF (eps == zero) THEN
eps = vsmall
END IF
IF (ABS(a-b) > eps) THEN
equals = .FALSE.
ELSE
equals = .TRUE.
ENDIF
END FUNCTION equals
SUBROUTINE error_est(etype)
USE parameters
USE constants
USE geometry
USE dep_vars
USE options
IMPLICIT NONE
INTEGER :: etype
INTEGER :: i,n1,n2,n3
REAL (KIND=mp) :: gradetax,gradetay,gradzx,gradzy,l
e_est = zero
SELECT CASE (etype)
CASE (1)
DO i = 1,n_elems
IF (h(i) < h_dry) CYCLE
n1 = grid(i)%vertex(1)
n2 = grid(i)%vertex(2)
n3 = grid(i)%vertex(3)
gradzx = (zvtx(n1)*(nodes(n2)%y - nodes(n3)%y) + &
zvtx(n2)*(nodes(n3)%y - nodes(n1)%y) + &
zvtx(n3)*(nodes(n1)%y - nodes(n2)%y))*half/grid(i)%area
gradzy = (zvtx(n1)*(nodes(n3)%x - nodes(n2)%x) + &
zvtx(n2)*(nodes(n1)%x - nodes(n3)%x) + &
zvtx(n3)*(nodes(n2)%x - nodes(n1)%x))*half/grid(i)%area
gradetax = gradzx + delh(i)%x
gradetay = gradzy + delh(i)%y
e_est(i) = SQRT(gradetax*gradetax + gradetay*gradetay)
END DO
CASE (2)
DO i = 1,n_elems
IF (h(i) < h_dry) CYCLE
n1 = grid(i)%vertex(1)
n2 = grid(i)%vertex(2)
n3 = grid(i)%vertex(3)
gradzx = (zvtx(n1)*(nodes(n2)%y - nodes(n3)%y) + &
zvtx(n2)*(nodes(n3)%y - nodes(n1)%y) + &
zvtx(n3)*(nodes(n1)%y - nodes(n2)%y))*half/grid(i)%area
gradzy = (zvtx(n1)*(nodes(n3)%x - nodes(n2)%x) + &
zvtx(n2)*(nodes(n1)%x - nodes(n3)%x) + &
zvtx(n3)*(nodes(n2)%x - nodes(n1)%x))*half/grid(i)%area
l = grid(i)%perim*one_third
gradetax = gradzx + delh(i)%x
gradetay = gradzy + delh(i)%y
e_est(i) = SQRT(gradetax*gradetax + gradetay*gradetay)*l/(z(i) + h(i))
END DO
CASE DEFAULT
PRINT *,''
PRINT *,'ERROR: invalid argument in subroutine error_est.'
CALL byebye('Program SToRM stopped.')
END SELECT
END SUBROUTINE error_est
SUBROUTINE flood_chan(t)
USE parameters
USE constants
USE geometry
USE dep_vars
USE memory
IMPLICIT NONE
REAL (KIND=mp) :: t
INTEGER :: i,j
REAL (KIND=mp) :: clength,hsink,hwater,s,scut,s_i(2),t0,vwave
t0 = zero
vwave = 1.0_mp
scut = 0.60_mp
clength = 28841.599_mp
s = vwave*(t + t0)/clength
s = MIN(s,scut)
IF (t > 86400) s = zero
s_i(1) = 0.90_mp; s_i(2) = one
IF (t > 93600_mp) s_i(1) = zero
hwater = 0.10_mp
hsink = 0.20_mp
s = zero
s_i(1) = zero
hsink = one
DO j = 1,n_chanbdr
DO i = 1,chantrigs(0,j)
IF (chant(i,j) < s) h(chantrigs(i,j)) = h(chantrigs(i,j)) + hwater
IF (chant(i,j) > s_i(1) .AND. chant(i,j) < s_i(2)) &
h(chantrigs(i,j)) = h(chantrigs(i,j))*(one - hsink)
END DO
END DO
END SUBROUTINE flood_chan
SUBROUTINE flux_ac(celledge,hlefty,uleft,vleft,hrighty,uright,vright,flux)
USE parameters
USE constants
IMPLICIT NONE
TYPE(edge), INTENT(IN) :: celledge
REAL (KIND=mp), INTENT(IN) :: hlefty,hrighty,uleft,uright,vleft,vright
REAL (KIND=mp), INTENT(OUT) :: flux(3)
INTEGER :: i
REAL (KIND=mp) :: cavg,denom,havg,hleft,hleftsqr,hright,hrightsqr,nx,ny, &
uavg,vavg
REAL (KIND=mp) :: a(3,3),aq(3),eigv(3),fileft(3),firight(3),gileft(3), &
giright(3),lev(3,3),qdiff(3),rev(3,3),temp1(3,3),temp2(3,3)
REAL (KIND=mp) :: delta = 0.10_mp
REAL (KIND=mp), PARAMETER :: small = 1.0e-9
nx = celledge%normal(1)
ny = celledge%normal(2)
hleft = hlefty
hright = hrighty
fileft(1) = hleft*uleft
fileft(2) = hleft*uleft*uleft + half*g*hleft*hleft
fileft(3) = hleft*uleft*vleft
gileft(1) = hleft*vleft
gileft(2) = hleft*uleft*vleft
gileft(3) = hleft*vleft*vleft + half*g*hleft*hleft
firight(1) = hright*uright
firight(2) = hright*uright*uright + half*g*hright*hright
firight(3) = hright*uright*vright
giright(1) = hright*vright
giright(2) = hright*uright*vright
giright(3) = hright*vright*vright + half*g*hright*hright
hleft = MAX(hleft,small)
hleftsqr=SQRT(hleft)
hright = MAX(hright,small)
hrightsqr = SQRT(hright)
denom = hrightsqr + hleftsqr
uavg = (uright*hrightsqr + uleft*hleftsqr)/denom
vavg = (vright*hrightsqr + vleft*hleftsqr)/denom
cavg = SQRT(half*g*(hleft + hright))
qdiff(1) = hright - hleft
qdiff(2) = hright*uright - hleft*uleft
qdiff(3) = hright*vright - hleft*vleft
rev(1,1) = zero; rev(1,2) = one; rev(1,3) = one
rev(2,1) = ny; rev(2,2) = uavg - cavg*nx; rev(2,3) = uavg + cavg*nx
rev(3,1) = -nx; rev(3,2) = vavg - cavg*ny; rev(3,3) = vavg + cavg*ny
lev(1,1) = -(uavg*ny - vavg*nx)
lev(1,2) = ny
lev(1,3) = -nx
lev(2,1) = half*(uavg*nx + vavg*ny)/cavg + half
lev(2,2) = -half*nx/cavg
lev(2,3) = -half*ny/cavg
lev(3,1) = -half*(uavg*nx + vavg*ny)/cavg + half
lev(3,2) = half*nx/cavg
lev(3,3) = half*ny/cavg
eigv(1) = MAX(ABS(uavg*nx + vavg*ny),delta)
eigv(2) = MAX(ABS(uavg*nx + vavg*ny - cavg),delta)
eigv(3) = MAX(ABS(uavg*nx + vavg*ny + cavg),delta)
a = zero
a(1,1) = eigv(1)
a(2,2) = eigv(2)
a(3,3) = eigv(3)
temp1 = MATMUL(rev,a)
temp2 = MATMUL(temp1,lev)
aq = MATMUL(temp2,qdiff)
DO i = 1,3
flux(i) = half*(firight(i)*nx + giright(i)*ny + fileft(i)*nx + &
gileft(i)*ny - aq(i))
END DO
END SUBROUTINE flux_ac
SUBROUTINE flux_bgnvc(celledge,hleft,uleft,vleft,hright,uright,vright,flux)
USE parameters
USE constants
IMPLICIT NONE
TYPE(edge), INTENT(IN) :: celledge
REAL (KIND=mp), INTENT(IN) :: hleft,hright,uleft,uright,vleft,vright
REAL (KIND=mp), INTENT(OUT) :: flux(3)
INTEGER :: i
REAL (KIND=mp) :: cavg,denom,hleftsqr,hrightsqr,nx,ny,uavg,vavg
REAL (KIND=mp) :: a(3,3),aq(3),fileft(3),firight(3),gileft(3),giright(3), &
lev(3,3),qdiff(3),rev(3,3),temp1(3,3),temp2(3,3)
nx = celledge%normal(1)
ny = celledge%normal(2)
fileft(1) = hleft*uleft*nx
fileft(2) = (hleft*uleft*uleft + half*g*hleft*hleft)*nx
fileft(3) = hleft*uleft*vleft*nx
gileft(1) = hleft*vleft*ny
gileft(2) = hleft*uleft*vleft*ny
gileft(3) = (hleft*vleft*vleft + half*g*hleft*hleft)*ny
firight(1) = hright*uright*nx
firight(2) = (hright*uright*uright + half*g*hright*hright)*nx
firight(3) = hright*uright*vright*nx
giright(1) = hright*vright*ny
giright(2) = hright*uright*vright*ny
giright(3) = (hright*vright*vright + half*g*hright*hright)*ny
hleftsqr=SQRT(hleft)
hrightsqr = SQRT(hright)
denom = hrightsqr + hleftsqr
uavg = (uright*hrightsqr + uleft*hleftsqr)/denom
vavg = (vright*hrightsqr + vleft*hleftsqr)/denom
cavg = SQRT(half*g*(hleft + hright))
qdiff(1) = hright - hleft
qdiff(2) = hright*uright - hleft*uleft
qdiff(3) = hright*vright - hleft*vleft
rev(1,1) = one; rev(1,2) = zero; rev(1,3) = one
rev(2,1) = uavg + cavg*nx; rev(2,2) = -cavg*ny; rev(2,3) = uavg - cavg*nx
rev(3,1) = vavg + cavg*ny; rev(3,2) = cavg*nx; rev(3,3) = vavg - cavg*ny
lev(1,1) = uavg*nx + vavg*ny + cavg
lev(1,2) = nx
lev(1,3) = ny
lev(2,1) = 2.0_mp*(uavg*ny - vavg*nx)
lev(2,2) = -2.0_mp*ny
lev(2,3) = 2.0_mp*nx
lev(3,1) = uavg*nx + vavg*ny + cavg
lev(3,2) = -nx
lev(3,3) = vavg - ny
lev = half*lev/cavg
a(1,1) = ABS(uavg*nx + vavg*ny + cavg); a(1,2) = zero; a(1,3) = zero
a(2,1) = zero; a(2,2) = ABS(uavg*nx + vavg*ny); a(2,3) = zero
a(3,1) = zero; a(3,2) = zero; a(3,3) = ABS(uavg*nx + vavg*ny - cavg)
temp1 = MATMUL(rev,a)
temp2 = MATMUL(temp1,lev)
aq = MATMUL(temp2,qdiff)
DO i = 1,3
flux(i) = half*(firight(i) + giright(i) + fileft(i) + gileft(i) - aq(i))
END DO
END SUBROUTINE flux_bgnvc
SUBROUTINE flux_HLL(celledge,hleft,uleft,vleft,hright,uright,vright,h0,flux)
USE parameters
USE constants
IMPLICIT NONE
TYPE(edge), INTENT(IN) :: celledge
REAL (KIND=mp), INTENT(IN) :: h0,hleft,hright,uleft,uright,vleft,vright
REAL (KIND=mp), INTENT(OUT) :: flux(3)
REAL (KIND=mp) :: aux1,aux2,cleft,cright,cstar,nx,ny,sleft,sright,ustar
REAL (KIND=mp) :: fleft(3),fright(3),gleft(3),gright(3)
nx = celledge%normal(1)
ny = celledge%normal(2)
cleft = SQRT(g*hleft)
cright = SQRT(g*hright)
ustar = half*((uleft + uright)*nx + (vleft + vright)*ny) + cleft - cright
cstar = half*(cleft + cright) + ((uleft - uright)*nx + &
(vleft - vright)*ny)/4.0_mp
IF (hleft < h0) THEN
sleft = uright*nx + vright*ny - 2.0_mp*cright
ELSE IF (hright < h0) THEN
sleft = uleft*nx + vleft*ny - cleft
ELSE
aux1 = uleft*nx + vleft*ny - cleft
aux2 = ustar - cstar
sleft = MIN(aux1,aux2)
END IF
IF (hleft < h0) THEN
sright = uright*nx + vright*ny + cright
ELSE IF (hright < h0) THEN
sright = uleft*nx + vleft*ny + 2.0_mp*cleft
ELSE
aux1 = uright*nx + vright*ny - cright
aux2 = ustar + cstar
sright = MAX(aux1,aux2)
END IF
fleft(1) = hleft*uleft*nx
fleft(2) = (hleft*uleft*uleft + half*g*hleft*hleft)*nx
fleft(3) = hleft*uleft*vleft*nx
gleft(1) = hleft*vleft*ny
gleft(2) = hleft*uleft*vleft*ny
gleft(3) = (hleft*vleft*vleft + half*g*hleft*hleft)*ny
fright(1) = hright*uright*nx
fright(2) = (hright*uright*uright + half*g*hright*hright)*nx
fright(3) = hright*uright*vright*nx
gright(1) = hright*vright*ny
gright(2) = hright*uright*vright*ny
gright(3) = (hright*vright*vright + half*g*hright*hright)*ny
IF (sleft > zero .AND. sright < zero) THEN
PRINT *,"Problem with HLL flux: debug."
CALL byebye("Program SToRM stopped.")
END IF
aux1 = one
IF (sleft > zero) THEN
flux(1) = fleft(1) + gleft(1)
flux(2) = fleft(2) + gleft(2)
flux(3) = fleft(3) + gleft(3)
ELSE IF (sright < zero) THEN
flux(1) = fright(1) + gright(1)
flux(2) = fright(2) + gright(2)
flux(3) = fright(3) + gright(3)
ELSE
aux1 = sright - sleft
flux(1) = sright*(fleft(1) + gleft(1)) - sleft*(fright(1) + gright(1)) + &
sleft*sright*(hright - hleft)/aux1
flux(2) = sright*(fleft(2) + gleft(2)) - sleft*(fright(2) + gright(2)) + &
sleft*sright*(hright*uright - hleft*uleft)/aux1
flux(2) = sright*(fleft(3) + gleft(3)) - sleft*(fright(3) + gright(3)) + &
sleft*sright*(hright*vright - hleft*vleft)/aux1
END IF
END SUBROUTINE flux_HLL
SUBROUTINE flux_rusanov(celledge,hleft,uleft,vleft,sqrtghleft,hright,uright, &
vright,sqrtghright,flux)
USE parameters
USE constants
IMPLICIT NONE
TYPE(edge), INTENT(IN) :: celledge
REAL (KIND=mp), INTENT(IN) :: hleft,hright,sqrtghleft,sqrtghright,uleft, &
uright,vleft,vright
REAL (KIND=mp), INTENT(OUT) :: flux(3)
INTEGER :: i
REAL (KIND=mp) :: nx,ny,splus,ul,ur
REAL (KIND=mp) :: fileft(3),firight(3),gileft(3),giright(3),qdiff(3)
nx = celledge%normal(1)
ny = celledge%normal(2)
fileft(1) = hleft*uleft
fileft(2) = hleft*uleft*uleft + half*g*hleft*hleft
fileft(3) = hleft*uleft*vleft
gileft(1) = hleft*vleft
gileft(2) = hleft*uleft*vleft
gileft(3) = hleft*vleft*vleft + half*g*hleft*hleft
firight(1) = hright*uright
firight(2) = hright*uright*uright + half*g*hright*hright
firight(3) = hright*uright*vright
giright(1) = hright*vright
giright(2) = hright*uright*vright
giright(3) = hright*vright*vright + half*g*hright*hright
ul = uleft*nx + vleft*ny
ur = uright*nx + vright*ny
splus = MAX(ABS(ul) + sqrtghleft,ABS(ur) + sqrtghright)
qdiff(1) = hright - hleft
qdiff(2) = hright*uright - hleft*uleft
qdiff(3) = hright*vright - hleft*vleft
DO i = 1,3
flux(i) = half*(firight(i)*nx + giright(i)*ny + fileft(i)*nx + &
gileft(i)*ny - splus*qdiff(i))
END DO
END SUBROUTINE flux_rusanov
SUBROUTINE fric_terms_impl(delta_t)
USE parameters
USE geometry
USE dep_vars
USE constants
USE options
USE RKparams
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: delta_t
INTEGER :: i
REAL (KIND=mp) :: ln,re
rough = zero
SELECT CASE (opt_friction)
CASE (0)
WHERE (.NOT. drycell) rough = g*cd*cd/h**one_third
CASE (1)
DO i = 1,n_elems
IF (cd(i) < fmachp) THEN
rough(i) = vlarge
ELSE
rough(i) = g/(cd(i)*cd(i))
END IF
END DO
CASE (2)
rough = cd
CASE (3)
DO i = 1,n_elems
IF (.NOT. drycell(i)) THEN
re = u_mag(i)*h(i)/visc
re = MAX(re,1000.0_mp)
re = MIN(re,REAL(2.5e7,mp))
rough(i) = cd(i)
IF (cd(i)/h(i) > 0.20_mp) rough(i) = 0.20_mp*h(i)
ln = LOG(1.725_mp/re + (rough(i)/(14.8_mp*h(i)))**1.11)
rough(i) = 0.204_mp/(ln*ln)
END IF
END DO
CASE (4)
DO i = 1,n_elems
IF (.NOT. drycell(i)) THEN
rough(i) = ((100.0_mp*cd(i))**0.16666666666666667_mp)/51.79_mp
rough(i) = g*rough(i)*rough(i)/h(i)**one_third
END IF
END DO
CASE (5)
DO i = 1,n_elems
IF (.NOT. drycell(i)) THEN
re = u_mag(i)*h(i)/visc
rough(i) = half*rcoef1(i)*re**(-rcoef2(i))
END IF
END DO
CASE DEFAULT
PRINT *,''
PRINT *,'ERROR: invalid value in opt_friction.'
CALL byebye('Program SToRM stopped.')
END SELECT
WHERE (.NOT. drycell) rough = h/(h + delta_t*rough*u_mag)
u = u*rough
v = v*rough
DO i = 1,n_elems
IF (.NOT. drycell(i)) THEN
rku(i,rkorder) = h(i)*u(i)
rkv(i,rkorder) = h(i)*v(i)
END IF
END DO
END SUBROUTINE fric_terms_impl
SUBROUTINE gaussgrad(phi,gradphi,ndim)
USE parameters
USE geometry
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: ndim
REAL (KIND=mp), INTENT(IN) :: phi(ndim)
TYPE(vector), INTENT(OUT) :: gradphi(ndim)
INTEGER :: i
DO i = 1,ndim
gwgrad(i)%x = gcoefx(i,1)*phi(gconn(i,1)) + gcoefx(i,2)*phi(gconn(i,2)) + &
gcoefx(i,3)*phi(gconn(i,3))
gwgrad(i)%y = gcoefy(i,1)*phi(gconn(i,1)) + gcoefy(i,2)*phi(gconn(i,2)) + &
gcoefy(i,3)*phi(gconn(i,3))
END DO
DO i = 1,ndim
gradphi(i)%x = half*(gwgrad(i)%x + gweight(i,1)*gwgrad(gconn(i,1))%x + &
gweight(i,2)*gwgrad(gconn(i,2))%x + &
gweight(i,3)*gwgrad(gconn(i,3))%x)
gradphi(i)%y = half*(gwgrad(i)%y + gweight(i,1)*gwgrad(gconn(i,1))%y + &
gweight(i,2)*gwgrad(gconn(i,2))%y + &
gweight(i,3)*gwgrad(gconn(i,3))%y)
END DO
END SUBROUTINE gaussgrad
SUBROUTINE gaussgrad_setup
USE parameters
USE constants
USE geometry
USE memory
IMPLICIT NONE
INTEGER :: i,ierror,j,k,l
REAL (KIND=mp) :: area,denom,x1,x2,xpto(3),y1,y2,ypto(3)
TYPE(vector) :: n(3)
ALLOCATE(gconn(n_elems,3),gcoefx(n_elems,3),gcoefy(n_elems,3), &
gweight(n_elems,3),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(4*n_elems*3 + 8*n_elems*9)
ALLOCATE(gwgrad(n_elems),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(v_size*n_elems)
DO i = 1,n_elems
DO j = 1,3
k = t2t3(i,j)
IF (k > 0) THEN
xpto(j) = grid(i)%xc
ypto(j) = grid(i)%yc
gconn(i,j) = k
ELSE
l = grid(i)%edge(j)
x1 = nodes(edges(l)%p(1))%x
y1 = nodes(edges(l)%p(1))%y
x2 = nodes(edges(l)%p(2))%x
y2 = nodes(edges(l)%p(2))%y
CALL ghost_cell(grid(i)%xc,grid(i)%yc,x1,y1,x2,y2,xpto(j),ypto(j))
gconn(i,j) = i
END IF
END DO
area = half*((xpto(2) - xpto(1))*(ypto(3) - ypto(1)) - &
(xpto(3) - xpto(1))*(ypto(2) - ypto(1)))
IF (area < zero) &
CALL byebye("Error is Gauss gradient set-up (gaussgrad_setup).")
n(1)%x = ypto(1) - ypto(2)  ;  n(1)%y = xpto(2) - xpto(1)
n(2)%x = ypto(2) - ypto(3)  ;  n(2)%y = xpto(3) - xpto(2)
n(3)%x = ypto(3) - ypto(1)  ;  n(3)%y = xpto(1) - xpto(3)
gcoefx(i,1) = half*(n(1)%x + n(3)%x)/area
gcoefx(i,2) = half*(n(1)%x + n(2)%x)/area
gcoefx(i,3) = half*(n(2)%x + n(3)%x)/area
gcoefy(i,1) = half*(n(1)%y + n(3)%y)/area
gcoefy(i,2) = half*(n(1)%y + n(2)%y)/area
gcoefy(i,3) = half*(n(2)%y + n(3)%y)/area
denom = (xpto(2) - xpto(1))*(ypto(3) - ypto(2)) - &
(xpto(3) - xpto(2))*(ypto(2) - ypto(1))
gweight(i,1) = ((xpto(2) - grid(i)%xc)*(ypto(3) - ypto(2)) - &
(xpto(3) - xpto(2))*(ypto(2) - grid(i)%yc))/denom
gweight(i,2) = ((xpto(3) - grid(i)%xc)*(ypto(1) - ypto(3)) - &
(xpto(1) - xpto(3))*(ypto(3) - grid(i)%yc))/denom
gweight(i,3) = ((xpto(1) - grid(i)%xc)*(ypto(2) - ypto(1)) - &
(xpto(2) - xpto(1))*(ypto(1) - grid(i)%yc))/denom
END DO
END SUBROUTINE gaussgrad_setup
SUBROUTINE ghost_cell(xc,yc,x1,y1,x2,y2,xg,yg)
USE parameters
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: xc,yc,x1,y1,x2,y2
REAL (KIND=mp), INTENT(OUT) :: xg,yg
REAL (KIND=mp) :: ax,ay,bx,by,cx,cy,denom
ax = xc - x1;  ay = yc - y1
bx = x2 - x1;  by = y2 - y1
denom = bx*bx + by*by
cx = (bx*(ax*bx + ay*by) + by*(ay*bx - ax*by))/denom
cy = (bx*(ax*by - ay*bx) + by*(ax*bx + ay*by))/denom
xg = cx + x1
yg = cy + y1
END SUBROUTINE ghost_cell
SUBROUTINE gradhordr(fegrad,gradphi)
USE parameters
USE geometry
USE constants
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: fegrad(n_elems,2)
REAL (KIND=mp), INTENT(OUT) :: gradphi(n_elems,2)
INTEGER :: i,ierror,j,k,n1,n2,n3
REAL (KIND=mp) :: tarea,weight
REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: phix,phiy
ALLOCATE(phix(n_pts),phiy(n_pts),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
DO i = 1,n_pts
tarea = zero
phix(i) = zero
phiy(i) = zero
DO j = 2,n2t(i,1) + 1
k = n2t(i,j)
weight = one/grid(k)%area
tarea = tarea + weight
phix(i) = phix(i) + weight*fegrad(k,1)
phiy(i) = phiy(i) + weight*fegrad(k,2)
END DO
phix(i) = phix(i)/tarea
phiy(i) = phiy(i)/tarea
END DO
DO i = 1,n_elems
n1 = grid(i)%vertex(1)
n2 = grid(i)%vertex(2)
n3 = grid(i)%vertex(3)
gradphi(i,1) = (phix(n1) + phix(n2) + phix(n3))/3.0_mp
gradphi(i,2) = (phiy(n1) + phiy(n2) + phiy(n3))/3.0_mp
END DO
DEALLOCATE(phix,phiy)
END SUBROUTINE gradhordr
SUBROUTINE hdlsq_grad_array(phi,gradphi,ndim)
USE parameters
USE geometry
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: ndim
REAL (KIND=mp), INTENT(IN) :: phi(ndim)
TYPE(vector), INTENT(OUT) :: gradphi(ndim)
INTEGER :: i,k,n
REAL (KIND=mp) :: di,dphi,ei
DO i = 1,n_elems
di = zero; ei = zero
n = t2tHD(i,1)
DO k = 1,n
dphi = phi(t2tHD(i,k+1)) - phi(i)
di = di + lsqweight(i,k)*dphi*dx(i,k)
ei = ei + lsqweight(i,k)*dphi*dy(i,k)
END DO
gradphi(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
gradphi(i)%y = (ai(i)*ei - di*bi(i))/det(i)
END DO
END SUBROUTINE hdlsq_grad_array
SUBROUTINE hdlsq_grad_array4h
USE parameters
USE geometry
USE constants
USE options
USE dep_vars
IMPLICIT NONE
INTEGER :: i,k,n
REAL (KIND=mp) :: di,dphi,ei,hnode
LOGICAL :: etac
DO i = 1,n_elems
IF (h(i) < h_dry) THEN
delh(i)%x = zero
delh(i)%y = zero
CYCLE
END IF
di = zero; ei = zero
etac = .FALSE.
n = t2tHD(i,1)
DO k = 1,n
hnode = h(t2tHD(i,k+1))
IF (hnode < h_dry .AND. z(t2tHD(i,k+1)) > h(i) + z(i)) THEN
etac = .TRUE.
EXIT
END IF
dphi = hnode - h(i)
di = di + lsqweight(i,k)*dphi*dx(i,k)
ei = ei + lsqweight(i,k)*dphi*dy(i,k)
END DO
IF (etac) THEN
delh(i)%x = -delz(i)%x
delh(i)%y = -delz(i)%y
ELSE
delh(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
delh(i)%y = (ai(i)*ei - di*bi(i))/det(i)
END IF
END DO
END SUBROUTINE hdlsq_grad_array4h
SUBROUTINE hdlsq_setup(wtype)
USE parameters
USE geometry
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: wtype
INTEGER :: i,ierror,j,k,l
k = 0
DO i = 1,n_elems
k = MAX(k,t2tHD(i,1))
END DO
ALLOCATE(ai(n_elems),bi(n_elems),ci(n_elems),det(n_elems),dx(n_elems,k), &
dy(n_elems,k),lsqweight(n_elems,k),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(8*(4 + 3*k)*n_elems)
DO i = 1,n_elems
j = t2tHD(i,1)
DO k = 2,j+1
l = t2tHD(i,k)
dx(i,k-1) = grid(l)%xc - grid(i)%xc
dy(i,k-1) = grid(l)%yc - grid(i)%yc
END DO
END DO
SELECT CASE (wtype)
CASE (1)
lsqweight = one
CASE (2)
DO i = 1,n_elems
j = t2tHD(i,1)
DO k = 1,j
lsqweight(i,k) = one/(dx(i,k)*dx(i,k) + dy(i,k)*dy(i,k))
END DO
END DO
CASE DEFAULT
PRINT *,"ERROR: least-squares procedure, invalid option."
CALL byebye('Program STORM stopped.')
END SELECT
ai = zero; bi = zero; ci = zero
DO i = 1,n_elems
j = t2tHD(i,1)
DO k = 1,j
ai(i) = ai(i) + lsqweight(i,k)*dx(i,k)*dx(i,k)
bi(i) = bi(i) + lsqweight(i,k)*dx(i,k)*dy(i,k)
ci(i) = ci(i) + lsqweight(i,k)*dy(i,k)*dy(i,k)
END DO
det(i) = ai(i)*ci(i) - bi(i)*bi(i)
END DO
END SUBROUTINE hdlsq_setup
FUNCTION heffective(cell_no)
USE parameters
USE constants
USE geometry
USE dep_vars
IMPLICIT NONE
REAL (KIND=mp) :: heffective
INTEGER, INTENT(IN) :: cell_no
REAL (KIND=mp) :: depth,eta
depth = h(cell_no)
eta = depth + z(cell_no)
IF (eta > zvtx(csortedz(cell_no)%p(3))) THEN
heffective = depth
RETURN
END IF
IF (eta > zvtx(csortedz(cell_no)%p(2))) THEN
heffective = one_third*(2.0_mp*eta - zvtx(csortedz(cell_no)%p(1)) - &
zvtx(csortedz(cell_no)%p(2)))
ELSE
heffective = one_third*(eta - zvtx(csortedz(cell_no)%p(1)))
END IF
END FUNCTION heffective
SUBROUTINE hequivalent(e,hleft,hright)
USE parameters
USE constants
USE geometry
USE dep_vars
USE options
IMPLICIT NONE
REAL (KIND=mp), INTENT(INOUT) :: hleft,hright
TYPE(edge), INTENT(IN) :: e
REAL (KIND=mp) :: h1,h2
IF (ABS(hleft - hright)/MAX(hleft,hright) < delta_hequiv) THEN
h1 = hvtx(e%p(1))
h2 = hvtx(e%p(2))
hright = SQRT((h1*h1 + h1*h2 + h2*h2)*one_third)
hleft = hright
END IF
END SUBROUTINE hequivalent
SUBROUTINE hgradfix(rkstep,solver)
USE geometry
USE constants
USE dep_vars
USE RKparams
USE options
IMPLICIT NONE
INTEGER, INTENT(IN) :: rkstep,solver
INTEGER :: i,k,l,total
REAL (KIND=mp) :: htotal,utotal,vtotal,ztotal
LOGICAL :: setgrad
INTEGER, EXTERNAL :: edge_in_element
bdrycell = .FALSE.
SELECT CASE (solver)
CASE (2)
DO i = 1,n_elems
IF (drycell(i)) CYCLE
setgrad = .FALSE.
DO k = 1,3
l = t2t3(k,i)
IF (l < 1) CYCLE
IF (drycell(l) .AND. z(l) > h(i) + z(i)) setgrad = .TRUE.
END DO
DO k = 1,3
l = t2t3(k,i)
IF (l < 1) CYCLE
IF (drycell(l) .AND. z(l) <= h(i) + z(i)) setgrad = .FALSE.
END DO
IF (setgrad) THEN
delh(i)%x = -delz(i)%x
delh(i)%y = -delz(i)%y
bdrycell(i) = .TRUE.
END IF
END DO
CASE (3)
DO i = 1,n_elems
IF (drycell(i)) CYCLE
setgrad = .FALSE.
DO k = 1,3
l = t2t3(k,i)
IF (l < 1) CYCLE
IF (drycell(l) .AND. z(l) > h(i) + z(i)) setgrad = .TRUE.
END DO
DO k = 1,3
l = t2t3(k,i)
IF (l < 1) CYCLE
IF (drycell(l) .AND. z(l) <= h(i) + z(i)) setgrad = .FALSE.
END DO
IF (setgrad) THEN
delh(i)%x = zero
delh(i)%y = zero
END IF
END DO
CASE DEFAULT
PRINT *,"ERROR: invalid option in hgradfix."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE hgradfix
SUBROUTINE hgradfix2(rkstep,solver)
USE geometry
USE constants
USE dep_vars
USE RKparams
USE options
IMPLICIT NONE
INTEGER, INTENT(IN) :: rkstep,solver
INTEGER :: i,k,l,total
REAL (KIND=mp) :: htotal,utotal,vtotal,ztotal
LOGICAL :: setgrad
INTEGER, EXTERNAL :: edge_in_element
bdrycell = .FALSE.
SELECT CASE (solver)
CASE (2)
DO i = 1,n_elems
IF (drycell(i)) CYCLE
setgrad = zeta(i) < zvtx(csortedz(i)%p(3))
DO k = 1,3
l = t2t3(k,i)
IF (l < 1) CYCLE
IF (drycell(l) .AND. z(l) < zeta(i)) setgrad = .FALSE.
END DO
IF (setgrad) THEN
delh(i)%x = -delz(i)%x
delh(i)%y = -delz(i)%y
bdrycell(i) = .TRUE.
u(i) = zero
v(i) = zero
rku(i,rkstep-1) = zero
rkv(i,rkstep-1) = zero
delu(i)%x = zero; delu(i)%y = zero
delv(i)%x = zero; delv(i)%y = zero
END IF
END DO
CASE (3)
PRINT *,"ERROR: invalid option in hgradfix2 (solver = 3)."
CALL byebye('Program STORM stopped.')
CASE DEFAULT
PRINT *,"ERROR: invalid option in hgradfix2."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE hgradfix2
SUBROUTINE ifname(v,vmax,fn,fnout)
IMPLICIT NONE
INTEGER, INTENT(IN) :: v,vmax
CHARACTER (LEN=*), INTENT(IN) :: fn
CHARACTER (LEN=*), INTENT(OUT) :: fnout
INTEGER :: i,j,k
CHARACTER (LEN=80) :: buffer
WRITE (buffer,*) vmax
buffer = ADJUSTL(buffer)
i = LEN_TRIM(buffer)
WRITE (buffer,*) v
buffer = ADJUSTL(buffer)
j = LEN_TRIM(buffer)
IF (j > i) THEN
fnout = ''
RETURN
END IF
fnout = TRIM(ADJUSTL(fn))
DO k = j,i - 1
fnout = TRIM(fnout) // "0"
END DO
fnout = TRIM(fnout) // TRIM(buffer) // ".dat"
END SUBROUTINE ifname
SUBROUTINE invisc_fluxes
USE parameters
USE geometry
USE dep_vars
USE vbc_arrays
USE options
USE constants
IMPLICIT NONE
INTEGER :: i,j,k,l,m,n
REAL (KIND=mp) :: atotal,dotprod,eps,f(3),h0,h1,h2,hleft,hright,incr, &
qfrac,sqrtghleft,sqrtghright,uedg,uleft,uright,vedg, &
vleft,vright,z0,zmiddle
LOGICAL :: set_wall
TYPE(edge) :: ework
INTEGER, EXTERNAL :: edge_in_element,local_edge
REAL (KIND=mp), EXTERNAL :: invert_edge_dirs
DO l = 1,flow_edges1
k = flowedg1(l)
i = edges(k)%e(1)
j = edges(k)%e(2)
CALL edge_copy(edges(k),ework)
eps = one
set_wall = .FALSE.
IF (drycell(i) .AND. drycell(j)) CYCLE
m = local_edge(grid(i),k)
IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)
IF (drycell(i)) THEN
uleft = zero
vleft = zero
hleft = zero
ELSE
hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
IF (hleft < zero) THEN
hleft = zero
partdry(i) = .TRUE.
END IF
uleft = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
vleft = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
END IF
n = local_edge(grid(j),k)
IF (drycell(j)) THEN
uright = zero
vright = zero
hright = zero
ELSE
hright = h(j) + delh(j)%x*rc(j,n)%x + delh(j)%y*rc(j,n)%y
IF (hright < zero) THEN
hright = zero
partdry(j) = .TRUE.
END IF
uright = u(j) + delu(j)%x*rc(j,n)%x + delu(j)%y*rc(j,n)%y
vright = v(j) + delv(j)%x*rc(j,n)%x + delv(j)%y*rc(j,n)%y
END IF
IF (hleft < h_dry .AND. hright < h_dry) CYCLE
IF (drycell(i) .OR. drycell(j)) THEN
h0 = MAX(hright,hleft)
z0 = half*(zvtx(ework%p(1)) + zvtx(ework%p(2)))
h0 = h0 + z0
IF (drycell(i)) THEN
IF (z(i) + h_wet > h0) THEN
set_wall = .TRUE.
hleft = hright
END IF
END IF
IF (drycell(j)) THEN
IF (z(j) + h_wet > h0) THEN
set_wall = .TRUE.
hright = hleft
END IF
END IF
END IF
IF (set_wall) THEN
f(1) = zero
hleft = MAX(hleft,hright)
f(2) = half*g*hleft*hleft*ework%normal(1)
f(3) = half*g*hleft*hleft*ework%normal(2)
ELSE
CALL flux_ac(ework,hleft,uleft,vleft,hright,uright,vright,f)
END IF
phi(i,1) = phi(i,1) - f(1)*ework%length
phi(i,2) = phi(i,2) - f(2)*ework%length
phi(i,3) = phi(i,3) - f(3)*ework%length
phi(j,1) = phi(j,1) + f(1)*ework%length
phi(j,2) = phi(j,2) + f(2)*ework%length
phi(j,3) = phi(j,3) + f(3)*ework%length
END DO
SELECT CASE (btype)
CASE (0)
CALL wenslp3_fluxes(h_dry)
CASE (1)
CALL weslp2_fluxes(h_dry)
END SELECT
CALL inflowQbyH
CALL inflowH
CALL outflowbyH
CALL crit_outflow
CALL free_flow
END SUBROUTINE invisc_fluxes
SUBROUTINE limiterBJ(beta,q,gradq,ndim)
USE parameters
USE constants
USE geometry
IMPLICIT NONE
INTEGER, INTENT(IN) :: ndim
REAL (KIND=mp), INTENT(IN) :: beta,q(ndim)
TYPE(vector), INTENT(INOUT) :: gradq(ndim)
INTEGER :: i,j,k
REAL (KIND=mp) :: a,b,delq,phi,q0,qj,qmax,qmin,qneighb,r(3)
LOGICAL, EXTERNAL :: equals
DO i = 1,n_elems
q0 = q(i)
DO j = 1,3
k = t2t3(j,i)
IF (k < 1) THEN
qneighb = q0
ELSE
qneighb = q(k)
END IF
qmin = MIN(q0,qneighb)
qmax = MAX(q0,qneighb)
qj = q0 + gradq(i)%x*rc(i,j)%x + gradq(i)%y*rc(i,j)%y
delq = qj - q0
IF (equals(delq,zero)) THEN
r(j) = one
ELSE IF (delq < 0) THEN
r(j) = (qmin - q0)/delq
ELSE
r(j) = (qmax - q0)/delq
END IF
a = MIN(beta*r(j),one)
b = MIN(r(j),beta)
r(j) = MAX(a,b)
END DO
phi = MIN(r(1),r(2),r(3))
gradq(i)%x = phi*gradq(i)%x
gradq(i)%y = phi*gradq(i)%y
END DO
END SUBROUTINE limiterBJ
SUBROUTINE lsq_grad_array(solver,phi,gradphi,ndim)
USE parameters
USE geometry
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: ndim,solver
REAL (KIND=mp), INTENT(IN) :: phi(ndim)
TYPE(vector), INTENT(OUT) :: gradphi(ndim)
INTEGER :: i,k,n
REAL (KIND=mp) :: di,dphi,ei
SELECT CASE (solver)
CASE (1)
DO i = 1,n_pts
di = zero; ei = zero
n = n2n(i,1)
DO k = 1,n
dphi = phi(n2n(i,k+1)) - phi(i)
di = di + lsqweight(i,k)*dphi*dx(i,k)
ei = ei + lsqweight(i,k)*dphi*dy(i,k)
END DO
gradphi(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
gradphi(i)%y = (ai(i)*ei - di*bi(i))/det(i)
END DO
CASE (2)
DO i = 1,n_elems
di = zero; ei = zero
n = t2t(i,1)
DO k = 1,n
dphi = phi(t2t(i,k+1)) - phi(i)
di = di + lsqweight(i,k)*dphi*dx(i,k)
ei = ei + lsqweight(i,k)*dphi*dy(i,k)
END DO
gradphi(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
gradphi(i)%y = (ai(i)*ei - di*bi(i))/det(i)
END DO
CASE DEFAULT
PRINT *,"ERROR: least-squares gradient computation, invalid option."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE lsq_grad_array
SUBROUTINE lsq_grad_array4h(solver)
USE parameters
USE geometry
USE constants
USE options
USE dep_vars
IMPLICIT NONE
INTEGER, INTENT(IN) :: solver
INTEGER :: i,k,n
REAL (KIND=mp) :: di,dphi,ei,hnode
LOGICAL :: etac
SELECT CASE (solver)
CASE (1)
DO i = 1,n_pts
di = zero; ei = zero
n = n2n(i,1)
DO k = 1,n
hnode = h(n2n(i,k+1))
IF (hnode < h_dry) hnode = z(i) + h(i) - z(n2n(i,k+1))
dphi = hnode - h(i)
di = di + lsqweight(i,k)*dphi*dx(i,k)
ei = ei + lsqweight(i,k)*dphi*dy(i,k)
END DO
delh(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
delh(i)%y = (ai(i)*ei - di*bi(i))/det(i)
END DO
CASE (2)
DO i = 1,n_elems
IF (drycell(i)) THEN
delh(i)%x = zero
delh(i)%y = zero
CYCLE
END IF
di = zero; ei = zero
etac = .FALSE.
n = t2t(i,1)
DO k = 1,n
hnode = h(t2t(i,k+1))
IF (drycell(t2t(i,k+1)) .AND. z(t2t(i,k+1)) > h(i) + z(i)) THEN
etac = .TRUE.
EXIT
END IF
dphi = hnode - h(i)
di = di + lsqweight(i,k)*dphi*dx(i,k)
ei = ei + lsqweight(i,k)*dphi*dy(i,k)
END DO
IF (etac) THEN
delh(i)%x = -delz(i)%x
delh(i)%y = -delz(i)%y
ELSE
delh(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
delh(i)%y = (ai(i)*ei - di*bi(i))/det(i)
END IF
END DO
CASE (3)
DO i = 1,n_elems
IF (drycell(i)) THEN
delh(i)%x = delz(i)%x
delh(i)%y = delz(i)%y
CYCLE
END IF
di = zero; ei = zero
n = t2t(i,1)
DO k = 1,n
IF (drycell(t2t(i,k+1))) THEN
di = zero
ei = zero
EXIT
END IF
dphi = zeta(t2t(i,k+1)) - zeta(i)
di = di + lsqweight(i,k)*dphi*dx(i,k)
ei = ei + lsqweight(i,k)*dphi*dy(i,k)
END DO
delh(i)%x = (ci(i)*di - bi(i)*ei)/det(i)
delh(i)%y = (ai(i)*ei - di*bi(i))/det(i)
END DO
CASE DEFAULT
PRINT *,"ERROR: least-squares gradient computation, invalid option."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE lsq_grad_array4h
SUBROUTINE lsq_gradient(solver,node,phi,ndim,phix,phiy)
USE parameters
USE geometry
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: ndim,node,solver
REAL (KIND=mp), INTENT(IN) :: phi(ndim)
REAL (KIND=mp), INTENT(OUT) :: phix,phiy
INTEGER :: i,k,n
REAL (KIND=mp) :: di,dphi,ei
SELECT CASE (solver)
CASE (1)
di = zero; ei = zero
n = n2n(node,1)
i = node
DO k = 1,n
dphi = phi(n2n(i,k+1)) - phi(i)
di = di + lsqweight(i,k)*dphi*dx(i,k)
ei = ei + lsqweight(i,k)*dphi*dy(i,k)
END DO
phix = (ci(i)*di - bi(i)*ei)/det(i)
phiy = (ai(i)*ei - di*bi(i))/det(i)
CASE (2)
di = zero; ei = zero
n = t2t(node,1)
i = node
DO k = 1,n
dphi = phi(t2t(i,k+1)) - phi(i)
di = di + lsqweight(i,k)*dphi*dx(i,k)
ei = ei + lsqweight(i,k)*dphi*dy(i,k)
END DO
phix = (ci(i)*di - bi(i)*ei)/det(i)
phiy = (ai(i)*ei - di*bi(i))/det(i)
CASE DEFAULT
PRINT *,"ERROR: least-squares gradient computation, invalid option."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE lsq_gradient
SUBROUTINE lsq_setup(solver,type)
USE parameters
USE geometry
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: type,solver
INTEGER :: i,ierror,j,k,l
SELECT CASE (solver)
CASE (1)
k = 0
DO i = 1,n_pts
k = MAX(k,n2n(i,1))
END DO
ALLOCATE(ai(n_pts),bi(n_pts),ci(n_pts),det(n_pts),dx(n_pts,k), &
dy(n_pts,k),lsqweight(n_pts,k),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(8*(4 + 3*k)*n_pts)
DO i = 1,n_pts
j = n2n(i,1)
DO k = 2,j+1
l = n2n(i,k)
dx(i,k-1) = nodes(l)%x - nodes(i)%x
dy(i,k-1) = nodes(l)%y - nodes(i)%y
END DO
END DO
SELECT CASE (type)
CASE (1)
lsqweight = one
CASE (2,3)
DO i = 1,n_pts
j = n2n(i,1)
DO k = 1,j
lsqweight(i,k) = one/(dx(i,k)*dx(i,k) + dy(i,k)*dy(i,k))
END DO
END DO
CASE DEFAULT
PRINT *,"ERROR: least-squares procedure, invalid option."
CALL byebye('Program STORM stopped.')
END SELECT
ai = zero; bi = zero; ci = zero
DO i = 1,n_pts
j = n2n(i,1)
DO k = 1,j
ai(i) = ai(i) + lsqweight(i,k)*dx(i,k)*dx(i,k)
bi(i) = bi(i) + lsqweight(i,k)*dx(i,k)*dy(i,k)
ci(i) = ci(i) + lsqweight(i,k)*dy(i,k)*dy(i,k)
END DO
det(i) = ai(i)*ci(i) - bi(i)*bi(i)
END DO
CASE (2)
k = 0
DO i = 1,n_elems
k = MAX(k,t2t(i,1))
END DO
ALLOCATE(ai(n_elems),bi(n_elems),ci(n_elems),det(n_elems),dx(n_elems,k), &
dy(n_elems,k),lsqweight(n_elems,k),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(8*(4 + 3*k)*n_elems)
DO i = 1,n_elems
j = t2t(i,1)
DO k = 2,j+1
l = t2t(i,k)
dx(i,k-1) = grid(l)%xc - grid(i)%xc
dy(i,k-1) = grid(l)%yc - grid(i)%yc
END DO
END DO
SELECT CASE (type)
CASE (1)
lsqweight = one
CASE (2)
DO i = 1,n_elems
j = t2t(i,1)
DO k = 1,j
lsqweight(i,k) = one/(dx(i,k)*dx(i,k) + dy(i,k)*dy(i,k))
END DO
END DO
CASE DEFAULT
PRINT *,"ERROR: least-squares procedure, invalid option."
CALL byebye('Program STORM stopped.')
END SELECT
ai = zero; bi = zero; ci = zero
DO i = 1,n_elems
j = t2t(i,1)
DO k = 1,j
ai(i) = ai(i) + lsqweight(i,k)*dx(i,k)*dx(i,k)
bi(i) = bi(i) + lsqweight(i,k)*dx(i,k)*dy(i,k)
ci(i) = ci(i) + lsqweight(i,k)*dy(i,k)*dy(i,k)
END DO
det(i) = ai(i)*ci(i) - bi(i)*bi(i)
END DO
CASE DEFAULT
PRINT *,"ERROR: least-squares procedure, invalid option."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE lsq_setup
SUBROUTINE margins
USE constants
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER :: count1,count2,i,j,n
DO i = 1,n_elems
IF (drycell(i)) CYCLE
IF (zeta(i) < zvtx(csortedz(i)%p(3))) THEN
DO j = 1,3
n = t2t3(j,i)
count1 = 0
count2 = 0
IF (n < 1) THEN
CYCLE
ELSE IF (drycell(n)) THEN
count1 = count1 + 1
ELSE IF (zeta(n) > zvtx(csortedz(n)%p(3))) THEN
count2 = count2 + 1
END IF
END DO
IF ((count1 > 0) .AND. (count2 > 0)) THEN
u(i) = zero
v(i) = zero
END IF
END IF
END DO
END SUBROUTINE margins
SUBROUTINE overdraft
USE parameters
USE constants
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER :: i,j,k
REAL (KIND=mp) :: vol
DO k = 1,3
DO i = 1,n_elems
IF (h(i) < zero) THEN
vol = h(i)*grid(i)%area
j = grid(i)%o_cell
h(j) = h(j) + vol/grid(j)%area
h(i) = zero
END IF
END DO
END DO
END SUBROUTINE overdraft
SUBROUTINE positivity_vtx
USE constants
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER :: i,j
DO i = 1,n_pts
IF (hvtx(i) < zero) THEN
hvtx(i) = zero
zetavtx(i) = zvtx(i)
END IF
END DO
END SUBROUTINE positivity_vtx
LOGICAL FUNCTION pto_in_triangle(e,x,y)
USE constants
USE geometry
IMPLICIT NONE
TYPE(triangle) :: e
REAL (KIND=mp), INTENT(IN) :: x,y
REAL(KIND=mp) :: s,t,tarea,x1,x2,x3,y1,y2,y3
pto_in_triangle = .FALSE.
tarea = e%area
x1 = nodes(e%vertex(1))%x
y1 = nodes(e%vertex(1))%y
x2 = nodes(e%vertex(2))%x
y2 = nodes(e%vertex(2))%y
x3 = nodes(e%vertex(3))%x
y3 = nodes(e%vertex(3))%y
s = half/tarea*(y1*x3 - x1*y3 + (y3 - y1)*x + (x1 - x3)*y)
t = half/tarea*(x1*y2 - y1*x2 + (y1 - y2)*x + (x2 - x1)*y)
IF (s < zero) RETURN
IF (t < zero) RETURN
IF (one-s-t < zero) RETURN
pto_in_triangle = .TRUE.
END FUNCTION pto_in_triangle
INTEGER FUNCTION query(q)
IMPLICIT NONE
LOGICAL :: q
query = 1
IF (q) query = 2
END FUNCTION query
LOGICAL FUNCTION readline(funit,line,lineno)
IMPLICIT NONE
INTEGER, INTENT(IN) :: funit
INTEGER, INTENT(INOUT) :: lineno
CHARACTER(LEN=40), INTENT(OUT) :: line
INTEGER :: i
CHARACTER, EXTERNAL :: to_upper
DO WHILE (.TRUE.)
line = ''
lineno = lineno + 1
READ (funit,'(A)',ERR = 1,END = 1) line
line = ADJUSTL(line)
IF (line(1:1) == '#') CYCLE
IF (LEN_TRIM(line) == 0) CYCLE
DO i = 1,LEN_TRIM(line)
line(i:i) = to_upper(line(i:i))
END DO
readline = .TRUE.
RETURN
END DO
1 readline = .FALSE.
END FUNCTION readline
INTEGER FUNCTION readp(string)
IMPLICIT NONE
CHARACTER (LEN=*), INTENT(INOUT) :: string
INTEGER :: j
IF (LEN_TRIM(string) == 0) THEN
readp = 0
ELSE
READ(string,*) readp
j = INDEX(string," ")
string = string(j:)
string = ADJUSTL(string)
END IF
END FUNCTION readp
SUBROUTINE rmv_array_pto(i,n,iarray)
IMPLICIT NONE
INTEGER, INTENT(IN) :: i
INTEGER, INTENT(INOUT) :: n
INTEGER, DIMENSION(n), INTENT(INOUT) :: iarray
INTEGER :: j
DO j = i,n-1
iarray(j) = iarray(j+1)
END DO
n = n - 1
END SUBROUTINE rmv_array_pto
FUNCTION skin_friction(ftype,ffactor,depth,vel,nu,vegc1,vegc2)
USE parameters
USE constants
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: depth,ffactor,nu,vegc1,vegc2,vel
INTEGER, INTENT(IN) :: ftype
REAL (KIND=mp) :: skin_friction
REAL (KIND=mp) :: ln,re,work
REAL (KIND=mp), PARAMETER :: small = 1.0E-9
SELECT CASE (ftype)
CASE (0)
IF (depth < small) THEN
skin_friction = one/small
ELSE
skin_friction = g*ffactor*ffactor/depth**one_third
END IF
CASE (1)
IF (ffactor < small) THEN
skin_friction = one/small
ELSE
skin_friction = g/(ffactor*ffactor)
END IF
CASE (2)
skin_friction = ffactor
CASE (3)
re = ABS(vel)*depth/nu
IF (depth < small) THEN
skin_friction = one/small
ELSE
re = MAX(re,1000.0_mp)
re = MIN(re,REAL(2.5e7,mp))
work = ffactor
IF (ffactor/depth > 0.20_mp) work = 0.20_mp*depth
ln = LOG(1.725_mp/re + (work/(14.8_mp*depth))**1.11)
skin_friction = 0.204_mp/(ln*ln)
END IF
CASE (4)
work = ((100.0_mp*ffactor)**0.16666666666666667_mp)/51.79_mp
IF (depth < small) THEN
skin_friction = one/small
ELSE
skin_friction = g*work*work/depth**one_third
END IF
CASE (5)
re = ABS(vel)*depth/nu
skin_friction = half*vegc1*re**(-vegc2)
CASE DEFAULT
PRINT *,''
PRINT *,'ERROR: invalid value in ftype, function skin_friction.'
CALL byebye('Program SToRM stopped.')
END SELECT
END FUNCTION skin_friction
CHARACTER FUNCTION to_upper(ch)
IMPLICIT NONE
CHARACTER, INTENT(IN) :: ch
INTEGER, PARAMETER :: strand = IACHAR("A") - IACHAR("a")
IF ("a" <= ch .AND. ch <= "z") THEN
to_upper = ACHAR(IACHAR(ch) + strand)
ELSE
to_upper = ch
END IF
END FUNCTION to_upper
SUBROUTINE vcrop(vclip,vceiling)
USE parameters
USE dep_vars
IMPLICIT NONE
LOGICAL, INTENT(IN) :: vclip
REAL(KIND=mp), INTENT(IN) :: vceiling
IF (.NOT. vclip) RETURN
WHERE (u_mag > vceiling)
u = vceiling*u/u_mag
v = vceiling*v/u_mag
u_mag = vceiling
END WHERE
END SUBROUTINE vcrop
SUBROUTINE visc_fluxes
USE geometry
USE dep_vars
USE options
USE constants
IMPLICIT NONE
INTEGER :: i,j,k,k1,k2,l,m,p1,p2,p3
REAL (KIND=mp) :: area,b(3),c(3),evisc,eps,fv(3),gv(3),h0,hedge,nx,ny,x3, &
xc,y3,yc,z0
TYPE(vector) :: gradu,gradv
TYPE(edge) :: ework
INTEGER, EXTERNAL :: edge_in_element,local_edge
REAL (KIND=mp), EXTERNAL :: invert_edge_dirs
DO i = 1,n_bpolygon
j = bpolygon(i)
k = edge_in_element(edges(j))
l = local_edge(grid(k),j)
CALL edge_copy(edges(j),ework)
IF (grid(k)%edge(l) < 0) THEN
p2 = edges(j)%p(2)
p3 = edges(j)%p(1)
eps = invert_edge_dirs(ework)
ELSE
p2 = edges(j)%p(1)
p3 = edges(j)%p(2)
END IF
hedge = half*(hvtx(p2) + hvtx(p3))
IF (hedge < h_dry) CYCLE
area = grid(k)%area*one_third*2.0_mp
b(1) = nodes(p2)%y - nodes(p3)%y ; c(1) = nodes(p3)%x - nodes(p2)%x
b(2) = nodes(p3)%y - grid(k)%yc  ; c(2) = grid(k)%xc - nodes(p3)%x
b(3) = grid(k)%yc - nodes(p2)%y  ; c(3) = nodes(p2)%x - grid(k)%xc
gradu%x = (b(1)*u(k) + b(2)*uvtx(p2) + b(3)*uvtx(p3))/area
gradu%y = (c(1)*u(k) + c(2)*uvtx(p2) + c(3)*uvtx(p3))/area
gradv%x = (b(1)*v(k) + b(2)*vvtx(p2) + b(3)*vvtx(p3))/area
gradv%y = (c(1)*v(k) + c(2)*vvtx(p2) + c(3)*vvtx(p3))/area
evisc = zero
IF (parabev) THEN
m = 0
IF (hvtx(p2) > h_dry) THEN
m = m + 1
evisc = evisc + ev(p2)
END IF
IF (hvtx(p3) > h_dry) THEN
m = m + 1
evisc = evisc + ev(p3)
END IF
evisc = evisc/REAL(m,mp)
END IF
evisc = evisc + visc
nx = ework%normal(1)*ework%length
ny = ework%normal(2)*ework%length
fv(2) = evisc*hedge*gradu%x*nx
fv(3) = evisc*hedge*gradv%x*nx
gv(2) = evisc*hedge*gradu%y*ny
gv(3) = evisc*hedge*gradv%y*ny
DO l = 2,3
phi(k,l) = phi(k,l) + fv(l) + gv(l)
END DO
END DO
DO i = 1,flow_edges1
j = flowedg1(i)
CALL edge_copy(edges(j),ework)
k1 = edges(j)%e(1)
l = local_edge(grid(k1),j)
k2 = edges(j)%e(2)
IF (drycell(k2)) THEN
hedge = h(k1) + delh(k1)%x*rc(k1,l)%x + delh(k1)%y*rc(k1,l)%y
IF (hedge < h_dry) CYCLE
IF (grid(k1)%edge(l) < 0) THEN
p2 = edges(j)%p(2)
p3 = edges(j)%p(1)
eps = invert_edge_dirs(ework)
ELSE
p2 = edges(j)%p(1)
p3 = edges(j)%p(2)
END IF
area = grid(k1)%area*one_third*2.0_mp
b(1) = nodes(p2)%y - nodes(p3)%y ; c(1) = nodes(p3)%x - nodes(p2)%x
b(2) = nodes(p3)%y - grid(k1)%yc ; c(2) = grid(k1)%xc - nodes(p3)%x
b(3) = grid(k1)%yc - nodes(p2)%y ; c(3) = nodes(p2)%x - grid(k1)%xc
gradu%x = (b(1)*u(k1) + b(2)*uvtx(p2) + b(3)*uvtx(p3))/area
gradu%y = (c(1)*u(k1) + c(2)*uvtx(p2) + c(3)*uvtx(p3))/area
gradv%x = (b(1)*v(k1) + b(2)*vvtx(p2) + b(3)*vvtx(p3))/area
gradv%y = (c(1)*v(k1) + c(2)*vvtx(p2) + c(3)*vvtx(p3))/area
ELSE IF (drycell(k1)) THEN
m = local_edge(grid(k2),j)
hedge = h(k2) + delh(k2)%x*rc(k2,m)%x + delh(k2)%y*rc(k2,m)%y
IF (hedge < h_dry) CYCLE
IF (grid(k1)%edge(m) < 0) THEN
p2 = edges(j)%p(2)
p3 = edges(j)%p(1)
ELSE
p2 = edges(j)%p(1)
p3 = edges(j)%p(2)
eps = invert_edge_dirs(ework)
END IF
area = grid(k2)%area*one_third*2.0_mp
b(1) = nodes(p2)%y - nodes(p3)%y ; c(1) = nodes(p3)%x - nodes(p2)%x
b(2) = nodes(p3)%y - grid(k2)%yc ; c(2) = grid(k2)%xc - nodes(p3)%x
b(3) = grid(k2)%yc - nodes(p2)%y ; c(3) = nodes(p2)%x - grid(k2)%xc
gradu%x = (b(1)*u(k2) + b(2)*uvtx(p2) + b(3)*uvtx(p3))/area
gradu%y = (c(1)*u(k2) + c(2)*uvtx(p2) + c(3)*uvtx(p3))/area
gradv%x = (b(1)*v(k2) + b(2)*vvtx(p2) + b(3)*vvtx(p3))/area
gradv%y = (c(1)*v(k2) + c(2)*vvtx(p2) + c(3)*vvtx(p3))/area
ELSE
IF (grid(k1)%edge(l) < 0) THEN
p1 = edges(j)%p(2)
p2 = edges(j)%p(1)
eps = invert_edge_dirs(ework)
ELSE
p1 = edges(j)%p(1)
p2 = edges(j)%p(2)
END IF
hedge = half*(hvtx(p1) + hvtx(p2))
area = (grid(k1)%area + grid(k2)%area)*one_third*2.0_mp
x3 = grid(k1)%xc  ;  y3 = grid(k1)%yc
xc = grid(k2)%xc  ;  yc = grid(k2)%yc
gradu%x = ((yc - y3)*(uvtx(p1) - uvtx(p2)) + (nodes(p1)%y - &
nodes(p2)%y)*(u(k1) - u(k2)))/area
gradu%y = ((x3 - xc)*(uvtx(p1) - uvtx(p2)) + (nodes(p2)%x - &
nodes(p1)%x)*(u(k1) - u(k2)))/area
gradv%x = ((yc - y3)*(vvtx(p1) - vvtx(p2)) + (nodes(p1)%y - &
nodes(p2)%y)*(v(k1) - v(k2)))/area
gradv%y = ((x3 - xc)*(vvtx(p1) - vvtx(p2)) + (nodes(p2)%x - &
nodes(p1)%x)*(v(k1) - v(k2)))/area
END IF
evisc = zero
IF (parabev) THEN
m = 0
IF (hvtx(p1) > h_dry) THEN
m = m + 1
evisc = evisc + ev(p1)
END IF
IF (hvtx(p2) > h_dry) THEN
m = m + 1
evisc = evisc + ev(p2)
END IF
evisc = evisc/REAL(m,mp)
END IF
evisc = evisc + visc
nx = ework%normal(1)*ework%length
ny = ework%normal(2)*ework%length
fv(2) = evisc*hedge*gradu%x*nx
fv(3) = evisc*hedge*gradv%x*nx
gv(2) = evisc*hedge*gradu%y*ny
gv(3) = evisc*hedge*gradv%y*ny
DO l = 2,3
phi(k1,l) = phi(k1,l) + fv(l) + gv(l)
phi(k2,l) = phi(k2,l) - fv(l) - gv(l)
END DO
END DO
END SUBROUTINE visc_fluxes
SUBROUTINE wenslp2_fluxes(hdry)
USE parameters
USE constants
USE geometry
USE dep_vars
USE vbc_arrays
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: hdry
INTEGER :: i,k,l,m
REAL (KIND=mp) :: eps,hedge
TYPE(edge) :: ework
INTEGER, EXTERNAL :: edge_in_element,local_edge
REAL (KIND=mp), EXTERNAL :: invert_edge_dirs
DO l = 1,wall_edges1
k = walledg1(l)
i = edge_in_element(edges(k))
CALL edge_copy(edges(k),ework)
eps = one
IF (drycell(i)) CYCLE
m = local_edge(grid(i),k)
hedge = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
IF (hedge < zero) THEN
hedge = zero
partdry(i) = .TRUE.
END IF
IF (hedge < hdry) CYCLE
IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)
phi(i,2) = phi(i,2) - half*g*hedge*hedge*ework%normal(1)*ework%length
phi(i,3) = phi(i,3) - half*g*hedge*hedge*ework%normal(2)*ework%length
END DO
END SUBROUTINE wenslp2_fluxes
SUBROUTINE wenslp3_fluxes(hdry)
USE parameters
USE constants
USE geometry
USE dep_vars
USE vbc_arrays
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: hdry
INTEGER :: i,k,l,m
REAL (KIND=mp) :: eps,h1,h2,hedge
TYPE(edge) :: ework
INTEGER, EXTERNAL :: edge_in_element,local_edge
REAL (KIND=mp), EXTERNAL :: invert_edge_dirs
DO l = 1,wall_edges1
k = walledg1(l)
i = edge_in_element(edges(k))
CALL edge_copy(edges(k),ework)
eps = one
IF (drycell(i)) CYCLE
h1 = hvtx(ework%p(1))
h2 = hvtx(ework%p(2))
hedge = SQRT((h1*h1 + h1*h2 + h2*h2)*one_third)
IF (hedge < zero) THEN
hedge = zero
partdry(i) = .TRUE.
END IF
IF (hedge < hdry) CYCLE
m = local_edge(grid(i),k)
IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)
phi(i,2) = phi(i,2) - half*g*hedge*hedge*ework%normal(1)*ework%length
phi(i,3) = phi(i,3) - half*g*hedge*hedge*ework%normal(2)*ework%length
END DO
END SUBROUTINE wenslp3_fluxes
SUBROUTINE wenslp_fluxes(hdry)
USE parameters
USE constants
USE geometry
USE dep_vars
USE vbc_arrays
USE options
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: hdry
INTEGER :: i,k,l,m
REAL (KIND=mp) :: eps,f(3),hleft,hright,uleft,uright,vleft,vright
TYPE(edge) :: ework
INTEGER, EXTERNAL :: edge_in_element,local_edge
REAL (KIND=mp), EXTERNAL :: invert_edge_dirs
uright = zero
vright = zero
DO l = 1,wall_edges1
k = walledg1(l)
i = edge_in_element(edges(k))
CALL edge_copy(edges(k),ework)
eps = one
IF (drycell(i)) CYCLE
m = local_edge(grid(i),k)
hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
IF (hleft < hdry) THEN
IF (hleft < zero) partdry(i) = .TRUE.
CYCLE
ELSE
uleft = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
vleft = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
END IF
IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)
hright = hleft
CALL flux_ac(ework,hleft,uleft,vleft,hright,uright,vright,f)
phi(i,2) = phi(i,2) - f(2)*ework%length
phi(i,3) = phi(i,3) - f(3)*ework%length
END DO
END SUBROUTINE wenslp_fluxes
SUBROUTINE weslp2_fluxes(hdry)
USE parameters
USE constants
USE geometry
USE dep_vars
USE vbc_arrays
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: hdry
INTEGER :: i,k,l,m
REAL (KIND=mp) :: dh,eps,f(3),hleft,hright,uleft,uright,vleft,vright
TYPE(edge) :: ework
INTEGER, EXTERNAL :: edge_in_element,local_edge
REAL (KIND=mp), EXTERNAL :: invert_edge_dirs
DO l = 1,wall_edges1
k = walledg1(l)
i = edge_in_element(edges(k))
CALL edge_copy(edges(k),ework)
eps = one
IF (drycell(i)) CYCLE
m = local_edge(grid(i),k)
hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
IF (hleft < hdry) THEN
IF (hleft < zero) partdry(i) = .TRUE.
CYCLE
ELSE
uleft = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
vleft = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
END IF
IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)
uright = uleft*ework%tang(1) + vleft*ework%tang(2)
uright = uright + 2.0_mp*SQRT(g*hleft)
hright = uright*uright/4.0_mp/g
f(1) = zero
f(2) = half*g*hright*hright*ework%normal(1)
f(3) = half*g*hright*hright*ework%normal(2)
dh = hvtx(ework%p(1)) - hvtx(ework%p(2))
dh = g*dh*dh/12.0_mp
dh = zero
phi(i,2) = phi(i,2) - (f(2) + dh*ework%normal(1))*ework%length
phi(i,3) = phi(i,3) - (f(3) + dh*ework%normal(2))*ework%length
END DO
END SUBROUTINE weslp2_fluxes
SUBROUTINE weslp_fluxes(hdry)
USE parameters
USE constants
USE geometry
USE dep_vars
USE vbc_arrays
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: hdry
INTEGER :: i,k,l,m
REAL (KIND=mp) :: eps,f(3),hleft,hright,uleft,uright,vleft,vright
TYPE(edge) :: ework
INTEGER, EXTERNAL :: edge_in_element,local_edge
REAL (KIND=mp), EXTERNAL :: invert_edge_dirs
DO l = 1,wall_edges1
k = walledg1(l)
i = edge_in_element(edges(k))
CALL edge_copy(edges(k),ework)
eps = one
IF (drycell(i)) CYCLE
m = local_edge(grid(i),k)
hleft = h(i) + delh(i)%x*rc(i,m)%x + delh(i)%y*rc(i,m)%y
IF (hleft < hdry) THEN
IF (hleft < zero) partdry(i) = .TRUE.
CYCLE
ELSE
uleft = u(i) + delu(i)%x*rc(i,m)%x + delu(i)%y*rc(i,m)%y
vleft = v(i) + delv(i)%x*rc(i,m)%x + delv(i)%y*rc(i,m)%y
END IF
IF (grid(i)%edge(m) < 0) eps = invert_edge_dirs(ework)
hright = hleft
vright = uleft*ework%tang(1) + vleft*ework%tang(2)
uright = vright*ework%tang(1)
vright = vright*ework%tang(2)
CALL flux_ac(ework,hleft,uleft,vleft,hright,uright,vright,f)
phi(i,2) = phi(i,2) - f(2)*ework%length
phi(i,3) = phi(i,3) - f(3)*ework%length
END DO
END SUBROUTINE weslp_fluxes
SUBROUTINE wetdryfix(itype,hv,hc,gradhc)
USE parameters
USE constants
USE geometry
USE options
IMPLICIT NONE
REAL (KIND=mp), INTENT(IN) :: hc(n_elems)
REAL (KIND=mp), INTENT(INOUT) :: hv(n_pts)
TYPE(vector), INTENT(IN) :: gradhc(n_elems)
INTEGER, INTENT(IN) :: itype
INTEGER :: counter,i,j,k
REAL (KIND=mp) :: newh
LOGICAL :: dryt,wett
SELECT CASE (itype)
CASE (1)
DO i = 1,n_pts
wett = .FALSE.
dryt = .FALSE.
DO j = 2,n2t(i,1) + 1
k = n2t(i,j)
IF (drycell(k)) THEN
dryt = .TRUE.
ELSE
wett = .TRUE.
END IF
END DO
IF (wett .AND. dryt) THEN
counter = 0
newh = zero
DO j = 2,n2t(i,1) + 1
k = n2t(i,j)
IF (drycell(k)) CYCLE
counter = counter + 1
newh = newh + hc(k) + (nodes(i)%x - grid(k)%xc)*gradhc(k)%x + &
(nodes(i)%y - grid(k)%yc)*gradhc(k)%y
END DO
hv(i) = newh/REAL(counter,mp)
END IF
END DO
CASE (2)
DO i = 1,n_pts
wett = .FALSE.
dryt = .FALSE.
DO j = 2,n2t2(i,1) + 1
k = n2t2(i,j)
IF (drycell(k)) THEN
dryt = .TRUE.
ELSE
wett = .TRUE.
END IF
END DO
IF (wett .AND. dryt) THEN
counter = 0
newh = zero
DO j = 2,n2t(i,1) + 1
k = n2t(i,j)
IF (drycell(k)) CYCLE
counter = counter + 1
newh = newh + hc(k) + (nodes(i)%x - grid(k)%xc)*gradhc(k)%x + &
(nodes(i)%y - grid(k)%yc)*gradhc(k)%y
END DO
hv(i) = newh/MAX(one,REAL(counter,mp))
END IF
END DO
CASE DEFAULT
PRINT *,"ERROR: invalid option in wetdryfix."
CALL byebye('Program STORM stopped.')
END SELECT
END SUBROUTINE wetdryfix
SUBROUTINE wetdryvfix
USE parameters
USE constants
USE geometry
USE options
USE dep_vars
IMPLICIT NONE
INTEGER :: i,j,k
LOGICAL :: dryt,wett
DO i = 1,n_pts
wett = .FALSE.
dryt = .FALSE.
DO j = 2,n2t(i,1) + 1
k = n2t(i,j)
IF (drycell(k)) THEN
dryt = .TRUE.
ELSE
wett = .TRUE.
END IF
END DO
IF (wett .AND. dryt) THEN
uvtx(i) = zero
vvtx(i) = zero
END IF
END DO
END SUBROUTINE wetdryvfix
SUBROUTINE windforcing
USE options
USE geometry
USE dep_vars
USE constants
IMPLICIT NONE
INTEGER :: f,i,ierror,lineno,n
REAL (KIND=mp) :: windmag
REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: tmp1,tmp2,tmp3,wcoef0
CHARACTER(LEN=40) :: buffer
LOGICAL, EXTERNAL :: get_iounit,readline
IF (.NOT. get_iounit(f)) THEN
WRITE (*,'(" ERROR: unable to open file ",A)') TRIM(wind_file)
CALL byebye('SToRM stopped.')
END IF
OPEN (f,FILE=wind_file,STATUS='OLD',IOSTAT=ierror)
IF (ierror /= 0) THEN
WRITE (*,'(" ERROR: unable to open file ",A)') TRIM(wind_file)
CALL byebye('SToRM stopped.')
END IF
lineno = 0
IF (readline(f,buffer,lineno)) THEN
READ (buffer,*) n
ELSE
WRITE (*,'(" Error reading file ",A)') TRIM(wind_file)
WRITE (*,'(" in line",I4)') lineno
CALL byebye('SToRM stopped.')
END IF
ALLOCATE(tmp1(n),tmp2(n),tmp3(n),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
DO i = 1,n
IF (readline(f,buffer,lineno)) THEN
READ (buffer,*) tmp1(i),tmp2(i),tmp3(i)
ELSE
WRITE (*,'(" Error reading file ",A)') TRIM(f_fric_facts)
WRITE (*,'(" in line",I4)') lineno
CALL byebye('SToRM stopped.')
END IF
END DO
SELECT CASE (opt_solver)
CASE (1)
ALLOCATE(wcoef0(n_pts),wcoef1(n_pts),wcoef2(n_pts),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(8*n_pts*2)
IF (n == n_elems) THEN
CALL ctr2vtx(nctr2vtx,tmp1,wcoef0)
CALL ctr2vtx(nctr2vtx,tmp2,wcoef1)
CALL ctr2vtx(nctr2vtx,tmp3,wcoef2)
ELSE IF (n == n_pts) THEN
wcoef0 = tmp1
wcoef1 = tmp2
wcoef2 = tmp3
ELSE
WRITE (*,'(" Error reading file ",A)') TRIM(f_fric_facts)
WRITE (*,'(" Invalid value in line",I4)') lineno
CALL byebye('SToRM stopped.')
END IF
DO i = 1,n_pts
windmag = SQRT(wcoef1(i)*wcoef1(i) + wcoef2(i)*wcoef2(i))
wcoef1(i) = wcoef0(i)*wcoef1(i)*windmag/rho
wcoef2(i) = wcoef0(i)*wcoef2(i)*windmag/rho
END DO
CASE (2,3)
ALLOCATE(wcoef0(n_elems),wcoef1(n_elems),wcoef2(n_elems),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(8*n_elems*2)
IF (n == n_elems) THEN
wcoef0 = tmp1
wcoef1 = tmp2
wcoef2 = tmp3
ELSE IF (n == n_pts) THEN
CALL vtx2ctr(tmp1,wcoef0)
CALL vtx2ctr(tmp2,wcoef1)
CALL vtx2ctr(tmp3,wcoef2)
ELSE
WRITE (*,'(" Error reading file ",A)') TRIM(f_fric_facts)
WRITE (*,'(" Invalid value in line",I4)') lineno
CALL byebye('SToRM stopped.')
END IF
DO i = 1,n_elems
windmag = SQRT(wcoef1(i)*wcoef1(i) + wcoef2(i)*wcoef2(i))
wcoef1(i) = wcoef0(i)*wcoef1(i)*windmag/rho
wcoef2(i) = wcoef0(i)*wcoef2(i)*windmag/rho
END DO
END SELECT
DEALLOCATE(tmp1,tmp2,tmp3,wcoef0)
CLOSE (f)
END SUBROUTINE windforcing
SUBROUTINE windforcing2
USE parameters
USE constants
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER :: i,ierror
REAL (KIND=mp) :: a
REAL (KIND=mp), ALLOCATABLE, DIMENSION(:) :: wcoef0
REAL (KIND=mp), EXTERNAL :: azimuth2angle
ALLOCATE(wcoef0(n_elems),wcoef1(n_elems),wcoef2(n_elems),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL mem_add(8*n_elems*2)
CALL vtx2ctr(w_fric,wcoef0)
CALL vtx2ctr(w_dir,wcoef1)
CALL vtx2ctr(w_mag,wcoef2)
DO i = 1,n_elems
a = azimuth2angle(wcoef1(i))
wcoef1(i) = wcoef0(i)*wcoef2(i)*COS(a)*wcoef2(i)/rho
wcoef2(i) = wcoef0(i)*wcoef2(i)*SIN(a)*wcoef2(i)/rho
END DO
DEALLOCATE(wcoef0,w_fric,w_mag,w_dir)
END SUBROUTINE windforcing2
SUBROUTINE zeta_bcs
USE parameters
USE dep_vars
USE vbc_arrays
IMPLICIT NONE
INTEGER :: i,i0,j
DO i0 = 1,n_outflowbdr
IF (htype(i0) /= 1) CYCLE
DO i = 1,n_hbc(i0)
j = hbc_nodes(i,i0)
zetavtx(j) = hvtx(j) + zvtx(j)
END DO
END DO
END SUBROUTINE zeta_bcs
