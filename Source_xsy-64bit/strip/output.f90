SUBROUTINE add_edge_to_array(edg,edges_array,n)
USE parameters
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: n
TYPE(edge), INTENT(IN) :: edg
TYPE(edge), INTENT(INOUT) :: edges_array(n+1)
n = n + 1
edges_array(n) = edg
END SUBROUTINE add_edge_to_array
SUBROUTINE alloc_err(i)
IMPLICIT NONE
INTEGER, INTENT(IN) :: i
WRITE (*,'("PROGRAM FAILURE: allocation of memory unsuccessful.",/, &
& "This is a computer systems failure mostly likely due to",/, &
& "insufficient memory or lack of other computing resources.",/, &
& "Reducing the size of the SToRM run may help.",/, &
& "System error no.",I5,/, &
& "Program SToRM stopped!
STOP
END SUBROUTINE alloc_err
SUBROUTINE bcs(funit)
IMPLICIT NONE
INTEGER, INTENT(IN) :: funit
CALL ioflow(funit)
CALL walls
END SUBROUTINE bcs
SUBROUTINE cell_velocity
USE parameters
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER :: i
DO i = 1,n_elems
u_avg(i,1) = (h(grid(i)%vertex(1)) + h(grid(i)%vertex(2)) + &
h(grid(i)%vertex(3)))/3.0_mp
u_avg(i,2) = (u(grid(i)%vertex(1)) + u(grid(i)%vertex(2)) + &
u(grid(i)%vertex(3)))/3.0_mp
u_avg(i,3) = (v(grid(i)%vertex(1)) + v(grid(i)%vertex(2)) + &
v(grid(i)%vertex(3)))/3.0_mp
END DO
END SUBROUTINE cell_velocity
MODULE constants
USE parameters
IMPLICIT NONE
SAVE
REAL(KIND=mp), PARAMETER :: fmachp = EPSILON(1.0_mp)
REAL(KIND=mp), PARAMETER :: vsmall = TINY(1.0_mp)
REAL(KIND=mp), PARAMETER :: half = 0.50_mp
REAL(KIND=mp), PARAMETER :: one = 1.0_mp
REAL(KIND=mp), PARAMETER :: one_third = one/3.0_mp
REAL(KIND=mp), PARAMETER :: zero = 0.0_mp
REAL(KIND=mp), PARAMETER :: g = 9.81_mp  !
END MODULE constants
MODULE dep_vars
USE parameters
IMPLICIT NONE
SAVE
REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: h,u,v,z
REAL(KIND=mp), ALLOCATABLE, DIMENSION(:) :: cd,rough
REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: u_avg
REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: flux
REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: source
REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: residual
REAL(KIND=mp), ALLOCATABLE, DIMENSION(:,:) :: phi
INTEGER :: n_h,n_qin
INTEGER, ALLOCATABLE, DIMENSION (:) :: h_nodes,qin_edges
REAL(KIND=mp) :: hbc,qin
END MODULE dep_vars
INTEGER FUNCTION distribution()
USE parameters
USE geometry
USE dep_vars
USE constants
USE options
IMPLICIT NONE
INTEGER :: counter,i,j,k,m,n
REAL(KIND=mp) :: hcell,nxi,nyi,ucell,vcell
REAL(KIND=mp) :: kpls(3,3),kmns(3,3),kplus(3,3,3),kminus(3,3,3),r(3,3)
LOGICAL :: flag
LOGICAL, EXTERNAL :: lda,narrow,nonlinearb,psi
phi = zero  !
counter = 0
DO j = 1,n_elems
hcell = u_avg(j,1)
ucell = u_avg(j,2)
vcell = u_avg(j,3)
IF (hcell .LT. h_dry) CYCLE  !
IF ((ucell*ucell + vcell*vcell) .LT. u_stag*u_stag) CYCLE  !
DO i = 1,3
k = MAX(1,MOD(i+1,4))  !
nxi = edges(ABS(grid(j)%edge(k)))%normal(1)
nyi = edges(ABS(grid(j)%edge(k)))%normal(2)
IF (grid(j)%edge(k) > 0) THEN
nxi = -nxi
nyi = -nyi
END IF
CALL kmatrices(kpls,kmns,nxi,nyi,hcell,ucell,vcell)
DO m = 1,3
DO n = 1,3
kplus(i,m,n) = kpls(m,n)
kminus(i,m,n) = kmns(m,n)
END DO
END DO
END DO
SELECT CASE (opt_rdscheme)
CASE (0)  !
flag = narrow(kplus,kminus,j,r)
IF (.NOT. flag) THEN  !
flag = lda(kplus,j,r)
IF (.NOT. flag) counter = counter + 1  !
END IF
CASE (1)  !
flag = lda(kplus,j,r)
IF (.NOT. flag) THEN  !
flag = narrow(kplus,kminus,j,r)
IF (.NOT. flag) counter = counter + 1  !
END IF
CASE (2)  !
flag = nonlinearb(kplus,kminus,j,r)
IF (.NOT. flag) counter = counter + 1
CASE (3)  !
PRINT *,''
PRINT *,'ERROR: the PSI scheme not yet implemented.'
PRINT *,'Program SToRM stopped.'
STOP
CASE DEFAULT
PRINT *,''
PRINT *,'ERROR: invalid value in opt_rdscheme.'
PRINT *,'Program SToRM stopped.'
STOP
END SELECT
DO m = 1,3
DO n = 1,3
phi(grid(j)%vertex(m),n) = phi(grid(j)%vertex(m),n) + r(m,n)
END DO
END DO
END DO  !
distribution = counter
END FUNCTION distribution
FUNCTION edge_in_array(edg,edges_array,n)
USE parameters
IMPLICIT NONE
LOGICAL :: edge_in_array
INTEGER, INTENT(IN) :: n
TYPE(edge), INTENT(IN) :: edg,edges_array(n+1)
INTEGER :: i
edge_in_array = .FALSE.
IF (n > 0) THEN
DO i = 1,n
IF ((edg%p(1) == edges_array(i)%p(1)  .AND. &
edg%p(2) == edges_array(i)%p(2)) .OR.  &
(edg%p(2) == edges_array(i)%p(1)  .AND. &
edg%p(1) == edges_array(i)%p(2))) THEN
edge_in_array = .TRUE.
RETURN
END IF
END DO
END IF
END FUNCTION edge_in_array
INTEGER FUNCTION edge_in_element(e)
USE geometry
IMPLICIT NONE
INTEGER, INTENT(IN) :: e
INTEGER :: i,j
DO i = 1,n_elems
DO j = 1,3
IF (e == ABS(grid(i)%edge(j))) THEN
edge_in_element = i  !
RETURN
END IF
END DO
END DO
edge_in_element = 0  !
END FUNCTION edge_in_element
SUBROUTINE edge_normals(side,nside,pto,npto)
USE parameters
IMPLICIT NONE
INTEGER, INTENT(IN) :: npto   !
INTEGER, INTENT(IN) :: nside  !
TYPE(edge), DIMENSION(nside) :: side  !
TYPE(point), DIMENSION(npto) :: pto  !
REAL(KIND=mp) :: x1,x2,y1,y2
INTEGER :: i
DO i = 1,nside
x1 = pto(side(i)%p(1))%x ; y1 = pto(side(i)%p(1))%y
x2 = pto(side(i)%p(2))%x ; y2 = pto(side(i)%p(2))%y
side(i)%normal(1) = y2 - y1  !
side(i)%normal(2) = x1 - x2  !
END DO
END SUBROUTINE edge_normals
FUNCTION equals(a,b) !
USE parameters  !
USE constants   !
IMPLICIT NONE
LOGICAL :: equals
REAL(KIND=mp), INTENT(IN) :: a,b
REAL(KIND=mp) :: eps
eps = ABS(a)*fmachp   !
IF (eps == zero) THEN
eps = vsmall        !
END IF
IF (ABS(a-b) > eps) THEN
equals = .FALSE.    !
ELSE
equals = .TRUE.     !
ENDIF
END FUNCTION equals
SUBROUTINE err_opts(string,lineno)
IMPLICIT NONE
INTEGER, INTENT(IN) :: lineno
CHARACTER (LEN=*) :: string
PRINT *,''
PRINT *,'ERROR: illegal value in option ',string
PRINT *,'in line',lineno,' of the options file.'
PRINT *,'Program SToRM stopped.'
STOP
END SUBROUTINE err_opts
SUBROUTINE fe_flux
USE parameters
USE geometry
USE dep_vars
USE constants
IMPLICIT NONE
INTEGER :: i,j,p1,p2
REAL(KIND=mp) :: h1,h2,q1,q2,u1,u2,v1,v2
TYPE(edge) :: side
DO j = 1,flow_edges
i = flowedg(j)
side = edges(i)
p1 = side%p(1)
p2 = side%p(2)
h1 = h(p1) ; u1 = u(p1) ; v1 = v(p1)
h2 = h(p2) ; u2 = u(p2) ; v2 = v(p2)
flux(i,1) = half*(h1*u1 + h2*u2)*side%normal(1) + &
half*(h1*v1 + h2*v2)*side%normal(2)
q1 = h1*u1*u1 + half*h1*h1*g
q2 = h2*u2*u2 + half*h2*h2*g
flux(i,2) = half*(q1 + q2)*side%normal(1)
q1 = h1*u1*v1
q2 = h2*u2*v2
flux(i,2) = flux(i,2) + half*(q1 + q2)*side%normal(2)
flux(i,3) = half*(q1 + q2)*side%normal(1)
q1 = h1*v1*v1 + half*h1*h1*g
q2 = h2*v2*v2 + half*h2*h2*g
flux(i,3) = flux(i,3) + half*(q1 + q2)*side%normal(2)
END DO
DO j = 1,wall_edges  !
i = walledg(j)
flux(i,1) = zero
flux(i,2) = zero
flux(i,3) = zero
END DO
END SUBROUTINE fe_flux
SUBROUTINE find_edges
USE parameters
USE geometry
IMPLICIT NONE
INTEGER :: i,ierror,j,k
REAL(KIND=mp) :: delx,dely
TYPE(edge) :: edg
TYPE(edge), ALLOCATABLE, DIMENSION (:) :: tmp_edges
LOGICAL, EXTERNAL :: edge_in_array
ALLOCATE(tmp_edges(n_elems*3 + 1),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
n_edges = 0
DO i = 1,n_elems
DO j = 1,3
edg%p(1) = grid(i)%vertex(j)
edg%p(2) = grid(i)%vertex(MAX(1,MOD(j+1,4)))
IF (.NOT. edge_in_array(edg,tmp_edges,n_edges)) &
CALL add_edge_to_array(edg,tmp_edges,n_edges)
END DO
END DO
ALLOCATE(edges(n_edges),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
DO i = 1,n_edges
edges(i) = tmp_edges(i)
END DO
DEALLOCATE(tmp_edges,STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
DO i = 1,n_elems
DO j = 1,3
edg%p(1) = grid(i)%vertex(j)
edg%p(2) = grid(i)%vertex(MAX(1,MOD(j+1,4)))
grid(i)%edge(j) = 0  !
DO k = 1,n_edges
IF (edg%p(1) == edges(k)%p(1) .AND. &
edg%p(2) == edges(k)%p(2)) THEN
grid(i)%edge(j) = k
EXIT
ELSE IF (edg%p(1) == edges(k)%p(2) .AND. &
edg%p(2) == edges(k)%p(1)) THEN
grid(i)%edge(j) = -k
EXIT
END IF
END DO
IF (grid(i)%edge(j) == 0) THEN  !
WRITE (*,'("PROGRAM ERROR: edge data structure incomplete.",/, &
& "This is a fatal programming error.  Please contact the",/, &
& "SToRM programming and maintenance team at the USGS.",/, &
& "Program stopped!
STOP
END IF
END DO
END DO
CALL edge_normals(edges,n_edges,nodes,n_pts)
DO k = 1,n_edges
delx = nodes(edges(k)%p(1))%x - nodes(edges(k)%p(2))%x
dely = nodes(edges(k)%p(1))%y - nodes(edges(k)%p(2))%y
edges(k)%length = SQRT(delx*delx + dely*dely)
END DO
END SUBROUTINE find_edges
LOGICAL FUNCTION fnames(funit,fname)
IMPLICIT NONE
INTEGER, INTENT(IN) :: funit
CHARACTER (LEN=80), INTENT(OUT) :: fname
DO WHILE (.TRUE.)
fname = ''
READ (funit,'(A)',ERR=1,END=1) fname
fname = ADJUSTL(fname)
IF (fname(1:1) == '#') CYCLE
IF (LEN_TRIM(fname) == 0) CYCLE
fnames = .TRUE.
RETURN
END DO
1 fnames = .FALSE.
END FUNCTION fnames
MODULE geometry
USE parameters
IMPLICIT NONE
SAVE
INTEGER :: n_pts  !
TYPE(point), ALLOCATABLE, DIMENSION(:) :: nodes
INTEGER :: n_elems  !
TYPE(triangle), ALLOCATABLE, DIMENSION(:) :: grid
INTEGER :: n_edges  !
TYPE(edge), ALLOCATABLE, DIMENSION(:) :: edges
INTEGER :: flow_edges,wall_edges,wall_elems
INTEGER, ALLOCATABLE, DIMENSION(:) :: flowedg,walledg,wallelm
END MODULE geometry
SUBROUTINE header(id,filename)
USE options
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER, INTENT(IN) :: id  !
CHARACTER (LEN=*), INTENT(IN) :: filename
INTEGER :: i,ierror,value(8)
CHARACTER :: hour*10,today*8,zone*5  !
CHARACTER (LEN=3) :: month(12),rds(4)
CHARACTER (LEN=9) :: req(5)
CHARACTER (LEN=80) :: buffer
month = (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct', &
'Nov','Dec'/)
req = (/"Manning's",". . Chezy","Drag Coef","Colebrook","Strickler"/)
rds = (/". N","LDA",". B","PSI"/)
IF (output_head) THEN
WRITE (id,'(A)')" ____  _____  ___   ____   __  __"
WRITE (id,'(A)')"/ ___||_   _|/ _ \ |  _ \ |  \/  |"
WRITE (id,'(A)')"\___ \  | | | | | || |_) || |\/| |"
WRITE (id,'(A)')" ___) | | | | |_| ||  _ < | |  | |"
WRITE (id,'(A)')"|____/  |_|  \___/ |_| \_\|_|  |_|"
WRITE (id,*)
WRITE (id,'(A)')"System for Transport and River Modeling [Version xo9]"
WRITE (id,'(A)')"U.S. Geological Survey"
WRITE (id,'(A)')"National Research Program"
WRITE (id,'(A)')"Lakewood, CO 80225 - U.S.A."
CALL DATE_AND_TIME(today,hour,zone,value)
WRITE (id,'(/A,I2,1X,A,1X,I4)')"Run started at "//hour(1:2)//":"// &
& hour(3:4)//":"//hour(5:6)//" on ",value(3),month(value(2)),value(1)
WRITE (id,'(A,A)')"Run data file is ",TRIM(filename)
END IF
IF (output_file) THEN
WRITE (id,'(/A/21("-"))')"Contents of run file:"
OPEN (10,FILE=filename,STATUS='OLD',IOSTAT=ierror)
IF (ierror /= 0) THEN
WRITE (id,'(" ERROR: unable to open file ",A)') TRIM(filename)
WRITE (id,'(" SToRM stopped.")')
STOP
END IF
DO
READ (10,'(A)',IOSTAT=ierror) buffer
IF (ierror /= 0) EXIT
WRITE (id,'(A)') TRIM(buffer)
END DO
CLOSE (10)
END IF
IF (output_opts) THEN
WRITE (id,'(/A/14("-"))')"Print options:"
WRITE (id,'(A,32(" ."),L)')"Print run file",output_file
WRITE (id,'(A,33(" ."),L)')"Print header",output_head
WRITE (id,'(A,32(" ."),L)')"Print options.",output_opts
WRITE (id,'(A,27(" ."),L)')"Print computational grid",output_geom
WRITE (id,'(A,27(" ."),L)')"Print initial conditions",output_icnd
WRITE (id,'(/A/16("-"))')"Problem options:"
WRITE (id,'(A,23(" ."),ES12.5)')"Water temperature (oC)",temp
WRITE (id,'(A,26(" ."),1X,A9)')"Roughness equation",req(opt_friction + 1)
WRITE (id,'(A,24(" ."),1X,A3)')"Residual distribution scheme", &
rds(opt_rdscheme + 1)
WRITE (id,'(A,25(" ."),ES12.5)')"Dry cell threshold",h_dry
WRITE (id,'(A,21(" ."),ES12.5)')"Stagnation cell threshold.",u_stag
WRITE (id,'(A,32(" ."),L)')"Wall friction.",wfrict
WRITE (id,'(A,24(" ."),ES12.5)')"Wall friction value.",wfrictv
END IF
IF (output_geom) THEN
WRITE (id,'(/10X,A)')"Node coordinates"
WRITE (id,'(7X,A)')"N        X             Y"
WRITE (id,'(37("-"))')
DO i = 1,n_pts
WRITE (id,'(I8,2(2X,ES12.5))') i,nodes(i)
END DO
WRITE (id,'(/22X,A)')"Control volume connectivity table"
WRITE (id,'(A)')"        N        I        J        K       E1       &
&E2       E3       AREA"
WRITE (id,'(78("-"))')
DO i = 1,n_elems
WRITE (id,'(7I9,2X,ES12.5)') i,grid(i)
END DO
WRITE (id,'(/21X,A)')"Edge data base"
WRITE (id,'(A)')"        N       P1       P2       NX            NY"
WRITE (id,'(56("-"))')
DO i = 1,n_edges
WRITE (id,'(3I9,2(2X,ES12.5))') i,edges(i)
END DO
END IF
IF (output_icnd) THEN
WRITE (id,'(/24X,A)')"Initial conditions"
WRITE (id,'(8X,"N",7X,"U",13X,"V",13X,"H",13X,"Z",13X,"CD")')
WRITE (id,'(80("-"))')
DO i = 1,n_pts
WRITE (id,'(I9,5(2X,ES12.5))') i,u(i),v(i),h(i),z(i),cd(i)
END DO
END IF
END SUBROUTINE header
SUBROUTINE initialize_options
USE constants
USE options
IMPLICIT NONE
REAL(mp), EXTERNAL :: viscosity
output_file = .TRUE.  !
output_geom = .TRUE.  !
output_head = .TRUE.  !
output_icnd = .TRUE.  !
output_opts = .TRUE.  !
opt_friction = 0  !
temp = REAL(18,mp)
visc = viscosity(temp)
opt_rdscheme = 0  !
h_dry = 0.001_mp  !
u_stag = 0.00001_mp  !
wfrict = .FALSE.  !
wfrictv = zero
END SUBROUTINE initialize_options
LOGICAL FUNCTION invert(a)
USE parameters
USE constants
IMPLICIT NONE
REAL(KIND=mp), INTENT(INOUT) :: a(3,3)
INTEGER :: i,ix,iy,j,k,m,mx(3),my(3),nx,ny
REAL(KIND=mp) :: b(3,3),dummy,factor,h,pivo,s1,s2
b = a
mx = 0
my = 0
DO i = 1,3
pivo = zero
DO ix = 1,3
IF (mx(ix) == 0) THEN
DO iy = 1,3
IF (my(iy) == 0) THEN
IF (ABS(b(ix,iy)) > ABS(pivo) ) THEN
pivo = b(ix,iy)
nx = ix
ny = iy
ENDIF
ENDIF
END DO
ENDIF
END DO
IF (ABS(pivo) <= 4.0_mp*fmachp) THEN
invert = .FALSE.
RETURN
ENDIF
mx(nx) = ny
my(ny) = nx
dummy = 1.0d0 / pivo
DO j = 1,3
IF (j /= nx) THEN
factor = b(j,ny)*dummy
DO k = 1,3
b(j,k) = b(j,k) - b(nx,k)*factor
END DO
b(j,ny) = -factor
ENDIF
END DO
DO k = 1,3
b(nx,k) = b(nx,k)*dummy
END DO
b(nx,ny) = dummy
END DO
DO i = 1,2
DO m = i,3
IF (mx(m) == i) EXIT
END DO
j = m
IF (j /= i) THEN
DO k = 1,3
h = b(i,k)
b(i,k) = b(j,k)
b(j,k) = h
END DO
mx(j) = mx(i)
mx(i) = i
ENDIF
DO m = i,3
IF (my(m) == i) EXIT
END DO
j = m
IF (j /= i) THEN
DO k = 1,3
h = b(k,i)
b(k,i) = b(k,j)
b(k,j) = h
END DO
my(j) = my(i)
my(i) = i
ENDIF
END DO
s1 = zero
s2 = zero
DO i = 1,3
DO j = 1,3
h = zero
DO k = 1,3
h = h + a(i,k)*b(k,j)
END DO
IF (i == j) THEN
s1 = s1 + ABS(h - one)
ELSE
s2 = s2 + ABS(h)
ENDIF
END DO
END DO
a = b
invert = .TRUE.
END FUNCTION invert
LOGICAL FUNCTION in_array(i,a)
IMPLICIT NONE
INTEGER, INTENT(IN) :: i
INTEGER, DIMENSION(:), INTENT(IN) :: a
INTEGER :: j,n
n = SIZE(a)
DO j = 1,n
IF (a(j) == i) THEN
in_array = .TRUE.
RETURN
END IF
END DO
in_array = .FALSE.
END FUNCTION in_array
MODULE io
IMPLICIT NONE
SAVE
INTEGER :: o_unit = 12
END MODULE io
SUBROUTINE ioflow(funit)
USE dep_vars
IMPLICIT NONE
INTEGER, INTENT(IN) :: funit
INTEGER :: counter,e,i,ierror,line,p,p1,p2
INTEGER, EXTERNAL :: is_edge,readp
CHARACTER (LEN=40) :: buffer
LOGICAL :: stats
LOGICAL, EXTERNAL :: readline
line = 0
stats = readline(funit,buffer,line)
READ(buffer,*) qin  !
stats = readline(funit,buffer,line)
READ(buffer,*) n_qin  !
ALLOCATE(qin_edges(n_qin),STAT=ierror)  !
IF (ierror /= 0) CALL alloc_err(ierror)
counter = 0
DO WHILE (counter < n_qin)
stats = readline(funit,buffer,line)  !
IF (.NOT.stats) EXIT
p = readp(buffer)
DO WHILE (p > 0)  !
counter = counter + 1
qin_edges(counter) = p
p = readp(buffer)
if (counter >= n_qin) EXIT
END DO
END DO
DO i = 1,n_qin - 1
p1 = qin_edges(i)
p2 = qin_edges(i+1)
e = is_edge(p1,p2)
IF (e < 1) THEN  !
PRINT *,"ERROR: inconsistent information, inflow boundary conditions."
PRINT *,"Please check the data in the boundary conditions input file."
PRINT *,"Program SToRM stopped."
STOP
ELSE
qin_edges(i) = e
END IF
END DO
n_qin = n_qin - 1
stats = readline(funit,buffer,line)
READ(buffer,*) hbc  !
stats = readline(funit,buffer,line)
READ(buffer,*) n_h  !
ALLOCATE(h_nodes(n_h),STAT=ierror)  !
IF (ierror /= 0) CALL alloc_err(ierror)
counter = 0
DO WHILE (counter < n_h)
stats = readline(funit,buffer,line)  !
IF (.NOT.stats) EXIT
p = readp(buffer)
DO WHILE (p > 0)  !
counter = counter + 1
h_nodes(counter) = p
p = readp(buffer)
if (counter >= n_h) EXIT
END DO
END DO
END SUBROUTINE ioflow
INTEGER FUNCTION is_edge(p1,p2)
USE geometry
IMPLICIT NONE
INTEGER, INTENT(IN) :: p1,p2
INTEGER :: e,e1,e2
IF (p1 > n_pts .OR. p1 < 1 .OR. p2 > n_pts .OR. p2 < 1) THEN
is_edge = -1  !
RETURN
END IF
DO e = 1,n_edges
e1 = edges(e)%p(1)
e2 = edges(e)%p(2)
IF ((e1 == p1 .AND. e2 == p2) .OR. (e1 == p2 .AND. e2 == p1)) THEN
is_edge = e  !
RETURN
END IF
END DO
is_edge = 0  !
END FUNCTION is_edge
SUBROUTINE kmatrices(kplus,kminus,nxi,nyi,h,u,v)
USE parameters
USE constants
IMPLICIT NONE
REAL(KIND=mp), INTENT(IN) :: nxi,nyi,h,u,v
REAL(KIND=mp), INTENT(OUT) :: kminus(3,3),kplus(3,3)
REAL(KIND=mp) :: c,f1,f2,lambda(3),ni,ni2,uni
REAL(KIND=mp) :: ki(3,3)
uni = u*nxi + v*nyi
ni2 = nxi*nxi + nyi*nyi
ni = SQRT(ni2)
c = SQRT(g*h)
lambda(1) = uni + c*ni
lambda(2) = uni
lambda(3) = uni - c*ni
ki(1,1) = zero;            ki(1,2) = nxi;         ki(1,3) = nyi
ki(2,1) = c*c*nxi - u*uni; ki(2,2) = u*nxi + uni; ki(2,3) = u*nyi
ki(3,1) = c*c*nyi - v*uni; ki(3,2) = v*nxi;       ki(3,3) = v*nyi + uni
ki = ki*half
IF (lambda(3) >= zero) THEN  !
kplus = ki
kminus = zero
ELSE IF (lambda(1) <= zero) THEN  !
kminus = ki
kplus = zero
ELSE IF (lambda(1) > zero .AND. lambda(2) <= zero) THEN
kplus(1,1) = c - uni/ni; kplus(1,2) = nxi/ni; kplus(1,3) = nyi/ni
f1 = (u*ni + c*nxi)/ni2
f2 = (v*ni + c*nyi)/ni2
kplus(2,1) = f1*(c*ni - uni)
kplus(2,2) = f1*nxi
kplus(2,3) = f1*nyi
kplus(3,1) = f2*(c*ni - uni)
kplus(3,2) = f2*nxi
kplus(3,3) = f2*nyi
kplus = kplus*lambda(1)/c/4.0_mp
kminus = ki - kplus
ELSE IF (lambda(3) < zero .AND. lambda(2) >= zero) THEN
kminus(1,1) = c + uni/ni; kminus(1,2)= -nxi/ni; kminus(1,3) = -nyi/ni
f1 = (c*nxi - u*ni)/ni2
f2 = (c*nyi - v*ni)/ni2
kminus(2,1) = -f1*(uni + c*ni)
kminus(2,2) = f1*nxi
kminus(2,3) = f1*nyi
kminus(3,1) = -f2*(uni + c*ni)
kminus(3,2) = f2*nxi
kminus(3,3) = f2*nyi
kminus = kminus*lambda(3)/c/4.0_mp
kplus = ki - kminus
ELSE
PRINT *,""
PRINT *,"FATAL ERROR: eigenvalue selection failure in Ki decomposition."
PRINT *,"Please report error to SToRM's development team."
PRINT *,"Program SToRM stopped."
STOP
END IF
END SUBROUTINE kmatrices
LOGICAL FUNCTION lda(kplus,n,r)
USE parameters
USE constants
USE dep_vars
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
REAL(KIND=mp), INTENT(IN) :: kplus(3,3,3)
REAL(KIND=mp), INTENT(OUT) :: r(3,3)
INTEGER :: i,j,l,m
REAL(KIND=mp) :: nhat(3,3),beta(3,3),k(3,3)
LOGICAL :: flag
LOGICAL, EXTERNAL :: invert
r = zero
nhat = zero
DO i = 1,3
DO j = 1,3
nhat(i,j) = kplus(1,i,j) + kplus(2,i,j) + kplus(3,i,j)
END DO
END DO
flag = invert(nhat)
IF (.NOT. flag) THEN  !
lda = .FALSE.       !
RETURN
END IF
DO i = 1,3
DO l = 1,3
DO m = 1,3
k(l,m) = kplus(i,l,m)
END DO
END DO
beta = MATMUL(k,nhat)  !
DO j = 1,3  !
r(i,1) = r(i,1) + beta(1,j)*residual(n,j)
r(i,2) = r(i,2) + beta(2,j)*residual(n,j)
r(i,3) = r(i,3) + beta(3,j)*residual(n,j)
END DO
END DO
lda = .TRUE.
END FUNCTION lda
LOGICAL FUNCTION narrow(kplus,kminus,n,r)
USE parameters
USE constants
USE dep_vars
USE geometry
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
REAL(KIND=mp), INTENT(IN) :: kplus(3,3,3),kminus(3,3,3)
REAL(KIND=mp), INTENT(OUT) :: r(3,3)
INTEGER :: i,j,k
REAL(KIND=mp) :: kiwi(3),nhat(3,3),w(3),wi(3),win(3)
LOGICAL :: flag
LOGICAL, EXTERNAL :: invert
r = zero
nhat = zero
DO i = 1,3
DO j = 1,3
nhat(i,j) = kminus(1,i,j) + kminus(2,i,j) + kminus(3,i,j)
END DO
END DO
flag = invert(nhat)
IF (.NOT. flag) THEN  !
narrow = .FALSE.    !
RETURN
END IF
kiwi = zero  !
DO i = 1,3   !
k = grid(n)%vertex(i) !
wi(1) = h(k)
wi(2) = u(k)
wi(3) = v(k)
DO j = 1,3
kiwi(1) = kiwi(1) + kplus(i,1,j)*wi(j)
kiwi(2) = kiwi(2) + kplus(i,2,j)*wi(j)
kiwi(3) = kiwi(3) + kplus(i,3,j)*wi(j)
END DO
END DO
w(1) = kiwi(1) + residual(n,1)
w(2) = kiwi(2) + residual(n,2)
w(3) = kiwi(3) + residual(n,3)
win = zero
DO j = 1,3
win(1) = win(1) - nhat(1,j)*w(j)
win(2) = win(2) - nhat(2,j)*w(j)
win(3) = win(3) - nhat(3,j)*w(j)
END DO
DO i = 1,3
k = grid(n)%vertex(i)
w(1) = h(k) - win(1)  !
w(2) = u(k) - win(2)
w(3) = v(k) - win(3)
DO j = 1,3
r(i,1) = r(i,1) - kplus(i,1,j)*w(j)  !
r(i,2) = r(i,2) - kplus(i,2,j)*w(j)
r(i,3) = r(i,3) - kplus(i,3,j)*w(j)
END DO
END DO
narrow = .TRUE.
END FUNCTION narrow
LOGICAL FUNCTION nonlinearb(kplus,kminus,n,r)
USE parameters
USE constants
USE dep_vars
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
REAL(KIND=mp), INTENT(IN) :: kplus(3,3,3),kminus(3,3,3)
REAL(KIND=mp), INTENT(OUT) :: r(3,3)
INTEGER :: i,j
REAL(KIND=mp) :: coef(3),rlda(3,3),rnarrow(3,3),sum
LOGICAL :: flag
LOGICAL, EXTERNAL :: equals,lda,narrow
r = zero
flag = narrow(kplus,kminus,n,rnarrow)
IF (.NOT. flag) THEN  !
flag = lda(kplus,n,rlda)
IF (.NOT. flag) THEN
nonlinearb = .FALSE.
RETURN
END IF
r = rlda
nonlinearb = .TRUE.
RETURN
END IF
flag = lda(kplus,n,rlda)
IF (.NOT. flag) THEN  !
flag = narrow(kplus,kminus,n,rnarrow)
IF (.NOT. flag) THEN
nonlinearb = .FALSE.
RETURN
END IF
r = rnarrow
nonlinearb = .TRUE.
RETURN
END IF
DO i = 1,3
sum = zero
DO j = 1,3
sum = sum + ABS(rnarrow(j,i))
END DO
IF (equals(sum,zero)) THEN
nonlinearb = .FALSE.
RETURN
END IF
coef(i) = residual(n,i)/sum
END DO
DO i = 1,3
DO j = 1,3
r(i,j) = coef(j)*rnarrow(i,j) + (one - coef(j))*rlda(i,j)
END DO
END DO
nonlinearb = .TRUE.
END FUNCTION nonlinearb
MODULE options
USE parameters
IMPLICIT NONE
SAVE
LOGICAL :: output_file,output_geom,output_head,output_icnd,output_opts
INTEGER :: opt_friction
REAL(KIND=mp) :: temp,visc
INTEGER :: opt_rdscheme
REAL(KIND=mp) :: h_dry
REAL(KIND=mp) :: u_stag
REAL(KIND=mp) :: wfrictv
LOGICAL :: wfrict
END MODULE options
SUBROUTINE output_data
USE options
USE geometry
USE dep_vars
USE io
IMPLICIT NONE
INTEGER :: value(8)
CHARACTER :: hour*10,today*8,zone*5  !
CHARACTER (LEN=3) :: month(12)
month = (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct', &
'Nov','Dec'/)
IF (output_head) THEN
CALL DATE_AND_TIME(today,hour,zone,value)
WRITE (o_unit,'(/A,I2,1X,A,1X,I4)')"Run ended at "//hour(1:2)//":"// &
& hour(3:4)//":"//hour(5:6)//" on ",value(3),month(value(2)),value(1)
END IF
CLOSE (o_unit)
END SUBROUTINE output_data
MODULE parameters
IMPLICIT NONE
SAVE
INTEGER, PARAMETER :: mp = KIND(1.0D0)  !
TYPE :: point
REAL(KIND=mp) :: x,y  !
END TYPE point
TYPE :: edge
INTEGER :: p(2)  !
REAL(KIND=mp) :: normal(2)  !
REAL(KIND=mp) :: length     !
END TYPE
TYPE :: triangle
INTEGER :: vertex(3)        !
INTEGER :: edge(3)          !
REAL(KIND=mp) :: area       !
END TYPE
END MODULE parameters
LOGICAL FUNCTION readline(funit,line,lineno)
IMPLICIT NONE
INTEGER, INTENT(IN) :: funit  !
INTEGER, INTENT(INOUT) :: lineno  !
CHARACTER(LEN=40), INTENT(OUT) :: line  !
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
RETURN  !
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
SUBROUTINE read_data
USE parameters
USE geometry
USE io
USE dflib !
IMPLICIT NONE
INTEGER :: ierror
CHARACTER (LEN=80) :: bcs_file,buffer,filename,grid_file,options_file, &
out_file
LOGICAL :: fnexist,options,stats
LOGICAL, EXTERNAL :: fnames
IF (NARGS() == 1) THEN
10  WRITE (*,'(" Enter the INPUT file name: ",$)') !
READ (*,'(A)') filename                        !
filename = ADJUSTL(filename)
INQUIRE (file=filename,EXIST=fnexist)
IF (.NOT.fnexist) THEN
WRITE (*,'(2A/)') 'ERROR: file not found, ', TRIM(filename)
WRITE (*,'("Press <Enter> to continue")'); READ (*,'(A)') buffer
GO TO 10
END IF
ELSE
CALL GETARG(1,filename)
filename = ADJUSTL(filename)
END IF
OPEN (10,FILE=filename,STATUS='OLD',IOSTAT=ierror)
IF (ierror /= 0) THEN
WRITE (*,'(" ERROR: unable to open file ",A)') TRIM(filename)
WRITE (*,'(" SToRM stopped.")')
STOP
END IF
stats = fnames(10,grid_file)
stats = fnames(10,bcs_file)
stats = fnames(10,out_file)
options = fnames(10,options_file)
CLOSE (10)
INQUIRE (file=grid_file,EXIST=fnexist)
IF (.NOT.fnexist) THEN
WRITE (*,'(A)') 'ERROR: file not found, ', TRIM(grid_file)
WRITE (*,'(A)') 'SToRM stopped.'
STOP
END IF
OPEN (11,FILE=grid_file,STATUS='OLD',IOSTAT=ierror)
IF (ierror /= 0) THEN
WRITE (*,'(" ERROR: unable to open grid file ",A)') TRIM(grid_file)
WRITE (*,'(" SToRM stopped.")')
STOP
END IF
CALL tecplot(11)
CLOSE (11)
CALL find_edges
CALL t_areas(grid,n_elems,nodes,n_pts)
CALL initialize_options
IF (options) THEN
OPEN (15,FILE=options_file,STATUS='OLD',IOSTAT=ierror)
IF (ierror /= 0) THEN
WRITE (*,'(" ERROR: unable to open options file ",A)') TRIM(options_file)
WRITE (*,'(" SToRM run will proceed with all default values.")')
ELSE
CALL read_options(15)
END IF
END IF
INQUIRE (file=bcs_file,EXIST=fnexist)
IF (.NOT.fnexist) THEN
WRITE (*,'(A)') 'ERROR: file not found, ', TRIM(bcs_file)
WRITE (*,'(A)') 'SToRM stopped.'
STOP
END IF
OPEN (12,FILE=grid_file,STATUS='OLD',IOSTAT=ierror)
IF (ierror /= 0) THEN
WRITE (*,'(" ERROR: unable to open boundary conditions file ",A)') &
TRIM(bcs_file)
WRITE (*,'(" SToRM stopped.")')
STOP
END IF
CALL bcs(12)
CLOSE (12)
OPEN (o_unit,FILE=out_file,IOSTAT=ierror)
IF (ierror /= 0) THEN
WRITE (*,'(" ERROR: unable to open output file ",A)') TRIM(out_file)
WRITE (*,'(" SToRM stopped.")')
STOP
END IF
CALL header(o_unit,filename)
END SUBROUTINE read_data
SUBROUTINE read_options(funit)
USE parameters
USE options
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: funit
INTEGER :: c,i,lineno,n
CHARACTER (LEN=40) :: line
LOGICAL :: flag
REAL(KIND=mp), EXTERNAL :: viscosity
CHARACTER, EXTERNAL :: to_upper
LOGICAL, EXTERNAL :: equals,readline
lineno = 0
DO WHILE (.TRUE.)
flag = readline(funit,line,lineno)
IF (.NOT. flag) EXIT  !
IF (line(1:12) == "TEMPERATURE") THEN
n = LEN_TRIM(line)
line(n:n) = to_upper(line(n:n))
c = IACHAR(line(n:n))  !
IF ((c >= 48 .AND. c <=57) .OR. c == 46) THEN  !
READ(line(13:40),*) temp
ELSE IF (c == 67) THEN  !
READ(line(13:n-1),*) temp
ELSE IF (c == 70) THEN  !
READ(line(13:n-1),*) temp
temp = REAL(5,mp)*(temp - REAL(32,mp))/REAL(9,mp)
ELSE
PRINT *,''
PRINT *,'ERROR: illegal value in option temperature (units)'
PRINT *,'in line',lineno,' of the options file.'
PRINT *,'Program SToRM stopped.'
STOP
END IF
visc = viscosity(temp)
IF (equals(visc,zero)) CALL err_opts('temperature',lineno)
ELSE IF (line(1:9) == "RD SCHEME") THEN
READ(line(10:40),*) opt_rdscheme
IF (opt_rdscheme < 0 .OR. opt_rdscheme > 3) &
CALL err_opts('rd scheme',lineno)
ELSE IF (line(1:9) == "ROUGHNESS") THEN
READ(line(10:40),*) opt_friction
IF (opt_friction < 0 .OR. opt_friction > 4) &
CALL err_opts('roughness',lineno)
ELSE IF (line(1:13) == "WALL FRICTION") THEN
READ(line(14:40),*) wfrictv
wfrict = .TRUE.
IF (wfrictv <= zero) CALL err_opts('wall friction',lineno)
ELSE IF (line(1:10) == "MASTERFILE") THEN
READ(line(11:40),*) i
IF (i <= 0) THEN
output_file = .FALSE.
ELSE
output_file = .TRUE.
END IF
ELSE IF (line(1:8) == "GEOMETRY") THEN
READ(line(9:40),*) i
IF (i <= 0) THEN
output_geom = .FALSE.
ELSE
output_geom = .TRUE.
END IF
ELSE IF (line(1:6) == "HEADER") THEN
READ(line(7:40),*) i
IF (i <= 0) THEN
output_head = .FALSE.
ELSE
output_head = .TRUE.
END IF
ELSE IF (line(1:14) == "INITINAL CONDS") THEN
READ(line(15:40),*) i
IF (i <= 0) THEN
output_icnd = .FALSE.
ELSE
output_icnd = .TRUE.
END IF
ELSE IF (line(1:7) == "OPTIONS") THEN
READ(line(8:40),*) i
IF (i <= 0) THEN
output_opts = .FALSE.
ELSE
output_opts = .TRUE.
END IF
ELSE IF (line(1:11) == "H THRESHOLD") THEN
READ(line(12:40),*) h_dry
h_dry = MIN(one,ABS(h_dry))
ELSE IF (line(1:11) == "U THRESHOLD") THEN
READ(line(12:40),*) u_stag
u_stag = MIN(one,ABS(u_stag))
END IF
END DO
END SUBROUTINE read_options
SUBROUTINE residuals
USE parameters
USE dep_vars
USE geometry
USE constants
IMPLICIT NONE
INTEGER :: i,j,k
CALL fe_flux  !
CALL sources_slope
DO i = 1,n_elems
DO k = 1,3
residual(i,k) = -flux(ABS(grid(i)%edge(1)),k)* &
SIGN(one,REAL(grid(i)%edge(1),mp))
END DO
DO j = 2,3
DO k = 1,3
residual(i,k) = residual(i,k) - flux(ABS(grid(i)%edge(j)),k)* &
SIGN(one,REAL(grid(i)%edge(j),mp))
END DO
END DO
DO k = 1,3
residual(i,k) = residual(i,k) + source(i,k)
END DO
END DO
END SUBROUTINE
SUBROUTINE solver
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER :: counter,ierror
INTEGER, EXTERNAL :: distribution
ALLOCATE(residual(n_elems,3),phi(n_pts,3),flux(n_edges,3), &
source(n_elems,3),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
ALLOCATE(u_avg(n_elems,3),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
CALL enforce_bcs
CALL cell_velocity
CALL residuals
counter = distribution()
END SUBROUTINE solver
SUBROUTINE sources_slope
USE parameters
USE geometry
USE dep_vars
USE constants
USE options
IMPLICIT NONE
INTEGER :: i
REAL(KIND=mp) :: f(3),h_cv,h_wall,ln,r,re,s0x,s0y,u_vel,v_vel,x1,x2,x3, &
y1,y2,y3
TYPE(triangle) :: cv
source = zero
DO i = 1,n_elems
IF (u_avg(i,1) < h_dry) CYCLE
IF ((u_avg(i,2)*u_avg(i,2) + u_avg(i,3)*u_avg(i,3)) < u_stag*u_stag) CYCLE
cv = grid(i)  !
x1 = nodes(cv%vertex(1))%x; y1 = nodes(cv%vertex(1))%y
x2 = nodes(cv%vertex(2))%x; y2 = nodes(cv%vertex(2))%y
x3 = nodes(cv%vertex(3))%x; y3 = nodes(cv%vertex(3))%y
s0x = -half*(z(cv%vertex(1))*(y2 - y3) + z(cv%vertex(2))*(y3 - y1) + &
z(cv%vertex(3))*(y1 - y2))
s0y = -half*(z(cv%vertex(1))*(x3 - x2) + z(cv%vertex(2))*(x1 - x3) + &
z(cv%vertex(3))*(x2 - x1))
h_cv = u_avg(i,1)
source(i,1) = zero        !
source(i,2) = g*h_cv*s0x  !
source(i,3) = g*h_cv*s0y  !
END DO
SELECT CASE (opt_friction)
CASE (0)  !
DO i = 1,n_pts
IF (h(i) < h_dry) THEN
rough(i) = zero
ELSE
rough(i) = g*cd(i)*cd(i)/h(i)**one_third
END IF
END DO
CASE (1)  !
DO i = 1,n_pts
rough(i) = g/(cd(i)*cd(i))
END DO
CASE (2)  !
rough = cd
CASE (3)  !
DO i = 1,n_pts
re = SQRT(u(i)*u(i) + v(i)*v(i))*h(i)/visc  !
IF (h(i) < h_dry) THEN
rough(i) = zero
ELSE
re = MAX(re,1000.0_mp)  !
re = MIN(re,REAL(2.5e7,mp))  !
rough(i) = cd(i)
IF (cd(i)/h(i) > 0.20_mp) rough(i) = 0.20_mp*h(i)
ln = LOG(1.725_mp/re + (rough(i)/(14.8_mp*h(i)))**1.11)
rough(i) = 0.204_mp/(ln*ln)
END IF
END DO
CASE (4)  !
rough = ((100.0_mp*cd)**0.16666666666666667_mp)/51.79_mp
DO i = 1,n_pts
IF (h(i) < h_dry) THEN
rough(i) = zero
ELSE
rough(i) = g*rough(i)*rough(i)/h(i)**one_third
END IF
END DO
CASE DEFAULT
PRINT *,''
PRINT *,'ERROR: invalid value in opt_friction.'
PRINT *,'Program SToRM stopped.'
STOP
END SELECT
DO i = 1,n_pts
rough(i) = -rough(i)*SQRT(u(i)*u(i) + v(i)*v(i))
END DO
DO i = 1,n_elems
IF (u_avg(i,1) < h_dry) CYCLE
IF ((u_avg(i,2)*u_avg(i,2) + u_avg(i,3)*u_avg(i,3)) < u_stag*u_stag) CYCLE
cv = grid(i)  !
source(i,1) = source(i,1) + zero  !
f(1) = rough(cv%vertex(1))*u(cv%vertex(1))
f(2) = rough(cv%vertex(2))*u(cv%vertex(2))
f(3) = rough(cv%vertex(3))*u(cv%vertex(3))
source(i,2) = source(i,2) + cv%area*(f(1) + f(2) + f(3))*one_third
f(1) = rough(cv%vertex(1))*v(cv%vertex(1))
f(2) = rough(cv%vertex(2))*v(cv%vertex(2))
f(3) = rough(cv%vertex(3))*v(cv%vertex(3))
source(i,3) = source(i,3) + cv%area*(f(1) + f(2) + f(3))*one_third
END DO
IF (wfrict) THEN
DO i = 1,wall_edges
h_wall = half*(h(edges(walledg(i))%p(1)) + h(edges(walledg(i))%p(2)))
u_vel = half*(u(edges(walledg(i))%p(1)) + u(edges(walledg(i))%p(2)))
v_vel = half*(v(edges(walledg(i))%p(1)) + v(edges(walledg(i))%p(2)))
IF (h_wall < h_dry) CYCLE  !
IF (u_vel*u_vel + v_vel*v_vel < u_stag*u_stag) CYCLE  !
SELECT CASE (opt_friction)
CASE (0)  !
r = -g*h_wall*wfrictv*wfrictv* &
(edges(walledg(i))%length/grid(wallelm(i))%area)**(4.0_mp/3.0_mp)
CASE (1)  !
r = -g*h_wall*edges(walledg(i))%length/grid(wallelm(i))%area/ &
(wfrictv*wfrictv)
CASE (2)  !
r = -wfrictv
CASE (3)  !
re = SQRT(u_vel*u_vel + v_vel*v_vel)*h_wall/visc
re = MAX(re,1000.0_mp)
re = MIN(re,REAL(2.5e7,mp))
r = wfrictv
IF (r/h_wall > 0.20_mp) r = 0.20_mp*h_wall
ln = LOG(1.725_mp/re + (r/(14.8_mp*h(i)))**1.11)
r = -0.204_mp/(ln*ln)
CASE (4)  !
r = (100.0_mp*wfrictv)**(one/6.0_mp)/51.79_mp  !
r = -g*h_wall*r*r* &
(edges(walledg(i))%length/grid(wallelm(i))%area)**(4.0_mp/3.0_mp)
END SELECT
r = r*SQRT(u_vel*u_vel + v_vel*v_vel)
source(wallelm(i),2) = source(wallelm(i),2) + r*u_vel
source(wallelm(i),3) = source(wallelm(i),3) + r*v_vel
END DO
END IF
END SUBROUTINE sources_slope
PROGRAM SToRM
IMPLICIT NONE
PRINT *,''
PRINT *,' _______ _______  _____   ______ _______'
PRINT *,' |______    |    |     | |_____/ |  |  |'
PRINT *,' ______|    |    |_____| |    \_ |  |  |'
PRINT *,''
PRINT *,' System for Transport and River Modeling'
PRINT *,' U.S. Geological Survey (NRP)'
PRINT *,' Lakewood, CO 80225 - USA'
PRINT *,' Version xo9 - 2004'
PRINT *,''
CALL read_data
CALL solver
CALL output_data
END PROGRAM SToRM
SUBROUTINE tecplot(id)
USE parameters
USE geometry
USE dep_vars
IMPLICIT NONE
INTEGER, INTENT(IN) :: id  !
INTEGER :: i,ierror,j
CHARACTER (LEN=80) :: buffer
DO i = 1,3
READ (id,'(A)') buffer
END DO
i = INDEX(buffer,'N=',.FALSE.)
j = INDEX(buffer,',',.FALSE.)
READ (buffer(i+2:j-1),*) n_pts
i = INDEX(buffer,'E=',.FALSE.)
j = INDEX(buffer(i:),',',.FALSE.)
READ (buffer(i+2:i+j-1),*) n_elems
ALLOCATE(nodes(n_pts),STAT=ierror)  !
IF (ierror /= 0) CALL alloc_err(ierror)
ALLOCATE(grid(n_elems),STAT=ierror)  !
IF (ierror /= 0) CALL alloc_err(ierror)
ALLOCATE(z(n_pts),STAT=ierror)  !
IF (ierror /= 0) CALL alloc_err(ierror)
ALLOCATE(h(n_pts),u(n_pts),v(n_pts),STAT=ierror)  !
IF (ierror /= 0) CALL alloc_err(ierror)
ALLOCATE(cd(n_pts),rough(n_pts),STAT=ierror)  !
IF (ierror /= 0) CALL alloc_err(ierror)
DO i = 1,n_pts
READ (id,*) nodes(i)%x,nodes(i)%y,z(i),h(i),u(i),v(i),cd(i)
END DO
DO i = 1,n_elems
READ (id,*) grid(i)%vertex(1),grid(i)%vertex(2),grid(i)%vertex(3)
END DO
END SUBROUTINE tecplot
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
SUBROUTINE t_areas(grid,ngrid,pto,npto)
USE parameters
USE constants
IMPLICIT NONE
INTEGER, INTENT(IN) :: ngrid  !
INTEGER, INTENT(IN) :: npto   !
TYPE(triangle), DIMENSION(ngrid) :: grid  !
TYPE(point), DIMENSION(npto) :: pto  !
TYPE(triangle) :: elem
REAL(KIND=mp) :: area,x1,x2,x3,y1,y2,y3
INTEGER :: i
DO i = 1,ngrid
elem = grid(i)  !
x1 = pto(elem%vertex(1))%x ; y1 = pto(elem%vertex(1))%y
x2 = pto(elem%vertex(2))%x ; y2 = pto(elem%vertex(2))%y
x3 = pto(elem%vertex(3))%x ; y3 = pto(elem%vertex(3))%y
area = half*ABS((x1 - x2)*(y1 + y2) + (x2 - x3)*(y2 + y3) + &
(x3 - x1)*(y3 + y1))
grid(i)%area = area
END DO
END SUBROUTINE t_areas
FUNCTION viscosity(t)
USE parameters
USE constants
IMPLICIT NONE
REAL(KIND=mp), INTENT(IN) :: t  !
REAL(KIND=mp) :: viscosity  !
viscosity = zero  !
IF (t > 0.0_mp .AND. t < 100.0_mp) THEN
viscosity = REAL(1.792e-6,mp)/(one + 0.0337_mp*t + 0.000221_mp*t*t)
END IF
END FUNCTION viscosity
SUBROUTINE walls
USE geometry
USE dep_vars
USE options
IMPLICIT NONE
INTEGER :: edge_count(n_edges),i,ierror,j,k,n
INTEGER, EXTERNAL :: edge_in_element
LOGICAL, EXTERNAL :: in_array
edge_count = 0
DO n = 1,n_elems
DO j = 1,3
i = ABS(grid(n)%edge(j))
edge_count(i) = edge_count(i) + 1
END DO
END DO
wall_edges = 0  !
flow_edges = 0  !
DO i = 1,n_edges
IF (edge_count(i) == 1) THEN
IF (.NOT.in_array(i,qin_edges)) wall_edges = wall_edges + 1
ELSE IF (edge_count(i) == 2) THEN
flow_edges = flow_edges + 1
ELSE
PRINT *,'ERROR: computational grid construction error, in edge'
PRINT *,'assignement (wall or inner edge). Please check the grid'
PRINT *,'for compliance with the definitions. SToRM stopped.'
STOP
END IF
END DO
IF (wall_edges + flow_edges /= n_edges) THEN  !
PRINT *,'ERROR: computational grid construction error, in edge'
PRINT *,'assignement (wall or inner edge). Please check the grid'
PRINT *,'for compliance with the definitions. SToRM stopped.'
STOP
END IF
wall_edges = MAX(1,wall_edges)
flow_edges = MAX(1,flow_edges)
ALLOCATE(walledg(wall_edges),flowedg(flow_edges),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
j = 0
k = 0
DO i = 1,n_edges
IF (edge_count(i) == 1) THEN
IF (.NOT.in_array(i,qin_edges)) THEN
j = j + 1
walledg(j) = i
END IF
ELSE IF (edge_count(i) == 2) THEN
k = k + 1
flowedg(k) = i
END IF
END DO
IF (j /= wall_edges .OR. k /= flow_edges) THEN
PRINT *,'ERROR: computational grid construction error, in edge'
PRINT *,'assignement (wall or inner edge). Please check the grid'
PRINT *,'for compliance with the definitions. SToRM stopped.'
STOP
END IF
IF (wfrict) THEN
wall_elems = wall_edges
ALLOCATE(wallelm(wall_elems),STAT=ierror)
IF (ierror /= 0) CALL alloc_err(ierror)
DO i = 1,wall_elems
wallelm(i) = edge_in_element(walledg(i))
IF (wallelm(i) == 0) THEN
PRINT *,'ERROR: computational grid construction error, in an edge'
PRINT *,'association to an element. Please check the grid for'
PRINT *,'compliance with the definitions. SToRM stopped.'
STOP
END IF
END DO
END IF
END SUBROUTINE walls
