subroutine read_matrix(fname, dim, M)
  implicit none
  integer, intent(in) :: dim
  character(*), intent(in) :: fname
  real ( kind(0Q0) ), dimension(dim,dim), intent(inout) :: M
  
  integer ( kind = 4 ) i,j
  integer u, io
  
  u = 15
  
  open(u, file=fname, iostat =io, action='read',status='old')
  if (io /= 0) then
     print *, '"',fname,'"' ,' does NOT exist'
     stop ',cannot open file'
  end if
  ! fortran is used row-major order
  read(u,*) ( (M(i,j), j=1,dim), i=1,dim)
  close(u)
  return
end subroutine read_matrix

subroutine eig(matz, dim, rfname, ifname, evals, evecs)
  !
  !! eig is an wrappe of subroutine of cg in (quad)-Eispack
  ! 
  !   Discussion:
  !		This subroutine calls cg in (quad)-EISPAC subrotines.
  !		To fine the eigenvalues and eigenvectors (if descrired)
  !		of a complex general matrix
  !
  !  Licensing:
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !    5 May 2013
  !
  !  Parameters:
  !    Input, integer ( kind = 4 ) matz, is 0 if only eigenvalues are desired, and
  !    nonzero if both eigenvalues and eigenvectors are to be computed.
  !
  !    Input, integer ( kind = 4 ) dim, the order of the matrix.
  !
  !    Input, character (*) rfname/ifname, the real/imagnary parts of complex matrix data file name
  !
  !    Output, complex ( kind(0Q0) ) evals(dim), the eigenvalues.
  !
  !    Output, complex ( kind(0Q0) ) evecs(dim,dim), the eigenvectors,
  !	 if matz is not zero.
  
  implicit none
  integer ( kind = 4 ), intent(in) :: matz
  integer ( kind = 4 ), intent(in) :: dim
  character (*),intent(in) :: rfname, ifname
  complex ( kind(0Q0) ), dimension(dim), intent(inout) :: evals
  complex ( kind(0Q0) ), dimension(dim,dim),intent(inout) :: evecs
  ! todo evecs(optional), raise Error
  ! Error: Dummy argument 'evecs' of procedure 'eig' at (1) has an attribute that requires an explicit interface for this procedure
  
  real ( kind(0Q0) ), allocatable, dimension(:,:) :: M_real, M_imag
  real ( kind(0Q0) ), dimension(dim,dim) :: xi,xr
  real ( kind(0Q0) ), dimension(dim) :: wi, wr
  integer ( kind = 4 ) i, ierr, j
  
  allocate( M_real(dim,dim) )
  allocate( M_imag(dim,dim) )
  
  call read_matrix(rfname, dim, M_real)
  call read_matrix(ifname, dim, M_imag)
  
  ! get eigen values (and eigen vectors)
  call cg ( dim, M_real, M_imag, wr, wi, matz, xr, xi, ierr )
  
  deallocate( M_real )
  deallocate( M_imag )	
  
  do i=1,dim
     evals(i)=cmplx(wr(i), wi(i), kind(0Q0))
  end do
  
  if ( ierr /= 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'eispack.cg - Warning!'
     write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
     return
  end if
  
  if ( matz /= 0 ) then
     do i = 1, dim
        do j = 1, dim
           ! fortran is used row-major order
           evecs(i,j) = cmplx(xr(j,i), xi(j,i), kind(0Q0))
        end do
     end do
  end if
  
  return
end subroutine eig

subroutine save_values(dim, evals, fname)
  implicit none
  integer, intent(in) :: dim
  complex ( kind(0Q0) ), intent(in),dimension(dim) :: evals
  character (*), intent(in) :: fname
  
  
  integer i, n, io
  n=15
  open(n,file=fname,action='write',status='replace', iostat=io)
  !	write(n, *) '#index, real, imag'
  do i=1,dim
     write (n, *) real(evals(i)),",", imag(evals(i)), ","
  end do
  close(n)
  return 
end subroutine save_values

subroutine file2call_eig(dim,length,rfname,ifname)
  implicit none
  integer, intent(in) :: dim, length
  character ( len = length ),intent(in) :: rfname, ifname
  
  complex ( kind(0Q0) ),allocatable, dimension(:) :: evals
  complex ( kind(0Q0) ),allocatable, dimension(:,:) :: evecs
  
  character (len=32) :: fname
  integer i, matz
  
  allocate ( evals(dim) )
  allocate ( evecs(dim, dim) )
  
  matz = 1
  call eig(matz,dim,rfname,ifname,evals,evecs)
  call save_values(dim,evals,'eigen_vals.dat')
  
  do i=1,dim
     write (fname, '(a,i0,a)'), 'eigen_vec_', i-1, '.dat'
     call save_values(dim,evecs(i,:), fname)
  end do
  
  deallocate(evals)
  deallocate(evecs)
end subroutine file2call_eig
