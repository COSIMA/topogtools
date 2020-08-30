program do_partial_cells
!
! Make sure that partial cells are not too thin.
!
! Usage: do_partial_cells topo_file min_thick min_frac
!
! Ensure that partial cells are greater than max(min_think,min_frac*dzt(kmax))
!

use iso_fortran_env
use netcdf
implicit none

real(real32), dimension(:,:),allocatable :: topog
integer, dimension(:,:),allocatable :: num_levels

real(real64), dimension(:),allocatable ::zeta
real(real32), dimension(:),allocatable ::zw,dz,mdz
real(real32) :: del

integer :: ierr, i,j,k,ni,nj,nzeta,nz,its,counter

integer :: ncid,vid
integer,dimension(2) :: dids

character*32 :: arg
character*128 :: topo_file
real(real32) :: min_thick, min_frac

if (command_argument_count() /= 3) then
   write(*,*) 'ERROR: Wrong number of arguments'
   write(*,*) 'Usage: do_partial_cells topo_file min_thick min_frac'
   stop 1
endif

call get_command_argument(1,topo_file)

call get_command_argument(2,arg)
read(arg,*) min_thick
call get_command_argument(3,arg)
read(arg,*) min_frac

write(*,*) 'Ensuring cell thicknesses are greater than ',min_thick,'m and ',min_frac,'times full thickness'



call handle_error(nf90_open('ocean_vgrid.nc',nf90_nowrite,ncid))
call handle_error(nf90_inq_varid(ncid,'zeta',vid))
call handle_error(nf90_inquire_variable(ncid,vid,dimids=dids))
call handle_error(nf90_inquire_dimension(ncid,dids(1),len=nzeta))

nz=nzeta/2

allocate(zeta(nzeta),zw(0:nz),dz(nz))

call handle_error(nf90_get_var(ncid,vid,zeta))
call handle_error(nf90_close(ncid))

zw(:)=zeta(1:nzeta:2)
dz=zw(1:nz)-zw(0:nz-1)
mdz=min(dz,max(min_frac*dz,min_thick))


call handle_error( nf90_open(trim(topo_file),nf90_write,ncid))
call handle_error(nf90_inq_varid(ncid,'depth',vid))
call handle_error(nf90_inquire_variable(ncid,vid,dimids=dids))
call handle_error(nf90_inquire_dimension(ncid,dids(1),len=ni))
call handle_error(nf90_inquire_dimension(ncid,dids(2),len=nj))

allocate(topog(ni,nj))
allocate(num_levels(ni,nj))
call handle_error(nf90_get_var(ncid,vid,topog))

! Find number of levels

num_levels=0

do j=1,nj
   do i=1,ni
      if(topog(i,j)>0.0) then
         kloop: do k=2,nz
            if(zw(k) >=topog(i,j)) then
               num_levels(i,j)=k
               exit kloop
            endif
          enddo kloop
       endif
   enddo
enddo
      
do j=1,nj
   do i=1,ni
      k=num_levels(i,j)
      if(k == 0 ) cycle
      if (topog(i,j) == zw(k)) cycle

      del=topog(i,j)-zw(k-1)

!If less than half mininimum allowable thickness shave

      if(del < 0.5*mdz(k)) then
         topog(i,j)=zw(k-1)
         cycle
      endif
! Deep if necessary and make sure there is no rounding problem

      topog(i,j)=min(max(topog(i,j),zw(k-1)+mdz(k)),zw(k))
   enddo
enddo

call handle_error(nf90_redef(ncid))
call handle_error(nf90_put_att(ncid,vid,'min_thick',min_thick))
call handle_error(nf90_put_att(ncid,vid,'min_frac',min_frac))
call handle_error(nf90_enddef(ncid))
call handle_error(nf90_put_var(ncid,vid,topog))
call handle_error(nf90_close(ncid))

contains

subroutine handle_error(error_flag,isfatal,err_string)
! Simple error handle for NetCDF
integer(int32),intent(in) :: error_flag
logical, intent(in),optional :: isfatal
character(*), intent(in),optional :: err_string
logical            :: fatal
fatal = .true.
if(present(isfatal)) fatal=isfatal
if ( error_flag  /= nf90_noerr ) then
   if ( fatal ) then
      write(*,*) 'FATAL ERROR:',nf90_strerror(error_flag)
      if (present(err_string)) write(*,*) trim(err_string)
      stop
   endif
endif
end subroutine handle_error

end program do_partial_cells
