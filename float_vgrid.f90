program float_vgrid
! Zeta is double precision. Convert to  single and rewrite. This stops small floating point errors.
! Write out 
use netcdf
use iso_fortran_env
implicit none
real(real32), dimension(:),allocatable :: zeta_float
real(real64), dimension(:),allocatable :: zeta_dp
integer :: ierr, nzeta
integer :: ncid,vid
integer,dimension(1) :: dids
ierr = nf90_open('ocean_vgrid.nc',nf90_write,ncid)
ierr= nf90_inq_varid(ncid,'zeta',vid)
ierr= nf90_inquire_variable(ncid,vid,dimids=dids)
ierr= nf90_inquire_dimension(ncid,dids(1),len=nzeta)
allocate(zeta_dp(nzeta),zeta_float(nzeta))
ierr= nf90_get_var(ncid,vid,zeta_dp)
zeta_float=zeta_dp
zeta_dp=zeta_float
ierr= nf90_put_var(ncid,vid,zeta_dp)
ierr= nf90_close(ncid)
end program float_vgrid
