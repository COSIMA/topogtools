program load_mosaic
!
! Create a file with all the wet points that we wish to map to.
!
! Mosaic input.
!
! Usage: load_mosaic
!
! Assumes there's a file 'mosaic.nc' lurking about....
!
! Output. A file called ...
!
!
use iso_fortran_env
use netcdf
implicit none

integer(int32) :: i,j,k
integer(int32) :: nx,ny                  ! Size of model grid
integer(int32) :: nxp,nyp                ! Size of model supergrid
integer(int32) :: nxt,nyt                ! Size of model T grid
integer(int32) :: nxc,nyc                ! Size of model corner grid

integer(int32) :: ncid,vid,did           ! NetCDF ids
integer(int32) :: ncid_topo, depth_id           ! NetCDF ids
integer(int32) :: xid,yid,hid          ! NetCDF ids
integer(int32) :: dids(2)           ! NetCDF ids
integer(int32) :: dids_topo(2)           ! NetCDF ids

real(real64),allocatable,dimension(:,:)   :: wrk,wrk_super
real(real64),allocatable,dimension(:,:)   :: x_c,y_c,x_t,y_t,area
real(real64),allocatable,dimension(:)     ::xtopo,ytopo,weight, x_rot
real(real32),allocatable,dimension(:,:)   :: topo_in, topo_out
integer(int16),allocatable,dimension(:,:)   :: itopo_in

character(len=128)  :: odir='',ofile=''  ! for ocean_mosaic.nc
character(len=128)  :: gdir='',gfile=''  ! for hgrid file
character(len=256)  :: dirfile=''        ! concatenation
character(len=256)  :: topo_file=''      ! 
character(len=16)   :: xname,yname

logical             :: fexist = .false., istripolar = .true.

integer(int32)      :: x_cyc, j_lat
real(real64)        :: x_shift

integer(int32)      :: xlen, ylen

integer(int32) :: istart,jstart
integer(int32) :: iend,jend
integer(int32) :: icount,jcount
integer(int32) :: iblock=100,jblock=100
integer(int32) :: imosaic,jmosaic
integer(int32) :: im_end,jm_end
integer(int32) :: ipoints,jpoints

real(real64) :: xstart,ystart
real(real64) :: xend,yend

! Get info on the grid from input

write(*,*) 'Getting model grid info'
! Get mosaic info
inquire(file=trim('mosaic.nc'),exist=fexist)
if ( .not. fexist ) then
   write(*,*) 'mosaic.nc does not exist. Bailing out' 
   stop 1
endif
call handle_error(nf90_open('mosaic.nc',nf90_nowrite,ncid))
call handle_error(nf90_inq_varid(ncid,'ocn_mosaic_dir',vid))
call handle_error(nf90_get_var(ncid,vid,odir))
call handle_error(nf90_inq_varid(ncid,'ocn_mosaic_file',vid))
call handle_error(nf90_get_var(ncid,vid,ofile))
call handle_error(nf90_close(ncid))
! Get horizontal grid
dirfile=odir(1:scan(odir,'/',back=.true.)) // ofile(1:scan(ofile,'c',back=.true.))
write(*,*) len_trim(dirfile),dirfile
inquire(file=trim(dirfile),exist=fexist)
if ( .not. fexist ) then
   write(*,*) 'ocn_mosaic_dir/ocn_mosaic_file =',trim(dirfile), ' does not exist. Bailing out' 
   stop 1
endif
call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
call handle_error(nf90_inq_varid(ncid,'gridlocation',vid))
call handle_error(nf90_get_var(ncid,vid,gdir))
call handle_error(nf90_inq_varid(ncid,'gridfiles',vid))
call handle_error(nf90_get_var(ncid,vid,gfile))
call handle_error(nf90_close(ncid))

!
! On mosaic "supergrid" we need to get every second point
!
write(*,*) 'Reading supergrid info'
! Read xt
dirfile=gdir(1:scan(gdir,'/',back=.true.)) // gfile(1:scan(gfile,'c',back=.true.))
inquire(file=trim(dirfile),exist=fexist)
if ( .not. fexist ) then
   write(*,*) 'gridlocation/gridfiles =',trim(dirfile), ' does not exist. Bailing out' 
   stop 1
endif
call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
call handle_error(nf90_inq_dimid(ncid,'nx',did))
call handle_error(nf90_inquire_dimension(ncid,did,len=nx))
nxp=nx+1
nxt=nx/2
nxc=nx/2+1
call handle_error(nf90_inq_dimid(ncid,'ny',did))
call handle_error(nf90_inquire_dimension(ncid,did,len=ny))
nyp=ny+1
nyt=ny/2
nyc=ny/2+1

allocate(x_c(nxc,nyc),y_c(nxc,nyc))
allocate(x_t(nxt,nyt),y_t(nxt,nyt))
allocate(wrk_super(nxp,nyp))

! Get x corners
call handle_error(nf90_inq_varid(ncid,'x',vid))
call handle_error(nf90_get_var(ncid,vid,wrk_super))
do j=1,nyc
   do i = 1,nxc
      x_c(i,j)= wrk_super(2*i-1,2*j-1)
   enddo
enddo
do j=1,nyt
   do i = 1,nxt
      x_t(i,j)= wrk_super(2*i,2*j)
   enddo
enddo

! Get y corners

call handle_error(nf90_inq_varid(ncid,'y',vid))
call handle_error(nf90_get_var(ncid,vid,wrk_super))
do j=1,nyc
   do i = 1,nxc
      y_c(i,j)= wrk_super(2*i-1,2*j-1)
   enddo
enddo
do j=1,nyt
   do i = 1,nxt
      y_t(i,j)= wrk_super(2*i,2*j)
   enddo
enddo

j_lat = nyc
if ( istripolar ) then
   j_lat =0
   do j = 1,nyc
      if( y_c(1,j) /= y_c(2,j) ) then
         j_lat = j-1
         exit
      endif
   enddo

   if ( j_lat == 0 ) then
      write(*,*) 'FATAL: unable to locate j_lat for tripolar grid'
      stop 1
   endif
endif

write(*,*) ' j_lat located at ', j_lat, 'latitude ', y_c(1,j_lat)

! Area, probably should do dxt*dyt correctly but I think this is ok.
deallocate(wrk_super)
!allocate(area(nxt,nyt))
!allocate(wrk_super(nx,ny))
!call handle_error(nf90_inq_varid(ncid,'area',vid))
!call handle_error(nf90_get_var(ncid,vid,wrk_super))
!call handle_error(nf90_close(ncid))
!write(*,*) nx,ny,shape(wrk_super),shape(area)
!do j = 1, nyt
!   do i = 1,nxt
!      area(i,j) = wrk_super(2*i-1,2*j-1)+wrk_super(2*i,2*j-1)+wrk_super(2*i-1,2*j)+wrk_super(2*i,2*j)
!   enddo
!enddo
!
!deallocate(wrk_super)

! Now load up spherical grid

!topo_file='/home/datalib/bathymetry/GEBCO_2008/gebco_08_2d.nc'
topo_file='gebco_08_2d_rot.nc'
call handle_error(nf90_open(trim(topo_file),nf90_nowrite,ncid))
call handle_error(nf90_inq_varid(ncid,'height',hid))
call handle_error(nf90_inquire_variable(ncid,hid,dimids=dids))
call handle_error(nf90_inquire_dimension(ncid,dids(1),xname,xlen))
call handle_error(nf90_inquire_dimension(ncid,dids(2),yname,ylen))
call handle_error(nf90_inq_varid(ncid,trim(xname),xid))
call handle_error(nf90_inq_varid(ncid,trim(yname),yid))

allocate(xtopo(xlen),ytopo(ylen),weight(ylen),x_rot(xlen))

call handle_error(nf90_get_var(ncid,xid,xtopo))
call handle_error(nf90_get_var(ncid,yid,ytopo))

! work in patches

! Get to southern edge

jstart=0
do j = 1,ylen
   jstart=jstart+1
   if(ytopo(jstart) >= y_c(1,1)) exit
enddo

print *, 'Mosaic grid starts at ', y_c(1,1), ' topography jstart = ',jstart,' lat = ',ytopo(jstart)

!
call handle_error(nf90_create('topog_new.nc',ior(nf90_netcdf4,nf90_clobber),ncid_topo))
call handle_error(nf90_def_dim(ncid_topo,'xx',nxt,dids_topo(1)))
call handle_error(nf90_def_dim(ncid_topo,'yy',nyt,dids_topo(2)))
call handle_error(nf90_def_var(ncid_topo,'depth',nf90_float,dids_topo,depth_id))
call handle_error(nf90_enddef(ncid_topo))



! Do 
do jmosaic=1,j_lat-1,jblock
write(*,*) 'jmosaic=',jmosaic

   jm_end=min(jmosaic+jblock-1,j_lat-1)
   ystart=y_c(1,jmosaic)
   yend=y_c(1,jm_end+1)
   
   do j=jstart+1,ylen-1
      jend=j
      if(ytopo(j+1) >= yend) exit
   enddo

   jcount=jend-jstart+1
   jpoints=jm_end-jmosaic+1
      

write(*,*) 'jstart,jend,jcount,ytopo(jstart),ytopo(jend)',jstart,jend,jcount,ytopo(jstart),ytopo(jend)
write(*,*) 'ystart,yend',ystart,yend
if(ytopo(jstart) < y_c(1,jmosaic)) then
  write(*,*) 'Err ytopo(jstart) < y_c(1,jmosaic)',jstart,jmosaic,ytopo(jstart),y_c(1,jmosaic)
  stop
endif
if(ytopo(jend) > y_c(1,jm_end+1)) then
  write(*,*) 'Err ytopo(jend) > y_c(1,jm_end+1)',jend,jmosaic,ytopo(jend),y_c(1,jm_end+1)
  stop
endif
   istart=1
   do imosaic=1,nxt,iblock
      im_end=min(imosaic+iblock-1,nxt)
      xstart=x_c(imosaic,1)
      xend=x_c(im_end+1,1)

      do i = istart+1, xlen-1
         iend=i
         if(xtopo(i+1) >=  xend ) exit
      enddo

      icount=iend-istart+1
      allocate(topo_in(icount,jcount))
      allocate(itopo_in(icount,jcount))
      call handle_error(nf90_get_var(ncid,hid,itopo_in,start=[istart,jstart],count=[icount,jcount]))
      topo_in=itopo_in
      ipoints=im_end-imosaic+1
!write(*,*) 'imosaic=',imosaic,im_end,icount,jcount,xtopo(istart),xtopo(iend),ytopo(jstart),ytopo(jend)

      allocate(topo_out(imosaic:im_end,jmosaic:jm_end))

      call make_topo_rect(topo_in,xtopo(istart:iend),ytopo(jstart:jend),            &
                     topo_out,x_c(imosaic:im_end+1,1),y_c(1,jmosaic:jm_end+1))
      call handle_error(nf90_put_var(ncid_topo,depth_id,topo_out,                   &
                     start=[imosaic,jmosaic],count=[ipoints,jpoints]))
      deallocate(topo_out,topo_in,itopo_in)

      istart=iend+1
   enddo
   jstart=jend+1
enddo
      

contains

subroutine make_topo_rect(topo_in,x_in,y_in,topo_out,x_out,y_out)
real(real32),dimension(:,:), intent(in) :: topo_in
real(real64),dimension(:), intent(in) :: x_in
real(real64),dimension(:), intent(in) :: y_in
real(real32),dimension(:,:), intent(out) :: topo_out
real(real64),dimension(:), intent(in) :: x_out
real(real64),dimension(:), intent(in) :: y_out

integer,dimension(:,:), allocatable :: npts

integer :: im,it,jm,jt,inext,jnext,itopo,jtopo
allocate(npts(size(topo_out,dim=1),size(topo_out,dim=2)))
npts=0
topo_out=0.0
jt=1
jnext=1
do jm=1,size(y_out)-1
   do jtopo=jt,size(y_in)
      if(y_in(jtopo) >= y_out(jm+1)) exit
      jnext=jtopo+1
      it=1
      inext=1
      do im=1,size(x_out)-1
         do itopo=it,size(x_in)
            if(x_in(itopo) >= x_out(im+1)) exit
            inext=itopo+1
            if(topo_in(itopo,jtopo) < 0.0 ) then
               topo_out(im,jm)=topo_out(im,jm)-topo_in(itopo,jtopo)
               npts(im,jm)= npts(im,jm)+1
            endif
         enddo
         it=inext
      enddo
   enddo
   jt=jnext
enddo
where(npts > 0) topo_out=topo_out/npts
deallocate(npts)
end subroutine make_topo_rect




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
    
end program load_mosaic
