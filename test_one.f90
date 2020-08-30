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
integer(int32) :: ncid_topo, depth_id,depth_all_id,frac_id          ! NetCDF ids
integer(int32) :: depth_med_id,depth_all_med_id          ! NetCDF ids
integer(int32) :: xid,yid,hid          ! NetCDF ids
integer(int32) :: dids(2)           ! NetCDF ids
integer(int32) :: dids_topo(2)           ! NetCDF ids

real(real64),allocatable,dimension(:,:)   :: wrk,wrk_super
real(real64),allocatable,dimension(:,:)   :: x_c,y_c,x_t,y_t,area
real(real64),allocatable,dimension(:)     ::xtopo,ytopo,weight, x_rot
real(real32),allocatable,dimension(:,:)   :: topo_in, topo_out,topo_all_out,frac
real(real32),allocatable,dimension(:,:)   :: topo_med_out,topo_all_med_out
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

real(real64) :: maxx,maxy,minx,miny
real(real64) :: xstart,ystart
real(real64) :: xend,yend

real(real64) ::xt_start,xt_delta
real(real64) ::yt_start,yt_delta

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
!topo_file='gebco_08_2d_rot.nc'
topo_file='gebco_2014_rot.nc'
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

yt_delta=(ytopo(ylen)-ytopo(1))/(ylen-1)
xt_delta=(xtopo(xlen)-xtopo(1))/(xlen-1)

jstart=nint((y_c(1,1)-ytopo(1))/yt_delta)+1
print *, 'Mosaic grid starts at ', y_c(1,1), ' topography jstart = ',jstart,' lat = ',ytopo(jstart)

!
call handle_error(nf90_create('topog_new.nc',ior(nf90_netcdf4,nf90_clobber),ncid_topo))
call handle_error(nf90_def_dim(ncid_topo,'nx',nxt,dids_topo(1)))
call handle_error(nf90_def_dim(ncid_topo,'ny',nyt,dids_topo(2)))
call handle_error(nf90_def_var(ncid_topo,'depth',nf90_float,dids_topo,depth_id))
call handle_error(nf90_def_var(ncid_topo,'depth_all',nf90_float,dids_topo,depth_all_id))
call handle_error(nf90_def_var(ncid_topo,'depth_med',nf90_float,dids_topo,depth_med_id))
call handle_error(nf90_def_var(ncid_topo,'depth_all_med',nf90_float,dids_topo,depth_all_med_id))
call handle_error(nf90_def_var(ncid_topo,'frac',nf90_float,dids_topo,frac_id))
call handle_error(nf90_enddef(ncid_topo))



! Do 
do jmosaic=1,j_lat-1,jblock

   jm_end=min(jmosaic+jblock-1,j_lat-1)
   ystart=y_c(1,jmosaic)
   yend=y_c(1,jm_end+1)
   
   jstart=nint((y_c(1,jmosaic)-ytopo(1))/yt_delta)+1
   jend=min(nint((y_c(1,jm_end+1)-ytopo(1))/yt_delta)+1,ylen)
   if(y_c(1,jmosaic)>ytopo(jstart)) jstart=jstart+1
   if(y_c(1,jm_end+1)<ytopo(jend)) jend=jend-1

   jcount=jend-jstart+1
   jpoints=jm_end-jmosaic+1
      
write(*,*) 'jmosaic=',jmosaic,ystart,yend,jstart,jend

   istart=1
   do imosaic=1,nxt,iblock
      im_end=min(imosaic+iblock-1,nxt)
      xstart=x_c(imosaic,1)
      xend=x_c(im_end+1,1)

      istart=nint((x_c(imosaic,1)-xtopo(1))/xt_delta)+1
      iend=min(nint((x_c(im_end+1,1)-xtopo(1))/xt_delta)+1,xlen)
      if(x_c(imosaic,1)>xtopo(istart)) istart=istart+1
      if(x_c(im_end+1,1)<xtopo(iend)) iend=iend-1
write(*,*) 'imosaic=',imosaic,xstart,xend,istart,iend

!      do i = istart+1, xlen-1
!         iend=i
!         if(xtopo(i+1) >=  xend ) exit
!      enddo

      icount=iend-istart+1
      allocate(topo_in(icount,jcount))
      allocate(itopo_in(icount,jcount))
      call handle_error(nf90_get_var(ncid,hid,itopo_in,start=[istart,jstart],count=[icount,jcount]))
      topo_in=itopo_in
      ipoints=im_end-imosaic+1
print *,ipoints

      allocate(topo_out(imosaic:im_end,jmosaic:jm_end))
      allocate(topo_all_out(imosaic:im_end,jmosaic:jm_end))
      allocate(topo_med_out(imosaic:im_end,jmosaic:jm_end))
      allocate(topo_all_med_out(imosaic:im_end,jmosaic:jm_end))
      allocate(frac(imosaic:im_end,jmosaic:jm_end))

      call make_topo_rect(topo_in,xtopo(istart:iend),ytopo(jstart:jend),            &
                     topo_out,x_c(imosaic:im_end+1,1),y_c(1,jmosaic:jm_end+1),topo_all_out,frac, &
                     topo_med_out, topo_all_med_out)
      call handle_error(nf90_put_var(ncid_topo,depth_id,topo_out,                   &
                     start=[imosaic,jmosaic],count=[ipoints,jpoints]))
      call handle_error(nf90_put_var(ncid_topo,depth_all_id,topo_all_out,                   &
                     start=[imosaic,jmosaic],count=[ipoints,jpoints]))
      call handle_error(nf90_put_var(ncid_topo,depth_med_id,topo_med_out,                   &
                     start=[imosaic,jmosaic],count=[ipoints,jpoints]))
      call handle_error(nf90_put_var(ncid_topo,depth_all_med_id,topo_all_med_out,                   &
                     start=[imosaic,jmosaic],count=[ipoints,jpoints]))
      call handle_error(nf90_put_var(ncid_topo,frac_id,frac,                   &
                     start=[imosaic,jmosaic],count=[ipoints,jpoints]))
      deallocate(topo_out,topo_in,itopo_in,topo_all_out,frac,topo_med_out, topo_all_med_out)

      istart=iend+1
   enddo
   jstart=jend+1
enddo

! Now we need to do the tripolar part

jblock=jblock/4
do jmosaic=j_lat,nyt,jblock
   jm_end=min(jmosaic+jblock-1,nyt)
   jpoints=jm_end-jmosaic+1
   do imosaic=1,nxt,iblock
      im_end=min(imosaic+iblock-1,nxt)
      ipoints=im_end-imosaic+1

      miny=min(minval(y_c(imosaic:im_end+1,jmosaic)),minval(y_c(imosaic:im_end+1,jm_end+1)), &
               minval(y_c(imosaic,jmosaic:jm_end+1)),minval(y_c(im_end+1,jmosaic:jm_end+1)))
      maxy=max(maxval(y_c(imosaic:im_end+1,jmosaic)),maxval(y_c(imosaic:im_end+1,jm_end+1)), &
               maxval(y_c(imosaic,jmosaic:jm_end+1)),maxval(y_c(im_end+1,jmosaic:jm_end+1)))
      minx=min(minval(x_c(imosaic:im_end+1,jmosaic)),minval(x_c(imosaic:im_end+1,jm_end+1)), &
               minval(x_c(imosaic,jmosaic:jm_end+1)),minval(x_c(im_end+1,jmosaic:jm_end+1)))
      maxx=max(maxval(x_c(imosaic:im_end+1,jmosaic)),maxval(x_c(imosaic:im_end+1,jm_end+1)), &
               maxval(x_c(imosaic,jmosaic:jm_end+1)),maxval(x_c(im_end+1,jmosaic:jm_end+1)))

      call get_range(xtopo,minx,maxx,istart,iend)
      call get_range(ytopo,miny,maxy,jstart,jend)
      icount=iend-istart+1
      jcount=jend-jstart+1
      allocate(topo_in(icount,jcount))
      allocate(itopo_in(icount,jcount))
      allocate(topo_out(imosaic:im_end,jmosaic:jm_end))
      allocate(topo_all_out(imosaic:im_end,jmosaic:jm_end))
      allocate(topo_med_out(imosaic:im_end,jmosaic:jm_end))
      allocate(topo_all_med_out(imosaic:im_end,jmosaic:jm_end))
      allocate(frac(imosaic:im_end,jmosaic:jm_end))
      call handle_error(nf90_get_var(ncid,hid,itopo_in,start=[istart,jstart],count=[icount,jcount]))
      topo_in=itopo_in
      call make_topo_gen(topo_in,xtopo(istart:iend),ytopo(jstart:jend),            &
                     topo_out,x_c(imosaic:im_end+1,jmosaic:jm_end+1),y_c(imosaic:im_end+1,jmosaic:jm_end+1), &
                     topo_all_out,frac,topo_med_out,topo_all_med_out)
      call handle_error(nf90_put_var(ncid_topo,depth_id,topo_out,                   &
                     start=[imosaic,jmosaic],count=[ipoints,jpoints]))
      call handle_error(nf90_put_var(ncid_topo,depth_all_id,topo_all_out,                   &
                     start=[imosaic,jmosaic],count=[ipoints,jpoints]))
      call handle_error(nf90_put_var(ncid_topo,depth_med_id,topo_med_out,                   &
                     start=[imosaic,jmosaic],count=[ipoints,jpoints]))
      call handle_error(nf90_put_var(ncid_topo,depth_all_med_id,topo_all_med_out,                   &
                     start=[imosaic,jmosaic],count=[ipoints,jpoints]))
      call handle_error(nf90_put_var(ncid_topo,frac_id,frac,                   &
                     start=[imosaic,jmosaic],count=[ipoints,jpoints]))
      deallocate(topo_out,topo_in,itopo_in,topo_all_out,frac,topo_med_out, topo_all_med_out)
   enddo
enddo




call handle_error(nf90_close(ncid_topo))
call handle_error(nf90_close(ncid))
      

contains

subroutine get_range1(vals,lower,upper,index_lo,index_hi)
! Get all values inside lower to upper

real(real64), dimension(:), intent(in) ::vals
real(real64), intent(in)               ::lower,upper
integer(int32), intent(out)             ::index_lo,index_hi

integer(int32) :: itmp,imx,imn,its,itsmax=20

imn=1
imx=size(vals)
its=1
if(vals(1)>=lower) then
   index_lo=1
else

   do
     its=its+1
     index_lo=(imx+imn)/2
     if(vals(index_lo) > lower) then
       imx=index_lo
     else
       imn=index_lo
     endif
     if(imx-imn<2) then
       index_lo=imx
       exit
     endif
     if(its>itsmax) then
        write(*,*) imn,imx,lower,vals(imn),vals(imx)
        stop 3
     endif
   enddo
endif
  
imn=1
imx=size(vals)
index_hi=imx/2
its=1
if(vals(imx)<=upper) then
   index_hi=imx
else
   do
     its=its+1
     index_hi=(imx+imn+1)/2
     if(vals(index_hi) > upper) then
       imx=index_hi
     else
       imn=index_hi
     endif
     if(imx-imn<2) then
       index_hi=imn
       exit
     endif
     if(its>itsmax) then
        write(*,*) imn,imx,upper,vals(imn),vals(imx)
        stop 4
     endif
   enddo
endif


end subroutine get_range1


subroutine get_range(vals,lower,upper,index_lo,index_hi)
! Get all values inside lower to upper

real(real64), dimension(:), intent(in) ::vals
real(real64), intent(in)               ::lower,upper
integer(int32), intent(out)             ::index_lo,index_hi

integer(int32) :: itmp,imx,imn,its,itsmax=20

imn=1
imx=size(vals)
its=1
if(vals(1)>=lower) then
   index_lo=1
else

   do
     its=its+1
     index_lo=(imx+imn)/2
     if(vals(index_lo) > lower) then
       imx=index_lo
     else
       imn=index_lo
     endif
     if(imx-imn<2) then
       index_lo=imx
       exit
     endif
     if(its>itsmax) then
        write(*,*) imn,imx,lower,vals(imn),vals(imx)
        stop 1
     endif
   enddo
endif
  
imn=1
imx=size(vals)
index_hi=imx/2
its=1
if(vals(imx)<=upper) then
   index_hi=imx
else
   do
     its=its+1
     index_hi=(imx+imn+1)/2
     if(vals(index_hi) > upper) then
       imx=index_hi
     else
       imn=index_hi
     endif
     if(imx-imn<2) then
       index_hi=imn
       exit
     endif
     if(its>itsmax) then
        write(*,*) imn,imx,upper,vals(imn),vals(imx)
        stop 2
     endif
   enddo
endif


end subroutine get_range


subroutine make_topo_gen(topo_in,x_in,y_in,topo_out,x_out,y_out,topo_all_out,frac,topo_med_out,topo_all_med_out)

use kdtree2_precision_module
use kdtree2_module

! Make topography for a general patch. 
! We know that y_out(i,:) is monotonic increasing
! We know that x_out(:,j) is monotonic increasing

real(real32),dimension(:,:), intent(in) :: topo_in
real(real64),dimension(:), intent(in) :: x_in
real(real64),dimension(:), intent(in) :: y_in
real(real32),dimension(:,:), intent(out) :: topo_out,topo_all_out,frac
real(real32),dimension(:,:), intent(out) :: topo_med_out,topo_all_med_out
real(real64),dimension(:,:), intent(in) :: x_out
real(real64),dimension(:,:), intent(in) :: y_out

logical,dimension(:,:),allocatable           :: mask

integer :: im,it,jm,jt,inext,jnext,itopo,jtopo

! tree stuff

real(kdkind), allocatable,dimension(:,:) :: possie
real(kdkind), allocatable,dimension(:)   :: csx,csy,six,siy
real(kdkind)                             :: cx(2,2),cy(2,2),cz(2,2)
real(kdkind)                             :: t_source(3),xt,yt, rad2
type(kdtree2),pointer    :: tree
type(kdtree2_result),allocatable,dimension(:) :: results
real(int32), allocatable,dimension(:) :: idx,jdx
real(real64), parameter :: DEG2RAD= asin(1.0_real64)/90.0_real64  ! PI/180
integer(int32)           :: num_found, n,ngd,num_max
real(real32),dimension(:),allocatable    :: t_s,t_s_all
integer(int32)           :: frst,lst


! bounding boxes etc.

real(real64) :: ys,yn
integer(int32) :: istart,iend, imask
real(real64) :: xc(2,2),yc(2,2)

real(real64) :: bndsx(4), bndsy(4)
integer(int32) :: wet_in_poly,in_poly,ib,jb,i,j

integer :: ifilt

im=size(topo_out,dim=1)
jm=size(topo_out,dim=2)
it=size(x_in)
jt=size(y_in)
topo_out=0.0
topo_all_out=0.0
topo_med_out=0.0
topo_all_med_out=0.0
frac=0.0

write(*,*) 'im,jm,it,jt',im,jm,it,jt
allocate(mask(it,jt))


! Cull points south of patch

mask=.true.

! Filter > 89

ifilt=1
do j=1,jt
   if(y_in(j) <89) cycle
   ifilt=2
   if ( y_in(j) > 89.5) ifilt=3
   if ( y_in(j) > 89.7) ifilt=5
   if ( y_in(j) > 89.8) ifilt=11
   if ( y_in(j) > 89.9) ifilt=17
   do i=1,it
      mask(i,j) = (mod(i,ifilt) == 0 ) .and. ( topo_in(i,j) <=0.0 )
   enddo
enddo

do i=1,im
   call get_range(x_in,x_out(i,1),x_out(i+1,1),istart,iend)
   ys=min(y_out(i,1),y_out(i+1,1))
   yn=max(y_out(i,1),y_out(i+1,1))
   do j= 1, jt
      if ( y_in(j) >= yn ) exit
      do imask=istart,iend
         if(y_in(j)<ys) then
            mask(imask,j)=.false.
         endif
      enddo
   enddo
enddo
 
! Cull points north of patch

!do i=1,im
!   call get_range(x_in,x_out(i,1),x_out(i+1,1),istart,iend)
!   ys=min(y_out(i,1),y_out(i+1,1))
!   yn=max(y_out(i,1),y_out(i+1,1))
!   do j= jt,1,-1
!      if ( y_in(j) <= ys ) exit
!      do imask=istart,iend
!         if(y_in(j)>yn) then
!            mask(imask,j)=.false.
!         endif
!      enddo
!   enddo
!enddo
!!write(*,*) 'Culled points N', count(.not.mask)
 
ngd=count(mask)

! Now pack points.

if(ngd > 0 ) then
   allocate(possie(3,ngd))
   allocate(idx(ngd),jdx(ngd))
   allocate(csy(jt),siy(jt))
   allocate(csx(it),six(it))
   do j=1,jt
      csy(j)=cos(y_in(j)*DEG2RAD)
      siy(j)=sin(y_in(j)*DEG2RAD)
   enddo
   do i=1,it
      csx(i)=cos(x_in(i)*DEG2RAD)
      six(i)=sin(x_in(i)*DEG2RAD)
   enddo
   imask=0
   do j=1,jt
      do i=1,it
         if(mask(i,j)) then
            imask=imask+1
            possie(1,imask)=csx(i)*csy(j)
            possie(2,imask)=six(i)*csy(j)
            possie(3,imask)=siy(j)
            idx(imask)=i
            jdx(imask)=j
         endif
      enddo
   enddo
   deallocate(mask)
      
! Create tree
   tree => kdtree2_create(possie,sort=.false.,rearrange=.true.)

   num_max=20000000
   allocate(t_s(num_max),t_s_all(num_max))
   
   allocate(results(num_max))
   do j=1,jm
      do i=1,im
   

         xc=x_out(i:i+1,j:j+1)*DEG2RAD
         yc=y_out(i:i+1,j:j+1)*DEG2RAD
         xt=sum(xc)/4.0
         yt=sum(yc)/4.0
         t_source(1)=cos(xt)*cos(yt)
         t_source(2)=sin(xt)*cos(yt)
         t_source(3)=sin(yt)
         cx=cos(xc)*cos(yc)
         cy=sin(xc)*cos(yc)
         cz=sin(yc)

         rad2=maxval((cx-t_source(1))**2+(cy-t_source(2))**2+(cz-t_source(3))**2)

         call kdtree2_r_nearest(tp=tree,qv=t_source,r2=rad2,nfound=num_found, &
                                             nalloc=num_max,results=results)

         if(num_found==0) then
!            write(*,*) 'Nothing found, i,j,x_t(i,j),y_t(i,j)',i,j,x_t(i,j),y_t(i,j)
!            write(*,*) x_c(i:i+1,j:j+1)
!            write(*,*) y_c(i:i+1,j:j+1)
!            write(*,*) rad2,minval(x_in),maxval(x_in),minval(y_in),maxval(y_in)
!            stop
            cycle
         endif

         wet_in_poly=0
         in_poly=0
         t_s=0
         t_s_all=0

         bndsx=[ x_out(i,j),x_out(i+1,j),x_out(i+1,j+1),x_out(i,j+1) ]
         bndsy=[ y_out(i,j),y_out(i+1,j),y_out(i+1,j+1),y_out(i,j+1) ]
         do n=1,min(num_found,num_max)
            ib=idx(results(n)%idx)
            jb=jdx(results(n)%idx)
!            if(topo_in(ib,jb) > 0.0) cycle
            if( pnt_in_quad( bndsx, bndsy, x_in(ib), y_in(jb))) then
                in_poly=in_poly+1
                topo_all_out(i,j)=topo_all_out(i,j)-topo_in(ib,jb)
                t_s_all(in_poly)=-topo_in(ib,jb)
                if(topo_in(ib,jb) < 0.0) then
                   wet_in_poly=wet_in_poly+1
                   topo_out(i,j)=topo_out(i,j)-topo_in(ib,jb)
                   t_s(wet_in_poly)=-topo_in(ib,jb)
                endif
            endif
         enddo
         if(in_poly /= 0 ) then
            topo_all_out(i,j)= topo_all_out(i,j)/in_poly
            frst=1;lst=in_poly
            call quicksort(t_s,frst,lst)
            topo_all_med_out(i,j)=t_s_all(max(in_poly/2,1))
         endif
         if(wet_in_poly /= 0 ) then
            topo_out(i,j)= topo_out(i,j)/wet_in_poly
            frac(i,j)=real(wet_in_poly)/in_poly
            frst=1;lst=wet_in_poly
            call quicksort(t_s,frst,lst)
            topo_med_out(i,j)=t_s(max(wet_in_poly/2,1))
         endif
      enddo
   enddo

   call kdtree2_destroy(tree)
endif
end subroutine make_topo_gen
             

   


      




subroutine make_topo_rect(topo_in,x_in,y_in,topo_out,x_out,y_out,topo_all_out,frac,topo_med_out,topo_all_med_out)
type med_type
     real(real32),dimension(:),allocatable :: topo
end type med_type
real(real32),dimension(:,:), intent(in) :: topo_in
real(real64),dimension(:), intent(in) :: x_in
real(real64),dimension(:), intent(in) :: y_in
real(real32),dimension(:,:), intent(out) :: topo_out,topo_all_out,frac
real(real32),dimension(:,:), intent(out) :: topo_med_out,topo_all_med_out
real(real64),dimension(:), intent(in) :: x_out
real(real64),dimension(:), intent(in) :: y_out

type(med_type),allocatable,dimension(:) :: t_s,t_s_all

integer,dimension(:,:), allocatable :: npts,npts_all

real(real32),dimension(:),allocatable    :: dummy
integer(int32)           :: frst,lst
integer :: im,it,jm,jt,inext,jnext,itopo,jtopo,n
allocate(npts(size(topo_out,dim=1),size(topo_out,dim=2)))
allocate(npts_all(size(topo_out,dim=1),size(topo_out,dim=2)))
allocate(t_s(size(topo_out,dim=1)))
allocate(t_s_all(size(topo_out,dim=1)))
npts=0
npts_all=0
topo_out=0.0
topo_all_out=0.0
topo_med_out=0.0
topo_all_med_out=0.0
frac=0.0
jt=1
jnext=1
do jm=1,size(y_out)-1
   do im=1,size(x_out)-1
      if(.not.allocated(t_s(im)%topo)) allocate(t_s(im)%topo(500))
      if(.not.allocated(t_s_all(im)%topo)) allocate(t_s_all(im)%topo(500))
   enddo
   do jtopo=jt,size(y_in)
      if(y_in(jtopo) >= y_out(jm+1)) exit
      jnext=jtopo+1
      it=1
      inext=1
      do im=1,size(x_out)-1
         do itopo=it,size(x_in)
            if(x_in(itopo) >= x_out(im+1)) exit
            inext=itopo+1
            topo_all_out(im,jm)=topo_all_out(im,jm)-topo_in(itopo,jtopo)
            npts_all(im,jm)= npts_all(im,jm)+1
            n=npts_all(im,jm)
            if(n > size(t_s_all(im)%topo) ) then
               allocate(dummy(n-1))
               dummy=t_s_all(im)%topo
               deallocate(t_s_all(im)%topo)
               allocate(t_s_all(im)%topo(2*(n-1)))
               t_s_all(im)%topo(1:n-1)=dummy
               deallocate(dummy)
            endif
            t_s_all(im)%topo(n)=-topo_in(itopo,jtopo)
            if(topo_in(itopo,jtopo) < 0.0 ) then
               topo_out(im,jm)=topo_out(im,jm)-topo_in(itopo,jtopo)
               npts(im,jm)= npts(im,jm)+1
               n=npts(im,jm)
               if(n > size(t_s(im)%topo) ) then
                  allocate(dummy(n-1))
                  dummy=t_s(im)%topo
                  deallocate(t_s(im)%topo)
                  allocate(t_s(im)%topo(2*(n-1)))
                  t_s(im)%topo(1:n-1)=dummy
                  deallocate(dummy)
               endif
               t_s(im)%topo(n)=-topo_in(itopo,jtopo)
            endif
         enddo
         it=inext
      enddo
   enddo
   jt=jnext
   do im=1,size(x_out)-1
      if(npts(im,jm)>0) then
         frst=1
         lst=npts(im,jm)
         call quicksort(t_s(im)%topo,frst,lst)
         topo_med_out(im,jm)=t_s(im)%topo((npts(im,jm)+1)/2)
      endif
      if(npts_all(im,jm)>0) then
         frst=1
         lst=npts_all(im,jm)
         call quicksort(t_s_all(im)%topo,frst,lst)
         topo_all_med_out(im,jm)=t_s_all(im)%topo((npts(im,jm)+1)/2)
      endif
   enddo
enddo
where(npts > 0) topo_out=topo_out/npts
where(npts_all > 0) topo_all_out=topo_all_out/npts_all
where(npts_all > 0) frac=real(npts)/npts_all
deallocate(npts,npts_all)
end subroutine make_topo_rect

logical function pnt_in_quad(x,y,PX,PY)

IMPLICIT NONE

real(real64), intent(in) :: x(4), y(4)
real(real64), intent(in) :: px, py
INTEGER I , J
real(real64) :: XI , YI , XJ , YJ
LOGICAL IX , IY , JX , JY

pnt_in_quad=.false.

do I = 1 , 4
  XI = X(I) - PX
  YI = Y(I) - PY
!        CHECK WHETHER THE POINT IN QUESTION IS AT THIS VERTEX.
  IF ( XI.EQ.0.0 .AND. YI.EQ.0.0 ) THEN
     pnt_in_quad=.true.
     RETURN
  ENDIF
!        J IS NEXT VERTEX NUMBER OF POLYGON.
  J = 1 + MOD(I,4)
  XJ = X(J) - PX
  YJ = Y(J) - PY
!  IS THIS LINE OF 0 LENGTH ?
  IF ( XI.EQ.XJ .AND. YI.EQ.YJ ) cycle
  IX = XI.GE.0.0
  IY = YI.GE.0.0
  JX = XJ.GE.0.0
  JY = YJ.GE.0.0
!  CHECK WHETHER (PX,PY) IS ON VERTICAL SIDE OF POLYGON.
  IF ( XI.EQ.0.0 .AND. XJ.EQ.0.0 .AND. (IY.neqv.JY) ) THEN
     pnt_in_quad=.true.
     RETURN
  ENDIF
!  CHECK WHETHER (PX,PY) IS ON HORIZONTAL SIDE OF POLYGON.
  IF ( YI.EQ.0.0 .AND. YJ.EQ.0.0 .AND. (IX.neqv.JX) ) THEN
     pnt_in_quad=.true.
     RETURN
  ENDIF
!       CHECK WHETHER BOTH ENDS OF THIS SIDE ARE COMPLETELY 1) TO RIGHT
!       OF, 2) TO LEFT OF, OR 3) BELOW (PX,PY).
   IF ( .NOT.((IY .OR. JY) .AND. (IX.neqv.JX)) ) cycle
!        DOES THIS SIDE OBVIOUSLY CROSS LINE RISING VERTICALLY FROM (PX,PY)
   IF ( .NOT.(IY .AND. JY .AND. (IX.neqv.JX)) ) THEN
      IF ( (YI*XJ-XI*YJ)/(XJ-XI).LT.0.0 ) THEN
         cycle
      ELSEIF ( (YI*XJ-XI*YJ)/(XJ-XI).EQ.0.0 ) THEN
         pnt_in_quad=.true.
         RETURN
      ELSE
         pnt_in_quad=.not.pnt_in_quad
      ENDIF
   ELSE
      pnt_in_quad=.not.pnt_in_quad
   ENDIF

ENDDO
end function pnt_in_quad



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
    
recursive subroutine quicksort(a,first,last)
  real(real32)  a(:), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort
end program load_mosaic
