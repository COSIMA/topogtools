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

integer(int32) :: ncid,vid,did           ! NetCDF ids

real(real64),allocatable,dimension(:,:)   :: wrk,wrk_super
real(real64),allocatable,dimension(:,:)   :: x_c,y_c,area

character(len=128)  :: odir='',ofile=''  ! for ocean_mosaic.nc
character(len=128)  :: gdir='',gfile=''  ! for hgrid file
character(len=256)  :: dirfile=''        ! concatenation

logical             :: fexist = .false.

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
call handle_error(nf90_inquire_dimension(ncid,did,nx)
nxp=nx+1
nxt=nx/2
nxc=nx/2+1
call handle_error(nf90_inq_dimid(ncid,'ny',did))
call handle_error(nf90_inquire_dimension(ncid,did,ny)
nyp=ny+1
nyt=ny/2
nxy=ny/2+1

allocate(x_c(nxc,nyc),y_c(nxc,nyc))
allocate(wrk_super(nxp,nyp)

! Get x corners
call handle_error(nf90_inq_varid(ncid,'x',vid))
call handle_error(nf90_get_var(ncid,vid,wrk_super))
do j=1,nyc
   do i = 1,nxc
      x_c(i,j)= wrk_super(2*i-1,2*j-1)
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

j_lat = nyc
if ( istripolar ) then
   j_lat =0
   do j = 1,1,nyc
      if( y_c(1,j) /= y_c(2,j) ) then
         jlat = j-1
         exit
      endif
   enddo

   if ( j_lat == 0 ) then
      write(*,*) 'FATAL: unable to locate j_lat for tripolar grid'
      stop 1
   endif
endif


! Area, probably should do dxt*dyt correctly but I think this is ok.
deallocate(wrk_super)
allocate(wrk_super(nx,ny))
call handle_error(nf90_inq_varid(ncid,'area',vid))
call handle_error(nf90_get_var(ncid,vid,wrk_super))
call handle_error(nf90_close(ncid))
do j = 1, ny
   do i = 1,nx
      area(i,j) = wrk_super(2*i-1,2*j-1)+wrk_super(2*i,2*j-1)+wrk_super(2*i-1,2*j)+wrk_super(2*i,2*j)
   enddo
enddo

deallocate(wrk_super)

! Now load up spherical grid

topo_file='/home/datalib/bathymetry/GEBCO_2008/gebco_08_2d.nc'
call handle_error(nf90_open(trim(topo_file),ncid))
call handle_error(nf90_inq_varid(ncid,'height',hid))
call handle_error(nf90_inquire_variable(ncid,hid,dimids=dids))
call handle_error(nf90_inq_dimension(ncid,dids(1),xname,xlen))
call handle_error(nf90_inq_dimension(ncid,dids(2),yname,ylen))
call handle_error(nf90_inq_varid(ncid,trim(xname),xid)
call handle_error(nf90_inq_varid(ncid,trim(yname),yid)
allocate(xtopo(xlen),ytopo(ylen),weight(ylen))
call handle_error(nf90_get_var(ncid,xid,xtopo))
call handle_error(nf90_get_var(ncid,yid,ytopo))

x_cyc = 0
x_shift=x_c(1,1)
if(x_shift < 0.0) then
   x_cyc=count(x_c(:,1) < 0.0)
else
   x_cyc=count(x_c(:,1) > 360.0)
endif

x_rot=mod(cshift(x_c(:,1),x_cyc),360.0)

write(*,*) 'Cycling ',x_cyc,' spaces x_c_start = ',x_c(1+x_cyc),'xtopo(1),xtopo(1),xtopo(2)

! Weights
ys=D2R*max(-90.0,1.5*ytopo(1)-0.5*ytopo(2))
yn=D2R*0.5*(ytopo(2)-ytopo(1)
weight(1)=sin(yn)-sin(ys)
do j = 2,ylen-1
   ys=yn
   yn=D2R*0.5*(ytopo(j)-ytopo(j-1)
   weight(i)=sin(yn)-sin(ys)
enddo
ys=y2
yn=D2R*min(90.0,1.5*ytopo(ylen)-0.5*ytopo(n-1))
weight(ylen)=sin(yn)-sin(ys)

! X weights

allocate(xw(xlen),dest(xlen))

ioff=1
xw=0
do i=2,nxp
   do itop = ioff,xlen
      if(xtopo(ioff) > x_rot(i) ) then
         exit
      endif
      xw(itop)=xw(itop)+1
      dest(itop)=i-1
   enddo
enddo
      



! Do by latitude
rows=0
do j = 1, jlat-1
   call get_rows(y_c(j:j+1),weight,ytopo,rows)
   call process_row(x_rot,xtopo,rows,topo_out)
enddo
   

do js = 1,ny,jblock
   je = min(js+jblock,ny)
   do is = 1,nx,iblock
      ie = min(is+iblock,nx)
      xmin=minval(x_c(is,js:je+1))
      xmax=maxval(x_c(ie+1,js:je+1))
      ymin=minval(y_c(is:ie+1,js))
      ymax=maxval(y_c(is:ie+1,je+1))
      ist=-1
      do i=1,xlen
         if(xtopo(i) >= xmin .and. xtopo(i) < xmax ) then
            ist=i
            exit
         endif
      enddo
      iet=ist
      do i=ist,xlen
         if(xtopo(i) >= xmax ) then
            exit
         endif
         iet=i
      enddo

      jst=-1
      do j=1,ylen
         if(ytopo(j) >= ymin .and. ytopo(j) < ymax ) then
            jst=j
            exit
         endif
      enddo
      jet=jst
      do j=jst,ylen
         if(ytopo(j) >= ymax ) then
            exit
         endif
         jet=j
      enddo
         
      call load_topo(topo_in,ist,iet,jst,jet)

      call process_2d(topo_in,xtopo(ist:iet),ytopo(jst:jet),



contains

subroutine get_rows(edges,lat,rows)
real(real32),intent(in) :: edges(2)
real(real32),dimension(:),intent(in) :: lat
integer(int32),intent(inout)           :: rows(2)

integer :: j
rows=0
do j=1,size(lat)
   if(lat(j)) >= edges(1) .and. lat(j)) < edges(2))
      rows = j
      exit
   else
      if(lat(j) >= edges(2)) exit
   endif
enddo

do j=rows(1)+1,size(lat)
   if(lat(j)) < edges(2))
      rows(2) = j
   else
      exit
   endif
enddo
end subroutine get_rows

subroutine process_rows(x,cyc,rows,top,xtopo)
real,dimension(:) :: x
integer           :: cyc
integer,dimension(2) :: rows
real,dimension(:) :: top

integer :: nrows
real,dimension(xlen) :: top_in

nrows=rows(2)-rows(1)+1
topx=0.0
sumj=0.0

do j=1,nrows
   jtop=j+rows(1)-1
   call handle_error(nf90_get_var(ncid,hid,top_in,start=[1,jtop],count=[xlen,1])
   top_in=-top_in
   wgt=weight(jtop)
   mask=top_in > 0.0_real64
   where(mask)
      sumj=sumj+wgt
      topx=topx+wgt*top_in
      cnt=cnt+1
      top_max=max(top_max,top_in)
      top_min=min(top_min,top_in)
   end where
enddo

is=0
ie=1
do it=1,xlen
   i=dest(it)
   sumi(i)=sum(i)+sumj(it)
   top(i)=top(i)+topx(it)
   cont(i)=cont(i)+cnt(it)
   topmax(i)=max(topmax(i),top_max(it))
   if(cnt(it) > 0 ) then topmin(i)=min(topmin(i),top_min(it))
enddo
where(cnt > 0 ) top=top/sumi
end subroutine process_rows
     


subroutine create_wet_file(wet)
   type(wet_type), intent(in) :: wet
   integer(int32) :: ncid, did_ic
   integer(int32) :: wet_i_id,wet_j_id
   integer(int32) :: wet_x_id,wet_y_id
   integer(int32) :: wet_area_id

   call handle_error(nf90_create('model_wet.nc',ior(NF90_CLOBBER,NF90_NETCDF4),ncid))
   call handle_error(nf90_def_dim(ncid,'iw',wet%npts,did_ic))
!
! These are for ALL the model wetal points
!
   call handle_error(nf90_def_var(ncid,'wet_i',nf90_int,did_ic,wet_i_id))
   call handle_error(nf90_put_att(ncid,wet_i_id,'long_name','model i index'))
   call handle_error(nf90_def_var(ncid,'wet_j',nf90_int,did_ic,wet_j_id))
   call handle_error(nf90_put_att(ncid,wet_j_id,'long_name','model j index'))
   call handle_error(nf90_def_var(ncid,'wet_x',nf90_double,did_ic,wet_x_id))
   call handle_error(nf90_put_att(ncid,wet_x_id,'long_name','model longitude'))
   call handle_error(nf90_put_att(ncid,wet_x_id,'units','degrees_E'))
   call handle_error(nf90_def_var(ncid,'wet_y',nf90_double,did_ic,wet_y_id))
   call handle_error(nf90_put_att(ncid,wet_y_id,'long_name','model latitude'))
   call handle_error(nf90_put_att(ncid,wet_y_id,'units','degrees_N'))
   call handle_error(nf90_def_var(ncid,'wet_area',nf90_double,did_ic,wet_area_id))
   call handle_error(nf90_put_att(ncid,wet_area_id,'long_name','model area'))
   call handle_error(nf90_put_att(ncid,wet_area_id,'units','m^2'))

   call handle_error(nf90_enddef(ncid,h_minfree=4096))

! Put it there
   call handle_error(nf90_put_var(ncid,wet_i_id,wet%i))
   call handle_error(nf90_put_var(ncid,wet_j_id,wet%j))
   call handle_error(nf90_put_var(ncid,wet_x_id,wet%x))
   call handle_error(nf90_put_var(ncid,wet_y_id,wet%y))
   call handle_error(nf90_put_var(ncid,wet_area_id,wet%area))

   call handle_error(nf90_close(ncid))

end subroutine create_wet_file

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
    
end program create_model_wet
