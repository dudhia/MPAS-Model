module module_ra_rrtmg_lw

use mpas_kind_types,only: kind_phys => RKIND
use mpas_log
! use mpas_atmphys_constants,only : cp,g=>gravity
! use ra_rrtmg_vinterp,only: vinterp_ozn

  use parrrtm, only : nbndlw, ngptlw
! use rrtmg_lw_init, only: rrtmg_lw_ini
! use rrtmg_lw_rad, only: rrtmg_lw
! use mcica_subcol_gen_lw, only: mcica_subcol_lw
  use rrtmg_lw_mpas,only: rrtmg_lw_mpas_timestep_init, &
                          rrtmg_lw_mpas_timestep_final
  use rrtmg_lw_rad,only : rrtmg_lw_run


!   ! For buffer layer adjustment.  Steven Cavallo, Dec 2010.
!   integer , save    :: nlayers    
!   real, PARAMETER :: deltap = 4.  ! Pressure interval for buffer layer in mb
    
contains

!=================================================================================================================
 subroutine rrtmg_lwrad( &
                       p3d,p8w,pi3d,t3d,t8w,dz8w,qv3d,qc3d,qr3d,  &
                       qi3d,qs3d,qg3d,cldfra3d,o33d,tsk,emiss,    &
                       xland,xice,snow,icloud,o3input,noznlevels, &
                       pin,o3clim,glw,olr,lwcf,rthratenlw,        &
                       has_reqc,has_reqi,has_reqs,re_cloud,       &
                       re_ice,re_snow,rre_cloud,rre_ice,rre_snow, &
                       lwupt,lwuptc,lwdnt,lwdntc,                 &
                       lwupb,lwupbc,lwdnb,lwdnbc,                 &
                       lwupflx, lwupflxc, lwdnflx, lwdnflxc,      &
                       errmsg, errflg,                            &
                       ids,ide, jds,jde, kds,kde,                 & 
                       ims,ime, jms,jme, kms,kme,                 &
                       its,ite, jts,jte, kts,kte                  &
                       )

 implicit none

!=================================================================================================================

!--- input arguments:
 integer,intent(in):: ids,ide, jds,jde, kds,kde, &
                      ims,ime, jms,jme, kms,kme, &
                      its,ite, jts,jte, kts,kte
 integer,intent(in):: icloud,has_reqc,has_reqi,has_reqs
 integer,intent(in),optional:: o3input

 real,intent(in),dimension(ims:ime,jms:jme):: emiss,tsk,snow,xice,xland
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme):: t3d,p3d,pi3d
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme):: dz8w,p8w,t8w

 real,intent(in),dimension(ims:ime,kms:kme,jms:jme):: re_cloud,re_ice,re_snow
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    cldfra3d,qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,o33d

!--- additional input arguments to use the CAM ozone climatology:
 integer,intent(in):: noznlevels
 real,intent(in),dimension(1:noznlevels),optional:: pin
 real,intent(in),dimension(ims:ime,1:noznlevels,jms:jme),optional:: o3clim

!--- inout arguments:
 real,intent(inout),dimension(ims:ime,jms:jme):: glw,olr,lwcf
 real,intent(inout),dimension(ims:ime,jms:jme),optional:: &
    lwupt,lwuptc,lwdnt,lwdntc,lwupb,lwupbc,lwdnb,lwdnbc

 real,intent(inout),dimension(ims:ime,kms:kme,jms:jme):: rthratenlw
 
!--- output arguments:
  real,intent(out),dimension(ims:ime,kms:kme,jms:jme),optional:: &
     rre_cloud,rre_ice,rre_snow
  real,intent(out),dimension(ims:ime,kms:kme+2,jms:jme ),optional:: &
     lwupflx,lwupflxc,lwdnflx,lwdnflxc
 character(len=*),intent(out):: &
    errmsg       ! output error message (-).

 integer,intent(out):: &
    errflg       ! output error flag (-).


!local variables only needed for conversion from ims->its,ime->ite,kms->kts,kme->kte,etc in call to
!subroutine rrtmg_lw_timestep_init:
 real(kind=kind_phys),dimension(its:ite):: &
    xland_hv,xice_hv,snow_hv,emiss_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: &
    t3d_hv,p3d_hv,pi3d_hv,qv3d_hv,re_cloud_hv,re_ice_hv,re_snow_hv,dz8w_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: &
    p8w_hv,t8w_hv

 real(kind=kind_phys),dimension(its:ite,kts:kte):: cldfra3d_hv,qc3d_hv,qr3d_hv,qi3d_hv,qs3d_hv,qg3d_hv,o33d_hv
 real(kind=kind_phys),dimension(its:ite,1:noznlevels):: o3clim_hv

!local variables only needed for call to rrtmg_lw_run (column arrays allocated in rrtmg_lw_init)
 integer:: ncol,nlayall,icld,inflglw,iceflglw,liqflglw
 integer,dimension(its:ite):: mpas_nlay
 real(kind=kind_phys),dimension(:,:),allocatable:: play,tlay,h2ovmr,o3vmr,co2vmr,o2vmr,ch4vmr,n2ovmr
 real(kind=kind_phys),dimension(:,:),allocatable:: cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr
 real(kind=kind_phys),dimension(:,:),allocatable:: relqmcl,reicmcl,resnmcl,rei,rel,res
 real(kind=kind_phys),dimension(:,:,:),allocatable:: cldfmcl,clwpmcl,ciwpmcl,cswpmcl,taucmcl
 real(kind=kind_phys),dimension(:,:),allocatable:: hr,hrc

 real(kind=kind_phys),dimension(:,:),allocatable:: plev,tlev
 real(kind=kind_phys),dimension(:,:),allocatable:: uflx,dflx,uflxc,dflxc

 real(kind=kind_phys),dimension(:,:,:),allocatable:: tauaer

 real(kind=kind_phys),dimension(its:ite)::tsfc 
 real(kind=kind_phys),dimension(its:ite,nbndlw)::emis

!local variables only needed for conversion from ims->its,ime->ite,kms->kts,kme->kte,etc in call to
!subroutine rrtmg_lw_timestep_final:
 real(kind=kind_phys),dimension(its:ite):: glw_hv,olr_hv
 real(kind=kind_phys),dimension(its:ite):: lwcf_hv,lwupt_hv,lwuptc_hv,lwdnt_hv,lwdntc_hv, &
                                           lwupb_hv,lwupbc_hv,lwdnb_hv,lwdnbc_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte)  :: rthratenlw_hv
 real(kind=kind_phys),dimension(its:ite,kts:kme+1):: lwupflx_hv,lwupflxc_hv,lwdnflx_hv,lwdnflxc_hv


 integer:: i,j,k,kk,n
 
!-----------------------------------------------------------------------------------------------------------------

 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine rrtmg_lwrad:')

 do j = jts,jte
    do i = its,ite
       xland_hv(i)  = xland(i,j)
       xice_hv(i)   = xice(i,j)
       snow_hv(i)   = snow(i,j)
       emiss_hv(i)  = emiss(i,j)
       tsfc(i)      = tsk(i,j)
    enddo

    do k = kts,kte
       do i = its,ite
          t3d_hv(i,k)  = t3d(i,k,j)
          p3d_hv(i,k)  = p3d(i,k,j)
          pi3d_hv(i,k) = pi3d(i,k,j)
          qv3d_hv(i,k) = qv3d(i,k,j)
          dz8w_hv(i,k) = dz8w(i,k,j)
          p8w_hv(i,k)  = p8w(i,k,j)
          t8w_hv(i,k)  = t8w(i,k,j)
          re_cloud_hv(i,k) = re_cloud(i,k,j)
          re_ice_hv(i,k)   = re_ice(i,k,j)
          re_snow_hv(i,k)  = re_snow(i,k,j)
       enddo
    enddo
    k = kte+1
    do i = its,ite
       p8w_hv(i,k)  = p8w(i,k,j)
       t8w_hv(i,k)  = t8w(i,k,j)
    enddo

    !--- optional arrays:
    if(present(qc3d)) then
       do k = kts,kte
          do i = its,ite
             qc3d_hv(i,k) = qc3d(i,k,j)
          enddo
       enddo
    endif
    if(present(qr3d)) then
       do k = kts,kte
          do i = its,ite
             qr3d_hv(i,k) = qr3d(i,k,j)
          enddo
       enddo
    endif
    if(present(qi3d)) then
       do k = kts,kte
          do i = its,ite
             qi3d_hv(i,k) = qi3d(i,k,j)
          enddo
       enddo
    endif
    if(present(qs3d)) then
       do k = kts,kte
          do i = its,ite
             qs3d_hv(i,k) = qs3d(i,k,j)
          enddo
       enddo
    endif
    if(present(qg3d)) then
       do k = kts,kte
          do i = its,ite
             qg3d_hv(i,k) = qg3d(i,k,j)
          enddo
       enddo
    endif
    if(present(cldfra3d)) then
       do k = kts,kte
          do i = its,ite
             cldfra3d_hv(i,k) = cldfra3d(i,k,j)
          enddo
       enddo
    endif
    if(present(o3clim)) then
       do k = 1,noznlevels
          do i = its,ite
             o3clim_hv(i,k) = o3clim(i,k,j)
          enddo
       enddo
    endif
    if(present(o33d)) then
       do k = kts,kte
          do i = its,ite
             o33d_hv(i,k) = o33d(i,k,j)
          enddo
       enddo
    endif
 enddo

       call rrtmg_lw_mpas_timestep_init( &
                       p3d_hv,p8w_hv,pi3d_hv,t3d_hv,t8w_hv,dz8w_hv,qv3d_hv,qc3d_hv,qr3d_hv,  &
                       qi3d_hv,qs3d_hv,qg3d_hv,cldfra3d_hv,o33d_hv,tsfc,emiss_hv,    &
                       xland_hv,xice_hv,snow_hv,icloud,o3input,noznlevels, &
                       pin,o3clim_hv,                                &
                       has_reqc,has_reqi,has_reqs,re_cloud_hv,       &
                       re_ice_hv,re_snow_hv,                            &
                       ncol,nlayall,mpas_nlay,icld,play,plev,tlay,  &
                       tlev,tsfc,h2ovmr,o3vmr,co2vmr,ch4vmr,  &
                       n2ovmr,o2vmr,cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr, &
                       emis,inflglw,iceflglw,liqflglw,cldfmcl,taucmcl, &
                       ciwpmcl,clwpmcl,cswpmcl,reicmcl,relqmcl,resnmcl, &
                       rei,rel,res,                            &
                       tauaer,uflx,dflx,hr,uflxc,dflxc,hrc, &
                       ids,ide, kds,kde,                          &
                       ims,ime, kms,kme,                          &
                       its,ite, kts,kte,                          &
                       errmsg, errflg )


       !--- CALL TO THE RRTMG LONG WAVE RADIATION MODEL:
       call rrtmg_lw_run &
                  (ncol    , nlayall , mpas_nlay, icld     , play     , plev     , tlay    , &
                   tlev    , tsfc    , h2ovmr   , o3vmr    , co2vmr   , ch4vmr  , &
                   n2ovmr  , o2vmr   , cfc11vmr , cfc12vmr , cfc22vmr , ccl4vmr , &
                   emis    , inflglw , iceflglw , liqflglw , cldfmcl  , taucmcl , & 
                   ciwpmcl , clwpmcl , cswpmcl  , reicmcl  , relqmcl  , resnmcl , &
                   tauaer  , uflx    , dflx     , hr       , uflxc    , dflxc   , &
                   hrc,  &
                   errmsg, errflg )

       call rrtmg_lw_mpas_timestep_final( &
                       pi3d_hv,icloud,cldfra3d_hv,                      &
                       dflx, dflxc, uflx, uflxc,                  &
                       hr, hrc, rel, rei, res,                         &
                       play, plev, tlay, &
                       tlev, h2ovmr,o3vmr,co2vmr,ch4vmr, &
                       n2ovmr  , o2vmr   , cfc11vmr , cfc12vmr , cfc22vmr , ccl4vmr , &
                       cldfmcl, taucmcl , &
                       ciwpmcl , clwpmcl , cswpmcl  , reicmcl  , relqmcl  , resnmcl , &
                       tauaer  , &
                       glw_hv,olr_hv,lwcf_hv,rthratenlw_hv,            &
                       lwupt_hv,lwuptc_hv,lwdnt_hv,lwdntc_hv,                 &
                       lwupb_hv,lwupbc_hv,lwdnb_hv,lwdnbc_hv,                 &
                       lwupflx_hv, lwupflxc_hv, lwdnflx_hv, lwdnflxc_hv,      &
                       rre_cloud,rre_ice,rre_snow,                &
                       ids,ide, kds,kde,                          &
                       ims,ime, kms,kme,                          &
                       its,ite, kts,kte,                          &
                       errmsg, errflg )

 do j = jts,jte
    do i = its,ite
       glw(i,j) = glw_hv(i)
       olr(i,j) = olr_hv(i)
       do k = kts,kte
          rthratenlw(i,k,j) = rthratenlw_hv(i,k)
       enddo

       !optional arguments:
       if(present(lwupt)) then
          lwcf(i,j)   = lwcf_hv(i)
          lwupt(i,j)  = lwupt_hv(i)
          lwuptc(i,j) = lwuptc_hv(i)
          lwdnt(i,j)  = lwdnt_hv(i)
          lwdntc(i,j) = lwdntc_hv(i)
          lwupb(i,j)  = lwupb_hv(i)
          lwupbc(i,j) = lwupbc_hv(i)
          lwdnb(i,j)  = lwdnb_hv(i)
          lwdnbc(i,j) = lwdnbc_hv(i)
       endif
       if(present(lwupflx)) then
          do k = kts, kme+1
             lwupflx(i,k,j)  = lwupflx_hv(i,k)
             lwupflxc(i,k,j) = lwupflxc_hv(i,k)
             lwdnflx(i,k,j)  = lwdnflx_hv(i,k)
             lwdnflxc(i,k,j) = lwdnflxc_hv(i,k)
          enddo
       endif
    enddo
 enddo


 call mpas_log_write('--- end subroutine rrtmg_lwrad:')
 call mpas_log_write(' ')

 end subroutine rrtmg_lwrad

end module module_ra_rrtmg_lw
