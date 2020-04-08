!=================================================================================================================
 module rrtmg_lw_mpas
 use mpas_kind_types,only: kind_phys => RKIND
 use mpas_log

 use ra_rrtmg_shared,only    : inirad,reicalc,relcalc,retab
 use ra_rrtmg_vinterp,only   : vinterp_ozn
 use mcica_subcol_gen_lw,only: mcica_subcol_lw
 use mpas_atmphys_constants,only : cp,g=>gravity
 use ra_rrtmg_vinterp,only: vinterp_ozn

 use parrrtm, only : nbndlw, ngptlw
 use rrtmg_lw_init, only: rrtmg_lw_ini
 use rrtmg_lw_rad,only : rrtmg_lw_run


 implicit none
 private
 public:: rrtmg_lw_mpas_run,           &
          rrtmg_lw_mpas_timestep_init, &
          rrtmg_lw_mpas_timestep_final

    ! For buffer layer adjustment.  Steven Cavallo, Dec 2010.
    integer , save    :: nlayers
    real, PARAMETER :: deltap = 4.  ! Pressure interval for buffer layer in mb

 contains

!=================================================================================================================
 subroutine rrtmg_lw_mpas_run(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out)::errflg

 errmsg = 'rrtmg_lw_mpas_run OK'
 errflg = 0

 end subroutine rrtmg_lw_mpas_run


!=================================================================================================================
 subroutine rrtmg_lw_mpas_timestep_init( &
                       p2d,p8w,pi2d,t2d,t8w,dz8w,qv2d,qc2d,qr2d,  &
                       qi2d,qs2d,qg2d,cldfra2d,o32d,tsk,emiss,    &
                       xland,xice,snow,icloud,o3input,noznlevels, &
                       pin,o3clim,                                &
                       has_reqc,has_reqi,has_reqs,re_cloud,       &
                       re_ice,re_snow,                            &
                       ncol,nlay,icld,play,plev,tlay,  &
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

 implicit none

!=================================================================================================================

!--- input arguments:
 integer,intent(in):: ids,ide, kds,kde, &
                      ims,ime, kms,kme, &
                      its,ite, kts,kte
 integer,intent(in):: icloud,has_reqc,has_reqi,has_reqs
 integer,intent(in),optional:: o3input

 real,intent(in),dimension(ims:ime):: emiss,tsk,snow,xice,xland
 real,intent(in),dimension(ims:ime,kms:kme):: t2d,p2d,pi2d
 real,intent(in),dimension(ims:ime,kms:kme):: dz8w,p8w,t8w

 real,intent(in),dimension(ims:ime,kms:kme):: re_cloud,re_ice,re_snow
 real,intent(in),dimension(ims:ime,kms:kme),optional:: &
    cldfra2d,qv2d,qc2d,qr2d,qi2d,qs2d,qg2d,o32d

!--- additional input arguments to use the CAM ozone climatology:
 integer,intent(in):: noznlevels
 real,intent(in),dimension(1:noznlevels),optional:: pin
 real,intent(in),dimension(ims:ime,1:noznlevels),optional:: o3clim

!--- inout arguments:
 
!--- output arguments:

 integer,intent(out):: ncol,nlay,icld,inflglw,iceflglw,liqflglw
 real,intent(out),dimension(:,:),allocatable:: &
    plev,tlev,play,tlay,h2ovmr,o3vmr,co2vmr,o2vmr,ch4vmr,n2ovmr,cfc11vmr,  &
    cfc12vmr,cfc22vmr,ccl4vmr,relqmcl,reicmcl,resnmcl,rel,rei,res
 real,intent(out),dimension(:,:),allocatable:: uflx,dflx,uflxc,dflxc,hr,hrc
 real,intent(out),dimension(:,:,:),allocatable:: cldfmcl,clwpmcl,ciwpmcl,cswpmcl,taucmcl,tauaer
 real,intent(out),dimension(1)::tsfc
 real,intent(out),dimension(1,nbndlw):: emis
 character(len=*),intent(out):: &
    errmsg       ! output error message (-).

 integer,intent(out):: &
    errflg       ! output error flag (-).

!local variables and arrays:
 integer:: nb
 integer:: iplon,irng,permuteseed
 integer:: pcols,pver
 integer:: idx_rei
 integer:: i,k,kk,n
 
 real:: corr
 real:: gliqwp,gicewp,gsnowp,gravmks
 real:: snow_mass_factor
 real,dimension(1):: landfrac,landm,snowh,icefrac
 real,dimension(1,1:kte-kts+1):: pdel,cliqwp,cicewp,csnowp,reliq,reice,resnow
 real,dimension(:,:,:),allocatable:: taucld
 real,dimension(:,:),allocatable:: clwpth,ciwpth,cswpth,cldfrac

 real,dimension(kts:kte):: &
    cldfra1d,dz1d,p1d,t1d,qv1d,qc1d,qr1d,qi1d,qs1d,qg1d,o31d
 real,dimension(kts:kte+1):: pw1d,tw1d

!--- additional local variables and arrays needed to include additional layers between the model top
!    and the top of the atmosphere:
 integer,dimension(its:ite):: mpas_nlay

 real,dimension(:),allocatable:: o3mmr,varint

!--- additional local variables and arrays needed for the CAM ozone climatologyL
 real,dimension(1:noznlevels):: o3clim1d

#if defined(mpas)
!MPAS specific (Dom Heinzeller):
 integer:: nlayers
#endif

!--- set trace gas volume mixing ratios, 2005 values, IPCC (2007):
!carbon dioxide (379 ppmv)
 real :: co2
 data co2 / 379.e-6 / 
!methane (1774 ppbv)
 real :: ch4
 data ch4 / 1774.e-9 / 
!nitrous oxide (319 ppbv)
 real :: n2o
 data n2o / 319.e-9 / 
!cfc-11 (251 ppt)
 real :: cfc11
 data cfc11 / 0.251e-9 / 
!cfc-12 (538 ppt)
 real :: cfc12
 data cfc12 / 0.538e-9 / 
!cfc-22 (169 ppt)
 real :: cfc22
 data cfc22 / 0.169e-9 / 
!ccl4 (93 ppt)
 real :: ccl4
 data ccl4 / 0.093e-9 / 

!--- set oxygen volume mixing ratio (for o2mmr=0.23143):
 real :: o2
 data o2 / 0.209488 /

!--- molecular weights and ratios for converting mmr to vmr units
! real :: amd       ! Effective molecular weight of dry air (g/mol)  
! real :: amw       ! Molecular weight of water vapor (g/mol)        
! real :: amo       ! Molecular weight of ozone (g/mol)              
! real :: amo2      ! Molecular weight of oxygen (g/mol)              
! Atomic weights for conversion from mass to volume mixing ratios                
! data amd   /  28.9660   /                                                  
! data amw   /  18.0160   /                                                  
! data amo   /  47.9998   /                                                  
! data amo2  /  31.9999   /
                                                                                 
  real :: amdw     ! Molecular weight of dry air / water vapor  
  real :: amdo     ! Molecular weight of dry air / ozone
  real :: amdo2    ! Molecular weight of dry air / oxygen
  data amdw /  1.607793 /                                                    
  data amdo /  0.603461 /
  data amdo2 / 0.905190 /
    
!--- added for top of model adjustment.  Steven Cavallo NCAR/MMM December 2010
 integer,parameter:: nproflevs = 60 ! Constant, from the table
 integer:: klev

 real:: wght,vark,vark1
 real,dimension(nproflevs):: pprof,tprof

!--- weighted mean pressure and temperature profiles from midlatitude summer (MLS),
!    midlatitude winter (MLW), sub-Arctic winter (SAW),sub-Arctic summer (SAS),
!    and tropical (TROP) standard atmospheres.
 data pprof /1000.00,855.47,731.82,626.05,535.57,458.16, &
              391.94,335.29,286.83,245.38,209.91,179.57, &
              153.62,131.41,112.42, 96.17, 82.27, 70.38, &
               60.21, 51.51, 44.06, 37.69, 32.25, 27.59, &
               23.60, 20.19, 17.27, 14.77, 12.64, 10.81, &
                9.25,  7.91, 6.77,   5.79,  4.95,  4.24, &
                3.63,  3.10, 2.65,   2.27,  1.94,  1.66, &
                1.42,  1.22, 1.04,   0.89,  0.76,  0.65, &
                0.56,  0.48, 0.41,   0.35,  0.30,  0.26, &
                0.22,  0.19, 0.16,   0.14,  0.12,  0.10/
 data tprof /286.96,281.07,275.16,268.11,260.56,253.02,  &
             245.62,238.41,231.57,225.91,221.72,217.79,  &
             215.06,212.74,210.25,210.16,210.69,212.14,  &
             213.74,215.37,216.82,217.94,219.03,220.18,  &
             221.37,222.64,224.16,225.88,227.63,229.51,  &
             231.50,233.73,236.18,238.78,241.60,244.44,  &
             247.35,250.33,253.32,256.30,259.22,262.12,  &
             264.80,266.50,267.59,268.44,268.69,267.76,  &
             266.13,263.96,261.54,258.93,256.15,253.23,  &
             249.89,246.67,243.48,240.25,236.66,233.86/	

!-----------------------------------------------------------------------------------------------------------------

!--- all fields are ordered vertically from bottom to top (pressures are in mb):
 ncol = 1

!--- initialize option for the calculation of the cloud optical properties:
 icld     = 2 ! with clouds using maximum/random cloud overlap in subroutine mcica_subcol_lw. 
 inflglw  = 2
 iceflglw = 3
 liqflglw = 1

!--- latitude loop:

!--- longitude loop:
    i_loop: do i = its,ite

       !--- set surface emissivity in each RRTMG longwave band:
       do nb = 1, nbndlw
          emis(ncol,nb) = emiss(i)
       enddo



       !--- INITIALIZE COLUMN SOUNDING (the call to the long wave radiation code is done one column at a time
       !    which is why we set ncol = 1 above):
       do k = kts, kte+1
          pw1d(k) = p8w(i,k) / 100.
          tw1d(k) = t8w(i,k)
       enddo

       do k = kts, kte
          p1d(k)      = p2d(i,k) / 100.
          dz1d(k)     = dz8w(i,k)
          t1d(k)      = t2d(i,k)
          qv1d(k)     = amax1(qv2d(i,k),1.e-12)
          o31d(k)     = 0.
          qc1d(k)     = 0.
          qr1d(k)     = 0.
          qi1d(k)     = 0.
          qs1d(k)     = 0.
          qg1d(k)     = 0.
          cldfra1d(k) = 0.
          if(present(o32d)) o31d(k) = o32d(i,k)
       enddo

       !--- initialize the local arrays containing the different cloud water and ice condenstates:
       if(icloud .gt. 0) then
          do k = kts,kte
             if(present(qc2d)) qc1d(k) = amax1(qc2d(i,k),0.)
             if(present(qr2d)) qr1d(k) = amax1(qr2d(i,k),0.)
             if(present(qi2d)) qi1d(k) = amax1(qi2d(i,k),0.)
             if(present(qs2d)) qs1d(k) = amax1(qs2d(i,k),0.)
             if((present(qc2d) .or. present(qi2d) .or. present(qs2d)) .and. present(cldfra2d)) &
                cldfra1d(k) = cldfra2d(i,k)
          enddo
       endif

       !--- add extra layers to include absorption and transmission between the top of the model and the top of
       !    the atmosphere: in contrast to WRF, MPAS columns have different model-top pressures since MPAS uses
       !    a height coordinate system. Therefore, we define nlayers for each individual column:
       nlayers = kte + max(nint(pw1d(kte+1)/deltap), 1)
       mpas_nlay(i) = nlayers-kte
!      write(0,101) i,kme,kte,nlayers,mpas_nlay(i),pw1d(kte+1),pw1d(kte+1)-mpas_nlay(i)*deltap
!      101 format(6i9,3(1x,f9.4))
       if(.not.allocated(o3mmr)   ) allocate(o3mmr(kts:nlayers)     )
       if(.not.allocated(varint)  ) allocate(varint(kts:nlayers+1)  )
       if(.not.allocated(plev)    ) allocate(plev(1,kts:nlayers+1)  )
       if(.not.allocated(tlev)    ) allocate(tlev(1,kts:nlayers+1)  )
       if(.not.allocated(play)    ) allocate(play(1,kts:nlayers)    )
       if(.not.allocated(tlay)    ) allocate(tlay(1,kts:nlayers)    )
       if(.not.allocated(h2ovmr)  ) allocate(h2ovmr(1,kts:nlayers)  )
       if(.not.allocated(o3vmr)   ) allocate(o3vmr(1,kts:nlayers)   )
       if(.not.allocated(co2vmr)  ) allocate(co2vmr(1,kts:nlayers)  )
       if(.not.allocated(o2vmr)   ) allocate(o2vmr(1,kts:nlayers)   )
       if(.not.allocated(ch4vmr)  ) allocate(ch4vmr(1,kts:nlayers)  )
       if(.not.allocated(n2ovmr)  ) allocate(n2ovmr(1,kts:nlayers)  )
       if(.not.allocated(cfc11vmr)) allocate(cfc11vmr(1,kts:nlayers))
       if(.not.allocated(cfc12vmr)) allocate(cfc12vmr(1,kts:nlayers))
       if(.not.allocated(cfc22vmr)) allocate(cfc22vmr(1,kts:nlayers))
       if(.not.allocated(ccl4vmr) ) allocate(ccl4vmr(1,kts:nlayers) )
       if(.not.allocated(clwpth)  ) allocate(clwpth(1,kts:nlayers)  )
       if(.not.allocated(ciwpth)  ) allocate(ciwpth(1,kts:nlayers)  )
       if(.not.allocated(cswpth)  ) allocate(cswpth(1,kts:nlayers)  )
       if(.not.allocated(rel)     ) allocate(rel(1,kts:nlayers)     )
       if(.not.allocated(rei)     ) allocate(rei(1,kts:nlayers)     )
       if(.not.allocated(res)     ) allocate(res(1,kts:nlayers)     )
       if(.not.allocated(cldfrac) ) allocate(cldfrac(1,kts:nlayers) )
       if(.not.allocated(relqmcl) ) allocate(relqmcl(1,kts:nlayers) )
       if(.not.allocated(reicmcl) ) allocate(reicmcl(1,kts:nlayers) )
       if(.not.allocated(resnmcl) ) allocate(resnmcl(1,kts:nlayers) )
       if(.not.allocated(uflx)    ) allocate(uflx(1,kts:nlayers+1)  )
       if(.not.allocated(dflx)    ) allocate(dflx(1,kts:nlayers+1)  )
       if(.not.allocated(uflxc)   ) allocate(uflxc(1,kts:nlayers+1) )
       if(.not.allocated(dflxc)   ) allocate(dflxc(1,kts:nlayers+1) )
       if(.not.allocated(hr)      ) allocate(hr(1,kts:nlayers)      )
       if(.not.allocated(hrc)     ) allocate(hrc(1,kts:nlayers)     )
       if(.not.allocated(taucld)  ) allocate(taucld(nbndlw,1,kts:nlayers) )
       if(.not.allocated(cldfmcl) ) allocate(cldfmcl(ngptlw,1,kts:nlayers))         
       if(.not.allocated(clwpmcl) ) allocate(clwpmcl(ngptlw,1,kts:nlayers))
       if(.not.allocated(ciwpmcl) ) allocate(ciwpmcl(ngptlw,1,kts:nlayers))
       if(.not.allocated(cswpmcl) ) allocate(cswpmcl(ngptlw,1,kts:nlayers))
       if(.not.allocated(taucmcl) ) allocate(taucmcl(ngptlw,1,kts:nlayers))
       if(.not.allocated(tauaer)  ) allocate(tauaer(1,kts:nlayers,nbndlw) )

       nlay = nlayers ! keep these indices the same.

       !--- initialize local arrays called in the calculation of the cloud optical properties and radiative
       !    fluxes:
       do n = 1, ncol
          do k = kts, kte
             reliq(ncol,k)  = 10.
             reice(ncol,k)  = 10.
             resnow(ncol,k) = 10.
             cliqwp(ncol,k) = 0.
             cicewp(ncol,k) = 0.
             csnowp(ncol,k) = 0.
          enddo

          do k = 1, nlayers
             clwpth(n,k)  = 0.
             ciwpth(n,k)  = 0.
             cswpth(n,k)  = 0.
             rel(n,k)     = 0.
             rei(n,k)     = 0.
             res(n,k)     = 0.
             cldfrac(n,k) = 0.
             relqmcl(n,k) = 0.
             reicmcl(n,k) = 0.
             resnmcl(n,k) = 0.
             uflx(n,k)    = 0.
             uflxc(n,k)   = 0.
             dflx(n,k)    = 0.
             dflxc(n,k)   = 0.
             hr(n,k)      = 0.
             hrc(n,k)     = 0.
             taucld(1:nbndlw,n,k)  = 0.
             tauaer(n,k,1:nbndlw)  = 0.
             cldfmcl(1:ngptlw,n,k) = 0.
             clwpmcl(1:ngptlw,n,k) = 0.
             ciwpmcl(1:ngptlw,n,k) = 0.
             cswpmcl(1:ngptlw,n,k) = 0.
             taucmcl(1:ngptlw,n,k) = 0.
          enddo
          uflx(n,nlayers+1)  = 0.
          uflxc(n,nlayers+1) = 0.
          dflx(n,nlayers+1)  = 0.
          dflxc(n,nlayers+1) = 0.
       enddo

       !--- fill local arrays with input sounding. convert water vapor mass mixing ratio to volume mixing ratio:
       plev(ncol,1) = pw1d(1)
       tlev(ncol,1) = tw1d(1)
       tsfc(ncol)   = tsk(i)
       do k = kts, kte
          plev(ncol,k+1)   = pw1d(k+1)
          tlev(ncol,k+1)   = tw1d(k+1)
          play(ncol,k)     = p1d(k)
          pdel(ncol,k)     = plev(ncol,k) - plev(ncol,k+1)
          tlay(ncol,k)     = t1d(k)
          h2ovmr(ncol,k)   = qv1d(k) * amdw
          co2vmr(ncol,k)   = co2
          o2vmr(ncol,k)    = o2
          ch4vmr(ncol,k)   = ch4
          n2ovmr(ncol,k)   = n2o
          cfc11vmr(ncol,k) = cfc11
          cfc12vmr(ncol,k) = cfc12
          cfc22vmr(ncol,k) = cfc22
          ccl4vmr(ncol,k)  = ccl4
       enddo

       !--- the sourcecode below follows Steven Cavallo's method to "fill" the atmospheric layers between the
       !    top of the model and the top of the atmosphere. check if the pressure at the top of the atmosphere
       !    is negative. if negative, set it to zero prior to the calculation of temperatures (tlev and tlay):
       do k=kte+1,nlayers,1
          plev(ncol,k+1) = max(plev(ncol,k) - deltap, 0.00)
          play(ncol,k) = 0.5*(plev(ncol,k) + plev(ncol,k+1))
       enddo          

       !--- add zero as top level. this gets the temperature max at the stratopause, reducing downward flux
       !    errors in the top levels.  If zero happened to be the top level already, this will add another
       !    level with zero, but will not affect the radiative transfer calculation.
       plev(ncol,nlayers+1) = 0.00
       play(ncol,nlayers) =  0.5*(plev(ncol,nlayers) + plev(ncol,nlayers+1))

       !--- interpolate the table temperatures to column pressure levels:
       do k = 1, nlayers+1, 1
          if(pprof(nproflevs) .lt. plev(ncol,k)) then
             do kk = 2, nproflevs, 1       
                if(pprof(kk) .lt. plev(ncol,k)) then           
                   klev = kk - 1
                   exit
                endif
             enddo
          else
             klev = nproflevs
          endif  
  
          if(klev .ne. nproflevs) then
             vark  = tprof(klev) 
             vark1 = tprof(klev+1)
             wght  = (plev(ncol,k)-pprof(klev) )/( pprof(klev+1)-pprof(klev))
          else
             vark  = tprof(klev) 
             vark1 = tprof(klev)
             wght  = 0.0
          endif
          varint(k) = wght*(vark1-vark)+vark
       enddo                   
       
       !--- match the interpolated table temperature profile:
       do k = kte+1, nlayers+1, 1
          tlev(ncol,k) = varint(k) + (tlev(ncol,kte) - varint(kte))
          !if(k .le. nlay) then
          tlay(ncol,k-1) = 0.5*(tlev(ncol,k) + tlev(ncol,k-1))
          !endif
       enddo

       !--- fill water vapor and chemical species volume mixing ratios, except ozone:
       do k = kte+1, nlayers, 1
          h2ovmr(ncol,k)   = h2ovmr(ncol,kte) 
          co2vmr(ncol,k)   = co2vmr(ncol,kte) 
          o2vmr(ncol,k)    = o2vmr(ncol,kte) 
          ch4vmr(ncol,k)   = ch4vmr(ncol,kte) 
          n2ovmr(ncol,k)   = n2ovmr(ncol,kte) 
          cfc11vmr(ncol,k) = cfc11vmr(ncol,kte) 
          cfc12vmr(ncol,k) = cfc12vmr(ncol,kte) 
          cfc22vmr(ncol,k) = cfc22vmr(ncol,kte) 
          ccl4vmr(ncol,k)  = ccl4vmr(ncol,kte) 
       enddo	 

       !--- initialize the ozone voume mixing ratio:
       call inirad(o3mmr,plev,kts,nlayers-1)
       if(o3input .eq. 2) then
          do k = 1, noznlevels
             o3clim1d(k) = o3clim(i,k)
          enddo
          call vinterp_ozn(1,ncol,ncol,nlayers,play,pin,noznlevels,o3clim1d,o3mmr)
          do k = kts,nlayers
             o3vmr(ncol,k) = o3mmr(k)
          enddo
       else
          do k = kts,nlayers
             o3vmr(ncol,k) = o3mmr(k)*amdo
          enddo
       endif



       !--- CALCULATE CLOUD OPTICAL PROPERTIES:
       if(inflglw .gt. 0) then 

          do k = kts, kte
             cldfrac(ncol,k) = cldfra1d(k)
          enddo

          !--- zero out cloud optical properties here (not used when passing physical properties to rrtmg_lw as
          !    taucld is calculated in radiation:
          do k = kts, kte
             do nb = 1, nbndlw
                taucld(nb,ncol,k) = 0.0
             enddo
          enddo

          pcols = ncol
          pver  = kte - kts + 1
          gravmks = g

          if(has_reqc .ne. 0) then
             !--- fill the effective radius for cloud water with that calculated in the Thompson cloud
             !    microphysics parameterization:
             inflglw = 3
             do k = kts, kte
                reliq(ncol,k) = max(2.5,re_cloud(i,k)*1.e6)
                if(reliq(ncol,k).le.2.5 .and. cldfrac(ncol,k).gt.0. .and. &
                   (xland(i)-1.5).gt.0.) then        !--- ocean.
                   reliq(ncol,k) = 10.5
                elseif(reliq(ncol,k).le.2.5 .and. cldfrac(ncol,k).gt.0. .and. &
                   (xland(i)-1.5).lt.0.) then        !--- land.
                   reliq(ncol,k) = 7.5
                endif
             enddo

          else

             !---  calculate the effective radius for cloud water in layers below the model top. we added the
             !     dimensions pcols, pver so that the calls to subroutines relcalc and reicalc remain the same
             !     as in WRF although the two dimensions are duplicate and not needed inside the subroutines:
             landfrac(ncol) = 2.-xland(i)
             landm(ncol)    = landfrac(ncol)
             snowh(ncol)    = 0.001*snow(i)
             icefrac(ncol)  = xice(i)

             call relcalc(ncol,pcols,pver,tlay,landfrac,landm,icefrac,reliq,snowh)
          endif

          if(has_reqi .ne. 0) then
             !--- fill the effective radius for cloud ice with that calculated in the Thompson cloud
             !    microphysics parameterization:
             inflglw  = 4
             iceflglw = 4
             do k = kts, kte
                reice(ncol,k) = max(5.,re_ice(i,k)*1.e6)
                if(reice(ncol,k).le.5. .and. cldfrac(ncol,k).gt.0.) then
                   idx_rei = int(t2d(i,k)-179.)
                   idx_rei = min(max(idx_rei,1),75)
                   corr = t2d(i,k) - int(t2d(i,k))
                   reice(ncol,k) = retab(idx_rei)*(1.-corr) + retab(idx_rei+1)*corr
                   reice(ncol,k) = max(reice(ncol,k),5.0)
                endif
             enddo

          else

             !--- calculate the effective radius for cloud ice in layers below the model top, following
             !    Kristjansson and Mitchell. limit upper bound of reice for Fu ice parameterization and
             !    convert effective radius to generalized effective size (*1.0315; Fu, 1996):

             call reicalc(ncol,pcols,pver,tlay,reice)

             if(iceflglw .eq. 3) then
                do k = kts, kte
                   reice(ncol,k) = reice(ncol,k) * 1.0315
                   reice(ncol,k) = min(140.0,reice(ncol,k))
                enddo
             endif
          endif

          if(has_reqs .ne. 0) then
             !--- fill the effective radius for snow with that calculated in the Thompson cloud
             !    microphysics parameterization:
             inflglw  = 5
             iceflglw = 5
             do k = kts, kte
                resnow(ncol,k) = max(10.,re_snow(i,k)*1.e6)
             enddo
          else
             do k = kts, kte
                resnow(ncol,k) = 10.
             enddo
          endif

          !--- calculate the cloud liquid water path in units of g/m2 below the model top:
          do k = kts, kte
             gliqwp = qc1d(k) * pdel(ncol,k)*100.0 / gravmks * 1000.0              !grid box liquid water path.
             cliqwp(ncol,k) = gliqwp / max(0.01,cldfrac(ncol,k))                   !in-cloud liquid water path.
          enddo

          !--- calculate the cloud ice path in units of g/m2 below the model top:
          if(iceflglw .eq. 3)then
             do k = kts, kte
                gicewp = (qi1d(k)+qs1d(k)) * pdel(ncol,k)*100.0 / gravmks * 1000.0 !grid box ice water path.
                cicewp(ncol,k) = gicewp / max(0.01,cldfrac(ncol,k))                !in-cloud ice water path.
             enddo
          elseif(iceflglw .ge. 4) then
             do k = kts, kte
                gicewp = qi1d(k) * pdel(ncol,k)*100.0 / gravmks * 1000.0           !grid box ice water path.
                cicewp(ncol,k) = gicewp / max(0.01,cldfrac(ncol,k))                !in-cloud ice water path.
             enddo
          endif

          !--- calculate the snow path in units of g/m2 below the model top:
          if(iceflglw.eq.5)then
              do k = kts, kte
                 snow_mass_factor = 1.0
                 if(resnow(ncol,k) .gt. 130.)then
                    snow_mass_factor = (130.0/resnow(ncol,k))*(130.0/resnow(ncol,k))
                    resnow(ncol,k)   = 130.0
                 endif
                 gsnowp = qs1d(k) * pdel(ncol,k)*100.0 / gravmks * 1000.0           ! grid box snow path.
                 csnowp(ncol,k) = snow_mass_factor * gsnowp / max(0.01,cldfrac(ncol,k))
              enddo
           endif

          !--- set the cloud physical property between the model top and the top of the atmosphere. do not allow
          !    clouds between the model top and the top of the atmosphere: 
          do k = kts, kte
             clwpth(ncol,k) = cliqwp(ncol,k)
             ciwpth(ncol,k) = cicewp(ncol,k)
             cswpth(ncol,k) = csnowp(ncol,k)
             rel(ncol,k) = reliq(ncol,k)
             rei(ncol,k) = reice(ncol,k)
             res(ncol,k) = resnow(ncol,k)
          enddo
          do k = kte+1, nlayers
             clwpth(ncol,k)  = 0.
             ciwpth(ncol,k)  = 0.
             cswpth(ncol,k)  = 0.
             rel(ncol,k)     = 10.
             rei(ncol,k)     = 10.
             res(ncol,k)     = 10.
             cldfrac(ncol,k) = 0.
             do nb = 1, nbndlw
                taucld(nb,ncol,k) = 0.
             enddo
          enddo	 

       endif

       !--- sub-column generator for McICA:
       iplon = 1
       irng  = 0
       permuteseed = 150

       call mcica_subcol_lw &
                  (iplon   , ncol    , nlay    , icld    , permuteseed , irng    , play    , &
                   cldfrac , ciwpth  , clwpth  , cswpth  , rei         , rel     , res     , &
                   taucld  , cldfmcl , ciwpmcl , clwpmcl , cswpmcl     , reicmcl , relqmcl , & 
                   resnmcl , taucmcl)

       !--- initialization of aerosol optical depths:
       do nb = 1, nbndlw
          do k = kts, nlayers
             tauaer(ncol,k,nb) = 0.
          enddo
       enddo
    end do i_loop !end longitude loop.

       if(allocated(o3mmr)   ) deallocate(o3mmr   )
       if(allocated(varint)  ) deallocate(varint  )
       if(allocated(clwpth)  ) deallocate(clwpth  )
       if(allocated(ciwpth)  ) deallocate(ciwpth  )
       if(allocated(cswpth)  ) deallocate(cswpth  )
       if(allocated(cldfrac) ) deallocate(cldfrac )
       if(allocated(taucld)  ) deallocate(taucld  )

 errmsg = 'rrtmg_lw_mpas_timestep_init OK'
 errflg = 0

 call mpas_log_write('--- end subroutine rrtmg_lw_mpas_timestep_init:')

 end subroutine rrtmg_lw_mpas_timestep_init



!      !--- CALL TO THE RRTMG LONG WAVE RADIATION MODEL:
!      call rrtmg_lw &
!                 (ncol    , nlay    , icld     , play     , plev     , tlay    , &
!                  tlev    , tsfc    , h2ovmr   , o3vmr    , co2vmr   , ch4vmr  , &
!                  n2ovmr  , o2vmr   , cfc11vmr , cfc12vmr , cfc22vmr , ccl4vmr , &
!                  emis    , inflglw , iceflglw , liqflglw , cldfmcl  , taucmcl , & 
!                  ciwpmcl , clwpmcl , cswpmcl  , reicmcl  , relqmcl  , resnmcl , &
!                  tauaer  , uflx    , dflx     , hr       , uflxc    , dflxc   , &
!                  hrc)
  


!=================================================================================================================
subroutine rrtmg_lw_mpas_timestep_final( &
                       pi2d,icloud,cldfra2d,                      &
                       dflx, dflxc, uflx, uflxc,                  &
                       hr, hrc, rel, rei, res,                         &
                       play, plev, tlay, &
                       tlev, h2ovmr,o3vmr,co2vmr,ch4vmr, &
                       n2ovmr  , o2vmr   , cfc11vmr , cfc12vmr , cfc22vmr , ccl4vmr , &
                       cldfmcl, taucmcl , &
                       ciwpmcl , clwpmcl , cswpmcl  , reicmcl  , relqmcl  , resnmcl , &
                       tauaer  , &
                       glw,olr,lwcf,rthratenlw,                   &
                       lwupt,lwuptc,lwdnt,lwdntc,                 &
                       lwupb,lwupbc,lwdnb,lwdnbc,                 &
                       lwupflx, lwupflxc, lwdnflx, lwdnflxc,      &
                       rre_cloud,rre_ice,rre_snow,                &
                       ids,ide, kds,kde,                          &
                       ims,ime, kms,kme,                          &
                       its,ite, kts,kte,                          &
                       errmsg, errflg )
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: ids,ide, kds,kde, &
                      ims,ime, kms,kme, &
                      its,ite, kts,kte
 integer,intent(in):: icloud

 real,intent(in),dimension(ims:ime,kms:kme):: pi2d,cldfra2d

!--- inout arguments:
 real,intent(inout),dimension(:,:),allocatable:: dflx, dflxc, uflx, uflxc
 real,intent(inout),dimension(:,:),allocatable:: hr, hrc, rel, rei, res

 real,intent(inout),dimension(:,:),allocatable:: &
    plev,tlev,play,tlay,h2ovmr,o3vmr,co2vmr,o2vmr,ch4vmr,n2ovmr,cfc11vmr,  &
    cfc12vmr,cfc22vmr,ccl4vmr,relqmcl,reicmcl,resnmcl
 real,intent(inout),dimension(:,:,:),allocatable:: cldfmcl,clwpmcl,ciwpmcl,cswpmcl,taucmcl,tauaer

 
!--- output arguments:
 real,intent(out),dimension(ims:ime,kms:kme),optional:: lwupflx, lwupflxc, lwdnflx, lwdnflxc
 real,intent(out),dimension(ims:ime,kms:kme),optional:: rre_cloud,rre_ice,rre_snow
 real,intent(out),dimension(ims:ime,kms:kme):: rthratenlw
 real,intent(out),dimension(ims:ime):: glw,olr,lwcf
 real,intent(out),dimension(ims:ime),optional:: lwupt,lwuptc,lwdnt,lwdntc
 real,intent(out),dimension(ims:ime),optional:: lwupb,lwupbc,lwdnb,lwdnbc
 character(len=*),intent(out):: &
    errmsg       ! output error message (-).

 integer,intent(out):: &
    errflg       ! output error flag (-).

!local variables and arrays:
 integer:: i,k,ncol
 real,dimension(kts:kte):: tten1d

    ncol = 1

    i_loop: do i = its,ite

       !--- OUTPUTS:
       glw(i)  = dflx(1,1)
       olr(i)  = uflx(1,nlayers+1)
       lwcf(i) = uflxc(1,nlayers+1) - uflx(1,nlayers+1)

       if(present(lwupt)) then
          !output up and down toa fluxes for total and clear sky:
          lwupt(i)  = uflx(1,nlayers+1)
          lwuptc(i) = uflxc(1,nlayers+1)
          lwdnt(i)  = dflx(1,nlayers+1)
          lwdntc(i) = dflxc(1,nlayers+1)
          !output up and down surface fluxes for total and clear sky:
          lwupb(i)  = uflx(1,1)
          lwupbc(i) = uflxc(1,1)
          lwdnb(i)  = dflx(1,1)
          lwdnbc(i) = dflxc(1,1)
       endif

       if(present(lwupflx)) then
          !output up and down fluxes:
          do k=kts,nlayers+1
             lwupflx(i,k)  = uflx(1,k)
             lwupflxc(i,k) = uflxc(1,k)
             lwdnflx(i,k)  = dflx(1,k)
             lwdnflxc(i,k) = dflxc(1,k)
          enddo
       endif

       !--- initialize the local radiative heating rate:
       do k = kts, kte
          tten1d(k) = 0.
       enddo

       !--- output heating rate tendency. convert heating rate from K per day to K per second:
       do k = kts, kte
          tten1d(k) = hr(ncol,k)/86400.
          rthratenlw(i,k) = tten1d(k)/pi2d(i,k)
       enddo

       !--- output the effective radii for cloud water, cloud ice, and snow:
       if(icloud.gt.0 .and. present(rre_cloud) .and. present(rre_ice) .and. present(rre_snow)) then
          do k = kts, kte
             if(cldfra2d(i,k) .gt. 0.) then
                rre_cloud(i,k) = rel(ncol,k)
                rre_ice(i,k)   = rei(ncol,k)
                rre_snow(i,k)  = res(ncol,k)
             endif
          enddo
       endif

       !--- deallocate local column arrays:
!      if(allocated(o3mmr)   ) deallocate(o3mmr   )
!      if(allocated(varint)  ) deallocate(varint  )
       if(allocated(plev)    ) deallocate(plev    )
       if(allocated(tlev)    ) deallocate(tlev    )
       if(allocated(play)    ) deallocate(play    )
       if(allocated(tlay)    ) deallocate(tlay    )
       if(allocated(h2ovmr)  ) deallocate(h2ovmr  )
       if(allocated(o3vmr)   ) deallocate(o3vmr   )
       if(allocated(co2vmr)  ) deallocate(co2vmr  )
       if(allocated(o2vmr)   ) deallocate(o2vmr   )
       if(allocated(ch4vmr)  ) deallocate(ch4vmr  )
       if(allocated(n2ovmr)  ) deallocate(n2ovmr  )
       if(allocated(cfc11vmr)) deallocate(cfc11vmr)
       if(allocated(cfc12vmr)) deallocate(cfc12vmr)
       if(allocated(cfc22vmr)) deallocate(cfc22vmr)
       if(allocated(ccl4vmr) ) deallocate(ccl4vmr )
!      if(allocated(clwpth)  ) deallocate(clwpth  )
!      if(allocated(ciwpth)  ) deallocate(ciwpth  )
!      if(allocated(cswpth)  ) deallocate(cswpth  )
       if(allocated(rel)     ) deallocate(rel     )
       if(allocated(rei)     ) deallocate(rei     )
       if(allocated(res)     ) deallocate(res     )
!      if(allocated(cldfrac) ) deallocate(cldfrac )
       if(allocated(relqmcl) ) deallocate(relqmcl )
       if(allocated(reicmcl) ) deallocate(reicmcl )
       if(allocated(resnmcl) ) deallocate(resnmcl )
       if(allocated(uflx)    ) deallocate(uflx    )
       if(allocated(dflx)    ) deallocate(dflx    )
       if(allocated(uflxc)   ) deallocate(uflxc   )
       if(allocated(dflxc)   ) deallocate(dflxc   )
       if(allocated(hr)      ) deallocate(hr      )
       if(allocated(hrc)     ) deallocate(hrc     )
!      if(allocated(taucld)  ) deallocate(taucld  )
       if(allocated(cldfmcl) ) deallocate(cldfmcl )         
       if(allocated(clwpmcl) ) deallocate(clwpmcl )
       if(allocated(ciwpmcl) ) deallocate(ciwpmcl )
       if(allocated(cswpmcl) ) deallocate(cswpmcl )
       if(allocated(taucmcl) ) deallocate(taucmcl )
       if(allocated(tauaer)  ) deallocate(tauaer  )

    end do i_loop !end longitude loop.

 errmsg = 'rrtmg_lw_mpas_timestep_final OK'
 errflg = 0

 call mpas_log_write('--- end subroutine rrtmg_lw_mpas_timestep_final:')


 end subroutine rrtmg_lw_mpas_timestep_final

!=================================================================================================================
 end module rrtmg_lw_mpas
!=================================================================================================================
