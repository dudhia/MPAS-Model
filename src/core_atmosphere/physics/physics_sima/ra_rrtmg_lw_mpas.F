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

 integer,intent(out):: ncol,icld,inflglw,iceflglw,liqflglw,nlayall
 integer,intent(out),dimension(its:ite)::mpas_nlay
 real,intent(out),dimension(:,:),allocatable:: &
    plev,tlev,play,tlay,h2ovmr,o3vmr,co2vmr,o2vmr,ch4vmr,n2ovmr,cfc11vmr,  &
    cfc12vmr,cfc22vmr,ccl4vmr,relqmcl,reicmcl,resnmcl,rel,rei,res
 real,intent(out),dimension(:,:),allocatable:: uflx,dflx,uflxc,dflxc,hr,hrc
 real,intent(out),dimension(:,:,:),allocatable:: cldfmcl,clwpmcl,ciwpmcl,cswpmcl,taucmcl,tauaer
 real,intent(out),dimension(its:ite)::tsfc
 real,intent(out),dimension(its:ite,nbndlw):: emis
 character(len=*),intent(out):: &
    errmsg       ! output error message (-).

 integer,intent(out):: &
    errflg       ! output error flag (-).

!local variables and arrays:
 integer:: nb,nlay
 integer:: iplon,irng,permuteseed
 integer:: pcols,pver
 integer:: idx_rei
 integer:: i,k,kk,n
 
 real:: corr
 real:: gliqwp,gicewp,gsnowp,gravmks
 real:: snow_mass_factor
 real,dimension(1):: landfrac,landm,snowh,icefrac
 real,dimension(its:ite,1:kte-kts+1):: pdel,cliqwp,cicewp,csnowp,reliq,reice,resnow
 real,dimension(:,:,:),allocatable:: taucld
 real,dimension(:,:),allocatable:: clwpth,ciwpth,cswpth,cldfrac

 real,dimension(kts:kte):: &
    cldfra1d,dz1d,p1d,t1d,qv1d,qc1d,qr1d,qi1d,qs1d,qg1d,o31d
 real,dimension(kts:kte+1):: pw1d,tw1d

!--- additional local variables and arrays needed to include additional layers between the model top
!    and the top of the atmosphere:

 real,dimension(:),allocatable:: o3mmr,varint,p1dlev,p1dlay,t1dlay,rel1d,rei1d

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

 call mpas_log_write('--- enter subroutine rrtmg_lw_mpas_timestep_init:')
!--- all fields are ordered vertically from bottom to top (pressures are in mb):
 ncol = ite-its+1

!--- initialize option for the calculation of the cloud optical properties:
 icld     = 2 ! with clouds using maximum/random cloud overlap in subroutine mcica_subcol_lw. 
 inflglw  = 2
 iceflglw = 3
 liqflglw = 1

       nlayers  = kte + 10 ! temporary value for allocation i-dependent
       nlayall = nlayers

       if(.not.allocated(o3mmr)   ) allocate(o3mmr(kts:nlayers)     )
       if(.not.allocated(varint)  ) allocate(varint(kts:nlayers+1)  )
       if(.not.allocated(p1dlev)  ) allocate(p1dlev(kts:nlayers+1)  )
       if(.not.allocated(p1dlay)  ) allocate(p1dlay(kts:nlayers)    )
       if(.not.allocated(t1dlay)  ) allocate(t1dlay(kts:nlayers)    )
       if(.not.allocated(rel1d)   ) allocate(rel1d(kts:nlayers)     )
       if(.not.allocated(rei1d)   ) allocate(rei1d(kts:nlayers)     )
       if(.not.allocated(plev)    ) allocate(plev(its:ite,kts:nlayers+1)  )
       if(.not.allocated(tlev)    ) allocate(tlev(its:ite,kts:nlayers+1)  )
       if(.not.allocated(play)    ) allocate(play(its:ite,kts:nlayers)    )
       if(.not.allocated(tlay)    ) allocate(tlay(its:ite,kts:nlayers)    )
       if(.not.allocated(h2ovmr)  ) allocate(h2ovmr(its:ite,kts:nlayers)  )
       if(.not.allocated(o3vmr)   ) allocate(o3vmr(its:ite,kts:nlayers)   )
       if(.not.allocated(co2vmr)  ) allocate(co2vmr(its:ite,kts:nlayers)  )
       if(.not.allocated(o2vmr)   ) allocate(o2vmr(its:ite,kts:nlayers)   )
       if(.not.allocated(ch4vmr)  ) allocate(ch4vmr(its:ite,kts:nlayers)  )
       if(.not.allocated(n2ovmr)  ) allocate(n2ovmr(its:ite,kts:nlayers)  )
       if(.not.allocated(cfc11vmr)) allocate(cfc11vmr(its:ite,kts:nlayers))
       if(.not.allocated(cfc12vmr)) allocate(cfc12vmr(its:ite,kts:nlayers))
       if(.not.allocated(cfc22vmr)) allocate(cfc22vmr(its:ite,kts:nlayers))
       if(.not.allocated(ccl4vmr) ) allocate(ccl4vmr(its:ite,kts:nlayers) )
       if(.not.allocated(clwpth)  ) allocate(clwpth(its:ite,kts:nlayers)  )
       if(.not.allocated(ciwpth)  ) allocate(ciwpth(its:ite,kts:nlayers)  )
       if(.not.allocated(cswpth)  ) allocate(cswpth(its:ite,kts:nlayers)  )
       if(.not.allocated(rel)     ) allocate(rel(its:ite,kts:nlayers)     )
       if(.not.allocated(rei)     ) allocate(rei(its:ite,kts:nlayers)     )
       if(.not.allocated(res)     ) allocate(res(its:ite,kts:nlayers)     )
       if(.not.allocated(cldfrac) ) allocate(cldfrac(its:ite,kts:nlayers) )
       if(.not.allocated(relqmcl) ) allocate(relqmcl(its:ite,kts:nlayers) )
       if(.not.allocated(reicmcl) ) allocate(reicmcl(its:ite,kts:nlayers) )
       if(.not.allocated(resnmcl) ) allocate(resnmcl(its:ite,kts:nlayers) )
       if(.not.allocated(uflx)    ) allocate(uflx(its:ite,kts:nlayers+1)  )
       if(.not.allocated(dflx)    ) allocate(dflx(its:ite,kts:nlayers+1)  )
       if(.not.allocated(uflxc)   ) allocate(uflxc(its:ite,kts:nlayers+1) )
       if(.not.allocated(dflxc)   ) allocate(dflxc(its:ite,kts:nlayers+1) )
       if(.not.allocated(hr)      ) allocate(hr(its:ite,kts:nlayers)      )
       if(.not.allocated(hrc)     ) allocate(hrc(its:ite,kts:nlayers)     )
       if(.not.allocated(taucld)  ) allocate(taucld(nbndlw,its:ite,kts:nlayers) )
       if(.not.allocated(cldfmcl) ) allocate(cldfmcl(ngptlw,its:ite,kts:nlayers))         
       if(.not.allocated(clwpmcl) ) allocate(clwpmcl(ngptlw,its:ite,kts:nlayers))
       if(.not.allocated(ciwpmcl) ) allocate(ciwpmcl(ngptlw,its:ite,kts:nlayers))
       if(.not.allocated(cswpmcl) ) allocate(cswpmcl(ngptlw,its:ite,kts:nlayers))
       if(.not.allocated(taucmcl) ) allocate(taucmcl(ngptlw,its:ite,kts:nlayers))
       if(.not.allocated(tauaer)  ) allocate(tauaer(its:ite,kts:nlayers,nbndlw) )

!--- latitude loop:

!--- longitude loop:
    i_loop: do i = its,ite

       !--- set surface emissivity in each RRTMG longwave band:
       do nb = 1, nbndlw
          emis(i,nb) = emiss(i)
       enddo


       !--- INITIALIZE COLUMN SOUNDING 

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

       !--- initialize the local arrays cont
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
!       mpas_nlay(i) = nlayers-kte
       mpas_nlay(i) = nlayers
!      write(0,*) i,kme,kte,nlayers,mpas_nlay(i),pw1d(kte+1),deltap
!      101 format(6i9,3(1x,f9.4))

       nlay = nlayall ! used for allocated vertical levels

       !--- initialize local arrays called in the calculation of the cloud optical properties and radiative
       !    fluxes:
       do k = kts, kte
          reliq(i,k)  = 10.
          reice(i,k)  = 10.
          resnow(i,k) = 10.
          cliqwp(i,k) = 0.
          cicewp(i,k) = 0.
          csnowp(i,k) = 0.
       enddo

       do k = 1, nlayers
          clwpth(i,k)  = 0.
          ciwpth(i,k)  = 0.
          cswpth(i,k)  = 0.
          rel(i,k)     = 0.
          rei(i,k)     = 0.
          res(i,k)     = 0.
          cldfrac(i,k) = 0.
          relqmcl(i,k) = 0.
          reicmcl(i,k) = 0.
          resnmcl(i,k) = 0.
          uflx(i,k)    = 0.
          uflxc(i,k)   = 0.
          dflx(i,k)    = 0.
          dflxc(i,k)   = 0.
          hr(i,k)      = 0.
          hrc(i,k)     = 0.
          taucld(1:nbndlw,i,k)  = 0.
          tauaer(i,k,1:nbndlw)  = 0.
          cldfmcl(1:ngptlw,i,k) = 0.
          clwpmcl(1:ngptlw,i,k) = 0.
          ciwpmcl(1:ngptlw,i,k) = 0.
          cswpmcl(1:ngptlw,i,k) = 0.
          taucmcl(1:ngptlw,i,k) = 0.
       enddo
       uflx(i,nlayers+1)  = 0.
       uflxc(i,nlayers+1) = 0.
       dflx(i,nlayers+1)  = 0.
       dflxc(i,nlayers+1) = 0.

    !--- fill local arrays with input sounding. convert water vapor mass mixing ratio to volume mixing ratio:
       plev(i,1) = pw1d(1)
       tlev(i,1) = tw1d(1)
       tsfc(i)   = tsk(i)
       do k = kts, kte
          plev(i,k+1)   = pw1d(k+1)
          tlev(i,k+1)   = tw1d(k+1)
          play(i,k)     = p1d(k)
          pdel(i,k)     = plev(i,k) - plev(i,k+1)
          tlay(i,k)     = t1d(k)
          h2ovmr(i,k)   = qv1d(k) * amdw
          co2vmr(i,k)   = co2
          o2vmr(i,k)    = o2
          ch4vmr(i,k)   = ch4
          n2ovmr(i,k)   = n2o
          cfc11vmr(i,k) = cfc11
          cfc12vmr(i,k) = cfc12
          cfc22vmr(i,k) = cfc22
          ccl4vmr(i,k)  = ccl4
       enddo

       !--- the sourcecode below follows Steven Cavallo's method to "fill" the atmospheric layers between the
       !    top of the model and the top of the atmosphere. check if the pressure at the top of the atmosphere
       !    is negative. if negative, set it to zero prior to the calculation of temperatures (tlev and tlay):
       do k=kte+1,nlayers,1
          plev(i,k+1) = max(plev(i,k) - deltap, 0.00)
          play(i,k) = 0.5*(plev(i,k) + plev(i,k+1))
       enddo          

       !--- add zero as top level. this gets the temperature max at the stratopause, reducing downward flux
       !    errors in the top levels.  If zero happened to be the top level already, this will add another
       !    level with zero, but will not affect the radiative transfer calculation.
       plev(i,nlayers+1) = 0.00
       play(i,nlayers) =  0.5*(plev(i,nlayers) + plev(i,nlayers+1))

       !--- interpolate the table temperatures to column pressure levels:
       do k = 1, nlayers+1, 1
          if(pprof(nproflevs) .lt. plev(i,k)) then
             do kk = 2, nproflevs, 1       
                if(pprof(kk) .lt. plev(i,k)) then           
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
             wght  = (plev(i,k)-pprof(klev) )/( pprof(klev+1)-pprof(klev))
          else
             vark  = tprof(klev) 
             vark1 = tprof(klev)
             wght  = 0.0
          endif
          varint(k) = wght*(vark1-vark)+vark
       enddo                   
       
       !--- match the interpolated table temperature profile:
       do k = kte+1, nlayers+1, 1
          tlev(i,k) = varint(k) + (tlev(i,kte) - varint(kte))
          !if(k .le. nlay) then
          tlay(i,k-1) = 0.5*(tlev(i,k) + tlev(i,k-1))
          !endif
       enddo

       !--- fill water vapor and chemical species volume mixing ratios, except ozone:
       do k = kte+1, nlayers, 1
          h2ovmr(i,k)   = h2ovmr(i,kte) 
          co2vmr(i,k)   = co2vmr(i,kte) 
          o2vmr(i,k)    = o2vmr(i,kte) 
          ch4vmr(i,k)   = ch4vmr(i,kte) 
          n2ovmr(i,k)   = n2ovmr(i,kte) 
          cfc11vmr(i,k) = cfc11vmr(i,kte) 
          cfc12vmr(i,k) = cfc12vmr(i,kte) 
          cfc22vmr(i,k) = cfc22vmr(i,kte) 
          ccl4vmr(i,k)  = ccl4vmr(i,kte) 
       enddo	 

       !--- initialize the ozone voume mixing ratio:
       pcols = 1
       do k = kts,nlayers
          p1dlay(k) = play(i,k)
          p1dlev(k) = plev(i,k)
       enddo
       p1dlev(nlayers+1) = plev(i,nlayers+1)
       call inirad(o3mmr,p1dlev,kts,nlayers-1)

       if(o3input .eq. 2) then
          do k = 1, noznlevels
             o3clim1d(k) = o3clim(i,k)
          enddo
          call vinterp_ozn(1,pcols,pcols,nlayers,p1dlay,pin,noznlevels,o3clim1d,o3mmr)
          do k = kts,nlayers
             o3vmr(i,k) = o3mmr(k)
          enddo
       else
          do k = kts,nlayers
             o3vmr(i,k) = o3mmr(k)*amdo
          enddo
       endif

       !--- initialization of aerosol optical depths:
       do nb = 1, nbndlw
          do k = kts, nlayers
             tauaer(i,k,nb) = 0.
          enddo
       enddo


       !--- CALCULATE CLOUD OPTICAL PROPERTIES:
       if(inflglw .gt. 0) then 

          do k = kts, kte
             cldfrac(i,k) = cldfra1d(k)
          enddo

          !--- zero out cloud optical properties here (not used when passing physical properties to rrtmg_lw as
          !    taucld is calculated in radiation:
          do k = kts, kte
             do nb = 1, nbndlw
                taucld(nb,i,k) = 0.0
             enddo
          enddo

!         pcols = 1
          pver  = kte - kts + 1
          gravmks = g

          if(has_reqc .ne. 0) then
             !--- fill the effective radius for cloud water with that calculated in the Thompson cloud
             !    microphysics parameterization:
             inflglw = 3
             do k = kts, kte
                reliq(i,k) = max(2.5,re_cloud(i,k)*1.e6)
                if(reliq(i,k).le.2.5 .and. cldfrac(i,k).gt.0. .and. &
                   (xland(i)-1.5).gt.0.) then        !--- ocean.
                   reliq(i,k) = 10.5
                elseif(reliq(i,k).le.2.5 .and. cldfrac(i,k).gt.0. .and. &
                   (xland(i)-1.5).lt.0.) then        !--- land.
                   reliq(i,k) = 7.5
                endif
             enddo

          else

             !---  calculate the effective radius for cloud water in layers below the model top. we added the
             !     dimensions pcols, pver so that the calls to subroutines relcalc and reicalc remain the same
             !     as in WRF although the two dimensions are duplicate and not needed inside the subroutines:
             landfrac(pcols) = 2.-xland(i)
             landm(pcols)    = landfrac(pcols)
             snowh(pcols)    = 0.001*snow(i)
             icefrac(pcols)  = xice(i)

             do k = kts,nlayers
                t1dlay(k) = tlay(i,k)
             enddo
             call relcalc(pcols,pcols,pver,t1dlay,landfrac,landm,icefrac,rel1d,snowh)
             do k = kts,nlayers
                rel(i,k) = rel1d(k)
             enddo

          endif

          if(has_reqi .ne. 0) then
             !--- fill the effective radius for cloud ice with that calculated in the Thompson cloud
             !    microphysics parameterization:
             inflglw  = 4
             iceflglw = 4
             do k = kts, kte
                reice(i,k) = max(5.,re_ice(i,k)*1.e6)
                if(reice(i,k).le.5. .and. cldfrac(i,k).gt.0.) then
                   idx_rei = int(tlay(i,k)-179.)
                   idx_rei = min(max(idx_rei,1),75)
                   corr = tlay(i,k) - int(tlay(i,k))
                   reice(i,k) = retab(idx_rei)*(1.-corr) + retab(idx_rei+1)*corr
                   reice(i,k) = max(reice(i,k),5.0)
                endif
             enddo

          else

             !--- calculate the effective radius for cloud ice in layers below the model top, following
             !    Kristjansson and Mitchell. limit upper bound of reice for Fu ice parameterization and
             !    convert effective radius to generalized effective size (*1.0315; Fu, 1996):

             do k = kts,nlayers
                t1dlay(k) = tlay(i,k)
             enddo
             call reicalc(pcols,pcols,pver,t1dlay,rei1d)
             do k = kts,nlayers
                rei(i,k) = rei1d(k)
             enddo

             if(iceflglw .eq. 3) then
                do k = kts, kte
                   reice(i,k) = reice(i,k) * 1.0315
                   reice(i,k) = min(140.0,reice(i,k))
                enddo
             endif
          endif

          if(has_reqs .ne. 0) then
             !--- fill the effective radius for snow with that calculated in the Thompson cloud
             !    microphysics parameterization:
             inflglw  = 5
             iceflglw = 5
             do k = kts, kte
                resnow(i,k) = max(10.,re_snow(i,k)*1.e6)
             enddo
          else
             do k = kts, kte
                resnow(i,k) = 10.
             enddo
          endif

          !--- calculate the cloud liquid water path in units of g/m2 below the model top:
          do k = kts, kte
             gliqwp = qc1d(k) * pdel(i,k)*100.0 / gravmks * 1000.0              !grid box liquid water path.
             cliqwp(i,k) = gliqwp / max(0.01,cldfrac(i,k))                   !in-cloud liquid water path.
          enddo

          !--- calculate the cloud ice path in units of g/m2 below the model top:
          if(iceflglw .eq. 3)then
             do k = kts, kte
                gicewp = (qi1d(k)+qs1d(k)) * pdel(i,k)*100.0 / gravmks * 1000.0 !grid box ice water path.
                cicewp(i,k) = gicewp / max(0.01,cldfrac(i,k))                !in-cloud ice water path.
             enddo
          elseif(iceflglw .ge. 4) then
             do k = kts, kte
                gicewp = qi1d(k) * pdel(i,k)*100.0 / gravmks * 1000.0           !grid box ice water path.
                cicewp(i,k) = gicewp / max(0.01,cldfrac(i,k))                !in-cloud ice water path.
             enddo
          endif

          !--- calculate the snow path in units of g/m2 below the model top:
          if(iceflglw.eq.5)then
              do k = kts, kte
                 snow_mass_factor = 1.0
                 if(resnow(i,k) .gt. 130.)then
                    snow_mass_factor = (130.0/resnow(i,k))*(130.0/resnow(i,k))
                    resnow(i,k)   = 130.0
                 endif
                 gsnowp = qs1d(k) * pdel(i,k)*100.0 / gravmks * 1000.0           ! grid box snow path.
                 csnowp(i,k) = snow_mass_factor * gsnowp / max(0.01,cldfrac(i,k))
              enddo
           endif

          !--- set the cloud physical property between the model top and the top of the atmosphere. do not allow
          !    clouds between the model top and the top of the atmosphere: 
          ! In sw code these 6 local variables are not used and this copy is not needed
          do k = kts, kte
             clwpth(i,k) = cliqwp(i,k)
             ciwpth(i,k) = cicewp(i,k)
             cswpth(i,k) = csnowp(i,k)
             rel(i,k) = reliq(i,k)
             rei(i,k) = reice(i,k)
             res(i,k) = resnow(i,k)
          enddo
          do k = kte+1, nlayers
             clwpth(i,k)  = 0.
             ciwpth(i,k)  = 0.
             cswpth(i,k)  = 0.
             rel(i,k)     = 10.
             rei(i,k)     = 10.
             res(i,k)     = 10.
             cldfrac(i,k) = 0.
             do nb = 1, nbndlw
                taucld(nb,i,k) = 0.
             enddo
          enddo	 

       endif
    end do i_loop !end longitude loop.

       !--- sub-column generator for McICA:
       iplon = 1
       irng  = 0
       permuteseed = 150

       call mcica_subcol_lw &
                  (iplon   , ncol    , nlay    , icld    , permuteseed , irng    , play    , &
                   cldfrac , ciwpth  , clwpth  , cswpth  , rei         , rel     , res     , &
                   taucld  , cldfmcl , ciwpmcl , clwpmcl , cswpmcl     , reicmcl , relqmcl , & 
                   resnmcl , taucmcl)


       if(allocated(o3mmr)   ) deallocate(o3mmr   )
       if(allocated(varint)  ) deallocate(varint  )
       if(allocated(p1dlev)  ) deallocate(p1dlev  )
       if(allocated(p1dlay)  ) deallocate(p1dlay  )
       if(allocated(t1dlay)  ) deallocate(t1dlay  )
       if(allocated(rel1d)   ) deallocate(rel1d   )
       if(allocated(rei1d)   ) deallocate(rei1d   )
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

 call mpas_log_write('--- enter subroutine rrtmg_lw_mpas_timestep_final:')
    ncol = ite-its+1

    i_loop: do i = its,ite

       !--- OUTPUTS:
       glw(i)  = dflx(i,1)
       olr(i)  = uflx(i,nlayers+1)
       lwcf(i) = uflxc(i,nlayers+1) - uflx(i,nlayers+1)

       if(present(lwupt)) then
          !output up and down toa fluxes for total and clear sky:
          lwupt(i)  = uflx(i,nlayers+1)
          lwuptc(i) = uflxc(i,nlayers+1)
          lwdnt(i)  = dflx(i,nlayers+1)
          lwdntc(i) = dflxc(i,nlayers+1)
          !output up and down surface fluxes for total and clear sky:
          lwupb(i)  = uflx(i,1)
          lwupbc(i) = uflxc(i,1)
          lwdnb(i)  = dflx(i,1)
          lwdnbc(i) = dflxc(i,1)
       endif

       if(present(lwupflx)) then
          !output up and down fluxes:
          do k=kts,nlayers+1
             lwupflx(i,k)  = uflx(i,k)
             lwupflxc(i,k) = uflxc(i,k)
             lwdnflx(i,k)  = dflx(i,k)
             lwdnflxc(i,k) = dflxc(i,k)
          enddo
       endif

       !--- initialize the local radiative heating rate:
       do k = kts, kte
          tten1d(k) = 0.
       enddo

       !--- output heating rate tendency. convert heating rate from K per day to K per second:
       do k = kts, kte
          tten1d(k) = hr(i,k)/86400.
          rthratenlw(i,k) = tten1d(k)/pi2d(i,k)
       enddo

       !--- output the effective radii for cloud water, cloud ice, and snow:
       if(icloud.gt.0 .and. present(rre_cloud) .and. present(rre_ice) .and. present(rre_snow)) then
          do k = kts, kte
             if(cldfra2d(i,k) .gt. 0.) then
                rre_cloud(i,k) = rel(i,k)
                rre_ice(i,k)   = rei(i,k)
                rre_snow(i,k)  = res(i,k)
             endif
          enddo
       endif

    end do i_loop !end longitude loop.

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


 errmsg = 'rrtmg_lw_mpas_timestep_final OK'
 errflg = 0

 call mpas_log_write('--- end subroutine rrtmg_lw_mpas_timestep_final:')


 end subroutine rrtmg_lw_mpas_timestep_final

!=================================================================================================================
 end module rrtmg_lw_mpas
!=================================================================================================================
