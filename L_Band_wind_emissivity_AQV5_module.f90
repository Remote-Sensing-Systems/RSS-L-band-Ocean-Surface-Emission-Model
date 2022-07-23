! RSS L-band wind emissivity
! based on Aquarius V5 release
! written by T. Meissner, RSS, September 29, 2018
! 
! References: 
! 1.	Meissner, T.; F. Wentz and D. Le Vine, Aquarius Salinity Retrieval Algorithm Theoretical Basis Document (ATBD), 
!       End of Mission Version; RSS Technical Report 120117; December 1, 2017; 
!       Available online at ftp://podaac-ftp.jpl.nasa.gov/allData/aquarius/docs/v5/AQ-014-PS-0017_Aquarius_ATBD-EndOfMission.pdf.
!
! 2.    Meissner, T, F. Wentz, and D, Le Vine, 2018, 
!       The Salinity Retrieval Algorithms for the NASA Aquarius Version 5 and SMAP Version 3 Releases, 
!       Remote Sensing 10, 1121, doi:10.3390/rs10071121. 
!
! 3.    Meissner, T. and F. Wentz, The complex dielectric constant of pure and sea water from microwave satellite observations, 
!       IEEE TGRS, 2004, 42(9), 1836 – 1849, doi:10.1109/TGRS.2004.831888.
!
! 4. 	Meissner, T. and F. Wentz, The emissivity of the ocean surface between 6 and 90 GHz 
!       over a large range of wind speeds and Earth incidence angles, 
!       IEEE TGRS, 2012, 50(8), 3004 – 3026, doi: 10.1109/TGRS.2011.2179662.   
!
! 5. 	Meissner, T., F. Wentz, F. and L. Ricciardulli, The emission and scattering of L-band microwave radiation 
!       from rough ocean surfaces and wind speed measurements from Aquarius, 
!       J. Geophys. Res. Oceans, 2014, 119, doi:10.1002/2014JC009837.
!


module L_Band_wind_emissivity_AQV5_module


implicit none
private
save


integer(4), parameter                       :: n_rad=3
real(4), parameter                          :: sst_ref0=20.0, sss_ref0=35.0, freq_aq=1.413, wcasp=11.0, wextrapol=17.0, teff=290.
real(4), parameter, dimension(n_rad)        :: tht_ref = (/29.36, 38.44, 46.29/)
real(4), parameter                          :: kappascale=1.4 ! empirical scale factor for delta form factor

! input tables
character(len=200), parameter               :: emiss_coeff_harm_file = &
                                            '\\Thomas\C\AQUARIUS\surface_roughness\tables\V_may_2013\coeffs\deW_harm_coeffs_V9A_MI.dat'

character(len=200), parameter               :: delta_file = '\\Thomas\C\AQUARIUS\V5.0\DTB_SST_WSPD\TABLES\delta_EW_V5_B.dat'


                      
integer(4), parameter                       :: iu=3

public                                      :: find_demiss_rough_LBAND

contains

subroutine find_demiss_rough_LBAND (wspd, sst, tht,    demiss_rough)
! roughness correction
! Aquarius V5.0
! interpolate/extrapolate EIA

implicit none

real(4), intent(in)                                 ::  wspd ! [m/s]
real(4), intent(in)                                 ::  sst  ! SST [Celsius]
real(4), intent(in)                                 ::  tht  ! Earth Incidence Angle [deg]

 
real(4), intent(out), dimension(2)                  ::  demiss_rough 
! isotropic part of wind induced surface emissivity [*290K]

real(4), dimension(2)                               ::  dew_0, dew_1
integer(4)                                          ::  irad_0,irad_1 
real(4)                                             ::  brief
real(4), dimension(2)                               ::  yy

real(4), dimension(0:n_rad)                         ::  thtfix


    thtfix(1:3) = tht_ref(1:3)
    thtfix(0)   = 0.0    

    if (tht<tht_ref(1)) then
        irad_0 = 0
        irad_1 = 1
        call find_demiss_rough_wspd_AQV5 (1, wspd, sst,    dew_1)   
        dew_0(1) = (dew_1(1)+dew_1(2))/2.0 ! (V+H)/2 at nadir
        dew_0(2) = (dew_1(1)+dew_1(2))/2.0 ! (V+H)/2 at nadir
    else if (tht>=tht_ref(1) .and. tht<tht_ref(2)) then
        irad_0 = 1
        irad_1 = 2
        call find_demiss_rough_wspd_AQV5 (1, wspd, sst,    dew_0)          
        call find_demiss_rough_wspd_AQV5 (2, wspd, sst,    dew_1)          
    else
        irad_0 = 2
        irad_1 = 3
        call find_demiss_rough_wspd_AQV5 (2, wspd, sst,    dew_0)          
        call find_demiss_rough_wspd_AQV5 (3, wspd, sst,    dew_1)                
    endif
    
     
    brief = (tht-thtfix(irad_0))/(thtfix(irad_1)-thtfix(irad_0)) 
    ! linear interpolation/extrapolation to/form Aquarius EIA
 
    yy = dew_0*(1.0-brief) + dew_1*brief   
    demiss_rough = yy                                                              

return
end subroutine find_demiss_rough_LBAND




! Aquarius V5.0
subroutine fd_TM_emiss_harmonics(irad,wspd,		aharm, daharm) 
implicit none

	integer(4), intent(in)					            :: irad
	real(4), intent(in)						            :: wspd
	
    real(4), dimension(0:2,2), intent(out)	            ::  aharm(0:2,2)  !1=V, 2=H
    real(4), dimension(0:2,2), intent(out), optional   	:: daharm(0:2,2)  !1=V, 2=H
 
    integer(4), parameter                               :: n_rad=3, npoly=5


	real(4), dimension(2)				                ::  A0,  A1,  A2  !1=V, 2=H
	real(4), dimension(2)				                :: dA0, dA1, dA2  !1=V, 2=H
	
	
	real(4)									            :: ww
	integer(4)								            :: ipol, iharm
	real(4)									            :: fval, dval

	real(4)         						            :: w0, w1, w2 ! linear extrapolation/cutoff points
	
    
	integer(4), save                                    :: istart=1	
    real(8), dimension(0:2,2,n_rad,npoly), save         :: acoef      ! harmonic coefficients for radiometer wind direction signal
    real(8), dimension(0:2,2,n_rad), save               :: wspd_max_a ! high wind speed for radiometer wind speed signal


	
	if (istart==1) then
	    istart=0
	    open(unit=3,file=emiss_coeff_harm_file, form='binary', action='read', status='old')
	    read(3) acoef
	    read(3) wspd_max_a
	    ! overwrite
	    wspd_max_a = wextrapol
	    close(3)
	endif 
	

	A0=0.0
	A1=0.0
	A2=0.0

	do ipol=1,2 

	! A0
	iharm=0
	w0 = wspd_max_a(iharm,ipol,irad)
	ww = wspd
	if (wspd >= w0) ww=w0 ! extrapolation at w0 
	fval = &
	ww*acoef(iharm,ipol,irad,1)  +  (ww**2)*acoef(iharm,ipol,irad,2) +       (ww**3)*acoef(iharm,ipol,irad,3) +       (ww**4)*acoef(iharm,ipol,irad,4) +       (ww**5)*acoef(iharm,ipol,irad,5) 		   
	dval = &
	   acoef(iharm,ipol,irad,1)  + (2.0*ww)*acoef(iharm,ipol,irad,2) + (3.0*(ww**2))*acoef(iharm,ipol,irad,3) + (4.0*(ww**3))*acoef(iharm,ipol,irad,4) + (5.0*(ww**4))*acoef(iharm,ipol,irad,5) 
	
	if (wspd<=w0) then
		A0(ipol) = fval
	else
		A0(ipol) = fval + dval*(wspd-w0)
	endif
    
    dA0(ipol) = dval
   
    
	! A1
	iharm=1
	w1 = wspd_max_a(iharm,ipol,irad)
	ww = wspd
	if (wspd >= w1) ww=w1 ! cutoff at w1 
	fval = ww*acoef(iharm,ipol,irad,1) + (ww**2) *acoef(iharm,ipol,irad,2) +       (ww**3)*acoef(iharm,ipol,irad,3)      +       (ww**4)*acoef(iharm,ipol,irad,4) +       (ww**5)*acoef(iharm,ipol,irad,5)   
	dval =    acoef(iharm,ipol,irad,1) + (2.0*ww)*acoef(iharm,ipol,irad,2) + (3.0*(ww**2))*acoef(iharm,ipol,irad,3)      + (4.0*(ww**3))*acoef(iharm,ipol,irad,4) + (5.0*(ww**4))*acoef(iharm,ipol,irad,5)  
	
	 A1(ipol)  = fval
	dA1(ipol)  = dval

	! A2
	iharm=2
	w2 = wspd_max_a(iharm,ipol,irad)
	ww = wspd
	if (wspd >= w2) ww=w2 ! cutoff at w2 
	fval = ww*acoef(iharm,ipol,irad,1) + (ww**2)      *acoef(iharm,ipol,irad,2) +        (ww**3)*acoef(iharm,ipol,irad,3) +       (ww**4)*acoef(iharm,ipol,irad,4) +       (ww**5)*acoef(iharm,ipol,irad,5)  
	dval =    acoef(iharm,ipol,irad,1) + (2.0*ww)     *acoef(iharm,ipol,irad,2) +  (3.0*(ww**2))*acoef(iharm,ipol,irad,3) + (4.0*(ww**3))*acoef(iharm,ipol,irad,4) + (5.0*(ww**4))*acoef(iharm,ipol,irad,5) 

	A2(ipol)  = fval
	dA2(ipol) = dval

	enddo !ipol


	do ipol=1,2
		aharm(0,ipol) = A0(ipol)
		aharm(1,ipol) = A1(ipol)
		aharm(2,ipol) = A2(ipol)
	enddo

	if (present(daharm)) then
	do ipol=1,2
		daharm(0,ipol) = dA0(ipol)
		daharm(1,ipol) = dA1(ipol)
		daharm(2,ipol) = dA2(ipol)
	enddo
	endif


return
end subroutine fd_TM_emiss_harmonics



subroutine find_demiss_rough_AQV5 (irad, wspd, phir, sst,    demiss_rough)
! full roughness correction
implicit none

integer(4), intent(in)                              ::  irad
real(4), intent(in)                                 ::  wspd ! [m/s]
real(4), intent(in)                                 ::  phir ! [deg] relative wind direction
real(4), intent(in)                                 ::  sst  ! Celsius

 
real(4), intent(out), dimension(2)                  ::  demiss_rough ! [*290K]

real(4)                                             ::  xsst, xwspd, xtht

real(4), dimension(2)                               ::  dew_1, dew_2, xem0, yem0

real(4), dimension(0:2,2)           	            ::  aharm  !1=V, 2=H
real(4), dimension(0:2,2)                       	:: daharm  !1=V, 2=H

real(4), dimension(2)                               ::  delta
integer(4)                                          ::  ipol

integer(4)                                          ::  ksst1, ksst2
real(4)                                             ::  x1, x2, brief, y1, y2

integer(4), parameter                               ::  msst=35, n_rad=3
integer(4), save                                    ::  istart=1



real(4), dimension(2,msst,n_rad), save              ::  dtab
real(4), parameter                                  ::  sst_step=1.0, sst0=0.0, sstmax=30.0, wmax=wcasp

if (istart==1) then
    istart=0
    open(unit=3,form='binary',file=delta_file,action='read',status='old')
    read(3) dtab
    close(3) 
endif

call fd_TM_emiss_harmonics(irad,wspd,		    aharm, daharm) 
if (phir > -900.) then
dew_1(1:2) = aharm(0,1:2) + aharm(1,1:2)*cosd(phir) + aharm(2,1:2)*cosd(2.0*phir) 	! *290K
else
dew_1(1:2) = aharm(0,1:2)
endif

xtht=tht_ref(irad)
call fdem0_meissner_wentz(freq_aq,xtht,sst,    sss_ref0, xem0) 
call fdem0_meissner_wentz(freq_aq,xtht,sst_ref0,sss_ref0, yem0) 

! sst adjustment
xsst=sst
if (xsst<sst0+sst_step/2)   xsst=sst0+sst_step/2
if (xsst>sstmax) xsst=sstmax

xwspd=wspd
if (xwspd<0.0)  xwspd=0.0
if (xwspd>wmax) xwspd=wmax

call fd_TM_emiss_harmonics(irad,xwspd,		    aharm, daharm) 
if (phir > -900.)  then
dew_2(1:2) = aharm(0,1:2) + aharm(1,1:2)*cosd(phir) + aharm(2,1:2)*cosd(2.0*phir) 	! *290K
else
dew_2(1:2) = aharm(0,1:2)
endif

! *290K dew_2 = dew_1 below wmax. above wmax it is kept constant


ksst1 = floor((xsst - (sst0+sst_step/2))/sst_step) + 1
if (ksst1< 1)           ksst1 = 1
if (ksst1> msst-1)      ksst1 = msst-1
ksst2 = ksst1 + 1

do ipol=1,2    
    x1 = sst0 + (ksst1-1)*sst_step + sst_step/2
    x2 = x1   + sst_step
    brief = (xsst-x1)/sst_step 
    y1 = dtab(ipol,ksst1,irad)
    y2 = dtab(ipol,ksst2,irad)
    delta(ipol) = y1*(1.0-brief) + y2*brief    
enddo ! ipol

! total roughness correction
do ipol=1,2
    demiss_rough(ipol) = dew_1(ipol)*(xem0(ipol)/yem0(ipol)) + kappascale*delta(ipol)*dew_2(ipol)
enddo

return
end subroutine find_demiss_rough_AQV5



subroutine find_demiss_rough_wspd_AQV5 (irad, wspd, sst,    demiss_rough)
! isotropic part of roughness correction
implicit none

integer(4), intent(in)                              ::  irad
real(4), intent(in)                                 ::  wspd ! [m/s]
real(4), intent(in)                                 ::  sst  ! Celsius

 
real(4), intent(out), dimension(2)                  ::  demiss_rough ! [*290K]

real(4)                                             ::  xsst, xwspd, xtht

real(4), dimension(2)                               ::  dew_1, dew_2, xem0, yem0

real(4), dimension(0:2,2)           	            ::  aharm  !1=V, 2=H
real(4), dimension(0:2,2)                       	:: daharm  !1=V, 2=H

real(4), dimension(2)                               ::  delta
integer(4)                                          ::  ipol

integer(4)                                          ::  ksst1, ksst2
real(4)                                             ::  x1, x2, brief, y1, y2

integer(4), parameter                               ::  msst=35, n_rad=3
integer(4), save                                    ::  istart=1

real(4), dimension(2,msst,n_rad), save              ::  dtab
real(4), parameter                                  ::  sst_step=1.0, sst0=0.0, sstmax=30.0, wmax=wcasp

if (istart==1) then
    istart=0
    open(unit=3,form='binary',file=delta_file,action='read',status='old')
    read(3) dtab
    close(3) 
endif

call fd_TM_emiss_harmonics(irad,wspd,		    aharm, daharm) 
dew_1(1:2) = aharm(0,1:2) ! *290K

xtht=tht_ref(irad)
call fdem0_meissner_wentz(freq_aq,xtht,sst,    sss_ref0, xem0) 
call fdem0_meissner_wentz(freq_aq,xtht,sst_ref0,sss_ref0, yem0) 

! sst adjustment
xsst=sst
if (xsst<sst0+sst_step/2)   xsst=sst0+sst_step/2
if (xsst>sstmax) xsst=sstmax

xwspd=wspd
if (xwspd<0.0)  xwspd=0.0
if (xwspd>wmax) xwspd=wmax

call fd_TM_emiss_harmonics(irad,xwspd,		    aharm, daharm) 
dew_2(1:2) = aharm(0,1:2) 	! *290K
! *290K dew_2 = dew_1 below wmax. above wmax it is kept constant


ksst1 = floor((xsst - (sst0+sst_step/2))/sst_step) + 1
if (ksst1< 1)           ksst1 = 1
if (ksst1> msst-1)      ksst1 = msst-1
ksst2 = ksst1 + 1

do ipol=1,2    
    x1 = sst0 + (ksst1-1)*sst_step + sst_step/2
    x2 = x1   + sst_step
    brief = (xsst-x1)/sst_step 
    y1 = dtab(ipol,ksst1,irad)
    y2 = dtab(ipol,ksst2,irad)
    delta(ipol) = y1*(1.0-brief) + y2*brief    
enddo ! ipol

! total roughness correction
do ipol=1,2
    demiss_rough(ipol) = dew_1(ipol)*(xem0(ipol)/yem0(ipol)) + kappascale*delta(ipol)*dew_2(ipol)
enddo

return
end subroutine find_demiss_rough_wspd_AQV5


subroutine fdem0_meissner_wentz(freq,tht,sst,salinity, em0) 
!    input:
!    name   parameter		unit		range
!      
!    freq     frequency		[GHz]		>0
!    tht      EIA           [deg]       [0, 90[
!    sst      SST			[C]			-25 c to 40 c for pure water
!										-2  c to 34 c for saline water
!    salinity salinity		[ppt]		0 to 40
!
!    output:
!    EM0     specular emissivity        [0,1]
!	         2-dimesnional vector, 1=v-pol, 2=h-pol 		 

implicit none

	real(4), intent(in)						  :: freq,tht,sst,salinity
	real(4), dimension(2), intent(out)		  :: em0 
	
	real(4), parameter						  :: f0=17.97510
 

	real(4)									  :: costht,sinsqtht
	real(4)									  :: e0s,e1s,e2s,n1s,n2s,sig
 
    complex(4)							      :: permit,esqrt,rh,rv
    complex(4), parameter					  :: j=(0.,1.)
	
 
	  call dielectric_meissner_wentz(sst,salinity,  e0s,e1s,e2s,n1s,n2s,sig)

	  costht=cosd(tht)
	  sinsqtht=1.-costht*costht


!     debye law (2 relaxation wavelengths)
      permit = (e0s - e1s)/(1.0 - j*(freq/n1s)) + (e1s - e2s)/(1.0 - j*(freq/n2s)) + e2s + j*sig*f0/freq
      permit = conjg(permit)
	
      esqrt=csqrt(permit-sinsqtht)
      rh=(costht-esqrt)/(costht+esqrt)
      rv=(permit*costht-esqrt)/(permit*costht+esqrt)
      em0(1)  =1.-rv*conjg(rv)
      em0(2)  =1.-rh*conjg(rh)
 
return
end subroutine fdem0_meissner_wentz


subroutine fdpermit_meissner_wentz(freq,sst,salinity, permit) 
!    input:
!    name   parameter		unit		range
!      
!    freq     frequency		[GHz]		>0
!    tht      EIA           [deg]       [0, 90[
!    sst      SST			[C]			-25 c to 40 c for pure water
!										-2  c to 34 c for saline water
!    salinity salinity		[ppt]		0 to 40
!
!    output:
!    permit  permittivity
!   (negative imaginary part to be consistent with wentz1 convention)

implicit none

	real(4), intent(in)						  :: freq,sst,salinity
	complex(4), intent(out)		              :: permit 
	
	real(4), parameter						  :: f0=17.97510
 

	real(4)									  :: e0s,e1s,e2s,n1s,n2s,sig
 
    complex(4), parameter					  :: j=(0.,1.)
	
 
	  call dielectric_meissner_wentz(sst,salinity,  e0s,e1s,e2s,n1s,n2s,sig)

!     debye law (2 relaxation wavelengths)
      permit = (e0s - e1s)/(1.0 - j*(freq/n1s)) + (e1s - e2s)/(1.0 - j*(freq/n2s)) + e2s + j*sig*f0/freq
      permit = conjg(permit)
	
return
end subroutine fdpermit_meissner_wentz
 

subroutine dielectric_meissner_wentz(sst_in,s,   e0s,e1s,e2s,n1s,n2s,sig)
!
!     complex dielectric constant: eps
!     [MW 2004, MW 2012].
!     References:
!     [MW 2004]:   T. Meissner and F. J. Wentz, 
!              "The complex dielectric constant of pure and sea water from microwave satellite observations," 
!              IEEE Trans. Geosci. Remote Sens., vol. 42, no.9, pp 1836 – 1849, 2004. 
!
!     [MW 2012]:   T. Meissner and F. J. Wentz, 
!              "The Emissivity of the Ocean Surface between 6 – 90 GHz over a Large Range of Wind Speeds and Earth Incidence Angles,"
!              IEEE Trans. Geosci. Remote Sens., vol. 50, no.8, pp 3004 - 3026, 2012.
!     
!     Changes from [MW 2012]:
!     1. Typo (sign) in the printed version of coefficient d3 in Table 7. Its value should be -0.35594E-06.
!     2. Changed SST behavior of coefficient b2 from:
!     b2 = 1.0 + s*(z(10) + z(11)*sst) to
!     b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30)) 
!
!!
!     input:
!     name   parameter  unit  range
!     sst      sst        [c]   -25 c to 40 c for pure water
!                               -2  c to 34 c for saline water
!     s      salinity   [ppt]  0 to 40
!
!     output:
!     Debye pparameters: e0s,e1s,e2s,n1s,n2s,sig
!     The permittivity can be calculated in the subroutine: fdpermit_meissner_wentz 
implicit none


      real(4), intent(in)  :: sst_in,s
      real(4), intent(out) :: e0s,e1s,e2s,n1s,n2s,sig
 
      real(4), dimension(11), parameter :: &
      x=(/ 5.7230e+00, 2.2379e-02, -7.1237e-04, 5.0478e+00, -7.0315e-02, 6.0059e-04, 3.6143e+00, &
           2.8841e-02, 1.3652e-01,  1.4825e-03, 2.4166e-04 /)
 
 
      real(4), dimension(13), parameter :: &
      z=(/ -3.56417e-03,  4.74868e-06,  1.15574e-05,  2.39357e-03, -3.13530e-05, &
            2.52477e-07, -6.28908e-03,  1.76032e-04, -9.22144e-05, -1.99723e-02, &
            1.81176e-04, -2.04265e-03,  1.57883e-04  /)  ! 2004

      real(4), dimension(3), parameter :: a0coef=(/ -0.33330E-02,  4.74868e-06,  0.0e+00/)
      real(4), dimension(5), parameter :: b1coef=(/0.23232E-02, -0.79208E-04, 0.36764E-05, -0.35594E-06, 0.89795E-08/)
 
      real(4) :: e0,e1,e2,n1,n2
      real(4) :: a0,a1,a2,b1,b2
      real(4) :: sig35,r15,rtr15,alpha0,alpha1

	  real(4) :: sst,sst2,sst3,sst4,s2

	  sst=sst_in
	  if(sst.lt.-30.16) sst=-30.16  !protects against n1 and n2 going zero for very cold water

	  sst2=sst*sst
	  sst3=sst2*sst
	  sst4=sst3*sst

	  s2=s*s
 
!     pure water
      e0    = (3.70886e4 - 8.2168e1*sst)/(4.21854e2 + sst) ! stogryn et al.
      e1    = x(1) + x(2)*sst + x(3)*sst2
      n1    = (45.00 + sst)/(x(4) + x(5)*sst + x(6)*sst2)
      e2    = x(7) + x(8)*sst
      n2    = (45.00 + sst)/(x(9) + x(10)*sst + x(11)*sst2)
 
!     saline water
!     conductivity [s/m] taken from stogryn et al.
      sig35 = 2.903602 + 8.60700e-2*sst + 4.738817e-4*sst2 - 2.9910e-6*sst3 + 4.3047e-9*sst4
      r15   = s*(37.5109+5.45216*s+1.4409e-2*s2)/(1004.75+182.283*s+s2)
 
      alpha0 = (6.9431+3.2841*s-9.9486e-2*s2)/(84.850+69.024*s+s2)
      alpha1 = 49.843 - 0.2276*s + 0.198e-2*s2
      rtr15 = 1.0 + (sst-15.0)*alpha0/(alpha1+sst)
 
      sig = sig35*r15*rtr15
 
!     permittivity
      a0 = exp(a0coef(1)*s + a0coef(2)*s2 + a0coef(3)*s*sst)  
      e0s = a0*e0

	  if(sst.le.30) then
		b1 = 1.0 + s*(b1coef(1) + b1coef(2)*sst + b1coef(3)*sst2 + b1coef(4)*sst3 + b1coef(5)*sst4)
	  else
		b1 = 1.0 + s*(9.1873715e-04 + 1.5012396e-04*(sst-30))
	  endif
      n1s = n1*b1


      a1  = exp(z(7)*s + z(8)*s2 + z(9)*s*sst)
      e1s = e1*a1

!     b2 = 1.0 + s*(z(10) + z(11)*sst)
      b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30))
      n2s = n2*b2


      a2 = 1.0  + s*(z(12) + z(13)*sst)
      e2s = e2*a2
 
return
end subroutine  dielectric_meissner_wentz




end module L_Band_wind_emissivity_AQV5_module



