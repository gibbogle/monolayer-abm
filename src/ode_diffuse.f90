!----------------------------------------------------------------------------------
! Note: The value of spcrad was first determined by writing out the value computed in rkc.
! Later it was just determined by trial, then made into a run parameter.
!----------------------------------------------------------------------------------
double precision function spcrad(neqn,t,y)
!DEC$ ATTRIBUTES DLLEXPORT :: spcrad
use global
integer :: neqn
double precision :: t, y(neqn)
spcrad = spcrad_value
end function

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
module ode_diffuse

use chemokine
use rkc_90

implicit none

!integer, parameter :: MAX_VARS = 2*max_nlist
!integer, parameter :: RKF45_SOLVE = 1
!integer, parameter :: RKSUITE_SOLVE = 2
!integer, parameter :: RKC_SOLVE = 3
!integer, parameter :: IRKC_SOLVE = 4
!logical, parameter :: EXPLICIT_INTRA = .false.
!real(REAL_KIND), allocatable :: allstate(:,:)
!real(REAL_KIND), allocatable :: allstatep(:,:)
!real(REAL_KIND), allocatable :: work_rkc(:)

integer :: ivdbug

real(REAL_KIND) :: work_rkc(8+5*2*MAX_CHEMO)
logical :: chemo_active(2*MAX_CHEMO)    ! flags necessity to solve for the constituent
contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine CheckDrugConcs
integer :: ndrugs_present, drug_present(3*MAX_DRUGTYPES), drug_number(3*MAX_DRUGTYPES)
integer :: idrug, iparent, im, kcell, ichemo, i
type(cell_type), pointer :: cp

ndrugs_present = 0
drug_present = 0
do idrug = 1,ndrugs_used
	iparent = TRACER + 1 + 3*(idrug-1)
	if (chemo(iparent)%present) then		! simulation with this drug has started
	    do im = 0,2
	        ichemo = iparent + im
	        ndrugs_present = ndrugs_present + 1
	        drug_present(ndrugs_present) = ichemo
	        drug_number(ndrugs_present) = idrug
	    enddo
	endif
enddo

do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
    cp => cell_list(kcell)
	do i = 1,ndrugs_present
	    ichemo = drug_present(i)
	    idrug = drug_number(i)
	    if (cp%conc(ichemo) > Cthreshold) drug_gt_cthreshold(idrug) = .true.
!	    if (cp%Cex(ichemo) > Cthreshold) drug_gt_cthreshold(idrug) = .true.
	enddo
enddo
    
end subroutine

!----------------------------------------------------------------------------------
! For constituent ichemo, the extracellular concentration is:
! Cex = chemo(ichemo)%conc
! In the case of oxygen this is determined from: 
!   depth, Kdiff, chemo(OXYGEN)%flux, chemo(OXYGEN)%bdry_conc
! where:
!   depth = depth of medium in the well
!   %flux = total O2 uptake rate by cells
!   %bdry_conc = specified O2 concentration at the medium-air boundary
! For other constituents %conc is the average concentration in the medium,
! i.e. the medium is considered to be fully mixed.  In this case:
!   dC/dt = -flux/V
! where:
!   V = total volume of medium
!   C = medium concentration %conc
!
! neqn = 2*ncvars = 2*number of constituents present
! ic > ncvars implies a medium concentration
! chemo_active(ic) = false means we do not solve for it (only medium variables)
!----------------------------------------------------------------------------------
subroutine f_rkc(neqn,t,y,dydt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, y(neqn), dydt(neqn)
integer :: ic, ichemo, idrug, im, ict, Ng, ncvars
real(REAL_KIND) :: dCsum, dCdiff, dCreact, vol_cm3, val, Cin(MAX_CHEMO), Cmedium(MAX_CHEMO), Cex
real(REAL_KIND) :: decay_rate, C, membrane_kin, membrane_kout, membrane_flux, area_factor, n_O2(0:2)
logical :: metabolised(MAX_CELLTYPES,0:2)
real(REAL_KIND) :: metab, dMdt, KmetC, vcell_actual, Kd(0:2), dC, C0
type(drug_type), pointer :: dp
real(REAL_KIND) :: average_volume = 1.2
logical :: use_average_volume = .true.

if (use_average_volume) then
    vol_cm3 = Vcell_cm3*average_volume	  ! not accounting for cell volume change
    area_factor = (average_volume)**(2./3.)
endif
!Vcell_actual = Vcell_cm3*cell_list(kcell)%volume
!vol_cm3 = Vcell_actual	            ! accounting for cell volume change
!Cin = cell_list(kcell)%conc
!ict = cell_list(kcell)%celltype
ncvars = neqn/2
do ic = 1,ncvars
	ichemo = chemomap(ic)
    Cin(ichemo) = y(ic)
    Cmedium(ichemo) = y(ncvars+ic)
enddo
!write(*,'(a,8e12.3)') 'Cin: ',Cin(1:ncvars)
!write(*,'(a,8e12.3)') 'Cmedium: ',Cmedium(1:ncvars)
ict = icase

do ic = 1,neqn
    if (ic <= ncvars) then
    	ichemo = chemomap(ic)
    else
    	ichemo = chemomap(ic-ncvars)
    endif
    if (ichemo == GLUCOSE) then
	    Ng = chemo(GLUCOSE)%Hill_N
    endif
    if (ichemo > TRACER) then
        idrug = (ichemo - TRACER - 1)/3 + 1
        im = ichemo - TRACER - 1 - 3*(idrug-1)		! 0 = drug, 1 = metab1, 2 = metab2
        dp => drug(idrug)
        metabolised(:,:) = (dp%Kmet0(:,:) > 0)
        if (idrug > 0) then
            n_O2(:) = dp%n_O2(ict,:)
        endif
    endif

    decay_rate = chemo(ichemo)%decay_rate
    membrane_kin = chemo(ichemo)%membrane_diff_in
    membrane_kout = chemo(ichemo)%membrane_diff_out

    Cex = Cmedium(ichemo)
	C = Cin(ichemo)     ! = y(ic)
	membrane_flux = area_factor*(membrane_kin*Cex - membrane_kout*C)
	dydt(ic) = 0
	if (ic <= ncvars .and. chemo_active(ic)) then      ! cell variable
	    dCreact = 0
	    if (ichemo == OXYGEN) then
		    metab = O2_metab(C)
		    dCreact = (-metab*chemo(ichemo)%max_cell_rate + membrane_flux)/vol_cm3	! convert mass rate (mol/s) to concentration rate (mM/s)
!		    write(*,'(a,6e12.3)') 'O2: ',C,metab,chemo(ichemo)%max_cell_rate,membrane_flux,vol_cm3,dCreact
	    elseif (ichemo == GLUCOSE) then
		    metab = C**Ng/(chemo(ichemo)%MM_C0**Ng + C**Ng)
		    dCreact = (-metab*chemo(ichemo)%max_cell_rate + membrane_flux)/vol_cm3	! convert mass rate (mol/s) to concentration rate (mM/s)
!		    write(*,'(a,6e11.3)') 'glucose: ',C,metab,chemo(ichemo)%max_cell_rate,membrane_flux,vol_cm3,dCreact
	    elseif (im == 0) then
	        if (metabolised(ict,0) .and. C > 0) then
			    KmetC = dp%Kmet0(ict,0)*C
			    if (dp%Vmax(ict,0) > 0) then
				    KmetC = KmetC + dp%Vmax(ict,0)*C/(dp%Km(ict,0) + C)
			    endif
			    dCreact = -(1 - dp%C2(ict,0) + dp%C2(ict,0)*dp%KO2(ict,0)**n_O2(0)/(dp%KO2(ict,0)**n_O2(0) + Cin(OXYGEN)**n_O2(0)))*KmetC
		    endif
		    dCreact = dCreact + membrane_flux/vol_cm3
	    elseif (im == 1) then	! ichemo-1 is the PARENT drug
		    if (metabolised(ict,0) .and. Cin(ichemo-1) > 0) then
			    dCreact = (1 - dp%C2(ict,0) + dp%C2(ict,0)*dp%KO2(ict,0)**n_O2(0)/(dp%KO2(ict,0)**n_O2(0) + Cin(OXYGEN)**n_O2(0)))*dp%Kmet0(ict,0)*Cin(ichemo-1)
		    endif
		    if (metabolised(ict,1) .and. C > 0) then
			    dCreact = dCreact - (1 - dp%C2(ict,1) + dp%C2(ict,1)*dp%KO2(ict,1)**n_O2(1)/(dp%KO2(ict,1)**n_O2(1) + Cin(OXYGEN)**n_O2(1)))*dp%Kmet0(ict,1)*C
		    endif
		    dCreact = dCreact + membrane_flux/vol_cm3
	    elseif (im == 2) then	! ichemo-1 is the METAB1
		    if (metabolised(ict,1) .and. Cin(ichemo-1) > 0) then
			    dCreact = (1 - dp%C2(ict,1) + dp%C2(ict,1)*dp%KO2(ict,1)**n_O2(1)/(dp%KO2(ict,1)**n_O2(1) + Cin(OXYGEN)**n_O2(1)))*dp%Kmet0(ict,1)*Cin(ichemo-1)
		    endif
		    if (metabolised(ict,2) .and. C > 0) then
			    dCreact = dCreact - (1 - dp%C2(ict,2) + dp%C2(ict,2)*dp%KO2(ict,2)**n_O2(2)/(dp%KO2(ict,2)**n_O2(2) + Cin(OXYGEN)**n_O2(2)))*dp%Kmet0(ict,2)*C
		    endif
		    dCreact = dCreact + membrane_flux/vol_cm3
	    endif
        dydt(ic) = dCreact - C*decay_rate
    else    ! medium variable 
        if (chemo_active(ic)) then
            dydt(ic) = -Ncells*membrane_flux/total_volume - C*decay_rate
        else
            dydt(ic) = 0
        endif
    endif
	if (isnan(dydt(ic))) then
		write(nflog,*) 'f_rkc: dydt isnan: ',ic,ichemo,dydt(ic)
		write(*,*) 'f_rkc: dydt isnan: ',ic,ichemo,dydt(ic)
		stop
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine Solver(it,tstart,dt,nc,ok)
integer :: it, nc
real(REAL_KIND) :: tstart, dt
logical :: ok
integer :: ichemo, ic, k, ict, ncvars, neqn, kcell
real(REAL_KIND) :: t, tend
real(REAL_KIND) :: C(2*MAX_CHEMO)
real(REAL_KIND) :: timer1, timer2
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1)
type(rkc_comm) :: comm_rkc(1)
logical :: solve_O2 = .true.

k = 0
do ic = 1,nchemo
	ichemo = chemomap(ic)
	k = k + 1
    chemo_active(k) = .not.chemo(ichemo)%constant
    if (ichemo == OXYGEN) then
        if (.not.solve_O2) then
            ! Suppress solving for cell oxygen
!	        Caverage(OXYGEN) = getCin(OXYGEN,Caverage(MAX_CHEMO+OXYGEN))
            chemo_active(k) = .false.
        endif
    endif
	C(k) = Caverage(ichemo)                ! average cell concentration
enddo
ncvars = k
! Note: ncvars = nchemo
do ic = 1,nchemo
	ichemo = chemomap(ic)
	k = k + 1
	C(k) = Caverage(MAX_CHEMO + ichemo)      ! average medium concentration
    chemo_active(k) = .not.chemo(ichemo)%constant
    if (ichemo == OXYGEN) then
        ! Suppress solving for medium oxygen
        chemo_active(k) = .false.
    endif
enddo
neqn = k
! Note: neqn = 2*ncvars

!write(*,*) 'solver: nchemo,neqn: ',nchemo,neqn
!write(*,'(10f7.3)') C(1:neqn)

ict = 1 ! for now just a single cell type

info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 1d-2
atol = rtol

idid = 0
t = tstart
tend = t + dt
call rkc(comm_rkc(1),neqn,f_rkc,C,t,tend,rtol,atol,info,work_rkc,idid,ict)
if (idid /= 1) then
	write(logmsg,*) 'Solver: Failed at t = ',t,' with idid = ',idid
	call logger(logmsg)
	ok = .false.
	return
endif

! This determines average cell concentrations, assumed the same for all cells
! Now put the concentrations into the cells

k = 0
do ic = 1,nchemo
    ichemo = chemomap(ic)
    k = k + 1
    Caverage(ichemo) = C(k)
    do kcell = 1,nlist
        if (cell_list(kcell)%state == DEAD) cycle
        cell_list(kcell)%conc(ichemo) = Caverage(ichemo)
    enddo
enddo
k = ncvars
do ic = 1,nchemo
    ichemo = chemomap(ic)
    k = k + 1
    Caverage(MAX_CHEMO + ichemo) = C(k)
enddo
! Note: medium oxygen is unchanged


end subroutine

!----------------------------------------------------------------------------------
! Note: This computes a rate of change of concentration! mM/s
! Currently only for O2!!!
! There are two options: use_Cex_Cin = true/false
!
! use_Cex_Cin = true
! ------------------
! The idea is that the speed of the intracellular reactions, compared with other
! processes, is so fast that effectively the intracellular concentration is always
! in equilibrium with the extracellular value.  This means that the rate of consumption
! in the cell matches the rate of transport across the cell membrane: both these rates 
! depend on Cin, therefore we can solve for Cin given Cex then deduce uptake rate
!
! use_Cex_Cin = false
! -------------------
! In this case we just use Cin = Cex to calculate the consumption rate - no
! dependence on chemo(OXYGEN)%membrane_diff
!----------------------------------------------------------------------------------
real(REAL_KIND) function UptakeRate(ichemo,Cex)
integer :: ichemo
real(REAL_KIND) :: Cex
real(REAL_KIND) :: vol, K1, Cin, flux
integer :: n, i

if (ichemo == OXYGEN) then
!    vol = Vsite_cm3
!    vol = Vsite_cm3 - Vextra_cm3	! this was used in the RKC solution
    vol = Vextra_cm3	! the current extracellular volume should be used I think !!!!!!!!!!!!!!!
	if (use_Cex_Cin) then
		Cin = getCin(ichemo,Cex)
!		flux = chemo(ichemo)%membrane_diff*(Cex - Cin)
		flux = (chemo(ichemo)%membrane_diff_in*Cex - chemo(ichemo)%membrane_diff_out*Cin)
	else	! 
		flux = O2_metab(Cex)*chemo(ichemo)%max_cell_rate
	endif
	if (dbug) write(nfout,'(a,2e12.4)') 'Cex, flux: ',Cex,flux
	UptakeRate = flux/vol	! concentration rate (mM/s)
else
	write(logmsg,*) 'ERROR: UptakeRate: currently only for OXYGEN'
	call logger(logmsg)
	stop
endif
end function

!----------------------------------------------------------------------------------
! Computes intracellular O2 concentration as a function of the extracellular level C,
! assuming equilibrium, i.e. rate of consumption = rate of membrane transport.
! Note that the cell's O2 uptake rate is taken to be independent of any other factors,
! e.g. independent of cell size.
! NOTE: Currently only for OXYGEN - OK because membrane_diff_in = membrane_diff_out
!----------------------------------------------------------------------------------
!real(REAL_KIND) function getCinO2(C)
real(REAL_KIND) function getCin(ichemo,C)
integer :: ichemo
real(REAL_KIND) :: C
real(REAL_KIND) :: K1, K2, K2K1, C0, a, b, cc, D, r(3), Cin
integer :: i, n

if (ichemo /= OXYGEN) then
	write(logmsg,*) 'ERROR: getCin: currently only for OXYGEN'
	call logger(logmsg)
	stop
endif
!ichemo = OXYGEN
!K1 = chemo(OXYGEN)%membrane_diff*(Vsite_cm3 - Vextra_cm3)
K1 = chemo(ichemo)%membrane_diff_in
K2 = chemo(ichemo)%max_cell_rate
K2K1 = K2/K1
C0 = chemo(ichemo)%MM_C0
if (chemo(ichemo)%Hill_N == 2) then
	a = K2K1 - C
	b = C0*C0
	cc = -b*C
	call cubic_roots(a,b,cc,r,n)
	if (n == 1) then
		Cin = r(1)
	else
		n = 0
		do i = 1,3
			if (r(i) > 0) then
				n = n+1
				Cin = r(i)
			endif
		enddo
		if (n > 1) then
			write(nflog,*) 'getCin: two roots > 0: ',r
			stop
		endif
	endif
elseif (chemo(ichemo)%Hill_N == 1) then
	b = K2K1 + C0 - C
	cc = -C0*C
	D = sqrt(b*b - 4*cc)
	Cin = (D - b)/2
endif
getCin = Cin
end function


!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!subroutine UpdateCbnd(dt)
!real(REAL_KIND) :: dt
!
!call UpdateCbnd_1D
!call UpdateChemomap
!end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine UpdateCbnd_1D
integer :: kpar = 0
integer :: i, ic, ichemo
real(REAL_KIND) :: tnow, alpha_Cbnd = 0.3
real(REAL_KIND) :: t_buffer = 3600	! one hour delay before applying smoothing to Cbnd
integer :: ndrugs_present, drug_present(3*MAX_DRUGTYPES), drug_number(3*MAX_DRUGTYPES)
integer :: idrug, iparent, im
logical :: present

tnow = istep*DELTA_T
ndrugs_present = 0
drug_present = 0
do idrug = 1,ndrugs_used
	iparent = TRACER + 1 + 3*(idrug-1)
	if (chemo(iparent)%present) then		! simulation with this drug has started
	    do im = 0,2
	        ichemo = iparent + im
	        ndrugs_present = ndrugs_present + 1
	        drug_present(ndrugs_present) = ichemo
	        drug_number(ndrugs_present) = idrug
	    enddo
	endif
enddo

!if ((tnow - t_lastmediumchange) > t_buffer) then
!	chemo(OXYGEN)%medium_Cbnd = alpha_Cbnd*chemo(OXYGEN)%medium_Cbnd + (1 - alpha_Cbnd)*chemo(OXYGEN)%medium_Cbnd_prev
!endif
!write(nflog,'(a,2e12.3)') 'medium_Cbnd: ',chemo(1:2)%medium_Cbnd
!chemo(OXYGEN)%medium_Cbnd_prev = chemo(OXYGEN)%medium_Cbnd

end subroutine

!---------------------------------------------------------------------------------- 
! The medium concentrations are updated explicitly, assuming a sphere with 
! known total uptake rate U.
! Compute U and update M, Cext, Cbnd.
! Note that for O2 Cext is fixed.
! Need to include decay
! Need to estimate U in a different way.
! Use the total cell uptake and decay within the blob.
! The change in total mass in the medium M is minus the sum of decay and U.dt
! The far-field concentration Cext is inferred from M.
! If the blob radius is R1 and the radius of the boundary layer is R2 = R1 + d_layer,
! then the concentration distribution in the layer at radius r is:
! C(r) = Cext + U/(4.pi.K).(1/R2 - 1/r)
! where K = diffusion coefficient in the unstirred layer
! By integrating C(r) and adding the portion at Cext we find the total mass is:
! M = Cext.Vmedium + (U/K)[(R2^3-R1^3)/(3R2) + (R2^2 - R1^2)/2)]
! which can be solved for Cext when U, M, R1 and R2 are known.
! Note that Vmedium = total_volume - Vblob = total_volume - (4/3)pi.R1^3
! From Cext, Cbnd = Cext + U/(4.pi.K).(1/R2 - 1/R1)
!----------------------------------------------------------------------------------
subroutine UpdateCbnd_mixed(dt)
real(REAL_KIND) :: dt
integer :: i, k, ichemo, ntvars
real(REAL_KIND) :: R1, R2, U(MAX_CHEMO)
real(REAL_KIND) :: a, Rlayer(MAX_CHEMO)
integer :: kcell, Nh, Nc
real(REAL_KIND) :: C, metab, dMdt, asum
real(REAL_KIND) :: Kin, Kout, decay_rate, vol_cm3, Cin, Cex
integer :: idrug, iparent, im

!write(*,*) 'UpdateCbnd_mixed'
! Start by looking at a conservative constituent (e.g. glucose)
! Contribution from cell uptake
!U = 0
!do ichemo = 1,MAX_CHEMO
!	if (.not.chemo(ichemo)%present) cycle
!	if (chemo(ichemo)%constant) then
!		chemo(ichemo)%medium_Cbnd = chemo(ichemo)%bdry_conc
!		cycle
!	endif
!	if (ichemo == OXYGEN .or. ichemo == GLUCOSE) then
!		Nh = chemo(ichemo)%Hill_N
!		asum = 0
!		Nc = 0
!		do kcell = 1,nlist
!			if (cell_list(kcell)%state == DEAD) cycle
!			Nc = Nc + 1
!			C = cell_list(kcell)%conc(ichemo)
!			metab = C**Nh/(chemo(ichemo)%MM_C0**Nh + C**Nh)
!			dMdt = metab*chemo(ichemo)%max_cell_rate 
!			asum = asum + dMdt
!		enddo
!		U(ichemo) = asum
!	else
!		! need to sum cell uptake and decay, and extracellular decay
!		decay_rate = chemo(ichemo)%decay_rate
!		Kin = chemo(ichemo)%membrane_diff_in
!		Kout = chemo(ichemo)%membrane_diff_out
!		asum = 0	
!		ntvars = ODEdiff%nextra + ODEdiff%nintra
!		do i = 1,ntvars
!			if (ODEdiff%vartype(i) == EXTRA) then
!				if (decay_rate == 0) cycle	! currently there is no distinction between intra- and extracellular decay rate
!				vol_cm3 = Vsite_cm3
!				Cex = allstate(i,ichemo)
!				if (i < ntvars) then
!					if (ODEdiff%vartype(i+1) == INTRA) then
!						kcell = ODEdiff%cell_index(i+1)		! for access to cell-specific parameters (note that the intra variable follows the extra variable)
!						vol_cm3 = Vsite_cm3 - Vcell_cm3*cell_list(kcell)%volume	! accounting for cell volume change
!					endif
!				endif
!				asum = asum + vol_cm3*Cex*decay_rate
!			else
!				Cin = allstate(i,ichemo)
!				kcell = ODEdiff%cell_index(i)
!				vol_cm3 = Vcell_cm3*cell_list(kcell)%volume
!				asum = asum + vol_cm3*Cin*decay_rate
!				! add cell uptake rate
!				Cex = allstate(i-1,ichemo)
!				asum = asum + Kin*Cex - Kout*Cin
!			endif
!		enddo
!		U(ichemo) = asum
!	endif
!enddo
!chemo(:)%medium_U = U(:)
!
!! First need the spheroid radius
!!call SetRadius(Nsites)
!R1 = blob_radius*DELTA_X		! cm
!Rlayer(:) = R1 + chemo(:)%medium_dlayer
!do ichemo = 1,MAX_CHEMO
!	if (.not.chemo(ichemo)%present) cycle
!	if (chemo(ichemo)%constant) cycle
!	R2 = Rlayer(ichemo)
!	if (ichemo /= OXYGEN) then	! update %medium_M, then %medium_Cext
!		chemo(ichemo)%medium_M = chemo(ichemo)%medium_M*(1 - chemo(ichemo)%decay_rate*dt) - U(ichemo)*dt
!		chemo(ichemo)%medium_M = max(0.0,chemo(ichemo)%medium_M)
!		chemo(ichemo)%medium_Cext = (chemo(ichemo)%medium_M - (U(ichemo)/(6*chemo(ichemo)%medium_diff_coef)*R2) &
!			*(R1*R1*(3*R2 - 2*R1) - R2*R2*R2))/(total_volume - 4*PI*R1*R1*R1/3.)
!	endif
!	a = (1/R2 - 1/R1)/(4*PI*chemo(ichemo)%medium_diff_coef)
!	chemo(ichemo)%medium_Cbnd = chemo(ichemo)%medium_Cext + a*chemo(ichemo)%medium_U
!	if (chemo(ichemo)%medium_Cbnd < 0) then
!		write(nflog,'(a,i2)') 'Setting negative medium_Cbnd to 0: Cbnd,M,Cext,a,U: ',ichemo
!		write(nflog,'(5e12.3)') chemo(ichemo)%medium_Cbnd,chemo(ichemo)%medium_M,chemo(ichemo)%medium_Cext,a,chemo(ichemo)%medium_U
!		chemo(ichemo)%medium_Cbnd = 0
!	endif
!enddo
!
!do idrug = 1,ndrugs_used
!	iparent = TRACER + 1 + 3*(idrug-1)
!	if (chemo(iparent)%present) then		! simulation with this drug has started
!	    do im = 0,2
!	        ichemo =iparent + im
!	        if (chemo(ichemo)%medium_Cext > Cthreshold) drug_gt_cthreshold(idrug) = .true.
!	        if (chemo(ichemo)%medium_Cbnd > Cthreshold) drug_gt_cthreshold(idrug) = .true.
!	    enddo
!	endif
!enddo

!write(nflog,'(a,10e12.3)') 'UpdateCbnd_mixed: ',chemo(:)%medium_Cbnd
end subroutine


end module

