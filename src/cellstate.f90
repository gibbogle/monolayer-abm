! Cancer cell state development

module cellstate
use global
use chemokine
use cycle_mod
implicit none

integer :: kcell_dividing = 0

contains

!-----------------------------------------------------------------------------------------
! Need to initialize site and cell concentrations when a cell divides and when there is
! cell death.
!-----------------------------------------------------------------------------------------
subroutine GrowCells(dose,dt,ok)
real(REAL_KIND) :: dose, dt
logical :: ok
logical :: changed
integer :: kcell, idrug, ichemo

!call logger('GrowCells: ')
tnow = istep*DELTA_T
ok = .true.
call grower(dt,changed,ok)
if (.not.ok) return
if (dose > 0) then
	call Irradiation(dose, ok)
	if (.not.ok) return
endif
if (use_death) then
	call CellDeath(dt,ok)
	if (.not.ok) return
endif
end subroutine

!-----------------------------------------------------------------------------------------
! The O2 concentration to use with cell kcell is either the intracellular concentration,
! or if use_extracellular_O2, the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine getO2conc(cp, C_O2)
type(cell_type), pointer :: cp
real(REAL_KIND) :: Cex, C_O2

if (use_extracellular_O2 .and. istep > 1) then		! fix 30/04/2015
    C_O2 = chemo(OXYGEN)%conc
else
    C_O2 = cp%Cin(OXYGEN)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! The glucose concentration to use with cell kcell is either the intracellular concentration,
! or if use_extracellular_O2 (!), the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine getGlucoseconc(cp, C_glucose)
type(cell_type), pointer :: cp
real(REAL_KIND) :: C_glucose

C_glucose = cp%Cin(GLUCOSE)
end subroutine

!-----------------------------------------------------------------------------------------
! Irradiate cells with dose (and duration tmin to be added)
!-----------------------------------------------------------------------------------------
subroutine Irradiation(dose,ok)
real(REAL_KIND) :: dose
logical :: ok
integer :: kcell, site(3), iv, ityp, idrug, im, ichemo, kpar=0
real(REAL_KIND) :: C_O2, SER, p_death, p_recovery, R, kill_prob, tmin
real(REAL_KIND) :: Cs							! concentration of radiosensitising drug
real(REAL_KIND) :: SER_max0, SER_Km, SER_KO2	! SER parameters of the drug
real(REAL_KIND) :: SERmax						! max sensitisation at the drug concentration
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

ok = .true.
call logger('Irradiation')
!tnow = istep*DELTA_T	! seconds
if (use_volume_method) then
    do kcell = 1,nlist
        if (colony_simulation) then
            cp => ccell_list(kcell)
        else
            cp => cell_list(kcell)
        endif
	    if (cp%state == DEAD) cycle
	    if (cp%radiation_tag) cycle	! we do not tag twice (yet)
	    ityp = cp%celltype
	    call getO2conc(cp,C_O2)
	    ! Compute sensitisation SER
	    SER = 1.0
	    do idrug = 1,Ndrugs_used
		    ichemo = 4 + 3*(idrug-1)
		    if (.not.chemo(ichemo)%present) cycle
		    do im = 0,2
			    ichemo = 4 + 3*(idrug-1) + im
			    if (drug(idrug)%sensitises(ityp,im)) then
				    Cs = cp%Cin(ichemo)	! concentration of drug/metabolite in the cell
				    SER_max0 = drug(idrug)%SER_max(ityp,im)
				    SER_Km = drug(idrug)%SER_Km(ityp,im)
				    SER_KO2 = drug(idrug)%SER_KO2(ityp,im)
				    SERmax = (Cs*SER_max0 + SER_Km)/(Cs + SER_Km)
				    SER = SER*(C_O2 + SER_KO2*SERmax)/(C_O2 + SER_KO2)
			    endif
		    enddo
	    enddo
	    call get_kill_probs(ityp,dose,C_O2,SER,p_recovery,p_death)
	    kill_prob = 1 - p_recovery
	    R = par_uni(kpar)
	    if (R < kill_prob) then
		    cp%radiation_tag = .true.
		    Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
		    cp%p_rad_death = p_death
		    if (LQ(ityp)%growth_delay_N > 0 .and. cp%Iphase) then
			    cp%growth_delay = .true.
			    cp%dt_delay = LQ(ityp)%growth_delay_factor*dose
			    cp%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
		    else
			    cp%growth_delay = .false.
		    endif
	    elseif (use_radiation_growth_delay_all .and. LQ(ityp)%growth_delay_N > 0) then
		    cp%growth_delay = .true.
		    cp%dt_delay = LQ(ityp)%growth_delay_factor*dose
		    cp%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
	    else
		    cp%growth_delay = .false.
	    endif
    enddo
else
    ccp => cc_parameters
    tmin = 1.0      ! for now...
    do kcell = 1,nlist
        if (colony_simulation) then
            cp => ccell_list(kcell)
        else
            cp => cell_list(kcell)
        endif
	    if (cp%state == DEAD) cycle
	    ityp = cp%celltype
	    call getO2conc(cp,C_O2)
	    ! Compute sensitisation SER
	    SER = 1.0
	    do idrug = 1,Ndrugs_used
		    ichemo = 4 + 3*(idrug-1)
		    if (.not.chemo(ichemo)%present) cycle
		    do im = 0,2
			    ichemo = 4 + 3*(idrug-1) + im
			    if (drug(idrug)%sensitises(ityp,im)) then
				    Cs = cp%Cin(ichemo)	! concentration of drug/metabolite in the cell
				    SER_max0 = drug(idrug)%SER_max(ityp,im)
				    SER_Km = drug(idrug)%SER_Km(ityp,im)
				    SER_KO2 = drug(idrug)%SER_KO2(ityp,im)
				    SERmax = (Cs*SER_max0 + SER_Km)/(Cs + SER_Km)
				    SER = SER*(C_O2 + SER_KO2*SERmax)/(C_O2 + SER_KO2)
			    endif
		    enddo
	    enddo
        call radiation_damage(cp, ccp, dose, SER, tmin)
    enddo
endif
call logger('did Irradiation')
end subroutine

!-----------------------------------------------------------------------------------------
! A cell that receives a dose of radiation either recovers completely before reaching 
! mitosis or retains damage that has a probability of causing cell death during mitosis.
! A damaged cell that does not die at this point passes the damage on to the progeny cells.
! The probability of complete recovery = p_recovery = p_r
! The probability of death for a damaged cell at mitosis = p_death = p_d
! To ensure that the short-term death probability is consistent with the previous
! LQ formulation, we require p_d(1-p_r) = kill_prob as previously calculated.
! If p_d is determined (currently it is fixed), then 1-p_r = kill_prob/p_d,
! therefore p_r = 1 - kill_prob/p_d
!-----------------------------------------------------------------------------------------
subroutine get_kill_probs(ityp,dose,C_O2,SER,p_recovery,p_death)
integer :: ityp
real(REAL_KIND) :: dose, C_O2, SER, p_recovery, p_death
real(REAL_KIND) :: OER_alpha_d, OER_beta_d, expon, kill_prob_orig

OER_alpha_d = dose*(LQ(ityp)%OER_am*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
OER_beta_d = dose*(LQ(ityp)%OER_bm*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)

OER_alpha_d = OER_alpha_d*SER
OER_beta_d = OER_beta_d*SER

expon = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_beta_d**2
p_recovery = exp(-expon)	! = SF
p_death = LQ(ityp)%death_prob
end subroutine

!-----------------------------------------------------------------------------------------
! This is the probability of death at time of division of cell that received a dose of 
! radiation and did not recover.
! In general one would expect this to depend of damage, i.e. on dose and C_O2, but
! for now it is a constant for a cell type
!-----------------------------------------------------------------------------------------
subroutine get_pdeath(ityp,dose,C_O2,p_death)
integer :: ityp
real(REAL_KIND) :: dose, C_O2, p_death

p_death = LQ(ityp)%death_prob
end subroutine


!-----------------------------------------------------------------------------------------
! Cells can be tagged to die, or finally die of anoxia or aglucosia, or they can be tagged 
! for death at division time if the drug is effective.
! Note: if simulating colony, no tagging, no death from anoxia, aglucosia
!-----------------------------------------------------------------------------------------
subroutine CellDeath(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, nlist0, site(3), i, ichemo, idrug, im, ityp, killmodel, kpar=0 
real(REAL_KIND) :: C_O2, C_glucose, Cdrug, n_O2, kmet, Kd, dMdt, kill_prob, dkill_prob, death_prob
logical :: anoxia_death, aglucosia_death
type(drug_type), pointer :: dp
type(cell_type), pointer :: cp

!call logger('CellDeath')
ok = .true.
if (colony_simulation) then
    return
endif

!tnow = istep*DELTA_T	! seconds
anoxia_death = chemo(OXYGEN)%controls_death
aglucosia_death = chemo(GLUCOSE)%controls_death
nlist0 = nlist
do kcell = 1,nlist
    cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ityp = cp%celltype
	call getO2conc(cp,C_O2)
	if (cp%anoxia_tag) then
		if (tnow >= cp%t_anoxia_die) then
			call CellDies(cp)
			Nanoxia_dead(ityp) = Nanoxia_dead(ityp) + 1
			cycle
		endif
	else
		if (anoxia_death .and. C_O2 < anoxia_threshold) then
			cp%t_anoxia = cp%t_anoxia + dt
			if (cp%t_anoxia > t_anoxia_limit) then
				cp%t_anoxia_die = tnow + anoxia_death_delay	! time that the cell will die
				if (.not.cp%anoxia_tag) then	! don't retag
					Nanoxia_tag(ityp) = Nanoxia_tag(ityp) + 1
				endif
				cp%anoxia_tag = .true.						! tagged to die later
			endif
		else
			cp%t_anoxia = 0
		endif
	endif
	
	call getGlucoseconc(cp,C_glucose)
	if (cp%aglucosia_tag) then
		if (tnow >= cp%t_aglucosia_die) then
			call CellDies(cp)
			Naglucosia_dead(ityp) = Naglucosia_dead(ityp) + 1
			cycle
		endif
	else
		if (aglucosia_death .and. C_glucose < aglucosia_threshold) then
			cp%t_aglucosia = cp%t_aglucosia + dt
			if (cp%t_aglucosia > t_aglucosia_limit) then
				cp%t_aglucosia_die = tnow + aglucosia_death_delay	! time that the cell will die
				if (.not.cp%aglucosia_tag) then	! don't retag
					Naglucosia_tag(ityp) = Naglucosia_tag(ityp) + 1
				endif
				cp%aglucosia_tag = .true.						! tagged to die later
			endif
		else
			cp%t_aglucosia = 0
		endif
	endif
	
	do idrug = 1,ndrugs_used	
		dp => drug(idrug)
		ichemo = TRACER + 1 + 3*(idrug-1)	
		kill_prob = 0
		death_prob = 0
		do im = 0,2
			if (.not.dp%kills(ityp,im)) cycle
			killmodel = dp%kill_model(ityp,im)		! could use %drugclass to separate kill modes
			Cdrug = cp%Cin(ichemo + im)
			Kd = dp%Kd(ityp,im)
			n_O2 = dp%n_O2(ityp,im)
			kmet = (1 - dp%C2(ityp,im) + dp%C2(ityp,im)*dp%KO2(ityp,im)**n_O2/(dp%KO2(ityp,im)**n_O2 + C_O2**n_O2))*dp%Kmet0(ityp,im)
			dMdt = kmet*Cdrug
			call getDrugKillProb(killmodel,Kd,dMdt,Cdrug,dt,dkill_prob)
			kill_prob = kill_prob + dkill_prob
			death_prob = max(death_prob,dp%death_prob(ityp,im))
		enddo
	    if (.not.cp%drug_tag(idrug) .and. par_uni(kpar) < kill_prob) then		! don't tag more than once
			cp%p_drug_death(idrug) = death_prob
			cp%drug_tag(idrug) = .true.
            Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) + 1
		endif
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getDrugKillProb(kill_model,Kd,dMdt,Cdrug,dt,dkill_prob)
integer :: kill_model
real(REAL_KIND) :: Kd, dMdt, Cdrug, dt, dkill_prob
real(REAL_KIND) :: SF, SF1, dtstep, kill_prob, c
integer :: Nsteps, istep

if (kill_model == 1) then
	c = Kd*dMdt
elseif (kill_model == 2) then
	c = Kd*dMdt*Cdrug
elseif (kill_model == 3) then
	c = Kd*dMdt**2
elseif (kill_model == 4) then
	c = Kd*Cdrug
elseif (kill_model == 5) then
	c = Kd*Cdrug**2
endif
SF = exp(-c*dt)
!write(nflog,'(a,2e12.3)') 'SF1, SF2: ',SF1, SF2
dkill_prob = 1 - SF
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_CellDies
integer :: kcell, i, kpar=0
type(cell_type), pointer :: cp

do i = 1,10
	kcell = random_int(1,nlist,kpar)
	cp => cell_list(kcell)
	call CellDies(cp)
enddo
stop
end subroutine

!-----------------------------------------------------------------------------------------
! If the dying cell site is less than a specified fraction f_migrate of the blob radius,
! the site migrates towards the blob centre.
! %indx -> 0
! If the site is on the boundary, it is removed from the boundary list, and %indx -> OUTSIDE_TAG
! The cell contents are released into the site.
!-----------------------------------------------------------------------------------------
subroutine CellDies(cp)
integer :: kcell
integer :: site(3), ityp, idrug
type(cell_type), pointer :: cp

cp%state = DEAD
ityp = cp%celltype
Ncells = Ncells - 1
Ncells_type(ityp) = Ncells_type(ityp) - 1
if (cp%anoxia_tag) then
	Nanoxia_tag(ityp) = Nanoxia_tag(ityp) - 1
endif
if (cp%aglucosia_tag) then
	Naglucosia_tag(ityp) = Naglucosia_tag(ityp) - 1
endif
do idrug = 1,ndrugs_used
	if (cp%drug_tag(idrug)) then
		Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) - 1
	endif
enddo
if (cp%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
endif

ngaps = ngaps + 1
if (ngaps > max_ngaps) then
    write(logmsg,'(a,i6,i6)') 'CellDies: ngaps > max_ngaps: ',ngaps,max_ngaps
    call logger(logmsg)
    stop
endif
gaplist(ngaps) = kcell
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine AddToMedium(cp,site)
integer :: kcell, site(3)
integer :: ichemo
real(REAL_KIND) :: V, Cex(MAX_CHEMO), Cin(MAX_CHEMO)
type(cell_type),pointer :: cp
return

!Cex = occupancy(site(1),site(2),site(3))%C
Cin = cp%Cin
V = cp%V
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	chemo(ichemo)%medium_M = chemo(ichemo)%medium_M + V*Cin(ichemo) + (Vsite_cm3 - V)*Cex(ichemo)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine RemoveFromMedium
integer :: ichemo

do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	if (chemo(ichemo)%constant) cycle
	chemo(ichemo)%medium_M = chemo(ichemo)%medium_M - Vsite_cm3*chemo(ichemo)%medium_Cbnd
enddo
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_CellDivision(ok)
logical :: ok
integer :: kcell, kpar=0
type(cell_type), pointer :: cp

kcell = random_int(1,nlist,kpar)
cp => cell_list(kcell)
call divider(cp,ok)
end subroutine

!-----------------------------------------------------------------------------------------
! Cell growth, death and division are handled here.  Division occurs when cell volume 
! exceeds the divide volume. 
! As the cell grows we need to adjust both Cin and Cex to maintain mass conservation.
! GROWTH DELAY
! When a cell has received a dose of radiation (or possibly drug - not yet considered)
! the cycle time is increased by an amount that depends on the dose.  The delay may be
! transmitted to progeny cells.
! 
! NOTE: now the medium concentrations are not affected by cell growth
!-----------------------------------------------------------------------------------------
subroutine grower(dt, changed, ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: k, kcell, nlist0, ityp, idrug, prev_phase, kpar=0
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: R
integer :: ndivide, divide_list(1000)
logical :: drugkilled
logical :: mitosis_entry, in_mitosis, divide

ok = .true.
changed = .false.
ccp => cc_parameters
nlist0 = nlist
ndivide = 0
!tnow = istep*DELTA_T !+ t_fmover
!if (colony_simulation) write(*,*) 'grower: ',nlist0,use_volume_method,tnow
do kcell = 1,nlist0
	kcell_now = kcell
	if (colony_simulation) then
	    cp => ccell_list(kcell)
	else
    	cp => cell_list(kcell)
    endif
	if (cp%state == DEAD) cycle
	ityp = cp%celltype
	divide = .false.
	mitosis_entry = .false.
	in_mitosis = .false.
	if (use_volume_method) then
!        if (colony_simulation) then
!            write(*,'(a,i6,L2,2e12.3)') 'kcell: ',kcell,cp%Iphase,cp%V,cp%divide_volume
!        endif
	    if (cp%Iphase) then
		    call growcell(cp,dt)
		    if (cp%V > cp%divide_volume) then	! time to enter mitosis
    	        mitosis_entry = .true.
	        endif
	    else
	        in_mitosis = .true.
	    endif
	else
	    prev_phase = cp%phase
        call timestep(cp, ccp, dt)
        if (cp%phase >= M_phase) then
            if (prev_phase == Checkpoint2) then
                mitosis_entry = .true.
            else
                in_mitosis = .true.
            endif
        endif
		if (cp%phase < Checkpoint2) then
		    call growcell(cp,dt)
		endif	
	endif
	if (mitosis_entry) then
		drugkilled = .false.
		do idrug = 1,ndrugs_used
			if (cp%drug_tag(idrug)) then
				call CellDies(cp)
				changed = .true.
				Ndrug_dead(idrug,ityp) = Ndrug_dead(idrug,ityp) + 1
				drugkilled = .true.
				exit
			endif
		enddo
		if (drugkilled) cycle
			
		if (use_volume_method) then
			if (cp%growth_delay) then
				if (cp%G2_M) then
					if (tnow > cp%t_growth_delay_end) then
						cp%G2_M = .false.
					else
						cycle
					endif
				else
					cp%t_growth_delay_end = tnow + cp%dt_delay
					cp%G2_M = .true.
					cycle
				endif
			endif
			! try moving death prob test to here
			if (cp%radiation_tag) then
				R = par_uni(kpar)
				if (R < cp%p_rad_death) then
					call CellDies(cp)
					changed = .true.
					Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
					cycle
				endif
			endif		
		else
		    ! Check for cell death by radiation lesions
		    ! For simplicity: 
		    !   cell death occurs only at mitosis entry
		    !   remaining L1 lesions and L2c misrepair (non-reciprocal translocation) are treated the same way
		    !   L2a and L2b are treated as non-fatal
		    if (cp%NL1 > 0 .or. cp%NL2(2) > 0) then
				call CellDies(cp)
				changed = .true.
				Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
				cycle
			endif		        
		endif
		
		cp%Iphase = .false.
		cp%mitosis = 0
		cp%t_start_mitosis = tnow
!		ncells_mphase = ncells_mphase + 1
	elseif (in_mitosis) then
		cp%mitosis = (tnow - cp%t_start_mitosis)/mitosis_duration
        if (cp%mitosis >= 1) then
			divide = .true.
		endif
	endif
	if (divide) then
		ndivide = ndivide + 1
		divide_list(ndivide) = kcell
	endif
enddo
do k = 1,ndivide
	changed = .true.
	kcell = divide_list(k)
	if (colony_simulation) then
	    cp => ccell_list(kcell)
	else
    	cp => cell_list(kcell)
    endif
	call divider(cp, ok)
	if (.not.ok) return
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine growcell(cp, dt)
type(cell_type), pointer :: cp
real(REAL_KIND) :: dt
real(REAL_KIND) :: Cin_0(NCONST), Cex_0(NCONST)		! at some point NCONST -> MAX_CHEMO
real(REAL_KIND) :: dVdt,  dV, metab_O2, metab_glucose, metab
logical :: oxygen_growth, glucose_growth

if (colony_simulation) then
    metab = 1
    dVdt = get_dVdt(cp,metab)
else
    oxygen_growth = chemo(OXYGEN)%controls_growth
    glucose_growth = chemo(GLUCOSE)%controls_growth
    Cin_0 = cp%Cin
    metab_O2 = O2_metab(Cin_0(OXYGEN))	! Note that currently growth depends only on O2
    metab_glucose = glucose_metab(Cin_0(GLUCOSE))
    if (oxygen_growth .and. glucose_growth) then
	    metab = metab_O2*metab_glucose
    elseif (oxygen_growth) then
	    metab = metab_O2
    elseif (glucose_growth) then
	    metab = metab_glucose
    endif
    dVdt = get_dVdt(cp,metab)
    if (suppress_growth) then	! for checking solvers
	    dVdt = 0
    endif
endif
cp%dVdt = dVdt
dV = dVdt*dt
if (use_cell_cycle .and. .not.(cp%phase == G1_phase .or. cp%phase == S_phase .or. cp%phase == G2_phase)) then
    dV = 0
endif
cp%V = cp%V + dV
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine SetInitialGrowthRate
integer :: kcell, ityp
real(REAL_KIND) :: C_O2, C_glucose, metab, metab_O2, metab_glucose, dVdt
logical :: oxygen_growth, glucose_growth
type(cell_type), pointer :: cp

oxygen_growth = chemo(OXYGEN)%controls_growth
glucose_growth = chemo(GLUCOSE)%controls_growth
do kcell = 1,nlist
    if (colony_simulation) then
        cp => ccell_list(kcell)
    else
        cp => cell_list(kcell)
    endif
	if (cp%state == DEAD) cycle
	C_O2 = chemo(OXYGEN)%bdry_conc
	C_glucose = cp%Cin(GLUCOSE)
	if (oxygen_growth .and. glucose_growth) then
	    metab_O2 = O2_metab(C_O2)
		metab_glucose = glucose_metab(C_glucose)
		metab = metab_O2*metab_glucose
	elseif (oxygen_growth) then
	    metab_O2 = O2_metab(C_O2)
		metab = metab_O2
	elseif (glucose_growth) then
		metab_glucose = glucose_metab(C_glucose)
		metab = metab_glucose
	endif
	dVdt = get_dVdt(cp,metab)
	if (suppress_growth) then	! for checking solvers
		dVdt = 0
	endif
	cp%dVdt = dVdt
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Need to account for time spent in mitosis.  Because there is no growth during mitosis,
! the effective divide_time must be reduced by mitosis_duration.
! Note that TERMINAL_MITOSIS is the only option.
!-----------------------------------------------------------------------------------------
function get_dVdt(cp, metab) result(dVdt)
type(cell_type), pointer :: cp
real(REAL_KIND) :: metab, dVdt
integer :: ityp
real(REAL_KIND) :: r_mean, c_rate

ityp = cp%celltype
if (use_cell_cycle) then
    if (use_constant_growthrate) then
        dVdt = max_growthrate(ityp)
    else
        dVdt = metab*max_growthrate(ityp)
    endif
else
    if (use_V_dependence) then
	    if (use_constant_divide_volume) then
		    dVdt = metab*log(2.0)*cp%V/(cp%divide_time - mitosis_duration)
	    else
		    c_rate = log(2.0)/(divide_time_mean(ityp) - mitosis_duration)
		    dVdt = c_rate*cp%V*metab
	    endif
    else
	    if (use_constant_divide_volume) then
		    dVdt = 0.5*metab*Vdivide0/(cp%divide_time  - mitosis_duration)
	    else
		    ityp = cp%celltype
		    r_mean = 0.5*Vdivide0/(divide_time_mean(ityp) - mitosis_duration)
		    dVdt = r_mean*metab
	    endif
    endif
endif
end function


!-----------------------------------------------------------------------------------------
! New version
!-----------------------------------------------------------------------------------------
subroutine divider(cp1, ok)
type(cell_type), pointer :: cp1
!integer :: kcell1
logical :: ok
integer :: kcell2, ityp, nbrs0
real(REAL_KIND) :: r(3), c(3), cfse0, cfse2, V0, Tdiv
type(cell_type), pointer :: cp2

!write(*,*) 'divider:'
!write(logmsg,*) 'divider: ',kcell1 
!call logger(logmsg)
ok = .true.
!tnow = istep*DELTA_T
!cp1 => cell_list(kcell1)
if (ngaps > 0) then
    kcell2 = gaplist(ngaps)
    ngaps = ngaps - 1
else
	nlist = nlist + 1
	if (nlist > MAX_NLIST) then
		ok = .false.
		call logger('Error: Maximum number of cells MAX_NLIST has been exceeded.  Increase and rebuild.')
		return
	endif
	kcell2 = nlist
endif
ncells = ncells + 1
if (colony_simulation) then
    cp2 => ccell_list(kcell2)
else
	cp2 => cell_list(kcell2)
endif
ityp = cp1%celltype
ncells_type(ityp) = ncells_type(ityp) + 1
ncells_mphase = ncells_mphase - 1
!cp2 => cell_list(kcell2)

cp1%state = ALIVE
cp1%generation = cp1%generation + 1
V0 = cp1%V/2
cp1%V = V0
cp1%birthtime = tnow
cp1%divide_volume = get_divide_volume(ityp,V0,Tdiv)
cp1%divide_time = Tdiv
cp1%mitosis = 0
cfse0 = cp1%CFSE
cp1%CFSE = generate_CFSE(cfse0/2)
cfse2 = cfse0 - cp1%CFSE

cp1%drug_tag = .false.
cp1%anoxia_tag = .false.
cp1%t_anoxia = 0
cp1%aglucosia_tag = .false.
cp1%t_aglucosia = 0

if (cp1%growth_delay) then
	cp1%N_delayed_cycles_left = cp1%N_delayed_cycles_left - 1
	cp1%growth_delay = (cp1%N_delayed_cycles_left > 0)
endif
cp1%G2_M = .false.
cp1%Iphase = .true.
cp1%phase = G1_phase

ndoublings = ndoublings + 1
doubling_time_sum = doubling_time_sum + tnow - cp1%t_divide_last
cp1%t_divide_last = tnow

! Second cell
!cell_list(kcell2) = cell_list(kcell1)
cp2 = cp1

! These are the variations from cp1
cp2%divide_volume = get_divide_volume(ityp,V0,Tdiv)
cp2%divide_time = Tdiv
cp2%CFSE = cfse2
if (cp2%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
endif
!if (colony_simulation) write(*,'(a,i6,2e12.3)') 'new cell: ',kcell2,cp2%V,cp2%divide_volume
end subroutine

!-----------------------------------------------------------------------------------------
! The daughter cell kcell1 is given the same characteristics as kcell0 and placed at site1.
! Random variation is introduced into %divide_volume.
!-----------------------------------------------------------------------------------------
subroutine CloneCell(kcell0,kcell1,site1,ok)
integer :: kcell0, kcell1, site1(3), ityp, idrug
logical :: ok
integer :: kpar = 0
real(REAL_KIND) :: V0, Tdiv, R, Vex, Cex(MAX_CHEMO)

ok = .true.
!tnow = istep*DELTA_T

if (use_gaplist .and. ngaps > 0) then
    kcell1 = gaplist(ngaps)
    ngaps = ngaps - 1
else
    nlist = nlist + 1
    if (nlist > max_nlist) then
		write(logmsg,*) 'Error: Dimension of cell_list() has been exceeded: increase max_nlist and rebuild'
		call logger(logmsg)
		ok = .false.
		return
	endif
    kcell1 = nlist
endif

ityp = cell_list(kcell0)%celltype
Ncells = Ncells + 1
Ncells_type(ityp) = Ncells_type(ityp) + 1
!if (cell_list(kcell0)%generation > 1 .and. cell_list(kcell0)%radiation_tag) then 
!	write(*,*) 'radiation-tagged cell gen > 1 divides: ',kcell0,cell_list(kcell0)%generation 
!	stop
!endif
cell_list(kcell0)%generation = cell_list(kcell0)%generation + 1
if (cell_list(kcell0)%growth_delay) then
	cell_list(kcell0)%N_delayed_cycles_left = cell_list(kcell0)%N_delayed_cycles_left - 1
	cell_list(kcell0)%growth_delay = (cell_list(kcell0)%N_delayed_cycles_left > 0)
endif
cell_list(kcell0)%G2_M = .false.
cell_list(kcell0)%t_divide_last = tnow
cell_list(kcell1)%celltype = cell_list(kcell0)%celltype
cell_list(kcell1)%state = cell_list(kcell0)%state
cell_list(kcell1)%generation = cell_list(kcell0)%generation
!cell_list(kcell1)%site = site1
!cell_list(kcell1)%ID = lastID
cell_list(kcell1)%ID = cell_list(kcell0)%ID
!if (cell_list(kcell0)%ID == 582) then
!	write(*,*) 'New cell with ID=582: ',kcell1
!endif
cell_list(kcell1)%p_rad_death = cell_list(kcell0)%p_rad_death
cell_list(kcell1)%p_drug_death = cell_list(kcell0)%p_drug_death
cell_list(kcell1)%radiation_tag = cell_list(kcell0)%radiation_tag
if (cell_list(kcell1)%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
endif
cell_list(kcell1)%drug_tag = cell_list(kcell0)%drug_tag
do idrug = 1,ndrugs_used
	if (cell_list(kcell1)%drug_tag(idrug)) then
		Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) + 1
	endif
enddo
cell_list(kcell1)%anoxia_tag = .false.
cell_list(kcell1)%aglucosia_tag = .false.
!cell_list(kcell1)%exists = .true.
cell_list(kcell1)%active = .true.
cell_list(kcell1)%growth_delay = cell_list(kcell0)%growth_delay
if (cell_list(kcell1)%growth_delay) then
	cell_list(kcell1)%dt_delay = cell_list(kcell0)%dt_delay
	cell_list(kcell1)%N_delayed_cycles_left = cell_list(kcell0)%N_delayed_cycles_left
endif
cell_list(kcell1)%G2_M = .false.
cell_list(kcell1)%t_divide_last = tnow
cell_list(kcell1)%dVdt = cell_list(kcell0)%dVdt
cell_list(kcell1)%V = cell_list(kcell0)%V
!R = par_uni(kpar)
!cell_list(kcell1)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
V0 = cell_list(kcell0)%V
cell_list(kcell1)%divide_volume = get_divide_volume(ityp, V0, Tdiv)
cell_list(kcell1)%divide_time = Tdiv
cell_list(kcell1)%t_anoxia = 0
cell_list(kcell1)%t_aglucosia = 0
cell_list(kcell1)%Cin = cell_list(kcell0)%Cin
!cell_list(kcell1)%Cex = cell_list(kcell0)%Cex
cell_list(kcell1)%dCdt = cell_list(kcell0)%dCdt
cell_list(kcell1)%dMdt = cell_list(kcell0)%dMdt
cell_list(kcell1)%M = cell_list(kcell0)%M
!occupancy(site1(1),site1(2),site1(3))%indx(1) = kcell1

!cell_list(kcell1)%Cex = occupancy(site1(1),site1(2),site1(3))%C
end subroutine

!----------------------------------------------------------------------------------
! Makes a slight modification to the Michaelis-Menten function to create a
! "soft landing" as C -> 0
! f1(C) = (C-d)/(C0 + C-d)  is a shifted version of the MM curve
! f0(C) = kC^2             is a function with derivative -> 0 as C -> 0
! At C=e, we want f0(e) = f1(e), and f0'(e) = f1'(e)
! =>
! ke^2 = (e-d)/(C0 + e-d)
! 2ke = C0/(Co + e-d)^2
! => e/2 = (e-d)(C0 + e-d)/C0
! Set x = e-d
! C0(x+d) = 2x(C0+x)
! => x = e-d = (sqrt(Co^2 + 8dC0) - C0)/4
! k = (e-d)/(e^2(C0 + e-d))
! We fix d (as small as possible, by trial and error) then deduce e, k.
! Note: This has really been superceded by the option of a Hill function with N = 2.
!----------------------------------------------------------------------------------
subroutine AdjustMM
real(REAL_KIND) :: deltaC, C0, C1

C0 = chemo(OXYGEN)%MM_C0
if (MM_THRESHOLD > 0) then
	deltaC = MM_THRESHOLD
	C1 = deltaC + (sqrt(C0*C0 + 8*C0*deltaC) - C0)/4
	ODEdiff%k_soft = (C1-deltaC)/(C1*C1*(C0+C1-deltaC))
	ODEdiff%C1_soft = C1
	ODEdiff%deltaC_soft = deltaC
else
	ODEdiff%k_soft = 0
	ODEdiff%C1_soft = 0
	ODEdiff%deltaC_soft = 0
endif
!write(logmsg,'(a,4e12.4)') 'AdjustMM: C0, deltaC, C1, k: ',C0, ODEdiff%deltaC_soft, ODEdiff%C1_soft, ODEdiff%k_soft
!call logger(logmsg)
end subroutine

end module
