! Cancer cell state development

module cellstate
use global
!use boundary
use chemokine
!use ode_diffuse
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
integer :: kcell, idrug, ichemo

!call logger('GrowCells')
ok = .true.
!if (use_radiation .and. dose > 0) then
if (dose > 0) then
	call Irradiation(dose, ok)
	if (.not.ok) return
endif
call CellGrowth(dt,ok)
if (.not.ok) return

if (use_death) then
	call CellDeath(dt,ok)
	if (.not.ok) return
endif
end subroutine

!-----------------------------------------------------------------------------------------
! The O2 concentration to use with cell kcell is either the intracellular concentration,
! or if use_extracellular_O2, the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine getO2conc(kcell, C_O2)
integer :: kcell
real(REAL_KIND) :: C_O2
integer :: iv, site(3)
real(REAL_KIND) :: tnow

if (use_extracellular_O2 .and. istep > 1) then		! fix 30/04/2015
!	iv = cell_list(kcell)%ivin
!	if (iv < 1) then
!		C_O2 = cell_list(kcell)%conc(OXYGEN)		! use this until %ivin is set in SetupODEDiff
!!		write(logmsg,*) 'getO2conc: ',istep,kcell,site,iv
!!		call logger(logmsg)
!!		stop
!!		tnow = istep*DELTA_T
!!		C_O2 = BdryConc(OXYGEN,tnow)	! assume that this is a site at the boundary
!	else
!		C_O2 = allstate(iv-1,OXYGEN)
!	endif
    C_O2 = chemo(OXYGEN)%conc
else
	C_O2 = cell_list(kcell)%conc(OXYGEN)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! The glucose concentration to use with cell kcell is either the intracellular concentration,
! or if use_extracellular_O2 (!), the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine getGlucoseconc(kcell, C_glucose)
integer :: kcell
real(REAL_KIND) :: C_glucose
integer :: iv, site(3)
real(REAL_KIND) :: tnow

if (use_extracellular_O2 .and. istep > 1) then		! fix 30/04/2015
!	iv = cell_list(kcell)%ivin
!	if (iv < 1) then
!		C_glucose = cell_list(kcell)%conc(GLUCOSE)		! use this until %ivin is set in SetupODEDiff
!	else
!		C_glucose = allstate(iv-1,GLUCOSE)
!	endif
    C_glucose = chemo(GLUCOSE)%conc
else
	C_glucose = cell_list(kcell)%conc(GLUCOSE)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Irradiate cells with dose.
!-----------------------------------------------------------------------------------------
subroutine Irradiation(dose,ok)
real(REAL_KIND) :: dose
logical :: ok
integer :: kcell, site(3), iv, ityp, idrug, im, ichemo, kpar=0
real(REAL_KIND) :: C_O2, SER, p_death, p_recovery, R, kill_prob, tnow
real(REAL_KIND) :: Cs							! concentration of radiosensitising drug
real(REAL_KIND) :: SER_max0, SER_Km, SER_KO2	! SER parameters of the drug
real(REAL_KIND) :: SERmax						! max sensitisation at the drug concentration

ok = .true.
call logger('Irradiation')
tnow = istep*DELTA_T	! seconds
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%radiation_tag) cycle	! we do not tag twice (yet)
	ityp = cell_list(kcell)%celltype
	call getO2conc(kcell,C_O2)
	! Compute sensitisation SER
	SER = 1.0
	do idrug = 1,Ndrugs_used
		ichemo = 4 + 3*(idrug-1)
		if (.not.chemo(ichemo)%present) cycle
		do im = 0,2
			ichemo = 4 + 3*(idrug-1) + im
			if (drug(idrug)%sensitises(ityp,im)) then
				Cs = cell_list(kcell)%conc(ichemo)	! concentration of drug/metabolite in the cell
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
		cell_list(kcell)%radiation_tag = .true.
		Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
		cell_list(kcell)%p_rad_death = p_death
		if (LQ(ityp)%growth_delay_N > 0) then
			cell_list(kcell)%growth_delay = .true.
			cell_list(kcell)%dt_delay = LQ(ityp)%growth_delay_factor*dose
			cell_list(kcell)%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
		else
			cell_list(kcell)%growth_delay = .false.
		endif
	elseif (use_radiation_growth_delay_all .and. LQ(ityp)%growth_delay_N > 0) then
		cell_list(kcell)%growth_delay = .true.
		cell_list(kcell)%dt_delay = LQ(ityp)%growth_delay_factor*dose
		cell_list(kcell)%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
	else
		cell_list(kcell)%growth_delay = .false.
	endif
enddo
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
!expon_err = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_alpha_d**2	! 07/08/2015

OER_alpha_d = OER_alpha_d*SER
OER_beta_d = OER_beta_d*SER

expon = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_beta_d**2
p_recovery = exp(-expon)	! = SF
p_death = LQ(ityp)%death_prob

!kill_prob_orig = 1 - exp(-expon)
!call get_pdeath(ityp,dose,C_O2,p_death)
!p_recovery = 1 - kill_prob_orig/p_death
!p_recovery = 1 - kill_prob_orig
!write(nflog,'(a,4e12.3)') 'get_kill_probs: OER_alpha_d,OER_beta_d,expon,expon_err: ',OER_alpha_d,OER_beta_d,expon,expon_err
!write(nflog,'(a,2e12.3)') 'kill_prob, kill_prob_err: ',kill_prob_orig,1 - exp(-expon_err)
!write(nflog,'(a,e12.3)') 'p_recovery: ',p_recovery
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
! Cells move to preferable nearby sites.
! For now this is turned off - need to formulate a sensible interpretation of "preferable"
!-----------------------------------------------------------------------------------------
!subroutine CellMigration(ok)
!logical :: ok
!integer :: kcell, j, indx, site0(3), site(3), jmax
!real(REAL_KIND) :: C0(MAX_CHEMO), C(MAX_CHEMO), v0, v, vmax, d0, d
!
!call logger('CellMigration is not yet implemented')
!ok = .false.
!return
!
!ok = .true.
!do kcell = 1,nlist
!	if (cell_list(kcell)%state == DEAD) cycle
!	site0 = cell_list(kcell)%site
!	C0 = occupancy(site0(1),site0(2),site0(3))%C(:)
!	v0 = SiteValue(C0)
!	d0 = cdistance(site0)
!	jmax = 0
!	vmax = -1.0e10
!	do j = 1,27
!		if (j == 14) cycle
!		site = site0 + jumpvec(:,j)
!		indx = occupancy(site(1),site(2),site(3))%indx(1)
!!		if (indx < -100) then	! necrotic site
!		if (indx == 0) then	!	vacant site
!			C = occupancy(site(1),site(2),site(3))%C(:)
!			v = SiteValue(C)
!			d = cdistance(site)
!			if (d > d0 .and. v > v0) then
!				vmax = v
!				jmax = j
!			endif
!		endif
!	enddo
!	if (jmax > 0) then
!		site = site0 + jumpvec(:,jmax)
!		indx = occupancy(site(1),site(2),site(3))%indx(1)
!		cell_list(kcell)%site = site
!		occupancy(site(1),site(2),site(3))%indx(1) = kcell
!		occupancy(site0(1),site0(2),site0(3))%indx(1) = indx
!	endif
!enddo
!end subroutine

!-----------------------------------------------------------------------------------------
! A measure of the attractiveness of a site with concentrations C(:)
!-----------------------------------------------------------------------------------------
!real(REAL_KIND) function SiteValue(C)
!real(REAL_KIND) :: C(:)
!
!SiteValue = C(OXYGEN)
!end function

!-----------------------------------------------------------------------------------------
! Cells can be tagged to die, or finally die of anoxia or aglucosia, or they can be tagged 
! for death at division time if the drug is effective.
!-----------------------------------------------------------------------------------------
subroutine CellDeath(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, nlist0, site(3), i, ichemo, idrug, im, ityp, killmodel, kpar=0 
real(REAL_KIND) :: C_O2, C_glucose, Cdrug, n_O2, kmet, Kd, dMdt, kill_prob, dkill_prob, death_prob, tnow
logical :: anoxia_death, aglucosia_death
type(drug_type), pointer :: dp

!call logger('CellDeath')
ok = .true.
tnow = istep*DELTA_T	! seconds
anoxia_death = chemo(OXYGEN)%controls_death
aglucosia_death = chemo(GLUCOSE)%controls_death
nlist0 = nlist
!write(logmsg,*) 'Nanoxia_tag: ',Nanoxia_tag(1),'  Nanoxia_dead: ',Nanoxia_dead(1)
!call logger(logmsg)
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ityp = cell_list(kcell)%celltype
	call getO2conc(kcell,C_O2)
	if (cell_list(kcell)%anoxia_tag) then
!		write(logmsg,*) 'anoxia_tag: ',kcell,cell_list(kcell)%state,tnow,cell_list(kcell)%t_anoxia_die
!		call logger(logmsg)
		if (tnow >= cell_list(kcell)%t_anoxia_die) then
			call CellDies(kcell)
			Nanoxia_dead(ityp) = Nanoxia_dead(ityp) + 1
!			write(logmsg,*) 'cell dies: ndead: ',Nanoxia_dead(ityp)
!			call logger(logmsg)
			cycle
		endif
	else
		if (anoxia_death .and. C_O2 < anoxia_threshold) then
			cell_list(kcell)%t_anoxia = cell_list(kcell)%t_anoxia + dt
			if (cell_list(kcell)%t_anoxia > t_anoxia_limit) then
				cell_list(kcell)%t_anoxia_die = tnow + anoxia_death_delay	! time that the cell will die
				if (.not.cell_list(kcell)%anoxia_tag) then	! don't retag
					Nanoxia_tag(ityp) = Nanoxia_tag(ityp) + 1
				endif
				cell_list(kcell)%anoxia_tag = .true.						! tagged to die later
			endif
		else
			cell_list(kcell)%t_anoxia = 0
		endif
	endif
	
	call getGlucoseconc(kcell,C_glucose)
	if (cell_list(kcell)%aglucosia_tag) then
!		write(logmsg,*) 'aglucosia_tag: ',kcell,cell_list(kcell)%state,tnow,cell_list(kcell)%t_aglucosia_die
!		call logger(logmsg)
		if (tnow >= cell_list(kcell)%t_aglucosia_die) then
			call CellDies(kcell)
			Naglucosia_dead(ityp) = Naglucosia_dead(ityp) + 1
!			write(logmsg,*) 'cell dies: ndead: ',Naglucosia_dead(ityp)
!			call logger(logmsg)
			cycle
		endif
	else
		if (aglucosia_death .and. C_glucose < aglucosia_threshold) then
			cell_list(kcell)%t_aglucosia = cell_list(kcell)%t_aglucosia + dt
			if (cell_list(kcell)%t_aglucosia > t_aglucosia_limit) then
				cell_list(kcell)%t_aglucosia_die = tnow + aglucosia_death_delay	! time that the cell will die
				if (.not.cell_list(kcell)%aglucosia_tag) then	! don't retag
					Naglucosia_tag(ityp) = Naglucosia_tag(ityp) + 1
				endif
				cell_list(kcell)%aglucosia_tag = .true.						! tagged to die later
			endif
		else
			cell_list(kcell)%t_aglucosia = 0
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
			Cdrug = cell_list(kcell)%conc(ichemo + im)
			Kd = dp%Kd(ityp,im)
			n_O2 = dp%n_O2(ityp,im)
			kmet = (1 - dp%C2(ityp,im) + dp%C2(ityp,im)*dp%KO2(ityp,im)**n_O2/(dp%KO2(ityp,im)**n_O2 + C_O2**n_O2))*dp%Kmet0(ityp,im)
			dMdt = kmet*Cdrug
			call getDrugKillProb(killmodel,Kd,dMdt,Cdrug,dt,dkill_prob)
			kill_prob = kill_prob + dkill_prob
			death_prob = max(death_prob,dp%death_prob(ityp,im))
		enddo
	    if (.not.cell_list(kcell)%drug_tag(idrug) .and. par_uni(kpar) < kill_prob) then		! don't tag more than once
			cell_list(kcell)%p_drug_death(idrug) = death_prob
			cell_list(kcell)%drug_tag(idrug) = .true.
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

!Nsteps = dt + 0.5
!dtstep = 1.0
!SF1 = 1.0
!do istep = 1,Nsteps
!	if (kill_model == 1) then
!		kill_prob = Kd*dMdt*dtstep
!	elseif (kill_model == 2) then
!		kill_prob = Kd*dMdt*Cdrug*dtstep
!	elseif (kill_model == 3) then
!		kill_prob = Kd*dMdt**2*dtstep
!	elseif (kill_model == 4) then
!		kill_prob = Kd*Cdrug*dtstep
!	elseif (kill_model == 5) then
!		kill_prob = Kd*(Cdrug**2)*dtstep
!	endif
!    kill_prob = min(kill_prob,1.0)
!    SF1 = SF1*(1 - kill_prob)
!enddo

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

do i = 1,10
	kcell = random_int(1,nlist,kpar)
	call CellDies(kcell)
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
subroutine CellDies(kcell)
integer :: kcell
integer :: site(3), ityp, idrug
real(REAL_KIND) :: V

cell_list(kcell)%state = DEAD
cell_list(kcell)%exists = .false.
ityp = cell_list(kcell)%celltype
Ncells = Ncells - 1
Ncells_type(ityp) = Ncells_type(ityp) - 1
if (cell_list(kcell)%anoxia_tag) then
	Nanoxia_tag(ityp) = Nanoxia_tag(ityp) - 1
endif
if (cell_list(kcell)%aglucosia_tag) then
	Naglucosia_tag(ityp) = Naglucosia_tag(ityp) - 1
endif
do idrug = 1,ndrugs_used
	if (cell_list(kcell)%drug_tag(idrug)) then
		Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) - 1
	endif
enddo
if (cell_list(kcell)%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
endif

!site = cell_list(kcell)%site
!occupancy(site(1),site(2),site(3))%indx(1) = 0
!if (associated(occupancy(site(1),site(2),site(3))%bdry)) then
!	call bdrylist_delete(site,bdrylist)
!    nullify(occupancy(site(1),site(2),site(3))%bdry)
!	occupancy(site(1),site(2),site(3))%indx(1) = OUTSIDE_TAG
!	Nsites = Nsites - 1
!	bdry_changed = .true.
!	call OutsideNeighbours(site)
!	call AddToMedium(kcell,site)
!else
!	V = cell_list(kcell)%volume*Vcell_cm3
!	occupancy(site(1),site(2),site(3))%C = ((Vsite_cm3 - V)*occupancy(site(1),site(2),site(3))%C + V*cell_list(kcell)%conc)/Vsite_cm3
!endif
!call NecroticMigration(site)
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
subroutine AddToMedium(kcell,site)
integer :: kcell, site(3)
integer :: ichemo
real(REAL_KIND) :: V, Cex(MAX_CHEMO), Cin(MAX_CHEMO)

return

!Cex = occupancy(site(1),site(2),site(3))%C
Cin = cell_list(kcell)%conc
V = cell_list(kcell)%volume*Vcell_cm3
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
! When a bdry cell dies, need to check the neighbours to see if a site is now made outside.
! To be made outside a site must have %indx(1) = 0 and at least one Neumann neighbour outside.
!-----------------------------------------------------------------------------------------
!subroutine OutsideNeighbours(site0)
!integer :: site0(3)
!integer :: k, j, x, y, z, xx, yy, zz
!logical :: isout, done
!
!done = .false.
!do while(.not.done)
!	done = .true.
!	do k = 1,27
!		if (k == 14) cycle
!		x = site0(1) + jumpvec(1,k)
!		y = site0(2) + jumpvec(2,k)
!		z = site0(3) + jumpvec(3,k)
!		if (outside_xyz(x,y,z)) cycle
!		if (occupancy(x,y,z)%indx(1) == 0) then
!			isout = .false.
!			do j = 1,6
!				xx = x + neumann(1,j)
!				yy = y + neumann(2,j)
!				zz = z + neumann(3,j)
!				if (outside_xyz(xx,yy,zz)) cycle
!				if (occupancy(xx,yy,zz)%indx(1) == OUTSIDE_TAG) then
!					isout = .true.
!					exit
!				endif
!			enddo
!			if (isout) then
!				done = .false.
!				occupancy(x,y,z)%indx(1) = OUTSIDE_TAG
!			endif
!		endif
!	enddo
!enddo
!end subroutine
!
!!-----------------------------------------------------------------------------------------
!! A necrotic site migrates towards the blob centre, stopping when another necrotic 
!! site is reached
!!-----------------------------------------------------------------------------------------
!subroutine NecroticMigration(site0)
!integer :: site0(3)
!integer :: site1(3), site2(3), site(3), j, jmin, kcell, tmp_indx
!real(REAL_KIND) :: d1, d2, dmin
!
!!write(logmsg,*) 'NecroticMigration: site0: ',site0
!!call logger(logmsg)
!site1 = site0
!do
!	d1 = cdistance(site1)
!	dmin = 1.0e10
!	jmin = 0
!	do j = 1,27
!		if (j == 14) cycle
!		site = site1 + jumpvec(:,j)
!!		if (occupancy(site(1),site(2),site(3))%indx(1) < 0) cycle	! do not swap with another necrotic site
!		if (occupancy(site(1),site(2),site(3))%indx(1) == 0) cycle	! do not swap with another necrotic site
!		d2 = cdistance(site)
!		if (d2 < dmin) then
!			dmin = d2
!			jmin = j
!		endif
!	enddo
!	if (dmin >= d1) then
!		site2 = site1
!		exit
!	endif
!	if (jmin <= 0) then
!		write(nflog,*) 'Error: jmin: ',jmin
!		stop
!	endif
!	site2 = site1 + jumpvec(:,jmin)
!	! Now swap site1 and site2
!	kcell = occupancy(site2(1),site2(2),site2(3))%indx(1)
!	tmp_indx = occupancy(site1(1),site1(2),site1(3))%indx(1)
!	occupancy(site1(1),site1(2),site1(3))%indx(1) = kcell
!	occupancy(site2(1),site2(2),site2(3))%indx(1) = tmp_indx
!	cell_list(kcell)%site = site1
!	site1 = site2
!enddo
!end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_CellDivision(ok)
logical :: ok
integer :: kcell, kpar=0
kcell = random_int(1,nlist,kpar)
call CellDivider(kcell,ok)
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
subroutine CellGrowth(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, nlist0, site(3), ityp, idrug, ichemo, kpar=0
integer :: divide_list(10000), ndivide, i
real(REAL_KIND) :: tnow, C_O2, C_glucose, metab, metab_O2, metab_glucose, dVdt, vol0, R		!r_mean(2), c_rate(2)
real(REAL_KIND) :: r_mean(MAX_CELLTYPES), c_rate(MAX_CELLTYPES)
real(REAL_KIND) :: Vin_0, Vex_0, dV, minVex
real(REAL_KIND) :: Cin_0(MAX_CHEMO), Cex_0(MAX_CHEMO)
character*(20) :: msg
logical :: drugkilled, oxygen_growth, glucose_growth, first_cycle
integer :: C_option = 1
type(cell_type), pointer :: cp

!call logger('CellGrowth')
ok = .true.
nlist0 = nlist
tnow = istep*DELTA_T
c_rate(1:2) = log(2.0)/divide_time_mean(1:2)		! Note: to randomise divide time need to use random number, not mean!
r_mean(1:2) = Vdivide0/(2*divide_time_mean(1:2))
oxygen_growth = chemo(OXYGEN)%controls_growth
glucose_growth = chemo(GLUCOSE)%controls_growth
ndivide = 0
minVex = 1.0e10
do kcell = 1,nlist0
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ityp = cp%celltype
	if (cp%volume < cp%divide_volume) then	! still growing - not delayed
		C_O2 = cp%conc(OXYGEN)
		C_glucose = cp%conc(GLUCOSE)
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
		Cin_0 = cp%conc
!		Cex_0 = occupancy(site(1),site(2),site(3))%C
        Cex_0 = Caverage(MAX_CHEMO+1:2*MAX_CHEMO)
		cp%dVdt = dVdt
		Vin_0 = cp%volume*Vcell_cm3	! cm^3
		Vex_0 = Vsite_cm3 - Vin_0					! cm^3
		dV = dVdt*dt*Vcell_cm3						! cm^3
		cp%volume = (Vin_0 + dV)/Vcell_cm3
		if (C_option == 1) then
			! Calculation based on transfer of an extracellular volume dV with constituents, i.e. holding extracellular concentrations constant
			cp%conc = (Vin_0*Cin_0 + dV*Cex_0)/(Vin_0 + dV)
!			occupancy(site(1),site(2),site(3))%C = (Vex_0*Cex_0 - dV*Cex_0)/(Vex_0 - dV)	! = Cex_0
		elseif (C_option == 2) then
			! Calculation based on change in volumes without mass transfer of constituents
			cp%conc = Vin_0*Cin_0/(Vin_0 + dV)
!			occupancy(site(1),site(2),site(3))%C = Vex_0*Cex_0/(Vex_0 - dV)
		endif
		minVex = min(minVex,Vex_0 - dV)
	endif
	
	if (cp%volume > cp%divide_volume) then	! time to divide
		drugkilled = .false.
		do idrug = 1,ndrugs_used
			if (cp%drug_tag(idrug)) then
				R = par_uni(kpar)
				if (R < cp%p_drug_death(idrug)) then
					call CellDies(kcell)
					Ndrug_dead(idrug,ityp) = Ndrug_dead(idrug,ityp) + 1
					drugkilled = .true.
					exit
				endif
			endif
		enddo
		if (drugkilled) cycle

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
				call CellDies(kcell)
				Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
				cycle
			endif
		endif
	    ndivide = ndivide + 1
	    divide_list(ndivide) = kcell
	endif
enddo
do i = 1,ndivide
    kcell = divide_list(i)
	kcell_dividing = kcell
	call CellDivider(kcell, ok)
	if (.not.ok) return
enddo
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
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	C_O2 = chemo(OXYGEN)%bdry_conc
	C_glucose = cp%conc(GLUCOSE)
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
!-----------------------------------------------------------------------------------------
function get_dVdt(cp, metab) result(dVdt)
type(cell_type), pointer :: cp
real(REAL_KIND) :: metab, dVdt
integer :: ityp
real(REAL_KIND) :: r_mean, c_rate

if (use_V_dependence) then
	if (use_constant_divide_volume) then
		dVdt = metab*log(2.0)*cp%volume/cp%divide_time
	else
		ityp = cp%celltype
		c_rate = log(2.0)/divide_time_mean(ityp)
		dVdt = c_rate*cp%volume*metab
	endif
else
	if (use_constant_divide_volume) then
		dVdt = metab*Vdivide0/(2*cp%divide_time)
	else
		ityp = cp%celltype
		r_mean = Vdivide0/(2*divide_time_mean(ityp))
		dVdt = r_mean*metab
	endif
endif
end function

!-----------------------------------------------------------------------------------------
! The dividing cell, kcell0, is at site0.
! The cases are:
!   (0) If there is a vacant interior neighbour site, this is used for the daughter cell, done.
! Otherwise, a neighbour site site01 is chosen randomly.  This becomes the site for the daughter cell.
! The cell occupying site01 must be moved to make way for the daughter cell.
! site01 is:
!   (1) outside the blob
!   (2) on the boundary
!   (3) inside the blob
!
! Revised simple mass balance
! ---------------------------
! In case (0), the effect on mass balance is insignificant, since the two sites are neighbours
! and therefore have very similar extracellular concentrations.  For now no adjustment is done.
! In case (1), same as case (0), no adjustment needed.
! In cases (2) and (3), cells are displaced on a line from site01 to the first (nearest) vacant site.
! The previously vacant site that receives a cell, at the end of the line of moved cells,
! is site2.  Mass conservation is (approximately) conserved by an adjustment to the 
! extracellular concentration at this site.  For this the extra conc at site0 is needed.
! The extracellular concentration at site2 is adjusted to ensure that mass is conserved
! approximately, on average.
! Let: 
! Vc = daughter cell volume = 0.8
! Vi = volume of displaced cells, on average = 1.2
! C0 = extracellular conc at site of division
! Ci = extracellular conc at sites along the line of displacement
! Cn = extracellular conc at unoccupied site that receives a cell (end of the line)
! Change in extracellular mass:
! + Vc.C0 at division site
! + (V1 - Vc).C1 at neighbour site that daughter cell moves into
! - Vn.Cn at the end site
! Note that on average there is no change in extra volume at intermediate sites, since on
! average the cells have the same volume = 1.2
! Therefore the total change (increase) in mass dM = Vc.C0 + (V1 - Vc).C1 - Vn.Cn
! Since Vc = 0.8 and V1 = Vn = 1.2 (on average), dM = 0.8C0 + 0.4C1 - 1.2Cn
! and since C1 is approximately equal to C0 (neighbour sites), approx dM = 1.2(C0 - Cn)
! The adjustment subtracts dM from the extra conc at the end site site2.
! Cn -> (Vn.Cn - 1.2(C0-Cn))/Vn
! where Vn is the normalised volume.
! Note: this "simple" adjustment does not work.  The required additional mass at site2
! is excessive, leads, for example, to O2 concentration that exceeds the boundary value.
!-----------------------------------------------------------------------------------------
subroutine CellDivider(kcell0, ok)
integer :: kcell0
logical :: ok
integer :: kpar=0
integer :: j, k, kcell1, ityp, site0(3), site1(3), site2(3), site01(3), site(3), ichemo, nfree, bestsite(3)
integer :: npath, path(3,200)
real(REAL_KIND) :: tnow, R, v, vmax, V0, Tdiv, Cex0(MAX_CHEMO), M0(MAX_CHEMO), M1(MAX_CHEMO), alpha(MAX_CHEMO)
real(REAL_KIND) :: cfse0, cfse1
logical :: use_simple_mass_balance = .false.

if (dbug) then
	write(logmsg,*) 'CellDivider: ',kcell0
	call logger(logmsg)
endif
ok = .true.
tnow = istep*DELTA_T
cell_list(kcell0)%t_divide_last = tnow
V0 = cell_list(kcell0)%volume/2
cell_list(kcell0)%volume = V0
cfse0 = cell_list(kcell0)%CFSE
cell_list(kcell0)%CFSE = generate_CFSE(cfse0/2)
cfse1 = cfse0 - cell_list(kcell0)%CFSE
cell_list(kcell0)%t_anoxia = 0
cell_list(kcell0)%t_aglucosia = 0

!R = par_uni(kpar)
!cell_list(kcell0)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
ityp = cell_list(kcell0)%celltype
cell_list(kcell0)%divide_volume = get_divide_volume(ityp,V0, Tdiv)
cell_list(kcell0)%divide_time = Tdiv
cell_list(kcell0)%M = cell_list(kcell0)%M/2
!write(nflog,'(a,i6,2f8.2)') 'divide: ',kcell0,cell_list(kcell0)%volume,cell_list(kcell0)%divide_volume
!write(logmsg,'(a,f6.1)') 'Divide time: ',tnow/3600
!call logger(logmsg)

site0 = [0,0,0]
call CloneCell(kcell0,kcell1,site0,ok)
cell_list(kcell1)%CFSE = cfse1

end subroutine

!-----------------------------------------------------------------------------------------
! The daughter cell kcell1 is given the same characteristics as kcell0 and placed at site1.
! Random variation is introduced into %divide_volume.
!-----------------------------------------------------------------------------------------
subroutine CloneCell(kcell0,kcell1,site1,ok)
integer :: kcell0, kcell1, site1(3), ityp, idrug
!real(REAL_KIND) :: Cex0(:)
logical :: ok
integer :: kpar = 0
real(REAL_KIND) :: tnow, V0, Tdiv, R, Vex, Cex(MAX_CHEMO)

ok = .true.
tnow = istep*DELTA_T

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
cell_list(kcell1)%exists = .true.
cell_list(kcell1)%active = .true.
cell_list(kcell1)%growth_delay = cell_list(kcell0)%growth_delay
if (cell_list(kcell1)%growth_delay) then
	cell_list(kcell1)%dt_delay = cell_list(kcell0)%dt_delay
	cell_list(kcell1)%N_delayed_cycles_left = cell_list(kcell0)%N_delayed_cycles_left
endif
cell_list(kcell1)%G2_M = .false.
cell_list(kcell1)%t_divide_last = tnow
cell_list(kcell1)%dVdt = cell_list(kcell0)%dVdt
cell_list(kcell1)%volume = cell_list(kcell0)%volume
!R = par_uni(kpar)
!cell_list(kcell1)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
V0 = cell_list(kcell0)%volume
cell_list(kcell1)%divide_volume = get_divide_volume(ityp, V0, Tdiv)
cell_list(kcell1)%divide_time = Tdiv
cell_list(kcell1)%t_anoxia = 0
cell_list(kcell1)%t_aglucosia = 0
cell_list(kcell1)%conc = cell_list(kcell0)%conc
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
