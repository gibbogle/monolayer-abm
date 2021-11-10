!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_summary1(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), i_hypoxia_cutoff,i_growth_cutoff
integer :: Nviable(MAX_CELLTYPES), Nlive(MAX_CELLTYPES), plate_eff_10(MAX_CELLTYPES)
integer :: nhypoxic(3), nclonohypoxic(3), ngrowth(3), &
    hypoxic_percent_10, clonohypoxic_percent_10, growth_percent_10, necrotic_percent_10, &
    EC_oxygen_1000, EC_glucose_1000, EC_drug_10000(2,0:2), &
    IC_oxygen_1000, IC_glucose_1000, IC_drug_10000(2,0:2), &
    medium_oxygen_1000, medium_glucose_1000, medium_drug_10000(2,0:2)
integer :: TNanoxia_dead, TNaglucosia_dead, TNradiation_dead, TNdrug_dead(2),  TNviable, &
           Ntagged_anoxia(MAX_CELLTYPES), Ntagged_aglucosia(MAX_CELLTYPES), Ntagged_radiation(MAX_CELLTYPES), &
           Ntagged_drug(2,MAX_CELLTYPES), &
           TNtagged_anoxia, TNtagged_aglucosia, TNtagged_radiation, TNtagged_drug(2)
integer :: Tplate_eff_10   
integer :: ityp, i, im, idrug, iparent
real(REAL_KIND) :: hour, plate_eff(MAX_CELLTYPES), doubling_time
real(REAL_KIND) :: EC(MAX_CHEMO), cmedium(MAX_CHEMO), clonohypoxic_fraction(3)

hour = istep*DELTA_T/3600.

Ntagged_anoxia(:) = Nanoxia_tag(:)			! number currently tagged by anoxia
Ntagged_aglucosia(:) = Naglucosia_tag(:)	! number currently tagged by aglucosia
Ntagged_radiation(:) = Nradiation_tag(:)	! number currently tagged by radiation
Ntagged_drug(1,:) = Ndrug_tag(1,:)			! number currently tagged by drugA
Ntagged_drug(2,:) = Ndrug_tag(2,:)			! number currently tagged by drugB

TNtagged_anoxia = sum(Ntagged_anoxia(1:Ncelltypes))
TNtagged_aglucosia = sum(Ntagged_aglucosia(1:Ncelltypes))
TNtagged_radiation = sum(Ntagged_radiation(1:Ncelltypes))
TNtagged_drug(1) = sum(Ntagged_drug(1,1:Ncelltypes))
TNtagged_drug(2) = sum(Ntagged_drug(2,1:Ncelltypes))

TNanoxia_dead = sum(Nanoxia_dead(1:Ncelltypes))
TNaglucosia_dead = sum(Naglucosia_dead(1:Ncelltypes))
TNradiation_dead = sum(Nradiation_dead(1:Ncelltypes))
TNdrug_dead(1) = sum(Ndrug_dead(1,1:Ncelltypes))
TNdrug_dead(2) = sum(Ndrug_dead(2,1:Ncelltypes))

call getNviable(Nviable, Nlive)
TNviable = sum(Nviable(1:Ncelltypes))

call getHypoxicCount(nhypoxic)
hypoxic_percent_10 = (1000.*nhypoxic(i_hypoxia_cutoff))/Ncells
call getClonoHypoxicCount(nclonohypoxic)
if (TNviable > 0) then
	clonohypoxic_percent_10 = (1000.*nclonohypoxic(i_hypoxia_cutoff))/TNviable
	clonohypoxic_fraction = nclonohypoxic(:)/real(TNviable)
else
	clonohypoxic_percent_10 = 0
	clonohypoxic_fraction = 0
endif		
call getGrowthCount(ngrowth)
growth_percent_10 = (1000.*ngrowth(i_growth_cutoff))/Ncells
do ityp = 1,Ncelltypes
	if (Nlive(ityp) > 0) then
		plate_eff(ityp) = real(Nviable(ityp))/Nlive(ityp)
	else
		plate_eff(ityp) = 0
	endif
enddo
plate_eff_10 = 1000.*plate_eff
Tplate_eff_10 = 0
do ityp = 1,Ncelltypes
	Tplate_eff_10 = Tplate_eff_10 + plate_eff_10(ityp)*celltype_fraction(ityp)
enddo

call getMediumConc(EC,cmedium)
EC_oxygen_1000 = EC(OXYGEN)*1000.
EC_glucose_1000 = EC(GLUCOSE)*1000.
do idrug = 1,2
	do im = 0,2
		iparent = DRUG_A + 3*(idrug-1)
		EC_drug_10000(idrug,im) = EC(iparent+im)*10000.
	enddo
enddo

medium_oxygen_1000 = cmedium(OXYGEN)*1000.
medium_glucose_1000 = cmedium(GLUCOSE)*1000.
do idrug = 1,2
	do im = 0,2
		iparent = DRUG_A + 3*(idrug-1)
		medium_drug_10000(idrug,im) = cmedium(iparent+im)*10000.
	enddo
enddo

IC_oxygen_1000 = caverage(OXYGEN)*1000.
IC_glucose_1000 = caverage(GLUCOSE)*1000.
do idrug = 1,2
	do im = 0,2
		iparent = DRUG_A + 3*(idrug-1)
		IC_drug_10000(idrug,im) = caverage(iparent+im)*10000.
	enddo
enddo

if (ndivided /= ndoublings) then
	write(*,*) 'ndivided /= ndoublings: ',ndivided,ndoublings
	stop
endif
if (ndoublings > 0) then
    doubling_time = doubling_time_sum/(3600*ndoublings)
else
    doubling_time = 0
endif

summaryData(1:42) = [ istep, Ncells, TNanoxia_dead, TNaglucosia_dead, TNdrug_dead(1), TNdrug_dead(2), TNradiation_dead, &
    TNtagged_anoxia, TNtagged_aglucosia, TNtagged_drug(1), TNtagged_drug(2), TNtagged_radiation, &
	hypoxic_percent_10, clonohypoxic_percent_10, growth_percent_10, Tplate_eff_10, &
	EC_oxygen_1000, EC_glucose_1000, EC_drug_10000(1,:), EC_drug_10000(2,:), &
	IC_oxygen_1000, IC_glucose_1000, IC_drug_10000(1,:), IC_drug_10000(2,:), &
	medium_oxygen_1000, medium_glucose_1000, medium_drug_10000(1,:), medium_drug_10000(2,:), &
	int(100*doubling_time), ndivided ]
write(nfres,'(a,a,2a12,i8,e12.4,22i7,36e12.4,i6)') trim(header),' ',gui_run_version, dll_run_version, &
	istep, hour, Ncells_type(1:2), &
    Nanoxia_dead(1:2), Naglucosia_dead(1:2), Ndrug_dead(1,1:2), &
    Ndrug_dead(2,1:2), Nradiation_dead(1:2), &
    Ntagged_anoxia(1:2), Ntagged_aglucosia(1:2), Ntagged_drug(1,1:2), &
    Ntagged_drug(2,1:2), Ntagged_radiation(1:2), &
	nhypoxic(:)/real(Ncells), clonohypoxic_fraction(:), ngrowth(:)/real(Ncells), &
	plate_eff(1:2), &
	EC(OXYGEN), EC(GLUCOSE), EC(DRUG_A:DRUG_A+2), EC(DRUG_B:DRUG_B+2), & 
	caverage(OXYGEN), caverage(GLUCOSE), caverage(DRUG_A:DRUG_A+2), caverage(DRUG_B:DRUG_B+2), &
	cmedium(OXYGEN), cmedium(GLUCOSE), cmedium(DRUG_A:DRUG_A+2), cmedium(DRUG_B:DRUG_B+2), & 
	doubling_time,  ndivided
	
!call sum_dMdt(GLUCOSE)
ndoublings = 0
doubling_time_sum = 0
ndivided = 0
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine Solver1(it,tstart,dt,nc,ok)
integer :: it, nc
real(REAL_KIND) :: tstart, dt
logical :: ok
integer :: ichemo, ic, k, ict, ncvars, neqn, kcell, i
real(REAL_KIND) :: t, tend
real(REAL_KIND) :: C(2*MAX_CHEMO), Csum
real(REAL_KIND) :: timer1, timer2
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1)
type(rkc_comm) :: comm_rkc(1)
logical :: solve_O2 = .true.
logical :: use_drugsolver = .true.

ok = .true.
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
	if (use_drugsolver .and. ichemo >= DRUG_A) then
        chemo_active(k) = .false.
	endif	
	C(k) = Caverage(ichemo)                ! average cell concentration
enddo
ncvars = k
! Note: ncvars = nchemo
do ic = 1,nchemo
	ichemo = chemomap(ic)
	k = k + 1
	C(k) = Caverage(MAX_CHEMO + ichemo)      ! average EC concentration
    chemo_active(k) = .not.chemo(ichemo)%constant
    if (ichemo == OXYGEN) then
        ! Suppress solving for medium oxygen
        chemo_active(k) = .false.
    endif
	if (use_drugsolver .and. ichemo >= DRUG_A) then
        chemo_active(k) = .false.
	endif	
enddo
neqn = k
! Note: neqn = 2*ncvars

!write(*,*) 'solver: nchemo,neqn: ',nchemo,neqn
!write(*,'(10f7.3)') C(1:neqn)
!write(*,'(a,3f8.5)') 'solver: metab1: ',Caverage(MAX_CHEMO+4:MAX_CHEMO+6)

ict = 1 ! for now just a single cell type

info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 1d-5
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
    if (.not.chemo_active(k)) cycle
    Caverage(ichemo) = C(k)
    do kcell = 1,nlist
        if (cell_list(kcell)%state == DEAD) cycle
        cell_list(kcell)%Cin(ichemo) = Caverage(ichemo)
    enddo
enddo
k = ncvars
do ic = 1,nchemo
    ichemo = chemomap(ic)
    k = k + 1
    if (.not.chemo_active(k)) cycle
    Caverage(MAX_CHEMO + ichemo) = C(k)
enddo

if (.not.use_drugsolver) return
if (chemo(DRUG_A)%present) then
	call DrugSolver(DRUG_A,tstart,dt,1,ok)
endif
if (chemo(DRUG_B)%present) then
	call DrugSolver(DRUG_B,tstart,dt,2,ok)
endif
end subroutine

!----------------------------------------------------------------------------------
! Radiation treatment is stored in protocol(0)
! NOT USED NOW
!----------------------------------------------------------------------------------
subroutine Treatment(radiation_dose)
real(REAL_KIND) :: radiation_dose
integer :: i, idrug, ichemo, nmetab, im	!, ichemo_metab

radiation_dose = 0
do i = 1,protocol(0)%n
	if (t_simulation >= protocol(0)%tstart(i) .and. .not.protocol(0)%started(i)) then
		radiation_dose = protocol(0)%dose(i)
		protocol(0)%started(i) = .true.
		protocol(0)%ended(i) = .true.
		write(nflog,*) 'Radiation started: dose: ',radiation_dose
		exit
	endif
enddo
do idrug = 1,Ndrugs_used
	ichemo = protocol(idrug)%ichemo
	nmetab = drug(idrug)%nmetabolites
	do i = 1,protocol(idrug)%n
		if (i == 1 .and. t_simulation < protocol(idrug)%tstart(i)) then
			chemo(ichemo)%bdry_conc = 0
			do im = 1,nmetab
				chemo(ichemo + im)%bdry_conc = 0
			enddo
			exit
		endif
		if (t_simulation >= protocol(idrug)%tstart(i) .and. .not.protocol(idrug)%started(i)) then
			chemo(ichemo)%bdry_conc = protocol(idrug)%conc(i)
			protocol(idrug)%started(i) = .true.
			protocol(idrug)%ended(i) = .false.
			chemo(ichemo)%present = .true.
			call InitConcs(ichemo)
			call SetupMedium(ichemo)
			do im = 1,nmetab
				if (chemo(ichemo + im)%used) then
					chemo(ichemo + im)%present = .true.
					call InitConcs(ichemo + im)
					call SetupMedium(ichemo + im)
				endif
			enddo
			write(nflog,*) 'Started DRUG: ',chemo(ichemo)%name,chemo(ichemo)%bdry_conc, i
			exit
		endif
	enddo
	do i = 1,protocol(idrug)%n
		if (t_simulation >= protocol(idrug)%tend(i) .and. .not.protocol(idrug)%ended(i)) then
			chemo(ichemo)%bdry_conc = 0
			protocol(idrug)%ended(i) = .true.
			call InitConcs(ichemo)
			call SetupMedium(ichemo)
			do im = 1,nmetab
				chemo(ichemo + im)%bdry_conc = 0
				if (chemo(ichemo + im)%used) then
					call InitConcs(ichemo + im)
					call SetupMedium(ichemo + im)
				endif
			enddo
			write(nflog,*) 'Ended DRUG: ',chemo(ichemo)%name,i
			exit
		endif
	enddo
enddo	
end subroutine

!call getBlobCentreRange(cntr,rng,radius) 
!
!if (axis == 1) then
!	x1 = cntr(1) - rng(1)/2 
!	x2 = cntr(1) + rng(1)/2
!	y0 = cntr(2)
!	z0 = cntr(3)
!	ix1 = x1/dx + 1
!	ix2 = x2/dx + 1
!	iy0 = y0/dx + 1
!	iz0 = z0/dx + 1
!	nxgpts = ix2 - ix1 + 2
!	nvars = 1 + MAX_CHEMO + N_EXTRA
!	allocate(ngc(NX))
!	allocate(ctemp(nxgpts,2,2,0:nvars-1))
!
!	ichemo = OXYGEN
!	write(nflog,'(a,40f8.5)') 'EC x O2: ',(Caverage(ix,iy0,iz0,ichemo),ix=ix1,ix2+1)
!
!	! First calculate averages at the enclosing grid pts
!	ctemp = 0
!	ngc = 0
!	do kx = 1,nxgpts
!		do ky = 1,2
!			do kz = 1,2
!				ix = ix1 + kx - 1
!				iy = iy0 + ky - 1
!				iz = iz0 + kz - 1
!				! Find the cells associated with this gridpt to compute cell-related variables
!				call get_gridptcells(ix,iy,iz,ngc(ix),gcell)
!				if (ngc(ix) > 0) then
!					do k = 1,ngc(ix)
!						kcell = gcell(k)
!						cp => cell_list(kcell)
!						ctemp(kx,ky,kz,0) = ctemp(kx,ky,kz,0) + cp%CFSE
!						ctemp(kx,ky,kz,GROWTH_RATE) = ctemp(kx,ky,kz,GROWTH_RATE) + cp%dVdt
!						ctemp(kx,ky,kz,CELL_VOLUME) = ctemp(kx,ky,kz,CELL_VOLUME) + cp%V
!						ctemp(kx,ky,kz,O2_BY_VOL) = ctemp(kx,ky,kz,O2_BY_VOL) + cp%Cin(OXYGEN)*cp%V
!					enddo
!					ctemp(kx,ky,kz,0:nvars-1) = ctemp(kx,ky,kz,0:nvars-1)/ngc(ix)
!				endif
!				! For constituents simply use the concentration fields
!				do ichemo = 1,MAX_CHEMO
!					ctemp(kx,ky,kz,ichemo) = Caverage(ix,iy,iz,ichemo)
!				enddo
!			enddo
!		enddo
!	enddo
!
!	dxc = (x2-x1)/(ns-1)
!
!	! Now need to interpolate values at ns points along the line from (x1,y0,z0) to (x2,y0,z0)
!	! alfax, alfay, alfaz each lies in (0,1)
!	y = y0
!	z = z0
!	iy = y/dx + 1
!	iz = z/dx + 1
!	alfay = (y - (iy-1)*dx)/dx
!	alfaz = (z - (iz-1)*dx)/dx
!	do ks = 1,ns
!		x = x1 + (ks-1)*dxc
!		ix = x/dx + 1
!		alfax = (x - (ix-1)*dx)/dx
!		kx = ix - ix1 + 1
!	!	if (kx >= nxgpts) then
!	!		write(*,*) 'bad kx: ',kx,nxgpts
!	!		write(*,'(a,i4,4f8.4,2i4)') 'ks,x,x1,x2,dxc,ix,kx: ',ks,x,x1,x2,dxc,ix,kx
!	!		stop
!	!	endif
!		do ichemo = 0, nvars-1
!			offset = ichemo*ns
!			k = offset - 1 + ks
!			ex_conc(k) = (1-alfax)*(1-alfay)*(1-alfaz)*ctemp(kx,1,1,ichemo) + &
!						 (1-alfax)*alfay*(1-alfaz)*ctemp(kx,2,1,ichemo) + &
!						 (1-alfax)*(1-alfay)*alfaz*ctemp(kx,1,2,ichemo) + &
!						 (1-alfax)*alfay*alfaz*ctemp(kx,2,2,ichemo) + &
!						 alfax*(1-alfay)*(1-alfaz)*ctemp(kx+1,1,1,ichemo) + &
!						 alfax*alfay*(1-alfaz)*ctemp(kx+1,2,1,ichemo) + &
!						 alfax*(1-alfay)*alfaz*ctemp(kx+1,1,2,ichemo) + &
!						 alfax*alfay*alfaz*ctemp(kx+1,2,2,ichemo)
!		enddo
!	enddo
!	!write(nflog,'(10f8.3)') ex_conc(0:nvars*ns-1)
!	deallocate(ctemp)
!	deallocate(ngc)
!elseif (axis == 2) then
!	y1 = cntr(2) - rng(2)/2
!	y2 = cntr(2) + rng(2)/2
!	x0 = cntr(1)
!	z0 = cntr(3)
!	iy1 = y1/dx + 1
!	iy2 = y2/dx + 1
!	ix0 = x0/dx + 1
!	iz0 = z0/dx + 1
!	nxgpts = iy2 - iy1 + 2
!	nvars = 1 + MAX_CHEMO + N_EXTRA
!	allocate(ngc(NX))
!	allocate(ctemp(2,nxgpts,2,0:nvars-1))
!
!	ichemo = OXYGEN
!	write(nflog,'(a,40f8.5)') 'EC y O2: ',(Caverage(ix0,iy,iz0,ichemo),iy=iy1,iy2+1)
!
!	! First calculate averages at the enclosing grid pts
!	ctemp = 0
!	ngc = 0
!	do ky = 1,nxgpts
!		do kx = 1,2
!			do kz = 1,2
!				iy = iy1 + ky - 1
!				ix = ix0 + kx - 1
!				iz = iz0 + kz - 1
!				! Find the cells associated with this gridpt to compute cell-related variables
!				call get_gridptcells(ix,iy,iz,ngc(iy),gcell)
!				if (ngc(ix) > 0) then
!					do k = 1,ngc(iy)
!						kcell = gcell(k)
!						cp => cell_list(kcell)
!						ctemp(kx,ky,kz,0) = ctemp(kx,ky,kz,0) + cp%CFSE
!						ctemp(kx,ky,kz,GROWTH_RATE) = ctemp(kx,ky,kz,GROWTH_RATE) + cp%dVdt
!						ctemp(kx,ky,kz,CELL_VOLUME) = ctemp(kx,ky,kz,CELL_VOLUME) + cp%V
!						ctemp(kx,ky,kz,O2_BY_VOL) = ctemp(kx,ky,kz,O2_BY_VOL) + cp%Cin(OXYGEN)*cp%V
!					enddo
!					ctemp(kx,ky,kz,0:nvars-1) = ctemp(kx,ky,kz,0:nvars-1)/ngc(ix)
!				endif
!				! For constituents simply use the concentration fields
!				do ichemo = 1,MAX_CHEMO
!					ctemp(kx,ky,kz,ichemo) = Caverage(ix,iy,iz,ichemo)
!				enddo
!			enddo
!		enddo
!	enddo
!
!	dxc = (y2-y1)/(ns-1)
!
!	! Now need to interpolate values at ns points along the line from (x0,y1,z0) to (x0,y2,z0)
!	! alfax, alfay, alfaz each lies in (0,1)
!	x = x0
!	z = z0
!	ix = x/dx + 1
!	iz = z/dx + 1
!	alfax = (x - (ix-1)*dx)/dx
!	alfaz = (z - (iz-1)*dx)/dx
!	do ks = 1,ns
!		y = y1 + (ks-1)*dxc
!		iy = y/dx + 1
!		alfay = (y - (iy-1)*dx)/dx
!		ky = iy - iy1 + 1
!		do ichemo = 0, nvars-1
!			offset = ichemo*ns
!			k = offset - 1 + ks
!			ex_conc(k) = (1-alfay)*(1-alfax)*(1-alfaz)*ctemp(1,ky,1,ichemo) + &
!						 (1-alfay)*alfax*(1-alfaz)*ctemp(2,ky,1,ichemo) + &
!						 (1-alfay)*(1-alfax)*alfaz*ctemp(1,ky,2,ichemo) + &
!						 (1-alfay)*alfax*alfaz*ctemp(2,ky,2,ichemo) + &
!						 alfay*(1-alfax)*(1-alfaz)*ctemp(1,ky+1,1,ichemo) + &
!						 alfay*alfax*(1-alfaz)*ctemp(2,ky+1,1,ichemo) + &
!						 alfay*(1-alfax)*alfaz*ctemp(1,ky+1,2,ichemo) + &
!						 alfay*alfax*alfaz*ctemp(2,ky+1,2,ichemo)
!		enddo
!	enddo
!	deallocate(ctemp)
!	deallocate(ngc)
!elseif (axis == 3) then
!	z1 = cntr(3) - rng(3)/2
!	z2 = cntr(3) + rng(3)/2
!	y0 = cntr(2)
!	x0 = cntr(1)
!	iz1 = z1/dx + 1
!	iz2 = z2/dx + 1
!	iy0 = y0/dx + 1
!	ix0 = x0/dx + 1
!	nxgpts = iz2 - iz1 + 2
!	nvars = 1 + MAX_CHEMO + N_EXTRA
!	allocate(ngc(NX))
!	allocate(ctemp(2,2,nxgpts,0:nvars-1))
!
!	ichemo = OXYGEN
!	write(nflog,'(a,40f8.5)') 'EC z O2: ',(Caverage(ix0,iy0,iz,ichemo),iz=iz1,iz2+1)
!
!	! First calculate averages at the enclosing grid pts
!	ctemp = 0
!	ngc = 0
!	do kz = 1,nxgpts
!		do ky = 1,2
!			do kx = 1,2
!				iz = iz1 + kz - 1
!				iy = iy0 + ky - 1
!				ix = ix0 + kx - 1
!				! Find the cells associated with this gridpt to compute cell-related variables
!				call get_gridptcells(ix,iy,iz,ngc(iz),gcell)
!				if (ngc(iz) > 0) then
!					do k = 1,ngc(iz)
!						kcell = gcell(k)
!						cp => cell_list(kcell)
!						ctemp(kx,ky,kz,0) = ctemp(kx,ky,kz,0) + cp%CFSE
!						ctemp(kx,ky,kz,GROWTH_RATE) = ctemp(kx,ky,kz,GROWTH_RATE) + cp%dVdt
!						ctemp(kx,ky,kz,CELL_VOLUME) = ctemp(kx,ky,kz,CELL_VOLUME) + cp%V
!						ctemp(kx,ky,kz,O2_BY_VOL) = ctemp(kx,ky,kz,O2_BY_VOL) + cp%Cin(OXYGEN)*cp%V
!					enddo
!					ctemp(kx,ky,kz,0:nvars-1) = ctemp(kx,ky,kz,0:nvars-1)/ngc(ix)
!				endif
!				! For constituents simply use the concentration fields
!				do ichemo = 1,MAX_CHEMO
!					ctemp(kx,ky,kz,ichemo) = Caverage(ix,iy,iz,ichemo)
!				enddo
!			enddo
!		enddo
!	enddo
!
!	dxc = (z2-z1)/(ns-1)
!
!	! Now need to interpolate values at ns points along the line from (x0,y0,z1) to (x0,y0,z2)
!	! alfax, alfay, alfaz each lies in (0,1)
!	y = y0
!	x = x0
!	iy = y/dx + 1
!	ix = x/dx + 1
!	alfay = (y - (iy-1)*dx)/dx
!	alfax = (x - (ix-1)*dx)/dx
!	do ks = 1,ns
!		z = z1 + (ks-1)*dxc
!		iz = z/dx + 1
!		alfaz = (z - (iz-1)*dx)/dx
!		kz = iz - iz1 + 1
!		do ichemo = 0, nvars-1
!			offset = ichemo*ns
!			k = offset - 1 + ks
!			ex_conc(k) = (1-alfaz)*(1-alfay)*(1-alfax)*ctemp(1,1,kz,ichemo) + &
!						 (1-alfaz)*alfay*(1-alfax)*ctemp(2,1,kz,ichemo) + &
!						 (1-alfaz)*(1-alfay)*alfax*ctemp(1,2,kz,ichemo) + &
!						 (1-alfaz)*alfay*alfax*ctemp(2,2,kz,ichemo) + &
!						 alfaz*(1-alfay)*(1-alfax)*ctemp(1,1,kz+1,ichemo) + &
!						 alfaz*alfay*(1-alfax)*ctemp(2,1,kz+1,ichemo) + &
!						 alfaz*(1-alfay)*alfax*ctemp(1,2,kz+1,ichemo) + &
!						 alfaz*alfay*alfax*ctemp(2,2,kz+1,ichemo)
!		enddo
!	enddo
!	deallocate(ctemp)
!	deallocate(ngc)
!endif