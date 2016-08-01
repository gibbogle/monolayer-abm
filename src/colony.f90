! To determine distribution of colony size

module colony

use global
use cellstate
implicit none

integer, parameter :: n_colony_days=10
integer :: nmax

contains

!---------------------------------------------------------------------------------------------------
! Simulate fate of cells grown with no nutrient constraints.
! The only determinants of the colony size for a cell are (considering radiation only):
! volume
! divide_time_mean(ityp)
! radiation_tag
! p_rad_death
! growth_delay
! G2_M
! dt_delay
! t_growth_delay_end
! N_delayed_cycles_left
! For simplicity, assume that each cell has the mean growth rate, and that for a cell that is not
! at the G2/M checkpoint, the fraction of time through the growth period is V0/divide_volume, where
! V0 = volume*Vcell_cm3
! and
! divide_volume = Vdivide0 + dVdivide*(2*R-1)
! The new method simply continues the simulation from where it ended, for 10 days.
!---------------------------------------------------------------------------------------------------
subroutine make_colony_distribution(dist,ddist,ndist) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: make_colony_distribution
use, intrinsic :: iso_c_binding
real(c_double) :: dist(*), ddist
integer(c_int) :: ndist
real(REAL_KIND) :: V0, dVdt, dt, t, tend
!integer, parameter :: ndist=40
!real(REAL_KIND) :: ddist, dist(ndist)
real(REAL_KIND) :: tnow_save
integer :: kcell, ityp, n, idist, ncycmax, ntot, nlist_save
type (cell_type), pointer :: cp

write(logmsg,*) 'make_colony_distribution: nlist: ',nlist
call logger(logmsg)
colony_simulation = .true.
nlist_save = nlist
tnow_save = tnow
ncycmax = 24*3600*n_colony_days/divide_time_mean(1) + 3
nmax = 2**ncycmax
allocate(ccell_list(nmax))
ddist = 50
dist(1:ndist) = 0
ntot = 0
tend = tnow + n_colony_days*24*3600    ! plate for 10 days
do kcell = 1, nlist_save
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ityp = cp%celltype
	
	! Now simulate colony growth from a single cell
	tnow = tnow_save
	call make_colony(kcell,tend,n)
	if (mod(kcell,100) == 0) then
	    write(logmsg,*) 'cell: n: ',kcell,n
	    call logger(logmsg)
	endif
	ntot = ntot + n
	idist = n/ddist + 1
	dist(idist) = dist(idist) + 1
enddo 
dist(1:ndist) = dist(1:ndist)/sum(dist(1:ndist))
write(logmsg,'(a,2i8,f8.1)') 'Colony size distribution: ', nlist_save,ntot,real(ntot)/nlist_save
call logger(logmsg)
write(nfout,'(a,2i8,f8.1)') 'Colony size distribution: ', nlist_save,ntot,real(ntot)/nlist_save
do idist = 1,ndist
	write(logmsg,'(i4,a,i4,f6.3)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
    call logger(logmsg)
	write(nfout,'(i4,a,i4,f6.3)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
enddo
deallocate(ccell_list)
colony_simulation = .false.
nlist = nlist_save
tnow = tnow_save
end subroutine

!!---------------------------------------------------------------------------------------------------
!!---------------------------------------------------------------------------------------------------
!subroutine simulate_cell_colony(kcell, n) bind(C)
!type (cell_type), pointer :: cp
!cp => cell_list(kcell)
!if (cp%state == DEAD) cycle
!ityp = cp%celltype
!
!! Now simulate colony growth from a single cell
!tnow = tnow_save
!call make_colony(kcell,tend,n)
!end subroutine
!
!!---------------------------------------------------------------------------------------------------
!subroutine close_make_colony_distribution() bind(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: close_make_colony_distribution
!
!end subroutine

!---------------------------------------------------------------------------------------------------
! The cell is at the point of division - possibly G2_M (arrested at G2/M checkpoint)
! For now only radiation tagging is handled
! Growth rate dVdt (mean) is used only to estimate the time of next division
!---------------------------------------------------------------------------------------------------
subroutine make_colony(kcell,tend,n)
integer :: kcell, n
real(REAL_KIND) :: tend, dt 
integer :: icell, ityp, nlist0, kpar=0
real(REAL_KIND) :: V0, Tdiv0, r_mean, c_rate, dVdt, Tmean, R
logical :: changed, ok
type (cell_type), pointer :: cp

!write(*,'(a,i6,2f8.0)') 'make_colony: ',kcell,tnow,tend
ccell_list(1) = cell_list(kcell)
ccell_list(1)%anoxia_tag = .false.
ccell_list(1)%aglucosia_tag = .false.
ccell_list(1)%drug_tag = .false.
!dt = 3600	! use big time steps, 1h
dt = DELTA_T
ityp = cell_list(kcell)%celltype
Tdiv0 = divide_time_mean(ityp)
r_mean = Vdivide0/Tdiv0
c_rate = log(2.0)/Tdiv0
!ccell_list(1)%V = cell_list(kcell)%divide_volume
!ccell_list(1)%t_divide_next = tnow

nlist = 1
do while (tnow < tend)
	tnow = tnow + dt
    call grower(dt,changed,ok)
enddo
n = 0
do icell = 1,nlist
	if (ccell_list(icell)%state /= DEAD) n = n+1
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine CloneColonyCell(kcell0,kcell1,ok)
integer :: kcell0, kcell1, ityp, idrug
logical :: ok
integer :: kpar = 0
real(REAL_KIND) :: R, Tdiv, Tdiv0

ok = .true.
ityp = ccell_list(kcell0)%celltype
Tdiv0 = divide_time_mean(ityp)
ccell_list(kcell0)%V = ccell_list(kcell0)%V/2
ccell_list(kcell0)%generation = ccell_list(kcell0)%generation + 1
if (ccell_list(kcell0)%growth_delay) then
	ccell_list(kcell0)%N_delayed_cycles_left = ccell_list(kcell0)%N_delayed_cycles_left - 1
	ccell_list(kcell0)%growth_delay = (ccell_list(kcell0)%N_delayed_cycles_left > 0)
endif
ccell_list(kcell0)%G2_M = .false.
ccell_list(kcell0)%t_divide_last = tnow
if (use_divide_time_distribution) then
	Tdiv = DivideTime(ityp)
	ccell_list(kcell0)%t_divide_next = tnow + Tdiv
else
	R = par_uni(kpar)
	ccell_list(kcell0)%t_divide_next = tnow + Tdiv0*(1 + (2*dVdivide/Vdivide0)*(2*R-1))
endif
!write(*,*) 'divide time: ',Tdiv0*(1 + (2*dVdivide/Vdivide0)*(2*R-1))/3600
ccell_list(kcell1)%celltype = ccell_list(kcell0)%celltype
ccell_list(kcell1)%state = ccell_list(kcell0)%state
ccell_list(kcell1)%V = ccell_list(kcell0)%V
ccell_list(kcell1)%generation = ccell_list(kcell0)%generation
ccell_list(kcell1)%ID = ccell_list(kcell0)%ID
ccell_list(kcell1)%p_rad_death = ccell_list(kcell0)%p_rad_death
ccell_list(kcell1)%p_drug_death = ccell_list(kcell0)%p_drug_death
ccell_list(kcell1)%radiation_tag = ccell_list(kcell0)%radiation_tag
ccell_list(kcell1)%anoxia_tag = .false.
ccell_list(kcell1)%aglucosia_tag = .false.
!ccell_list(kcell1)%exists = .true.
ccell_list(kcell1)%active = .true.
ccell_list(kcell1)%growth_delay = ccell_list(kcell0)%growth_delay
if (ccell_list(kcell1)%growth_delay) then
	ccell_list(kcell1)%dt_delay = ccell_list(kcell0)%dt_delay
	ccell_list(kcell1)%N_delayed_cycles_left = ccell_list(kcell0)%N_delayed_cycles_left
endif
ccell_list(kcell1)%G2_M = .false.
ccell_list(kcell1)%t_divide_last = tnow
if (use_divide_time_distribution) then
	Tdiv = DivideTime(ityp)
	ccell_list(kcell1)%t_divide_next = tnow + Tdiv
else
	R = par_uni(kpar)
	ccell_list(kcell1)%t_divide_next = tnow + Tdiv0*(1 + (2*dVdivide/Vdivide0)*(2*R-1))
endif
ccell_list(kcell1)%t_anoxia = 0
ccell_list(kcell1)%t_aglucosia = 0
ccell_list(kcell1)%Cin = ccell_list(kcell0)%Cin
!ccell_list(kcell1)%Cex = ccell_list(kcell0)%Cex
ccell_list(kcell1)%dCdt = ccell_list(kcell0)%dCdt
ccell_list(kcell1)%dMdt = ccell_list(kcell0)%dMdt
ccell_list(kcell1)%M = ccell_list(kcell0)%M
end subroutine


end module
