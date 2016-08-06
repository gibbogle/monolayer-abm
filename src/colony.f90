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
! The new method simply continues the simulation from where it ended, for 10 days.
!---------------------------------------------------------------------------------------------------
subroutine make_colony_distribution(dist,ddist,ndist) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: make_colony_distribution
use, intrinsic :: iso_c_binding
real(c_double) :: dist(*), ddist
integer(c_int) :: ndist
real(REAL_KIND) :: V0, dVdt, dt, t, tend
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
dt = DELTA_T

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

end module
