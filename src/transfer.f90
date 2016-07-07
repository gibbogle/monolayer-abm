! Transferring results to the GUI

module transfer

use global
use chemokine
!use envelope
use ode_diffuse
use, intrinsic :: iso_c_binding

#include "../src/version.h"

implicit none

!integer, parameter :: X_AXIS = 1
!integer, parameter :: Y_AXIS = 2
!integer, parameter :: Z_AXIS = 3

type, bind(C) :: celldata_type
	integer(c_int) :: tag
	real(c_double) :: radius
	real(c_double) :: centre(3)
	integer(c_int) :: celltype
	integer(c_int) :: status
end type

type, bind(C) :: fielddata_type
    integer(c_int) :: NX, NY, NZ, NCONST
    real(c_double) :: DX
    type(c_ptr) :: Conc_ptr   ! Cslice(NX,NY,NZ,NCONST)
    integer(c_int) :: ncells
    type(c_ptr) :: cell_ptr
end type

!type(celldata_type) :: cdata(4000)

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_DLL_build_version(version_array,array_len) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: get_dll_build_version
use, intrinsic :: iso_c_binding
character(c_char) :: version_array(12)
integer(c_int) :: array_len
integer :: k

dll_version = DLL_BUILD_VERSION
gui_version = GUI_BUILD_VERSION
!write(nflog,*) 'get_DLL_build_version: ',dll_version
do k = 1,12
	version_array(k) = dll_version(k:k)
!	write(nflog,'(i2,a,a)') k,' ',version_array(k)
	if (version_array(k) == ' ') then
		version_array(k) = char(0)
		array_len = k
		exit
	endif
enddo
!write(nflog,*) 'array_len: ',array_len
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(nsteps_dim, deltat, maxchemo, nextra, cused) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: nsteps_dim, maxchemo, nextra
real(c_double) :: deltat
logical(c_bool) :: cused(*)
integer :: ichemo

write(nflog,*) 'get_dimensions: MAX_CHEMO: ',MAX_CHEMO
nsteps_dim = nsteps
deltat = DELTA_T
maxchemo = MAX_CHEMO
nextra = N_EXTRA
do ichemo = 1,MAX_CHEMO
	cused(ichemo+1) = chemo(ichemo)%used
enddo
cused(1) = .true.			! CFSE
cused(MAX_CHEMO+2) = .true.	! Growth rate
end subroutine


!-----------------------------------------------------------------------------------------
! Rendered cell colour may depend on stage, state, receptor expression level.
! col(:) = (r,g,b)
!-----------------------------------------------------------------------------------------
subroutine CellColour(kcell,highlight,col)
integer :: kcell, highlight, col(3)
integer :: stage, status
integer, parameter :: WHITE(3) = (/255,255,255/)
integer, parameter :: RED(3) = (/255,0,0/)
integer, parameter :: GREEN(3) = (/0,255,0/)
integer, parameter :: BLUE(3) = (/0,0,255/)
integer, parameter :: DEEPRED(3) = (/200,0,0/)
integer, parameter :: DEEPBLUE(3) = (/30,20,255/)
integer, parameter :: DEEPGREEN(3) = (/0,150,0/)
integer, parameter :: LIGHTRED(3) = (/255,70,90/)
integer, parameter :: LIGHTBLUE(3) = (/0,200,255/)
integer, parameter :: LIGHTGREEN(3) = (/50,255,150/)
integer, parameter :: DEEPORANGE(3) = (/240,70,0/)
integer, parameter :: LIGHTORANGE(3) = (/255,130,0/)
integer, parameter :: YELLOW(3) = (/255,255,0/)
integer, parameter :: DEEPPURPLE(3) = (/180,180,30/)
integer, parameter :: LIGHTPURPLE(3) = (/230,230,100/)
integer, parameter :: DEEPBROWN(3) = (/130,70,0/)
integer, parameter :: LIGHTBROWN(3) = (/200,100,0/)
integer, parameter :: GRAY(3) = (/128,128,128/)

integer, parameter :: Qt_white = 3
integer, parameter :: Qt_black = 2
integer, parameter :: Qt_red = 7
integer, parameter :: Qt_darkRed = 13
integer, parameter :: Qt_green = 8
integer, parameter :: Qt_darkGreen = 14
integer, parameter :: Qt_blue = 9
integer, parameter :: Qt_darkBlue = 15
integer, parameter :: Qt_cyan = 10
integer, parameter :: Qt_darkCyan = 16
integer, parameter :: Qt_magenta = 11
integer, parameter :: Qt_darkMagenta = 17
integer, parameter :: Qt_yellow = 12
integer, parameter :: Qt_darkYellow = 18
integer, parameter :: Qt_gray = 5
integer, parameter :: Qt_darkGray = 4
integer, parameter :: Qt_lightGray = 6

if (highlight == 0) then
    col = LIGHTORANGE
else
    col = LIGHTRED
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Pack the colours (r,g,b) into an integer.
!-----------------------------------------------------------------------------------------
integer function rgb(col)
integer :: col(3)

rgb = ishft(col(1),16) + ishft(col(2),8) + col(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getMediumConc(cmedium)
real(REAL_KIND) :: cmedium(:)

cmedium = Caverage(MAX_CHEMO+1:2*MAX_CHEMO)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getNecroticFraction(necrotic_fraction,totvol_cm3)
real(REAL_KIND) :: necrotic_fraction, totvol_cm3	! vol_cm3 not used here, needed in scell
real(REAL_KIND) :: cellvol_cm3, dvol
!necrotic_fraction = (Nsites-Ncells)/real(Nsites)
cellvol_cm3 = Ncells*DELTA_X**3
dvol = totvol_cm3-cellvol_cm3
necrotic_fraction = dvol/totvol_cm3
if (necrotic_fraction < 0.005) necrotic_fraction = 0
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), i_hypoxia_cutoff,i_growth_cutoff
integer :: Nviable(MAX_CELLTYPES), Nlive(MAX_CELLTYPES), plate_eff_10(MAX_CELLTYPES)
integer :: nhypoxic(3), nclonohypoxic(3), ngrowth(3), &
    hypoxic_percent_10, clonohypoxic_percent_10, growth_percent_10, necrotic_percent_10, &
    medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(2)
integer :: TNanoxia_dead, TNaglucosia_dead, TNradiation_dead, TNdrug_dead(2),  TNviable, &
           Ntagged_anoxia(MAX_CELLTYPES), Ntagged_aglucosia(MAX_CELLTYPES), Ntagged_radiation(MAX_CELLTYPES), &
           Ntagged_drug(2,MAX_CELLTYPES), &
           TNtagged_anoxia, TNtagged_aglucosia, TNtagged_radiation, TNtagged_drug(2)
integer :: Tplate_eff_10   
integer :: ityp
real(REAL_KIND) :: hour, plate_eff(MAX_CELLTYPES)
real(REAL_KIND) :: cmedium(MAX_CHEMO)

hour = istep*DELTA_T/3600.

Ntagged_anoxia(:) = Nanoxia_tag(:)			! number currently tagged by anoxia
Ntagged_aglucosia(:) = Naglucosia_tag(:)	! number currently tagged by aglucosia
Ntagged_radiation(:) = Nradiation_tag(:)	! number currently tagged by radiation
Ntagged_drug(1,:) = Ndrug_tag(1,:)			! number currently tagged by drugA
Ntagged_drug(2,:) = Ndrug_tag(2,:)			! number currently tagged by drugA

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
clonohypoxic_percent_10 = (1000.*nclonohypoxic(i_hypoxia_cutoff))/TNviable
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
call getMediumConc(cmedium)
medium_oxygen_1000 = cmedium(OXYGEN)*1000.
medium_glucose_1000 = cmedium(GLUCOSE)*1000.
medium_drug_1000(1) = cmedium(DRUG_A)*1000.
medium_drug_1000(2) = cmedium(DRUG_B)*1000.

summaryData(1:20) = [ istep, Ncells, TNanoxia_dead, TNaglucosia_dead, TNdrug_dead(1), TNdrug_dead(2), TNradiation_dead, &
    TNtagged_anoxia, TNtagged_aglucosia, TNtagged_drug(1), TNtagged_drug(2), TNtagged_radiation, &
	hypoxic_percent_10, clonohypoxic_percent_10, growth_percent_10, Tplate_eff_10, &
	medium_oxygen_1000, medium_glucose_1000, medium_drug_1000(1), medium_drug_1000(2) ]
write(nfres,'(a,a,2a12,i8,e12.4,22i7,15e12.4)') trim(header),' ',gui_run_version, dll_run_version, &
	istep, hour, Ncells_type(1:2), &
    Nanoxia_dead(1:2), Naglucosia_dead(1:2), Ndrug_dead(1,1:2), &
    Ndrug_dead(2,1:2), Nradiation_dead(1:2), &
    Ntagged_anoxia(1:2), Ntagged_aglucosia(1:2), Ntagged_drug(1,1:2), &
    Ntagged_drug(2,1:2), Ntagged_radiation(1:2), &
	nhypoxic(:)/real(Ncells), nclonohypoxic(:)/real(TNviable), ngrowth(:)/real(Ncells), &
	plate_eff(1:2), &
	cmedium(OXYGEN), cmedium(GLUCOSE), cmedium(DRUG_A), cmedium(DRUG_B)
		
!call sum_dMdt(GLUCOSE)

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getHypoxicCount(nhypoxic)
integer :: nhypoxic(3)
integer :: kcell, i

nhypoxic = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	do i = 1,3
		if (cell_list(kcell)%conc(OXYGEN) < O2cutoff(i)) nhypoxic(i) = nhypoxic(i) + 1
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getClonoHypoxicCount(nclonohypoxic)
integer :: nclonohypoxic(3)
integer :: kcell, i, idrug
logical :: tagged

nclonohypoxic = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%anoxia_tag) cycle
	if (cell_list(kcell)%aglucosia_tag) cycle
	if (cell_list(kcell)%radiation_tag) cycle
	tagged = .false.
	do idrug = 1,MAX_DRUGTYPES
	    if (cell_list(kcell)%drug_tag(idrug)) tagged = .true.
	enddo
	if (tagged) cycle
	do i = 1,3
		if (cell_list(kcell)%conc(OXYGEN) < O2cutoff(i)) then
			nclonohypoxic(i) = nclonohypoxic(i) + 1
		endif
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Need to compare growth rate with a fraction of average growth rate
!--------------------------------------------------------------------------------
subroutine getGrowthCount(ngrowth)
integer :: ngrowth(3)
integer :: kcell, i, ityp
real(REAL_KIND) :: r_mean(2)

r_mean(1:2) = Vdivide0/(2*divide_time_mean(1:2))
ngrowth = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	ityp = cell_list(kcell)%celltype
	do i = 1,3
		if (cell_list(kcell)%dVdt < growthcutoff(i)*r_mean(ityp)) ngrowth(i) = ngrowth(i) + 1
	enddo
enddo

end subroutine

!--------------------------------------------------------------------------------
! Compute total uptake rate for a constituent
!--------------------------------------------------------------------------------
subroutine sum_dMdt(ichemo)
integer :: ichemo
integer :: kcell, Nh, Nc
real(REAL_KIND) :: C, metab, dMdt, asum, msum, Csum

if (ichemo > GLUCOSE) then
	write(*,*) 'Error: sum_dMdt: only for oxygen and glucose'
	stop
endif
Nh = chemo(ichemo)%Hill_N
asum = 0
Csum = 0
msum = 0
Nc = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	Nc = Nc + 1
	C = cell_list(kcell)%conc(ichemo)
	Csum = Csum + C
	metab = C**Nh/(chemo(ichemo)%MM_C0**Nh + C**Nh)
	msum = msum + metab
	dMdt = metab*chemo(ichemo)%max_cell_rate 
	asum = asum + dMdt
enddo
total_dMdt = total_dMdt + asum
!write(*,'(a,2i6,2e12.3)') 'sum_dMdt: ',ichemo,Nc,asum,total_dMdt*3600 
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getNviable(Nviable, Nlive)
integer :: Nviable(:), Nlive(:)
integer :: kcell, ityp, idrug
logical :: tag

Nviable = 0
Nlive = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
    ityp = cell_list(kcell)%celltype
    Nlive(ityp) = Nlive(ityp) + 1
	if (cell_list(kcell)%anoxia_tag .or. cell_list(kcell)%aglucosia_tag .or. cell_list(kcell)%radiation_tag) cycle
    tag = .false.
    do idrug = 1,ndrugs_used
		if (cell_list(kcell)%drug_tag(idrug)) tag = .true.
	enddo
	if (tag) cycle
	Nviable(ityp) = Nviable(ityp) + 1
enddo	
end subroutine


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
function getstatus(cp) result(status)
type(cell_type), pointer :: cp
integer :: status

if (cp%anoxia_tag) then
	status = 2	! tagged to die of anoxia
elseif (cp%aglucosia_tag) then
	status = 4	! tagged to die of aglucosia
elseif (cp%radiation_tag) then
	status = 10
elseif (cp%drug_tag(1)) then
	status = 11
elseif (cp%drug_tag(1)) then
	status = 12
elseif (cp%conc(OXYGEN) < hypoxia_threshold) then
	status = 1	! radiobiological hypoxia
!elseif (cp%mitosis > 0) then
elseif (cp%volume > 0.9*cp%divide_volume) then  ! just a surrogate for mitosis
	status = 3	! in mitosis
else
	status = 0
endif
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_constituents(nvars,cvar_index,nvarlen,name_array,narraylen) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_constituents
use, intrinsic :: iso_c_binding
character(c_char) :: name_array(0:*)
integer(c_int) :: nvars, cvar_index(0:*), nvarlen, narraylen
integer :: ivar, k, ichemo
character*(24) :: name
character(c_char) :: c

write(nflog,*) 'get_constituents'
nvarlen = 24
ivar = 0
k = ivar*nvarlen
cvar_index(ivar) = 0	! CFSE
name = 'CFSE'
call copyname(name,name_array(k),nvarlen)
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	ivar = ivar + 1
	k = ivar*nvarlen
	cvar_index(ivar) = ichemo
	name = chemo(ichemo)%name
	write(nflog,*) 'get_constituents: ',ichemo,name
	call copyname(name,name_array(k),nvarlen)
enddo
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = GROWTH_RATE
name = 'Growth rate'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = CELL_VOLUME
name = 'Cell volume'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = O2_BY_VOL
name = 'Cell O2xVol'
call copyname(name,name_array(k),nvarlen)
nvars = ivar + 1
write(nflog,*) 'did get_constituents'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine copyname(name,name_array,n)
character*(*) :: name
character :: name_array(*)
integer :: n
integer :: k

do k = 1,n
	name_array(k) = name(k:k)
enddo
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_string(bufptr) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_string
use, intrinsic :: iso_c_binding
type(c_ptr) :: bufptr
character(c_char) :: buf(1024)
integer :: buflen
character*(1024), save :: string

string = 'A test string'
buflen = len(trim(string))
!write(*,*) 'buflen: ',buflen
string(buflen+1:buflen+1) = char(0)
bufptr = c_loc(string)
end subroutine


end module
