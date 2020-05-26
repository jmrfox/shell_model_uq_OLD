!
!  Program STANDARDORDERINT
!
!  reorders BIGSTICK interaction files into a standard order
!  a <= b <= c <= c, increasing J, increasing T
!
!  by CWJ @ SDSU Dec 2015
!
!  VERSION 2 writes out/reads in a vector
!
!  added in July 2017 @ INT: option to add in phase from l-s ordering
!  that is, a phase we pick up for (l x s )j  vs (s x l) j

module information
    implicit none
	integer :: norbiti,norbitf  ! number of orbits in initial and final spaces
	                            ! must have norbiti >= norbitf
!------------ CREATE A DEFINED TYPE-----------------------------------------
	type orb
		integer :: nr            ! radial quantum number
		integer :: j             ! 2 x j
		integer :: l             ! L
		integer :: par           ! parity = +/- 1
		integer :: map           ! maps to destination
		integer :: w             ! weighting; needed for adding in single particle energies
		logical :: core          ! if in "core"
    end type orb
	type (orb),allocatable,target :: orbqni(:), orbqnf(:) ! initial, final orbqn(iorb); assume single species

	real, allocatable :: spe(:)
	integer :: ntbmei,ntbmef  ! initial, final # of TBMEs
    integer, allocatable :: qlist(:,:)   ! labels for each TBME
	real, allocatable :: vtbme(:)
	logical :: inducedspes  ! whether or not to created induced spes from sum over "core"
	real :: ecore
	logical :: lsphase
	integer, allocatable :: lsphaser(:)

end module information

!================= MAIN PROGRAM ==============================================
use information
implicit none
character menu_char
integer nme

print*,' '
print*,' Program to reorder .int file into standard order   '
print*,' '
print*,' Do you want to... '
print*,' (o) order .int file into another .int file? '
print*,' (w) write .int file as a vector? '
print*,' (r) read vector and save as a .int file?'
read(5,'(a)')menu_char

call open_sps_file('i')
call standardize(.true.,nme)
call standardize(.false.,nme)

!call mapspaces
select case (menu_char)
case('o','O')
call makelsphase(.true.)
call openhamiltonian('b')
allocate(vtbme(nme))
vtbme = 0.0
call readhamiltonian(nme)
ntbmei = nme
call writehamiltonian

case('w','W')
call makelsphase(.false.)

call openhamiltonian('i')
allocate(vtbme(nme))
vtbme = 0.0
call readhamiltonian(nme)
ntbmei = nme
call writestdvector

case('r','R')
call makelsphase(.false.)

ntbmei = nme
call openhamiltonian('f')
allocate(vtbme(nme))
allocate(spe(norbiti))
call readstdvector
call writehamiltonian

end select


end

!================= END MAIN PROGRAM ==========================================

subroutine open_sps_file(ichar)
	
	use information
	implicit none
	character :: ichar  ! determines if initial or final space
	
	type (orb),pointer :: curorbqn(:)
	integer :: iorb,norb
    character*20 :: filename
	character*3 :: dummy
	integer  :: ilast
	logical  :: success
	integer :: ierr
	real    :: xn,xl,xj
	integer :: w
	
	
	success = .false.
	
	do while(.not.success)
	
  	   if(ichar=='i')then
	       print*,' Enter name of INITIAL .sps file (leave off extension )'	 
	   else
	       print*,' Enter name of FINAL .sps file (leave off extension )'	 
  	   end if
	   read(5,'(a)')filename
	   ilast = index(filename,' ')-1
	   open(unit=1,file=filename(1:ilast)//".sps",status='old',iostat=ierr)
	   if(ierr/=0)then
		   print*,filename(1:ilast)//".sps",' does not exist '
	   else
		   success=.true.
	   end if
    end do
   	read(1,*)dummy   ! SHOULD BE 'iso'
	read(1,*)norb
	if(ichar=='i')then
		allocate(orbqni(norb))
		curorbqn=> orbqni
		norbiti = norb
	else
		allocate(orbqnf(norb))
		curorbqn=> orbqnf
		norbitf = norb
	end if
	do iorb =1,norb
		read(1,*)xn,xl,xj,w
		curorbqn(iorb)%nr = nint(xn)
		curorbqn(iorb)%l  =  nint(xl)
		curorbqn(iorb)%j = nint(xj*2)
        curorbqn(iorb)%map = 0         ! initialize
		curorbqn(iorb)%w   = w        ! used to determine core
		curorbqn(iorb)%core = .false.   ! not always used
	end do
	close(1)	

	return
	
end subroutine open_sps_file


!=============================================================================
!
! create standard ordering of TBMES
!
! CALLED BY: main routine
! 
subroutine standardize(countup,nme)
	use information
	implicit none
	logical :: countup
	
	integer :: a,b,c,d,J,T
	integer :: nme,ime
	integer :: dstart
	integer :: jmin,jmax
	logical :: same
	integer :: tstart,tstop
	
	
	ime = 0
	print*,norbiti
	do a = 1,norbiti
		do b = a,norbiti
			do c = a,norbiti
				if(a == c)then
					dstart = b
				else
					dstart = c
				end if
				do d = dstart,norbiti
					same = .false.
					if(a==b .or. c==d)same =.true.
					jmin = abs(orbqni(a)%j - orbqni(b)%j)/2
					jmin = max(jmin, abs(orbqni(c)%j - orbqni(d)%j)/2    )
					jmax= (orbqni(a)%j + orbqni(b)%j)/2
					jmax = min(jmax, (orbqni(c)%j + orbqni(d)%j)/2    )
					do j = jmin,jmax
						if(.not.same)then
							tstart = 0
							tstop  = 1
						else
							tstart = mod(J+1,2)
							tstop = tstart
						
						end if
						do t = tstart,tstop
							ime = ime + 1
							
							if(.not.countup)then
								qlist(ime,1)=a
								qlist(ime,2)=b
								qlist(ime,3)=c
								qlist(ime,4)=d
								qlist(ime,5)=j
								qlist(ime,6)=t
								
								
							end if
							
						end do
					end do
					
					
				end do
				
				
			end do ! c
				
		
		
		end do ! b
	end do ! a
	if(countup)then
		nme = ime

		print*,' there are ',nme,' matrix elements'
		allocate(qlist(nme,6))
	end if
	return
	
	
end subroutine standardize

!=============================================================================
subroutine find_index_and_phase(nme, nlist,j,t,indx,phase)
	use information
	implicit none
	integer :: nme
	integer:: nlist(4)
	integer:: j,t,indx,phase
	integer :: a,b,c,d,x,y
	integer :: i
	logical :: success
	
	phase = 1
	if (nlist(1) > nlist(2))then
		a = nlist(2)
		b = nlist(1)
		phase = phase * (-1)**(j + t + (orbqni(a)%j + orbqni(b)%j)/2)
		
	else
		a = nlist(1)
		b = nlist(2)
		
	end if
	if (nlist(3) > nlist(4))then
		c = nlist(4)
		d = nlist(3)
		phase = phase * (-1)**(j + t + (orbqni(c)%j + orbqni(d)%j)/2)
		
	else
		c = nlist(3)
		d = nlist(4)
		
	end if
	
	if( c < a .or. (c == a .and. d < b))then
		x = a
		y = b
		a = c
		b = d
		c = x
		d = y
		
	end if
	
	success = .false.
	do i = 1,nme
		if(a /= qlist(i,1))cycle
		if(b /= qlist(i,2))cycle
		if(c /= qlist(i,3))cycle
		if(d /= qlist(i,4))cycle
		if(j /= qlist(i,5))cycle
		if(t==qlist(i,6))then
			success = .true.
			indx = i
			exit
		end if
	end do
	
!....... ADD IN LS PHASE...........	
	
    phase=phase*lsphaser(a)*lsphaser(b)*lsphaser(c)*lsphaser(d)
	if(.not.success)then
		print*,' could not find it ! '
		print*,a,b,c,d,j,t
		stop
	end if
	return
	
	
end subroutine find_index_and_phase

!=============================================================================
!
!  open Hamiltonian file; check for compatibility
!

subroutine openhamiltonian(filetype)
	implicit none
	character    :: filetype
    character*60 :: filename
	character*3 :: dummy
	integer  :: ilast
	logical  :: success
	integer :: ierr
	character :: ychar
	
	success=.false.
	if(filetype == 'b' .or. filetype == 'i')then
	do while(.not.success)
		print*,' Enter name of initial .int interaction file '
		read(5,'(a)')filename
		ilast = index(filename,' ')-1
		open(unit=3,file=filename(1:ilast)//".int",status='old',iostat=ierr)
		if(ierr/=0)then
			print*,filename(1:ilast)//".int",' does not exist '
		else
			success=.true.
		end if
	end do
        print*,' successfully opened '
	end if

	if(filetype == 'b' .or. filetype == 'f')then
		
		success=.false.
		do while(.not.success)
			print*,' Enter name of final .int interaction file '
			read(5,'(a)')filename
			ilast = index(filename,' ')-1
			open(unit=33,file=filename(1:ilast)//".int",status='new',iostat=ierr)
			if(ierr/=0)then
				print*,filename(1:ilast)//".int",' already exists; overwrite (y/n)? '
				read(5,'(a)')ychar
				if(ychar=='y' .or. ychar=='Y')then
					open(unit=33,file=filename(1:ilast)//".int",status='old',iostat=ierr)
					success=.true.
				end if
			else
				success=.true.
			end if
		end do	
	end if
		
	return
	
end subroutine openhamiltonian

!=============================================================================
!
!  count up how many TBMES are in the new Hamiltonian
!
subroutine readhamiltonian(nme)
	use information
	implicit none
	character*70 :: dummy
	logical :: endofheader
	integer :: i,n,jj,tt
	real    :: v
	integer(4) :: nlab(4)
	logical :: keep
	integer :: delta
	integer :: k
	integer :: iphase,indx
	integer:: nme
    print*,' '
!--------- READ PAST HEADER-------	
	endofheader=.false.
	do while(.not.endofheader)
		read(3,'(a)')dummy
		if(dummy(1:1)/='!' .and. dummy(1:1)/='#')then
			backspace(3)
			endofheader=.true.
		else
			print*,dummy
			write(33,'(a70)')dummy
		end if
	end do
	allocate(spe(norbiti))
	read(3,*)ntbmei,(spe(i),i=1,min(norbiti,10))
	
	if(norbiti>10)then
		k = 10
		do while(k < norbiti)
			read(3,*)(spe(k+i),i=1,min(norbiti-k,10))
			k = k+10
			
		end do
		
	end if
	
	ntbmef = 0

!.... LOOP OVER MATRIX ELEMENTS
	do i = 1,ntbmei
		read(3,*)(nlab(n),n=1,4),jj,tt,v
!----------- CHECK CONSISTENCY ----------------		
        call find_index_and_phase(nme,nlab,jj,tt,indx,iphase)
		vtbme(indx)=v*iphase
 

	end do

	close(3)
	return
end subroutine readhamiltonian

!=============================================================================

subroutine writehamiltonian
	use information
	implicit none
    character*70 :: filename
	character*1 :: ychar
	integer  :: ilast
	logical  :: success
	integer :: ierr
	integer :: iorb
	integer :: i,j
    integer nline,nres
	
	
	write(33,'("!SP ORB  N  L  J" )')
	do iorb = 1, norbiti
		write(33,'("!SP ",4i3,"/2")')iorb,orbqni(iorb)%nr,orbqni(iorb)%l,orbqni(iorb)%j
	end do  ! iorbf 

	write(33,'(i8,10f9.4)')ntbmei,(spe(j),j=1,min(10,norbiti))
	if(norbiti > 10)then
        nline = norbiti/10
        nres = norbiti-nline*10
		print*,nline,nres
        do j = 2,nline
               write(33,'(10f9.4)')(spe((i+(j-1)*10)),i=1,10)
        end do
        if(nres > 0) write(33,'(10f9.4)')(spe((i+(nline)*10)),i=1,nres)
		print*, ' Testing sum rule ', nres+nline*10,norbiti
		
	end if
	do i =1,ntbmei
		write(33,'(4i4,2x,2i4,2x,f10.5)')(qlist(i,j),j=1,6),vtbme(i)
	end do
			
	close(33)
	
	return
end subroutine writehamiltonian
!=============================================================================
subroutine writestdvector
	use information
	implicit none
    character*60 :: filename
	integer  :: ilast
	logical  :: success,foundit
	integer :: ierr
	integer :: i
	
	success=.false.
	do while(.not.success)
		print*,' Enter entire name of output vector file '
		read(5,'(a)')filename
		ilast = index(filename,' ')
		inquire(file=filename(1:ilast),exist=foundit)
		if(foundit)then
			print*,' That file already exists ',filename(1:ilast)
			cycle
		else
			success = .true.
		end if
		open(unit=44,file=filename(1:ilast),status='new',iostat=ierr)
	end do
    print*,' successfully opened '
	write(44,*)ntbmei+norbiti
	do i = 1,ntbmei
		write(44,*)vtbme(i)
	end do
	do i = 1,norbiti
		write(44,*)spe(i)
	end do
	close(44)
	return		
	
	
end subroutine writestdvector
!=============================================================================
subroutine readstdvector
	use information
	implicit none
    character*60 :: filename
	integer  :: ilast
	logical  :: success
	integer :: ierr
	integer :: i
	integer :: n
	
	success=.false.
	do while(.not.success)
		print*,' Enter entire name of input vector file '
		read(5,'(a)')filename
		ilast = index(filename,' ')
		inquire(file=filename(1:ilast),exist=success)
		if(.not.success)then
			print*,' That file does not exists ',filename(1:ilast)
			cycle
		end if
		open(unit=55,file=filename(1:ilast),status='old',iostat=ierr)
	end do
    print*,' successfully opened '
	read(55,*)n
	if(n/=ntbmei+norbiti)then
		print*,' mismatch in dimensions ',n,ntbmei+norbiti
		stop
	end if
	do i = 1,ntbmei
		read(55,*)vtbme(i)
	end do
	do i = 1,norbiti
		read(55,*)spe(i)
	end do
	close(55)
	return		
	
	
end subroutine readstdvector

!..........................
!

subroutine makelsphase(ask)
	use information
	implicit none
	character :: ychar
	integer :: iorb
	logical ask
	
	allocate(lsphaser(norbiti))
	lsphaser=1
	if(.not.ask)return
	
	print*,' Do you want to have a phase from L + S = J ordering (y/n)'
	read(5,'(a)')ychar
	if(ychar=='y' .or. ychar=='Y')then
		do  iorb = 1,norbiti
			lsphaser(iorb) = (-1)**( (orbqni(iorb)%l*2 + 1 - orbqni(iorb)%j )/2)
			
			
		end do
		
		
	end if
	
	return
end subroutine makelsphase
