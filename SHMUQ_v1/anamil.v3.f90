!
!  code to analyze MILCOMs -- 
!  Most Important Linear Combinations Of Matrix elements
!
! initiated late April 2017 by CWJ @ SDSU
!
!
! version 2, March 2019 
!
! version 3, April 2019, allows for more general format
!
module milcom_data
	
	implicit none
	
	integer :: num_nuclides
	integer, allocatable :: list_nuclides(:,:)  ! list_nuclides(i,j) i = which nuclides, j = 1,2 for Z, N
	integer, allocatable :: list_levels(:,:)  ! list_nuclides(i,j) i = level, j = 1,2 for Z, N

	real, allocatable    :: list_energies(:,:)    ! list of energies(i,j) i = which nuclides, j = which energy
	real, allocatable    :: list_expect(:,:,:)   ! list of < H_k >   i,j,k  i  = which nuclide, j = which energy, k = which operator
	
	integer :: nlevels ! # of level analyzed
	integer :: nlevels_assumed   ! assumed # of levels per nuclide
	
	integer :: nops    ! # of operators, e.g., 66 in sd shell
	
	integer :: milfile  =55 ! unit number for file containing milcoms
	
	
	integer :: myzmin,myzmax,mynmin,mynmax
	
	integer :: nlevels_chosen
	integer :: ncases
	
	logical :: excited_only = .false.
	logical :: scale_me = .true.
	integer :: A0_scale = 18
	real    :: X_scale = 0.3
	integer :: nspes = 3
	
	
	real(8), allocatable :: B(:,:),A(:,:),lambda2(:),vec(:,:)
	
	
	
contains

!===============================================================================	
	
	subroutine open_milcom_file
		implicit none
		
		character*70 :: filename
		integer      :: ilast
		logical      :: success
		
		success = .false.
		do while (.not. success)
			write(6,*)' '
			write(6,*)' Enter full name of file with MILCOM data '
			read(5,'(a)')filename
			ilast = index(filename,' ')-1
			inquire(file=filename(1:ilast),exist=success)
			if(.not.success)then
				write(6,*)filename(1:ilast),' does not seem to exist '
				cycle
			end if
			open(unit=milfile,file=filename(1:ilast),status='old')
			
		end do
		print*,' unit opened '
		return
	end subroutine open_milcom_file
	
!===============================================================================	
	
	subroutine survey_milcom_file(fill)
		implicit none
		
		integer :: maxN,maxZ,minN,minZ
		integer :: Nnukes
		integer :: maxops,minops
		integer :: z,n,iop
		integer :: oldz,oldn
		integer :: ilevel,i
		real    :: e,xj,xt,xval
		logical :: endofroad,endoflocal
		character*50 :: lineread
		integer :: maxlevels,minlevels
		logical fill
		
		endofroad = .false.
		maxN = 0
		maxZ = 0
		minZ = 10000
		minN = 10000
		maxops = 0000
		minops = 100000
		maxlevels = 00000
		minlevels = 1000000
		Nnukes = 0
		oldz= -1
		oldn = -1
		
		if(.not.fill)then
			print*,' Enter the number of levels per nuclide (assuming uniform!)'
		    read*,nlevels_assumed
		end if
		do while (.not.endofroad)
			read(milfile,*,end=201,err=201)z,n,iop
			if(z/=oldz .or. n /=oldn)then
				nnukes = nnukes+1
				oldz = z
				oldn = n
				if(fill)then
					list_nuclides(nnukes,1)=z
					list_nuclides(nnukes,2)=n
					
				end if
			end if
			maxZ = max(Z,maxZ)
			minZ = min(Z,minZ)
			maxN = max(N,maxN)
			minN = min(N,minN)
			maxops = max(iop,maxops)
			minops = min(iop,minops)
			read(milfile,'(a50)')lineread


          do i = 1,nlevels_assumed
				read(milfile,*,err=101,end=201)ilevel,e,xj,xt,xval
				if(fill .and. iop == 1)then
					list_energies(nnukes,ilevel)=e
					
				end if
				if(fill)then
					list_expect(nnukes,ilevel,iop)=xval
				end if
				if(ilevel/=i)then
					print*,' mismatch '
					stop
				end if
				minlevels = min(minlevels,ilevel)
				maxlevels = max(maxlevels,ilevel)
				cycle
101             continue
                endoflocal = .true.
                exit
				
			end do
			
		end do
201     continue	


if(.not.fill)then	
        print*,' There appear to be ',nnukes,' nuclides in the files'
		print*,' min/max Z = ',minZ,maxZ
		print*,' min/max N = ',minN,maxN
		print*,' min/max of levels = ',minlevels,maxlevels
		print*,' min/max of ops = ',minops,maxops
		nops = maxops
end if
		num_nuclides= nnukes
		return
		
		
	end subroutine survey_milcom_file

!===============================================================================	
!
! revised April 2019 for version 3, to use JF's new format for more generalized input
!	
	subroutine new_survey_milcom_file(fill)
		implicit none
		
		integer :: maxN,maxZ,minN,minZ
		integer :: Nnukes
		integer :: maxops,minops
		integer :: z,n,iop
		integer :: oldz,oldn
		integer :: ilevel,i,jlevel
		real    :: e,xj,xt2,xval,xi
		logical :: endofroad,endoflocal
		character*50 :: lineread
		integer :: maxlevels,minlevels
		integer :: nvals
		real :: scale
		logical fill
		integer :: iline
		
		endofroad = .false.
		maxN = 0
		maxZ = 0
		minZ = 10000
		minN = 10000
		maxops = 0000
		minops = 100000
		maxlevels = 00000
		minlevels = 1000000
		Nnukes = 0
		oldz= -1
		oldn = -1
		nvals = 0
		
		
		if(fill)then
			list_nuclides(:,:)=-999
			num_nuclides = 0
		end if
		ilevel = 0
        iline = 0
		do while (.not.endofroad)
			read(milfile,*,end=201,err=201)z,n,xj,xt2,xi,e,iop,xval
			iline = iline+1
			

			maxZ = max(Z,maxZ)
			minZ = min(Z,minZ)
			maxN = max(N,maxN)
			minN = min(N,minN)
			maxops = max(iop,maxops)
			minops = min(iop,minops)
			nvals = nvals+1
			
			if(.not.fill .and. iop==1)nlevels = nlevels+1
			
			if(fill)then
!............... CHECK TO SEE IF NUCLIDE ALREADY FOUND....	
                if(num_nuclides > 0)then	
				    do i = 1,num_nuclides
					
					    if(Z==list_nuclides(i,1) .and. n==list_nuclides(i,2))go to 3
			   	    end do ! i
			    end if
				num_nuclides = num_nuclides + 1
				list_nuclides(num_nuclides,1)=Z
				list_nuclides(num_nuclides,2)=N
3               continue			
				
				if(iop==1)then
					ilevel = ilevel + 1
					list_levels(ilevel,1)=Z
					list_levels(ilevel,2)=N
					list_energies(1,ilevel)=e
!					if(z==8 .and. n==16)print*,ilevel,' found ',Z,N,e
					jlevel = ilevel
				else   ! search for which level it is
!					print*,' repeat level ',ilevel
					do jlevel = 1,ilevel
						if(list_levels(jlevel,1)/=Z .or. list_levels(jlevel,2)/=N )cycle
						
!						if(Z==8 .and. n==16)print*,e, list_energies(1,jlevel)
						if(abs(list_energies(1,jlevel)-e) > 0.001)cycle
						
						go to 4
						
					end do
					print*,' Missing level? line ',iline
					print*,Z,N,e,iop
					stop

					
				end if
4				continue				
				
				scale =1.0
				if(scale_me .and. iop <= 63)then
					scale = (float(A0_scale)/float(Z+N))**X_scale
				end if
				list_expect(1,jlevel,iop)=xval*scale

			end if
			
			
			

			
		end do
201     continue	


if(.not.fill)then	
		print*,' min/max Z = ',minZ,maxZ
		print*,' min/max N = ',minN,maxN
		print*,' min/max of ops = ',minops,maxops
		nops = maxops
		print*,nvals,' vals to analyze'
		print*,' Have a total of ',nlevels,' levels counted ',nlevels*66,iline
		print*,' check ',nvals/nlevels,' ops? '
		num_nuclides = (maxZ-minZ+1)*(maxN-minZ+1)  ! just temp
		nlevels_chosen = nlevels
!		num_nuclides = 1  ! just for default
end if
if(fill)then
	print*,' Have a total of ',num_nuclides,' unique nuclides '
	num_nuclides = 1  ! just for default
end if


		return
		
		
	end subroutine new_survey_milcom_file
	
!===============================================================================	
subroutine choose_limits
	implicit none
	integer zi,ni
	
	print*,' Enter # of levels per nuclide desired '
	read*,nlevels_chosen
	
	
	print*,' Enter min,max of Z '
	read*,myzmin,myzmax
	print*,' Enter min,max of N '
	read*,mynmin,mynmax
	
	ncases = 0
	
	do zi = myzmin,myzmax
		do ni = mynmin,mynmax
			if(ni < zi )cycle
			ncases = ncases + 1
			
			
		end do
	end do
	print*,' a total of ',ncases,' nuclides kept, '
	print*,' or a total of ',ncases*nlevels_chosen,' levels included '
	
	return
	
	
end subroutine choose_limits
!===============================================================================	


subroutine setup_milcoms
	
	implicit none
		
	integer :: i,j,k,ii,iselect,icase
	integer :: myz,myn
	integer :: levelstart, levelstop
	
	allocate(B(ncases*nlevels_chosen,nops))
	
	iselect = 0
	
	if(excited_only)then
		levelstart = 2
		levelstop = nlevels_chosen+1
	else
		levelstart = 1
		levelstop = nlevels_chosen
	end if
	
	do i = 1,num_nuclides
		myz = list_nuclides(i,1)
		myn = list_nuclides(i,2)
		
		if(myz < myzmin.or. myz > myzmax )cycle
		if(myn < mynmin.or. myn > mynmax )cycle
		
		do j = levelstart,levelstop
			iselect = iselect +1
			
			do k = 1,nops

				if(excited_only)then
					B(iselect,k) = list_expect(i,j,k)-list_expect(i,1,k)
				else
					B(iselect,k) = list_expect(i,j,k)		
!					write(87,*)iselect,k,		list_expect(i,j,k)		
				end if

			end do
		end do
	end do
	print*,iselect,' total cases considered '
	
	return
	
	
end subroutine setup_milcoms
!===============================================================================	


subroutine setup_milcoms_new
	
	implicit none
		
	integer :: i,j,k,ii,iselect,icase
	integer :: myz,myn
	integer :: levelstart, levelstop
	
	allocate(B(nlevels_chosen,nops))


		
		do j = 1,nlevels_chosen
			
			do k = 1,nops

					B(j,k) = list_expect(1,j,k)		
!					write(87,*)iselect,k,		list_expect(i,j,k)		

			end do
		end do
		ncases =1 
	
	return
	
	
end subroutine setup_milcoms_new
!===============================================================================	
subroutine solve_milcoms
	implicit none
	
	integer i,k,kk
	real(8) :: dtmp,dtrace
	real(8), allocatable :: work(:)
	integer info
	
	allocate(A(nops,nops))
	A = 0.d0
	dtrace= 0.d0
	print*,ncases,' cases '
	do k = 1,nops
		do kk = k,nops
			dtmp = 0.d0
			do i = 1,ncases*nlevels_chosen
				dtmp = dtmp + B(i,k)*B(i,kk)
	!			if(k==65 .and. kk==65)print*,i,B(i,k)
				
			end do
!.... MAKE MATRIX NEGATIVE TO GET INTO DESCENDING ORDER....			
			A(k,kk)=-dtmp/real(ncases*nlevels_chosen-nops,8)
			A(kk,k)=-dtmp/real(ncases*nlevels_chosen-nops,8)
			
		end do
		dtrace = dtrace  -A(k,k)
!		print*,k,A(k,k)
	end do
	allocate(lambda2(nops),work(3*nops))
	call dsyev('V','U',nops,A,nops,lambda2,work,3*nops,info)
	print*,nops,nlevels_chosen,info

	do k = 1,nops
		lambda2(k)=-lambda2(k)
		write(88,*)k,lambda2(k)
	end do
	print*,' trace = ',dtrace
	return
	
	
end subroutine solve_milcoms

subroutine write_svd_evals_to_file(n,eigval)
	implicit none
	integer :: n
	real(8) :: eigval(n)
	character*70 :: filename
	integer      :: ilast
	logical      :: success
	integer :: i
	
	print*,' '
	print*,' Enter name of file to write SVD eigenvalues to '
	read(5,'(a)')filename
	ilast = index(filename,' ')-1
	open(unit=31,file=filename(1:ilast)//".dat",status='unknown')
	print*,' file ',filename(1:ilast),'.dat opened '
	do i = 1,n
		write(31,*)i,eigval(i)
	end do

	close(31)
	return
	
end subroutine write_svd_evals_to_file

subroutine write_milcoms_to_file(n,eigval,eigvec)
	implicit none
	integer :: n
	real(8) :: eigval(n),eigvec(n,n)
	character*70 :: filename
	integer      :: ilast
	logical      :: success
	integer :: i,j
	
	print*,' '
	print*,' Enter name of file to write MILCOMs to '
	read(5,'(a)')filename
	ilast = index(filename,' ')-1
	open(unit=33,file=filename(1:ilast)//".milcom",status='unknown')
	print*,' file ',filename(1:ilast),'.milcom opened '
	write(33,*)n
	do i = 1,n
	write(33,*)eigval(i)
end do
do i = 1,n
	do j = 1,n
	write(33,*)eigvec(i,j)
end do
end do
	close(33)
	return
	
end subroutine write_milcoms_to_file


subroutine read_milcoms_from_file(n,eigval,eigvec)
	implicit none
	integer :: n
	real(8) :: eigval(n),eigvec(n,n)
	character*70 :: filename,tmpline
	integer      :: ilast
	logical      :: success
	integer      :: ntmp,i,j
	real(8) :: dnorm
	
	print*,' '
	success =.false.
	do while(.not.success)
	   print*,' Enter name of file to read MILCOMs from '
	   read(5,'(a)')filename
	   ilast = index(filename,' ')-1
	   inquire(file=filename(1:ilast)//".milcom",exist=success)
	   if(.not.success)then
		   write(6,*)' file ',filename(1:ilast),'.milcom does not appear to exist '
	   else
		   
	   open(unit=34,file=filename(1:ilast)//".milcom",status='old')
	   print*,' opened '
      end if
    end do
! READ PAST ANY HEADERS
!    success = .false.
!	print*,'(A)'

!    do while(.not.success)
!		read(34,'(a70)')tmpline
!		if(tmpline(1:1)=='!' .or. tmpline(1:1)=='#')then
!			write(6,'(a70)')tmpline
!		else
!			success =.true.
!			backspace(34)
!		end if
		
!	end do		
       
	read(34,*)ntmp
	print*,ntmp
	if(ntmp /= n)then
		print*,' mismatch in dimensions ',n,ntmp
		stop
	end if
	do i = 1,n
	read(34,*)eigval(i)
!	print*,i,eigval(i)
end do
do i = 1,n
	do j = 1,n
!	read(34,*)eigvec(i,j)
!   READ COLUMN-WISE (Fox)
read(34,*)eigvec(j,i)
!	print*,i,j,eigvec(i,j)
end do
end do

	close(34)
	print*,' closed '
	return
	
end subroutine read_milcoms_from_file


subroutine readin_int_vec(localvec)
	implicit none
	character*70 :: filename,tmpline
	integer      :: ilast
	logical      :: success
	integer      :: ntmp,i,j
	real(8)      :: localvec(nops)
	print*,' '
	success =.false.
	do while(.not.success)
	   print*,' Enter full name of file to read interaction vector from '
	   read(5,'(a)')filename
	   ilast = index(filename,' ')-1
	   inquire(file=filename(1:ilast),exist=success)
	   if(.not.success)then
		   write(6,*)' file ',filename(1:ilast),' does not appear to exist '
	   else
		   
	   open(unit=55,file=filename(1:ilast),status='old')
	   print*,' opened '
      end if
    end do	
	read(55,*)ntmp
	if(ntmp /= nops)then
		print*,' Mismatch in # of operators (dimension of vector ) ',ntmp,nops
		stop
	end if
	do i = 1,nops
		read(55,*)localvec(i)
	end do 
	close(55)

	
	return
	
end subroutine readin_int_vec

subroutine writeout_int_vec(localvec)
	implicit none
	character*70 :: filename,tmpline
	integer      :: ilast
	logical      :: success,whoops
	integer      :: ntmp,i,j
	real(8)      :: localvec(nops)
	print*,' '
	success =.false.
	do while(.not.success)
	   print*,' Enter name of file to read interaction vector to '
	   read(5,'(a)')filename
	   ilast = index(filename,' ')-1
	   inquire(file=filename(1:ilast),exist=whoops)
	   if(whoops)then
		   write(6,*)' file ',filename(1:ilast),'.milcom already exists!'
		   cycle
	   else
		   
	   open(unit=56,file=filename(1:ilast),status='new')
	   success = .true.
	   print*,' opened '
      end if
    end do	
	write(56,*)nops

	do i = 1,nops
		write(56,*)localvec(i)
	end do 
	close(56)

	
	return
	
end subroutine writeout_int_vec


!----- 
!  FORWARD = .true. transforms interaction vector INTO milcom-space
!          =.false. transforms from milcom-space INTO 'standard' space
!
subroutine transform_int_vec(umat,vecin,vecout,forward)
	implicit none
    logical      :: forward
	integer      :: ntmp,i,j
	real(8)      :: umat(nops,nops)
	real(8)      :: vecin(nops),vecout(nops),dtmp

!.............. NOW TRANSFORM..........................................	
!
!  MILCOMS SVD vectors stored in U  so take
!  vectrans = U^T*vec
!

    do i = 1,nops
		dtmp = 0.d0
		if(forward)then
		   do j = 1,nops
			   dtmp = dtmp + umat(j,i)*vecin(j)
		   end do
	    else
 		   do j = 1,nops
 			   dtmp = dtmp + umat(i,j)*vecin(j)
 		   end do
		end if
		vecout(i) = dtmp
	
	end do
			
	
	return
	
end subroutine transform_int_vec
		
end module milcom_data

!=============================== MAIN PROGRAM ====================

use milcom_data
implicit none

character*1 :: menu_char

logical     :: finished

print*,' Welcome to creation and analysis of MILCOMs '
print*,' Version 3, April 2019'

finished = .false.

do while (.not.finished)
	print*,'  Enter menu choice : '
	print*,' (m) *M*ake MILCOM (old format)'
	print*,' (n) Make MILCOM (*n*ew - 4/2019 format)'
	print*,' (c) *C*ompare MILCOMs '
	print*,' (t) *T*ransform vectors via MILCOM space'
	print*,' (p) *P*eturb a interaction vector'
	print*,' (x) E*x*it code '
	
	
	read(5,'(a)')menu_char
	
	select case(menu_char)
	
	   case ('m','M')
	      call basic_milcom_routine
   	   case ('n','N')
   	      call new_milcom_routine	  
	   case('c','C')
	      call compare_milcoms
		  
	   case('t','T')
	      call transform_wrt_milcoms
		  
	   case ('p','P')	  
	     call perturb_intvec
		  
	   case ('x','X')
	      finished = .true.
	      cycle
	
    end select
	
end do ! while not finished

end

! original data format 

subroutine basic_milcom_routine
	use milcom_data
	implicit none
	character :: ychar
	
	print*,' '
	print*,' Some information on this run '
	if(excited_only)print*,' Using only excitation energies '
	if(scale_me)then
		print*,' Not yet set up for scaling '
!		print*,' Scaling TBMEs by (',A0_scale,' / A)^',X_scale
		
	end if
	
	print*,' '

	call open_milcom_file
	call survey_milcom_file(.false.)
	print*,num_nuclides, nlevels_assumed,nops
	allocate(list_nuclides(num_nuclides,2))
	allocate(list_energies(num_nuclides,nlevels_assumed))
	allocate(list_expect(num_nuclides,nlevels_assumed,nops))
	rewind(milfile)
	call survey_milcom_file(.true.)

	call choose_limits
	call setup_milcoms
	call solve_milcoms
	print*,' Do you want to write SVD spectrum to file ?'
	read(5,'(a)')ychar
	if(ychar=='y' .or. ychar=='Y')then
		call write_svd_evals_to_file(nops,lambda2)
	end if
	print*,' Do you want to save MILCOMs to file? '
	read(5,'(a)')ychar
	if(ychar=='y' .or. ychar=='Y')then
		call write_milcoms_to_file(nops,lambda2,A)
		
	end if
	
!	print*,' Do you want to transform an interaction vector? (y/n)'
!	read(5,'(a)')ychar
!	if(ychar=='y' .or. ychar=='Y')then
!   end if	
	
	return
end subroutine basic_milcom_routine

subroutine new_milcom_routine
	use milcom_data
	implicit none
	character :: ychar

	
	print*,' '
	print*,' Some information on this run '
	if(excited_only)print*,' Using only excitation energies '
	if(scale_me)then
		print*,' Scaling TBMEs by (',A0_scale,' / A)^',X_scale
		
	end if
	
	print*,' '

	call open_milcom_file
	call new_survey_milcom_file(.false.)
!	print*,num_nuclides, nlevels_assumed,nops
	allocate(list_nuclides(num_nuclides,2))
	allocate(list_energies(1,nlevels),list_levels(nlevels,2))
	allocate(list_expect(1,nlevels,nops))
	rewind(milfile)
	print*,' Rereading '
	call new_survey_milcom_file(.true.)
	print*,' Done rereading'

!	call choose_limits
	call setup_milcoms_new
	call solve_milcoms
	print*,' Do you want to write SVD spectrum to file ?'
	read(5,'(a)')ychar
	if(ychar=='y' .or. ychar=='Y')then
		call write_svd_evals_to_file(nops,lambda2)
	end if
	print*,' Do you want to save MILCOMs to file? '
	read(5,'(a)')ychar
	if(ychar=='y' .or. ychar=='Y')then
		call write_milcoms_to_file(nops,lambda2,A)
		
	end if

	
	return
end subroutine new_milcom_routine


subroutine compare_milcoms
	use milcom_data
	implicit none
	integer nn,n1,n2
	real(8), allocatable :: milcom1(:,:),milcom2(:,:),eval(:)
	
	print*,' Enter # of operators '
	read*,nn
	allocate(milcom1(nn,nn),milcom2(nn,nn),eval(nn))
	print*,size(eval),size(milcom1),size(milcom2)
	call read_milcoms_from_file(nn,eval,milcom1)
	print*,' How many MILCOMs do you want to keep from this?'
	read*,n1
	call read_milcoms_from_file(nn,eval,milcom2)
	print*,' How many MILCOMs do you want to keep from this?'
	read*,n2
	
	call subspace_overlaps(nn,milcom1,n1,milcom2,n2)
	
	
end subroutine compare_milcoms
	
!
! subroutine to determine overlap between two subspaces
!
!
subroutine subspace_overlaps(n,vecs1,n1,vecs2,n2)
	implicit none
	integer :: n
	real(8)  :: vecs1(n,n),vecs2(n,n)
	integer :: n1, n2,nn
	real(8), allocatable :: myb(:,:),mya(:,:)
	integer i1,i2,i
	real(8) :: dtmp
	real(8), allocatable :: work(:),ee(:)
	integer info
	
	
	allocate(myb(n1,n2))
	
	do i1 = 1,n1
		do i2 = 1,n2
			dtmp = 0.d0
			do i = 1,n
			    dtmp = dtmp + vecs1(i,n-i1+1)*vecs2(i,n-i2+1)	
				
			end do ! i
			myb(i1,i2)=dtmp
!			print*,i1,i2,dtmp
 
		end do ! i2
		
	end do ! i1
	
	nn = min(n1,n2)
	print*,nn
	allocate(mya(nn,nn))
	
		do i1 = 1,nn
			
			do i2 = 1,nn
				dtmp = 0.d0
				
				if(n1 < n2)then
				
				   do i = 1,n2
					  dtmp = dtmp + myb(i1,i)*myb(i2,i)
					
				   end do
			   else
				
				   do i = 1,n1
					  dtmp = dtmp + myb(i,i1)*myb(i,i2)
					
				   end do
			    end if
!				print*,i1,i2,dtmp
                mya(i1,i2)=dtmp
				
			end do
		end do

		allocate(ee(nn),work(3*nn))
!		print*,nn
!        do i1 = 1,nn
!			write(6,*)(mya(i1,i2),i2=1,nn)
!		end do

		call dsyev('V','U',nn,mya,nn,ee,work,3*nn,info)
		print*,' check  ',info
		do i = 1,nn
			print*,i,ee(i)
		end do
	
	
	return
	
end subroutine subspace_overlaps

!---------------------------------------------------------------------------------
!.......... ROUTINE TO TRANSFORM INTERACTION VECTORS INTO/FROM MILCOM SPACE........


subroutine transform_wrt_milcoms
	use milcom_data
	implicit none
	integer nn,n1,n2
	real(8), allocatable :: milcom1(:,:),eval(:),vecin(:),vecout(:)
	character :: transchar
	
	print*,' Transforming interaction vectors to/from MILCOM space '
	print*,' Enter # of operators '
	read*,nn
	nops = nn
	allocate(milcom1(nn,nn),eval(nn))
	call read_milcoms_from_file(nn,eval,milcom1)	
	allocate(vecin(nops),vecout(nops))
	call readin_int_vec(vecin)
	
	print*,' Do you want to transform INTO MILCOM space (i) or FROM MILCOM space (f)?'
	read(5,'(a)')transchar
	
	select case (transchar)
	case('i')
	   call  transform_int_vec(milcom1,vecin,vecout,.true.)
	case('f')
	   call  transform_int_vec(milcom1,vecin,vecout,.false.)
	case default
	   print*,transchar,' is not a recognized response '
	   stop
    end select

	call writeout_int_vec(vecout)
	
	return
end subroutine transform_wrt_milcoms
	
!===============================================================================	
subroutine perturb_intvec
	use milcom_data
	implicit none
	integer nn,n1,n2
	real(8), allocatable ::vecin(:),vecout(:)
	logical :: done
	real(8) :: dv
	
	print*,' Perturbing an interaction vector  '
	print*,' Enter # of operators '
	read*,nn
	nops= nn

	allocate(vecin(nops),vecout(nops))
	call readin_int_vec(vecin)
	

    do n1 = 1,nn
		vecout(n1)=vecin(n1)
	end do	
	done = .false.
	do while (.not.done)
		print*,' Enter index, dv for perturbation (enter index=0 to stop)'
		read*,n1,dv
		if(n1 < 1)then
			done=.true.
			exit
		end if
		vecout(n1)=vecout(n1)+dv
		
	end do
	
	print*,' All done '

	call writeout_int_vec(vecout)
	
	return
end subroutine perturb_intvec
	
!===============================================================================		
	
