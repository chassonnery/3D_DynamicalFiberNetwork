module lecture_ecriture

implicit none

save

contains

subroutine LectureParam(Pmod,Psim,Pdom)
!- Read the parameters for the simulation from the file PARAMETERS.txt
use definitions
	implicit none
	TYPE(ParametersSim), intent(out)  :: Psim
	TYPE(ParametersDom), intent(out)  :: Pdom
	TYPE(ParametersMod), intent(out)  :: Pmod
	character(8)                      :: temp
	Double Precision, PARAMETER       :: pi = 3.14159265358979323846

	! Read file parameters.txt ------------------------------------------------
	open(unit=15,file='PARAMETERS.txt',status='old',action='read')
	read(15,*) temp
	read(15,*) temp
	read(15,*) Pmod%Ncell
	read(15,*) Pmod%Nmax
	read(15,*) Pmod%Rmin
	read(15,*) Pmod%Rmax
	read(15,*) Pmod%Nfib
	read(15,*) Pmod%Lf
	read(15,*) Pmod%Rf
	read(15,*) temp
	read(15,*) Pmod%alpha_repCC
	read(15,*) Pmod%alpha_repCF
	read(15,*) Pmod%alpha_repFF
	read(15,*) Pmod%alpha_align
	read(15,*) Pmod%alpha_rappel
	read(15,*) Pmod%cumul_replink
	read(15,*) temp
	read(15,*) Pmod%freq_ens
	read(15,*) Pmod%Kcroiss
	read(15,*) temp
	read(15,*) Pmod%freq_link
	read(15,*) Pmod%freq_dlink
	read(15,*) Pmod%LinkType
	read(15,*) Pmod%eps
	read(15,*) Pmod%d0
	read(15,*) temp
	read(15,*) Pdom%xmaxR
	read(15,*) Pdom%ymaxR
	read(15,*) Pdom%zmaxR
	read(15,*) Pmod%DirFillingRate
	read(15,*) Pmod%DirLim(1:6)
	read(15,*) temp
	read(15,*) Psim%TF
	read(15,*) Psim%dt
	read(15,*) temp
	read(15,*) Psim%InitType
	read(15,*) temp
	read(15,*) Psim%Tenr
	read(15,*) Psim%Nomdossier
	close(unit=15)

	! Check validity of the parameters ----------------------------------------
	if (Pmod%Ncell.lt.0) then
		write(*,*) "The initial number of cells must be a positive or null integer"
		STOP
	elseif (Pmod%Nmax.lt.Pmod%Ncell) then
		write(*,*) "The maximal number of cells must be an integer greater than the initial number of cells"
		STOP
	elseif (Pmod%Rmin.le.0) then
		write(*,*) "The cells' minimal radius must be positive"
		STOP
	elseif (Pmod%Rmax.lt.Pmod%Rmin) then
		write(*,*) "The cells' maximal radius must be greater than the cells' minimal radius"
		STOP
	elseif (Pmod%Nfib.lt.0) then
		write(*,*) "The number of fibers must be a positive or null integer"
		STOP
	elseif (Pmod%Lf.lt.0) then
		write(*,*) "The fibers' length must be positive or null"
		STOP
	elseif (Pmod%Rf.le.0) then
		write(*,*) "The fibers' radius must be positive"
		STOP
	elseif ((Pmod%cumul_replink.ne.0).and.(Pmod%cumul_replink.ne.1)) then
		write(*,*) "Parameter cumuk_replink can only be equal to 0 or 1"
		STOP
	elseif (Pmod%freq_ens.lt.0) then
		write(*,*) "The insemination frequency must be positive or null"
		STOP
	elseif (Pmod%Kcroiss.lt.0) then
		write(*,*) "The cell growth's speed must be positive or null"
		STOP
	elseif (Pmod%freq_link.lt.0) then
		write(*,*) "The fiber's linking frequency must be positive or null"
		STOP
	elseif (Pmod%freq_dlink.lt.0) then
		write(*,*) "The fiber's unlinking frequency must be positive or null"
		STOP
	elseif ((Pmod%LinkType.ne.0).and.(Pmod%LinkType.ne.1).and.(Pmod%LinkType.ne.2)) then
		write(*,*) "Parameter LinkType can only be equal to 0, 1 or 2"
		STOP
	elseif (Pmod%eps.lt.0) then
		write(*,*) "The fiber's linking distance must be positive or null"
		STOP
	elseif (Pmod%d0.lt.0) then
		write(*,*) "The fiber link's equilibrium distance must be positive or null"
		STOP
	elseif ((Pdom%xmaxR.le.0).or.(Pdom%ymaxR.le.0).or.(Pdom%zmaxR.le.0)) then
		write(*,*) "The domain's length in the three dimensions must be positive"
		STOP
	elseif (Pmod%DirFillingRate.lt.0) then
		write(*,*) "The filling rate of the dirichlet layer must be positive or null"
		STOP
	elseif (Psim%Tf.le.0) then
		write(*,*) "The final simulated time must be positive"
		STOP
	elseif (Psim%dt.le.0) then
		write(*,*) "The maximal time-step allowed must be positive"
		STOP
	elseif ((Psim%InitType.ne.0).and.(Psim%InitType.ne.1).and.(Psim%InitType.ne.2)) then
		write(*,*) "Parameter InitType can only be equal to 0, 1 or 2"
		STOP
	elseif (Psim%Tenr.le.0) then
		write(*,*) "The data's saving frequency must be positive"
		STOP
	end if
	
end subroutine LectureParam



subroutine LectureInput(Xcell,Rcell,Xfib,omega,link,distpoints,Pmod,Pdom,BoxCell,BoxFib)
use definitions
use fonctions
	implicit none
	double precision, dimension(:,:), intent(inout) :: Xcell, Xfib, omega, link, distpoints
	double precision, dimension(:), intent(inout)   :: Rcell
	TYPE(ParametersMod), intent(in)                 :: Pmod
	type(ParametersDom), intent(in)                 :: Pdom
	type(ParametersBox), intent(inout)              :: BoxCell,BoxFib
	integer :: i,k,ios
	double precision :: lik,lki
	
	! Read cells' data
	if (Pmod%Ncell.ne.0) then
		open(unit=15,file='INPUTCELLS.dat',status='old',action='read')
		do i=1,Pmod%Ncell
			read(15,*) Xcell(i,:), Rcell(i)
		end do
		close(unit=15)
	end if
	
	! Read system fibers' data
	if (Pmod%Nfib.ne.0) then
		open(unit=15,file='INPUTFIBERS.dat',status='old',action='read')
		do i=1,Pmod%Nfib
			read(15,*) Xfib(i,:), omega(i,:)
		end do
		close(unit=15)
	end if
	
	! Read dirichlet-fibers' data
	if (Pmod%Ndirtot.ne.0) then
		open(unit=15,file='INPUTDIRICHLET.dat',status='old',action='read')
		do i=Pmod%Nfib+1,Pmod%Nfib+Pmod%Ndirtot
			read(15,*) Xfib(i,:), omega(i,:)
		end do
		close(unit=15)
	end if
	
	! Read links' data
	link = 0
	distpoints = 0.0
	if (Pmod%LinkType.ne.2) then
		open(unit=15,file='INPUTLINKS.dat',status='old',action='read')
		read(15,*,iostat=ios) i,k,lik,lki
		do while(ios==0)
			link(i,k) = 1
			if(k<=Pmod%Nfib) then
				link(k,i) = 1
			end if
			distpoints(i,k) = lik
			distpoints(k,i) = lki
			read(15,*,iostat=ios) i,k,lik,lki
		end do
		close(unit=15)
	end if

	! Allocate agents to corresponding boxes
	do i = 1,Pmod%Ncell
		call UpdateVerlet(i,Xcell(i,:),Pdom,BoxCell)
	end do
	call FillGhostBoxes(Pmod,Pdom,BoxCell)
	
	do i = 1,Pmod%Nfib+Pmod%Ndirtot
		call UpdateVerlet(i,Xfib(i,:),Pdom,BoxFib)
	end do
	call FillGhostBoxes(Pmod,Pdom,BoxFib)
end subroutine LectureInput



subroutine LecturePreviousSim(Xcell,Rcell,Xfib,omega,link,distpoints,NbStepInPreviousSim,&
	nomdossier_old,Pmod,Pdom,BoxCell,BoxFib)
use definitions
use fonctions
	implicit none
	double precision, dimension(:,:), intent(inout) :: Xcell, Xfib, omega, link, distpoints
	double precision, dimension(:), intent(inout)   :: Rcell
	integer, intent(inout)                          :: NbStepInPreviousSim
	character(len=80), intent(in)                   :: nomdossier_old
	TYPE(ParametersMod), intent(inout)              :: Pmod
	type(ParametersDom), intent(in)                 :: Pdom
	type(ParametersBox), intent(inout)              :: BoxCell,BoxFib
	integer :: i,k,ios
	character(len=30) :: ch_iter
	double precision :: x,y,z,r,lik,lki

	! Find the number of step done in the previous simulation
	call system('find '//trim(adjustl(nomdossier_old))//' -maxdepth 1 -name "fibers*.dat" |wc -l >NbStepInPreviousSim.dat')
	open(unit=15, file='NbStepInPreviousSim.dat', status='old', action='read')
	read(15,*) NbStepInPreviousSim
	close(15)
	call system('rm NbStepInPreviousSim.dat')
	
	NbStepInPreviousSim = NbStepInPreviousSim !- 1 to account for the file fibers0.dat
	write(ch_iter,*) NbStepInPreviousSim
	
	
	! Read cells' data
	Pmod%Ncell = 0
	if (access(trim(adjustl(nomdossier_old))//'/cells'//trim(adjustl(ch_iter))//'.dat', ' ' ) .eq. 0) then ! check if the file exists
		open(unit=15,file=trim(adjustl(nomdossier_old))//'/cells'//trim(adjustl(ch_iter))//'.dat',status='old',action='read')
		read(15,*,iostat=ios) x, y, z, r
		do while (ios==0)
			Pmod%Ncell = Pmod%Ncell + 1
			Xcell(Pmod%Ncell,1) = x
			Xcell(Pmod%Ncell,2) = y
			Xcell(Pmod%Ncell,3) = z
			Rcell(Pmod%Ncell) = r
			read(15,*,iostat=ios) x, y, z, r 
		end do
		close(unit=15)
	end if
	
	! Read system fibers' data
	if (Pmod%Nfib.ne.0) then
		open(unit=15,file=trim(adjustl(nomdossier_old))//'/fibers'//trim(adjustl(ch_iter))//'.dat',status='old',action='read')
		do i=1,Pmod%Nfib
			read(15,*) Xfib(i,:), omega(i,:)
		end do
		close(unit=15)
	end if
	
	! Read dirichlet-fibers' data
	if (Pmod%Ndirtot.ne.0) then
		open(unit=15,file=trim(adjustl(nomdossier_old))//'/dirichlet_layer.dat',status='old',action='read')
		do i=Pmod%Nfib+1,Pmod%Nfib+Pmod%Ndirtot
			read(15,*) Xfib(i,:), omega(i,:)
		end do
		close(unit=15)
	end if

	! Read links' data
	link = 0
	distpoints = 0.0
	if (Pmod%LinkType.eq.0) then
		open(unit=15,file=trim(adjustl(nomdossier_old))//'/links'//trim(adjustl(ch_iter))//'.dat',status='old',action='read')
		read(15,*,iostat=ios) i,k,lik,lki
		do while(ios==0)
			link(i,k) = 1
			if(k<=Pmod%Nfib) then
				link(k,i) = 1
			end if
			distpoints(i,k) = lik
			distpoints(k,i) = lki
			read(15,*,iostat=ios) i,k,lik,lki
		end do
		close(unit=15)
	elseif (Pmod%LinkType.eq.1) then
		open(unit=15,file=trim(adjustl(nomdossier_old))//'/links0.dat',status='old',action='read')
		read(15,*,iostat=ios) i,k,lik,lki
		do while(ios==0)
			link(i,k) = 1
			if(k<=Pmod%Nfib) then
				link(k,i) = 1
			end if
			distpoints(i,k) = lik
			distpoints(k,i) = lki
			read(15,*,iostat=ios) i,k,lik,lki
		end do
		close(unit=15)
	end if

	! Allocate agents to corresponding box
	do i = 1,Pmod%Ncell
		call UpdateVerlet(i,Xcell(i,:),Pdom,BoxCell)
	end do
	call FillGhostBoxes(Pmod,Pdom,BoxCell)
	
	do i = 1,Pmod%Nfib+Pmod%Ndirtot
		call UpdateVerlet(i,Xfib(i,:),Pdom,BoxFib)
	end do
	call FillGhostBoxes(Pmod,Pdom,BoxFib)
end subroutine LecturePreviousSim



subroutine WriteOutput_onefile(ksave0,ksavef,Saving,Pmod,Psim)
use definitions
	implicit none
	integer, intent(in)             :: ksave0, ksavef
	TYPE(SaveState),     intent(in) :: Saving
	TYPE(ParametersMod), intent(in) :: Pmod
	TYPE(ParametersSim), intent(in) :: Psim
	integer :: i,k
	
	open(101,file=trim(adjustl(Psim%nomdossier))//'/cells.dat',action='write')
	open(102,file=trim(adjustl(Psim%nomdossier))//'/fibers.dat',action='write')
	open(103,file=trim(adjustl(Psim%nomdossier))//'/links.dat',action='write')
		
	do k = ksave0,ksavef
		! Cells --------------------------------------------------------------------------------------
		do i=1,Saving%VariousInt(k,3)
			write(101,2000) k*Psim%Tenr, Saving%CellState(k,i,:)
		end do
		
		! Fibers -------------------------------------------------------------------------------------
		do i=1,Pmod%Nfib
			write(102,2000) k*Psim%Tenr, Saving%FibState(k,i,:)
		end do
		
		! Links --------------------------------------------------------------------------------------
		do i=1,Saving%VariousInt(k,1)+Saving%VariousInt(k,2) ! total number of links
			write(103,3000) k*Psim%Tenr, int(Saving%LinkState(k,i,1)), int(Saving%LinkState(k,i,2)), Saving%LinkState(k,i,3:4)
		end do
	end do
	
	close(101)
	close(102)
	close(103)

2000 FORMAT(f10.4,6(1x,e23.17e2))
3000 FORMAT(f10.4,2(1x,i15),2(1x,e23.17e2))

end subroutine WriteOutput_onefile



subroutine WriteOutput_manyfiles(Nsave0,Nsavef,offset,Saving,Pmod,Psim)
use definitions
	implicit none
	integer, intent(in)             :: Nsave0, Nsavef, offset
	TYPE(SaveState),     intent(in) :: Saving
	TYPE(ParametersMod), intent(in) :: Pmod
	TYPE(ParametersSim), intent(in) :: Psim
	character(len=30)               :: ch_iter
	integer :: i,k
	
	do k = Nsave0,Nsavef
		write(ch_iter,*) k+offset

		! Cells --------------------------------------------------------------------------------------
		if (Pmod%Nmax.ne.0) then 
			open(101,file=trim(adjustl(Psim%nomdossier))//'/cells'//trim(adjustl(ch_iter))//'.dat',action='write')
			do i=1,Saving%VariousInt(k,3)
				write(101,2000) Saving%CellState(k,i,:)
			end do
			close(101)
		end if
		
		! Fibers -------------------------------------------------------------------------------------
		if (Pmod%Nfib.ne.0) then
			open(102,file=trim(adjustl(Psim%nomdossier))//'/fibers'//trim(adjustl(ch_iter))//'.dat',action='write')
			do i=1,Pmod%Nfib
				write(102,2000) Saving%FibState(k,i,:)
			end do
			close(102)
		end if 
		
		! Links --------------------------------------------------------------------------------------
		if (Pmod%LinkType.eq.0) then 
			open(103,file=trim(adjustl(Psim%nomdossier))//'/links'//trim(adjustl(ch_iter))//'.dat',action='write')
			do i=1,Saving%VariousInt(k,1)+Saving%VariousInt(k,2) ! total number of links
				write(103,3000) int(Saving%LinkState(k,i,1)), int(Saving%LinkState(k,i,2)), Saving%LinkState(k,i,3:4)
			end do
			close(103)
		end if
	end do

2000 FORMAT(6(e23.17e2,1x))
3000 FORMAT(i15,1x,i15,1x,e23.17e2,1x,e23.17e2)

end subroutine WriteOutput_manyfiles



subroutine WriteInput(Xcell,Rcell,Xfib,omega,link,distpoints,Pmod,Pdom,Psim)
use definitions
	implicit none
	double precision, dimension(:,:), intent(inout) :: Xcell, Xfib, omega, link, distpoints
	double precision, dimension(:), intent(inout)   :: Rcell
	TYPE(ParametersMod), intent(in)                 :: Pmod
	TYPE(ParametersDom), intent(in)                 :: Pdom
	TYPE(ParametersSim), intent(in)                 :: Psim
	integer :: i,j
		
	! Cells --------------------------------------------------------------------------------------
	if (Pmod%Ncell.ne.0) then
		open(2,file=trim(adjustl(Psim%nomdossier))//'/INPUTCELLS.dat',action='write')
		do i=1,Pmod%Ncell
			write(2,2000) Xcell(i,1),Xcell(i,2),Xcell(i,3),Rcell(i)
		end do
		close(2)
	end if

	! Fibers -------------------------------------------------------------------------------------
	if (Pmod%Nfib.ne.0) then
		open(2,file=trim(adjustl(Psim%nomdossier))//'/INPUTFIBERS.dat',action='write')
		do i=1,Pmod%Nfib
			write(2,2000) Xfib(i,1), Xfib(i,2), Xfib(i,3), omega(i,1), omega(i,2), omega(i,3)
		end do
		close(2)
	end if

	! Links --------------------------------------------------------------------------------------
	if (Pmod%LinkType.ne.2) then 
		open(2,file=trim(adjustl(Psim%nomdossier))//'/INPUTLINKS.dat',action='write')
		do i=1,Pmod%Nfib
			do j=i+1,Pmod%Nfib+Pmod%Ndirtot
				if (link(i,j)==1) then
					write(2,3000) i,j,distpoints(i,j),distpoints(j,i)
				endif
			end do
		end do
		close(2)
	end if

2000 FORMAT(6(e23.17e2,' '))
3000 FORMAT(i15,' ',i15,' ',e23.17e2,' ',e23.17e2)

end subroutine WriteInput





end module lecture_ecriture