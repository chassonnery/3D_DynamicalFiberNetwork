program main

use definitions
use lecture_ecriture
use fonctions
use omp_lib

implicit none

!------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------------------DECLARATION PARAMETRES---------------------------------------------------!
!----------------------------------------------------UTILITAIRES---------------------------------------------------------!
double precision, dimension(:,:), allocatable :: Xcell, Xfib, omega, link, distpoints, gradC, gradF
double precision, dimension(:), allocatable   :: Rcell
double precision                              :: alea, dik, lik, lki, Xtemp(3), gradtemp(6), dist
integer, dimension(:), allocatable            :: IsInDirichlet, IsANeighbor
integer                                       :: i,j,k,l,l1,l2,l3,ltemp
integer                                       :: Nlinks_int, Nlinks_ext
!-----------------------------------------------MODELE et DOMAINE--------------------------------------------------------!
type(ParametersMod)                           :: Pmod
type(ParametersDom)                           :: Pdom
type(ParametersBox)                           :: BoxCell, BoxFib
!----------------------------------------------------NUMERIQUES----------------------------------------------------------!
type(ParametersSim)                           :: Psim
double precision, parameter                   :: pi=3.141592653589793238
double precision                              :: tcpu0,tcpuf,tcpuglob0,tcpuglobf,tcpu_days,tcpu_hours,tcpu_secs,tcpu_write
double precision                              :: t_end_ins
!----------------------------------------------------REDEMARRAGE SIM-----------------------------------------------------!
character(len=80)                             :: nomdossier_old
integer                                       :: NbStepInPreviousSim
!----------------------------------------------------ENREGISTREMENT DONNEES----------------------------------------------!
type(SaveState)                               :: Saving
integer, parameter                            :: Nmemorymax = 1000
integer                                       :: Nsavestep0, Nsavestep, Nsaveblock, Nbetween
double precision                              :: time, tnext, dt, dt_min, dt_max
character(len=8)                              :: realdate
character(len=10)                             :: realtime, timezone
integer(kind=4)                               :: timevalues(8)



call omp_set_num_threads(10)
tcpuglob0 = omp_get_wtime()
!------------------------------------------------------------------------------------------------------------------
!---------------------------- INITIALISATION ----------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------

! Initialize the seed for random processes
call init_random_seed(0) ! 0 for a random run and 1 for a repeatable run

! Lecture des paramètres
call LectureParam(Pmod,Psim,Pdom)




!------- Création du domaine --------------------------------------------------------------------------------------
dist = max(Pmod%Lf + 2*Pmod%Rf + Pmod%eps, Pmod%Lf + 1.5*Pmod%d0, 2*Pmod%Rmax) ! minimal size of a box = maximal length of interaction

Pdom%Nx = int(Pdom%xmaxR*1.0/dist) ! number of boxes that can fit in the dimension x
Pdom%Ny = int(Pdom%ymaxR*1.0/dist)
Pdom%Nz = int(Pdom%zmaxR*1.0/dist)

if (Pdom%Nx.le.1) then
	write(*,*) "ATTENTION : Nx < 2 !"
	write(*,*) "  Le nombre de boîtes numériques est insuffisant pour assurer le bon fonctionnement des conditions aux limites."
	write(*,*) "  Veuilliez augmenter la taille Lx du domaine ou réduire la distance maximale d'interaction."
	STOP
else if (Pdom%Ny.le.1) then
	write(*,*) "ATTENTION : Ny < 2 !"
	write(*,*) "  Le nombre de boîtes numériques est insuffisant pour assurer le bon fonctionnement des conditions aux limites."
	write(*,*) "  Veuilliez augmenter la taille Ly du domaine ou réduire la distance maximale d'interaction."
	STOP
else if (Pdom%Nz.le.1) then
	write(*,*) "ATTENTION : Nz < 2 !"
	write(*,*) "  Le nombre de boîtes numériques est insuffisant pour assurer le bon fonctionnement des conditions aux limites."
	write(*,*) "  Veuilliez augmenter la taille Lz du domaine ou réduire la distance maximale d'interaction."
	STOP
end if 

Pdom%dx = Pdom%xmaxR*1.0/Pdom%Nx ! size of the boxes in the dimension x
Pdom%dy = Pdom%ymaxR*1.0/Pdom%Ny
Pdom%dz = Pdom%zmaxR*1.0/Pdom%Nz

Pdom%xmax = (Pdom%Nx+1)*Pdom%dx ! limit of the domain in x including ghost boxes
Pdom%ymax = (Pdom%Ny+1)*Pdom%dy ! limit of the domain in y including ghost boxes
Pdom%zmax = (Pdom%Nz+1)*Pdom%dz ! limit of the domain in z including ghost boxes




!------- Nombre maximal de cellules -------------------------------------------------------------------------------
! Pmod%Nmax = int(8*Pdom%xmaxR*Pdom%ymaxR*Pdom%zmaxR /( 2.0*(4.0/3.0)*pi*Pmod%Rmax**3 )) 
! volume du domaine divisé par volume maximal d'une cellule, le tout sur 2




!------- Couche de Dirichlet --------------------------------------------------------------------------------------
allocate(BoxCell%BoxType(1:8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1),4))
allocate(BoxFib%BoxType(1:8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1),4))

call computeBoxType(Pdom,Pmod,BoxFib)
BoxCell%BoxType(:,:) = BoxFib%BoxType(:,:)

Pmod%Ndirtot = 0
do i = 1,8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)
	if (BoxFib%BoxType(i,1).eq.3) then
		Pmod%Ndirtot = Pmod%Ndirtot + 1
	end if
end do
Pmod%Ndirtot = Pmod%Ndirtot * int(Pmod%DirFillingRate * Pdom%dx*Pdom%dy*Pdom%dz)




!------- Allocation des variables ---------------------------------------------------------------------------------
allocate(Xcell(1:Pmod%Nmax,3))
allocate(Rcell(1:Pmod%Nmax))
allocate(gradC(1:Pmod%Nmax,3))

allocate(Xfib(1:Pmod%Nfib+Pmod%Ndirtot,3))
allocate(omega(1:Pmod%Nfib+Pmod%Ndirtot,3))
allocate(gradF(1:Pmod%Nfib,6))
allocate(IsInDirichlet(Pmod%Nfib))

allocate(link(Pmod%Nfib,Pmod%Nfib+Pmod%Ndirtot))
allocate(distpoints(Pmod%Nfib+Pmod%Ndirtot,Pmod%Nfib+Pmod%Ndirtot))
allocate(IsANeighbor(Pmod%Nfib+Pmod%Ndirtot))

allocate(BoxCell%FirstAgent(1:8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)))
allocate(BoxCell%LastAgent(1:8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)))
allocate(BoxCell%AgentBox(1:Pmod%Nmax))
allocate(BoxCell%VerletList(1:Pmod%Nmax))

allocate(BoxFib%FirstAgent(1:8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)))
allocate(BoxFib%LastAgent(1:8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)))
allocate(BoxFib%AgentBox(1:Pmod%Nfib+Pmod%Ndirtot))
allocate(BoxFib%VerletList(1:Pmod%Nfib+Pmod%Ndirtot))

BoxCell%FirstAgent(:) = 0
BoxCell%LastAgent(:)  = 0
BoxCell%AgentBox(:)   = 0
BoxCell%VerletList(:) = 0
BoxFib%FirstAgent(:) = 0
BoxFib%LastAgent(:)  = 0
BoxFib%AgentBox(:)   = 0
BoxFib%VerletList(:) = 0
IsInDirichlet(:) = 0

allocate(Saving%CellState(Nmemorymax,Pmod%Nmax,4))
allocate(Saving%FibState(Nmemorymax,Pmod%Nfib,6))
allocate(Saving%LinkState(Nmemorymax,Pmod%Nfib*(Pmod%Nfib-1)/2 + Pmod%Nfib*Pmod%Ndirtot, 4))
allocate(Saving%VariousReal(Nmemorymax,2))
allocate(Saving%VariousInt(Nmemorymax,4))




!-------- Création du dossier de sauvegarde -----------------------------------------------------------------------
call date_and_time(realdate,realtime,timezone,timevalues)

if (Psim%InitType.eq.2) then
	nomdossier_old = Psim%nomdossier
	Psim%nomdossier = trim(adjustl(Psim%nomdossier))//'_continued_'//realdate(3:4)&
		//'.'//realdate(5:6)//'.'//realdate(7:8)//'_'//realtime(1:2)//'.'//realtime(3:4)
else
	Psim%nomdossier = trim(adjustl(Psim%nomdossier))//'_'//realdate(3:4)&
		//'.'//realdate(5:6)//'.'//realdate(7:8)//'_'//realtime(1:2)//'.'//realtime(3:4)
end if

call system('mkdir '//trim(adjustl(Psim%nomdossier)))
call system('cp PARAMETERS.txt '//trim(adjustl(Psim%nomdossier))//'/')
	
open(101,file=trim(adjustl(Psim%nomdossier))//'/Quantif.dat',action='write')
write(101,*) "# time, real elapsed time, internal links, external links, Ncell, number of step since previous saving"
close(101)




!-------- Initialisation des agents -------------------------------------------------------------------------------
if (Psim%InitType==0) then
	! Initialize randomly
	call RandomInitialization(Xcell,Rcell,Xfib,omega,link,distpoints,Pmod,Pdom,BoxCell,BoxFib)
	! Save the initial state thus generated
	call WriteInput(Xcell,Rcell,Xfib,omega,link,distpoints,Pmod,Pdom,Psim)
	
elseif (Psim%InitType==1) then
	! Initialize by reading data from pre-existing files
	call LectureInput(Xcell,Rcell,Xfib,omega,link,distpoints,Pmod,Pdom,BoxCell,BoxFib)
	! Save the initial state used
	call WriteInput(Xcell,Rcell,Xfib,omega,link,distpoints,Pmod,Pdom,Psim)

elseif (Psim%InitType==2) then
	! Initialize by reading data from end of another simulation
	call LecturePreviousSim(Xcell,Rcell,Xfib,omega,link,distpoints,NbStepInPreviousSim,nomdossier_old,Pmod,Pdom,BoxCell,BoxFib)

end if

if (Pmod%Ndirtot.ne.0) then
	open(2,file=trim(adjustl(Psim%nomdossier))//'/dirichlet_layer.dat',action='write')
	write(2,*) "# Number of Dirichlet's fibers per box =", int(Pmod%DirFillingRate * Pdom%dx*Pdom%dy*Pdom%dz)
	write(2,*) "# Total number of Dirichlet's fibers =", Pmod%Ndirtot
	do i=Pmod%Nfib+1,Pmod%Nfib+Pmod%Ndirtot
		write(2,'(6(e23.17e2,1x))') Xfib(i,1), Xfib(i,2), Xfib(i,3), omega(i,1), omega(i,2), omega(i,3)
	end do
	close(2)
	
	write(*,*) "Total number of Dirichlet's fibers = ", Pmod%Ndirtot
end if




!-------- Initialisation du pas de temps --------------------------------------------------------------------------
dt_min = 0.001

dt_max = Psim%dt
if (Pmod%Ncell.lt.Pmod%Nmax) then
	dt_max = min(dt_max,0.5/Pmod%freq_ens)
end if
if (Pmod%LinkType.eq.0) then
	dt_max = min(dt_max,0.5/Pmod%freq_link,0.5/Pmod%freq_dlink)
end if




!------------------------------------------------------------------------------------------------------------------------
!--------------------------------- LOOP IN TIME -------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------
if (Psim%InitType.eq.2) then
	time = NbStepInPreviousSim*Psim%Tenr
	Nsaveblock = int(NbStepInPreviousSim/Nmemorymax)
	Nsavestep0 = NbStepInPreviousSim + 1 - Nsaveblock*Nmemorymax
else
	time = 0
	Nsaveblock = 0
	Nsavestep0 = 1
end if
Nsavestep = Nsavestep0
Nbetween = 0



tcpu0 = omp_get_wtime()
do while (time<=Psim%Tf)
!---------------------------- REINITIALIZE --------------------------------------------------------------	
	gradC(:,:) = 0.0
	gradF(:,:) = 0.0

!---------------------------- COMPUTE FORCES EXERTED ON THE CELLS ---------------------------------------
!$OMP PARALLEL DO private(l,l1,l2,l3,ltemp,k,Xtemp,gradtemp)
	do i=1,Pmod%Ncell
	
		l = BoxCell%AgentBox(i)
		do l1 = -1,1,1 ! for each neighbor k of cell i, compute interaction
			do l2 = -1,1,1
				do l3 = -1,1,1
					ltemp = l + l1 + 2*l2*(Pdom%Nx+1) + 4*l3*(Pdom%Nx+1)*(Pdom%Ny+1) ! depends if we have ghost agents or not

					if ( (ltemp.lt.1).or.(ltemp.gt.8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)) ) then
						write(*,*) "BoxCell number is out of range at motion-step : time = ", time
						write(*,*) "Cell ", i, " is in box ", l, " (with l1 =", l1, " l2 =", l2, " l3 =", l3, " and ltemp = ", ltemp, ")."
						write(*,*) "Its coordinates are ", Xcell(i,:)
						STOP
					end if
				
					k = BoxCell%FirstAgent(ltemp)
					do while (k /= 0)  ! k is a ghost/real cell different from i
						if (i/=k) then
							! Take into account periodic boundary condition 
							if (BoxCell%BoxType(ltemp,1).eq.1) then
								Xtemp(:) = Xcell(k,:) + BoxCell%BoxType(ltemp,2:4)
							else
								Xtemp(:) = Xcell(k,:)
							end if
												
							! Compute interaction force
							call grad_cellcell(gradC(i,:),Xcell(i,:),Xtemp,Rcell(i),Rcell(k),Pmod)
						end if
						k = BoxCell%VerletList(k)
					end do
					
					k = BoxFib%FirstAgent(ltemp)
					do while (k /= 0)
						if (k<=Pmod%Nfib) then ! k is a moving fiber (ghost or real)
							! Take into account periodic boundary condition 
							if (BoxFib%BoxType(ltemp,1).eq.1) then
								Xtemp(:) = Xfib(k,:) + BoxFib%BoxType(ltemp,2:4)
							else
								Xtemp(:) = Xfib(k,:)
							end if
												
							! Compute interaction force
							gradtemp(:) = 0.0
							call grad_cellfib(gradtemp,Xcell(i,:),Xtemp,Rcell(i),omega(k,:),Pmod)
							gradC(i,:) = gradC(i,:) - gradtemp(1:3)
							
						elseif (k>Pmod%Nfib) then ! k is a fixed fiber from the Dirichlet layer
							gradtemp(:) = 0.0
							call grad_cellfib(gradtemp,Xcell(i,:),Xfib(k,:),Rcell(i),omega(k,:),Pmod)
							gradC(i,:) = gradC(i,:) - gradtemp(1:3)
						end if
						
						k = BoxFib%VerletList(k)
					end do
				
				end do
			end do
		end do
	end do
!$OMP END PARALLEL DO
	
	
	
!---------------------------- COMPUTE FORCES EXERTED ON THE FIBERS --------------------------------------
!$OMP PARALLEL DO private(l,l1,l2,l3,ltemp,k,Xtemp)
	do i=1,Pmod%Nfib
		if (IsInDirichlet(i).ne.0) goto 999 ! fiber i has entered the Dirichlet layer and is now motionless
		
		l = BoxFib%AgentBox(i)
		do l1 = -1,1,1 ! for each neighbor k of fiber i, compute fiber-fiber interaction
			do l2 = -1,1,1
				do l3 = -1,1,1
					ltemp = l + l1 + 2*l2*(Pdom%Nx+1) + 4*l3*(Pdom%Nx+1)*(Pdom%Ny+1) ! depends if we have ghost agents or not

					if ( (ltemp.lt.1).or.(ltemp.gt.8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)) ) then
						write(*,*) "BoxFib number is out of range at motion-step : time = ", time
						write(*,*) "Fiber ", i, " is in box ", l, " (with l1 =", l1, " l2 =", l2, " l3 =", l3, " and ltemp = ", ltemp, ")."
						write(*,*) "Its coordinates are ", Xfib(i,:)
						STOP
					end if
					
					k = BoxFib%FirstAgent(ltemp)
					do while (k /= 0)
						if ( (k<=Pmod%Nfib).and.(i/=k) ) then ! k is a moving fiber (ghost or real) different from i
							! Take into account periodic boundary condition 
							if (BoxFib%BoxType(ltemp,1).eq.1) then
								Xtemp(:) = Xfib(k,:) + BoxFib%BoxType(ltemp,2:4)
							else
								Xtemp(:) = Xfib(k,:)
							end if
												
							! Compute interaction force
							call grad_fibfib(gradF(i,:),Xfib(i,:),Xtemp,omega(i,:),omega(k,:),i,k,Pmod,link,distpoints)
														
						elseif (k > Pmod%Nfib) then ! k is a fixed fiber from the Dirichlet layer
							call grad_fibfib(gradF(i,:),Xfib(i,:),Xfib(k,:),omega(i,:),omega(k,:),i,k,Pmod,link,distpoints)
						end if
						
						k = BoxFib%VerletList(k)
					end do
					
					k = BoxCell%FirstAgent(ltemp)
					do while (k /= 0)  ! k is a ghost/real cell
						! Take into account periodic boundary condition 
						if (BoxCell%BoxType(ltemp,1).eq.1) then
							Xtemp(:) = Xcell(k,:) + BoxCell%BoxType(ltemp,2:4)
						else
							Xtemp(:) = Xcell(k,:)
						end if
						! Compute interaction force
						call grad_cellfib(gradF(i,:),Xtemp,Xfib(i,:),Rcell(k),omega(i,:),Pmod)						
						k = BoxCell%VerletList(k)
					end do
					
				end do
			end do
		end do
999		CONTINUE
	end do
!$OMP END PARALLEL DO
	
	
	
!---------------------------- ADAPTATIVE TIME-STEP ----------------------------------------------------
	dt = dt_max
	do i=1,Pmod%Ncell
		dt = min(dt, 0.5*Rcell(i)**2/norm(gradC(i,:)) )
	end do
	do i=1,Pmod%Nfib
		dt = min(dt, 0.5*Pmod%Rf*Pmod%Lf/norm(gradF(i,1:3)), 0.1*Pmod%Lf**3/norm(gradF(i,4:6)))
	end do
	dt = max(dt,dt_min)
	
	tnext = (Nsavestep + Nsaveblock*Nmemorymax)*Psim%Tenr
	if ((time < tnext).and.(time + dt > tnext)) then
		dt = tnext - time
	end if
	time = time + dt
	Nbetween = Nbetween + 1
	
	
			
!---------------------------------- MOTION ------------------------------------------------------------
	do i=1,Pmod%Ncell
		Xtemp(1) = Xcell(i,1) + dt*gradC(i,1)/Rcell(i)
		Xtemp(2) = Xcell(i,2) + dt*gradC(i,2)/Rcell(i)
		Xtemp(3) = Xcell(i,3) + dt*gradC(i,3)/Rcell(i)
		
		call HandleCellBoundaryCondition(Xtemp,Pmod,Pdom,BoxCell)
		
		Xcell(i,:) = Xtemp(:)
	end do
	
	do i=1,Pmod%Nfib
		Xtemp(1) = Xfib(i,1) + dt*gradF(i,1)/Pmod%Lf
		Xtemp(2) = Xfib(i,2) + dt*gradF(i,2)/Pmod%Lf
		Xtemp(3) = Xfib(i,3) + dt*gradF(i,3)/Pmod%Lf
		
		call HandleFibBoundaryCondition(Xtemp,IsInDirichlet(i),Pmod,Pdom,BoxFib)
		
		Xfib(i,:) = Xtemp(:)
		
		omega(i,1) = omega(i,1) + dt*gradF(i,4)/Pmod%Lf**3
		omega(i,2) = omega(i,2) + dt*gradF(i,5)/Pmod%Lf**3
		omega(i,3) = omega(i,3) + dt*gradF(i,6)/Pmod%Lf**3
		omega(i,:) = omega(i,:)/norm(omega(i,:))
	end do
	
	
	
!------------------------------- SAVE CURRENT STATE TO RAM --------------------------------------------
	if (time.ge.(Nsavestep + Nsaveblock*Nmemorymax)*Psim%Tenr) then
!		call SaveCurrentState(Saving,Xcell,Rcell,Xfib,omega,link,distpoints,Pmod)
		do i = 1,Pmod%Ncell
			Saving%CellState(Nsavestep,i,:) = (/ Xcell(i,:), Rcell(i) /)
		end do
		
		do i = 1,Pmod%Nfib
			Saving%FibState(Nsavestep,i,:) = (/ Xfib(i,:), omega(i,:) /)
		end do
		
		k = 0
		do i = 1,Pmod%Nfib
			do j = i+1,Pmod%Nfib+Pmod%Ndirtot
				if (link(i,j)==1) then
					k = k+1
					Saving%LinkState(Nsavestep,k,:) = (/ i*1d0,j*1d0,distpoints(i,j),distpoints(j,i) /)
				endif
			end do
		end do
		
		Nlinks_int = sum(link(:,1:Pmod%Nfib))/2.0
		Nlinks_ext = sum(link(:,Pmod%Nfib+1:))
		tcpuf = omp_get_wtime()
		Saving%VariousReal(Nsavestep,:) = (/time, tcpuf - tcpu0 /)
		Saving%VariousInt(Nsavestep,:) = (/ Nlinks_int, Nlinks_ext, Pmod%Ncell, Nbetween /)

		tcpu0 = omp_get_wtime()
		Nsavestep = Nsavestep + 1
		Nbetween = 0
	endif
	
	
	
!--------------------------------- DUMP RAM INTO SAVE FILES -------------------------------------------
	if (Nsavestep.gt.Nmemorymax) then
		call WriteOutput_manyfiles(Nsavestep0,Nmemorymax,Nsaveblock*Nmemorymax,Saving,Pmod,Psim)
	
		open(101,file=trim(adjustl(Psim%nomdossier))//'/Quantif.dat',action='write',position='append')
		do i = Nsavestep0,Nmemorymax
			write(101,5000) Saving%VariousReal(i,:), Saving%VariousInt(i,:)
		end do
		close(101)
		
		Nsavestep0 = 1
		Nsavestep  = 1
		Nsaveblock = Nsaveblock+1

		Saving%CellState(:,:,:) = 0
		Saving%FibState(:,:,:)  = 0
		Saving%LinkState(:,:,:) = 0
		Saving%VariousReal(:,:) = 0
		Saving%VariousInt(:,:)  = 0
	end if
	
	
	
	
!---------------------------------- CELL GROWTH -------------------------------------------------------
	do i = 1,Pmod%Ncell
		Rcell(i) = (Rcell(i)**3 + dt*Pmod%Kcroiss)**(1./3.)
		if (Rcell(i).gt.Pmod%Rmax) Rcell(i) = Pmod%Rmax
	end do
	
	
	
!----------------------------- CELL INSEMINATION ------------------------------------------------------
	call random_number(alea)
	if ((alea < 1-exp(-Pmod%freq_ens*dt)).and.(Pmod%Ncell<Pmod%Nmax)) then
		call insemination(Xcell,Rcell,Pmod,Pdom)
		write(*,2000) time, Pmod%Ncell

		if (Pmod%Ncell.eq.Pmod%Nmax) then
			t_end_ins = time
			if (Pmod%LinkType.eq.0) then
				dt_max = min(Psim%dt,0.5/Pmod%freq_link,0.5/Pmod%freq_dlink)
			else
				dt_max = Psim%dt
			end if
		end if
	endif
	
	
			
!---------------------------- CLEAN BOXES AND REATTRIBUTE ---------------------------------------------
	BoxCell%FirstAgent(:) = 0
	BoxCell%LastAgent(:)  = 0
	BoxCell%AgentBox(:)   = 0
	BoxCell%VerletList(:) = 0
	
	do i = 1,Pmod%Ncell
		call UpdateVerlet(i,Xcell(i,:),Pdom,BoxCell)
	end do
	call FillGhostBoxes(Pmod,Pdom,BoxCell)



	call cleanInsideDomain(Pmod,Pdom,BoxFib)
	do i = 1,Pmod%Nfib
		if (IsInDirichlet(i).ne.2) then
			BoxFib%AgentBox(i)   = 0
			BoxFib%VerletList(i) = 0
		end if
	end do
	
	do i = 1,Pmod%Nfib
		if (IsInDirichlet(i).ne.2) then
			call UpdateVerlet(i,Xfib(i,:),Pdom,BoxFib)
		end if
	end do
	call FillGhostBoxes(Pmod,Pdom,BoxFib)
	
	
	
!----------------------------- LINKS' UPDATE ---------------------------------------------------------
	if (Pmod%LinkType.ne.0) goto 881 ! skip link actualization for all fibers if LinkType /= 0	

!$OMP PARALLEL DO private(IsANeighbor,l,l1,l2,l3,ltemp,k,Xtemp,dik,lik,lki,alea)
	do i=1,Pmod%Nfib
		if (IsInDirichlet(i).eq.0) then ! 1- it is meaningless to actualize links between two motionless fibers and 2- fibers in a dirichlet box do not have neighbor in all directions.	
			IsANeighbor(:) = 0
			l = BoxFib%AgentBox(i)
			do l1 = -1,1,1 ! search all neighbors k of fiber i
				do l2 = -1,1,1
					do l3 = -1,1,1
						ltemp = l + l1 + 2*l2*(Pdom%Nx+1) + 4*l3*(Pdom%Nx+1)*(Pdom%Ny+1) ! depends if we have ghost agents or not
        	
						if ( (ltemp.lt.1).or.(ltemp.gt.8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)) ) then
							write(*,*) "BoxFib number is out of range at linking-step : time = ", time
							write(*,*) "Fiber ", i, " is in box ", l, " ( with l1 =", l1, " l2 =", l2, " l3 =", l3, " and ltemp = ", ltemp, ")"
							write(*,*) "its coordinates are ", Xfib(i,:)
							STOP
						end if
					
						k = BoxFib%FirstAgent(ltemp)
						do while (k /= 0)
							IsANeighbor(k) = 1
							if ( (k>i).or.(IsInDirichlet(k).ne.0) ) then ! k is a moving fiber > i or a dirichlet fiber
								! Take into account periodic boundary condition 
								if (BoxFib%BoxType(ltemp,1).eq.1) then
									Xtemp(:) = Xfib(k,:) + BoxFib%BoxType(ltemp,2:4)
								else
									Xtemp(:) = Xfib(k,:)
								end if
						
								call find_overlap_fibfib(Xfib(i,:),Xtemp,omega(i,:),omega(k,:),Pmod,dik,lik,lki)
								call random_number(alea)
								if (( link(i,k)==0 ).and.( dik.lt.2*Pmod%Rf+Pmod%eps ).and.( alea<1-exp(-Pmod%freq_link*dt) )) then
									link(i,k) = 1
									link(k,i) = 1
									distpoints(i,k) = lik
									distpoints(k,i) = lki
								elseif (( link(i,k)==1 ).and.( alea<1-exp(-Pmod%freq_dlink*dt) )) then
									link(i,k) = 0
									link(k,i) = 0
									distpoints(i,k) = 0
									distpoints(k,i) = 0
								end if
							end if
							k = BoxFib%VerletList(k)
						end do
					end do
				end do
			end do
			
			! Erase links between fibers that are not in neighboring boxes anymore
			where (IsANeighbor(1:Pmod%Nfib)==0)	
				link(i,1:Pmod%Nfib) = 0
				link(:,i) = 0
				distpoints(i,1:Pmod%Nfib) = 0
				distpoints(:,i) = 0
			end where
			where (IsANeighbor(Pmod%Nfib+1:)==0)
				link(i,Pmod%Nfib+1:) = 0
				distpoints(i,Pmod%Nfib+1:) = 0
			end where	
			
		end if
	end do
!$OMP END PARALLEL DO
		
881	CONTINUE
	
	
	
	
end do





!------------------------------------------------------------------------------------------------------------------
!------------------------------ FINALISATION ----------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
tcpuglobf  = omp_get_wtime()
tcpu_days  = aint( (tcpuglobf - tcpuglob0)/86400 )
tcpu_hours = aint( (tcpuglobf - tcpuglob0)/3600 - tcpu_days*24 )
tcpu_secs  = tcpuglobf - tcpuglob0 - tcpu_days*86400 - tcpu_hours*3600

call WriteOutput_manyfiles(Nsavestep0,Nsavestep-1,Nsaveblock*Nmemorymax,Saving,Pmod,Psim)
	
open(101,file=trim(adjustl(Psim%nomdossier))//'/Quantif.dat',action='write',position='append')
do i = Nsavestep0,Nsavestep-1,1
	write(101,5000) Saving%VariousReal(i,:), Saving%VariousInt(i,:)
end do
tcpu_write  = omp_get_wtime()
write(101,*)
write(101,*)
write(101,*) "# Total time of computation : ",tcpu_days," days, ",tcpu_hours," hours and ",tcpu_secs," seconds."
write(101,*) "# Total time of computation (in seconds) : ",tcpuglobf - tcpuglob0
write(101,*) "# Time needed to write data (in seconds) : ",tcpu_write - tcpuglobf
write(101,*) "# Internal time of the simulation at the end of the insemination proccess : ",t_end_ins
close(101)



2000 FORMAT("New insemination at time ",f10.4,", number of cells is now : ",i5)
5000 FORMAT(f10.4,1x,e23.17e2,4(1x,i15))


deallocate(Xcell)
deallocate(Rcell)
deallocate(gradC)

deallocate(Xfib)
deallocate(omega)
deallocate(gradF)
deallocate(IsInDirichlet)

deallocate(link)
deallocate(distpoints)
deallocate(IsANeighbor)

deallocate(BoxCell%VerletList)
deallocate(BoxCell%AgentBox)
deallocate(BoxCell%FirstAgent)
deallocate(BoxCell%LastAgent)
deallocate(BoxFib%VerletList)
deallocate(BoxFib%AgentBox)
deallocate(BoxFib%FirstAgent)
deallocate(BoxFib%LastAgent)

deallocate(Saving%CellState)
deallocate(Saving%FibState)
deallocate(Saving%LinkState)
deallocate(Saving%VariousReal)
deallocate(Saving%VariousInt)

end program main