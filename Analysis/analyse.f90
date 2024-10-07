program main

use definitions
use lecture_ecriture
use fonctions
use jacobi
use omp_lib

implicit none

!------------------------------------------------------------------------------------------------------------------------!
!---DECLARATION VARIABLES------------------------------------------------------------------------------------------------!
!---------UTILITAIRES----------------------------------------------------------------------------------------------------!
double precision, dimension(:,:), allocatable :: Xfib, omega
double precision                              :: dij, lij, lji, Xtemp(3), dist
integer                                       :: i,j,k,l,l1,l2,l3,ltemp,lreal
double precision                              :: Rneigh, chi_theo, chi_real
integer                                       :: Nneigh_threshold = 40
integer, dimension(:), allocatable            :: Nlinks, Nrep, Ncross, Nneigh
double precision                              :: Nlinks_mean, Nlinks_std, Nrep_mean, Nrep_std, Ncross_mean, Ncross_std
double precision                              :: Nneigh_mean, Nneigh_std
integer                                       :: Nlinks_min, Nlinks_max, Nrep_min, Nrep_max, Ncross_min, Ncross_max
integer                                       :: Nneigh_min, Nneigh_max
integer                                       :: Nlinks_tot, Nrep_tot, Ncross_tot, Nfat ! Nfat = Nb of fibers above threshold (i.e. with more neighbours than Nneigh_threshold)
double precision, dimension(:), allocatable   :: Al
double precision                              :: Al_mean, Al_std, Al_min, Al_max
double precision                              :: Mat_proj(3,3), Pval(3), Pvec(3,3), Pval_mean  ! Projection matrix, its eigenvalues and its eigenvectors
double precision                              :: Intermed(1,3) ! Intermediate for matrix product
integer                                       :: it_num, rot_num ! For the subroutine Jacobi_eigenvalue
!--------MODELE et DOMAINE-----------------------------------------------------------------------------------------------!
type(ParametersMod)                           :: Pmod
type(ParametersDom)                           :: Pdom
type(ParametersBox)                           :: BoxFib
type(ParametersSim)                           :: Psim
!--------ENREGISTREMENT DONNEES------------------------------------------------------------------------------------------!
integer                                       :: index0, indexmax
character(len=30)                             :: ch_iter
character(len=200)                            :: name_fibers, name_links
character(8)                                  :: skip
integer                                       :: ios = 0



!------------------------------------------------------------------------------------------------------------------
!---------------------------- INITIALISATION ----------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
call omp_set_num_threads(10)

!------- Lecture des paramètres -----------------------------------------------------------------------------------
write(*,*) "Entrer le chemin d'accès au dossier contenant les données :"
read(*,'(a)') Psim%nomdossier

call LectureParam(Pmod,Psim,Pdom,trim(adjustl(Psim%nomdossier))//"/PARAMETERS.txt")
index0 = 0
indexmax = int(Psim%Tf/Psim%Tenr)
Rneigh = 2*Pmod%Rf + max(Pmod%eps, 2*Pmod%Rf, 2*Pmod%Rmax)




!------- Création du domaine --------------------------------------------------------------------------------------
if (Pmod%Nmax.gt.0) then
	dist = max(Pmod%Lf + 2*Pmod%Rf + Pmod%eps, Pmod%Lf + 1.5*Pmod%d0, 2*Pmod%Rmax) ! minimal size of a box = maximal length of interaction
else
	dist = max(Pmod%Lf + 2*Pmod%Rf + Pmod%eps, Pmod%Lf + 1.5*Pmod%d0)
end if

Pdom%Nx = int(Pdom%xmaxR*1.0/dist) ! number of boxes that can fit in the dimension x
Pdom%Ny = int(Pdom%ymaxR*1.0/dist)
Pdom%Nz = int(Pdom%zmaxR*1.0/dist)

Pdom%dx = Pdom%xmaxR*1.0/Pdom%Nx ! size of the boxes in the dimension x
Pdom%dy = Pdom%ymaxR*1.0/Pdom%Ny
Pdom%dz = Pdom%zmaxR*1.0/Pdom%Nz

Pdom%xmax = (Pdom%Nx+1)*Pdom%dx ! limit of the domain in x including ghost boxes
Pdom%ymax = (Pdom%Ny+1)*Pdom%dy ! limit of the domain in y including ghost boxes
Pdom%zmax = (Pdom%Nz+1)*Pdom%dz ! limit of the domain in z including ghost boxes




!------- Couche de Dirichlet --------------------------------------------------------------------------------------
allocate(BoxFib%BoxType(1:8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1),4))
call computeBoxType(Pdom,Pmod,BoxFib)

Pmod%Ndirtot = 0
do i = 1,8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)
	if (BoxFib%BoxType(i,1).eq.3) then
		Pmod%Ndirtot = Pmod%Ndirtot + 1
	end if
end do
Pmod%Ndirtot = Pmod%Ndirtot * int(Pmod%DirFillingRate * Pdom%dx*Pdom%dy*Pdom%dz)




!------- Allocation des variables ---------------------------------------------------------------------------------
allocate(Nlinks(1:Pmod%Nfib))
allocate(Nrep(1:Pmod%Nfib))
allocate(Ncross(1:Pmod%Nfib))
allocate(Nneigh(1:Pmod%Nfib))
allocate(Al(1:Pmod%Nfib))

allocate(Xfib(1:Pmod%Nfib+Pmod%Ndirtot,3))
allocate(omega(1:Pmod%Nfib+Pmod%Ndirtot,3))

allocate(BoxFib%FirstAgent(1:8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)))
allocate(BoxFib%LastAgent(1:8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)))
allocate(BoxFib%AgentBox(1:Pmod%Nfib+Pmod%Ndirtot))
allocate(BoxFib%VerletList(1:Pmod%Nfib+Pmod%Ndirtot))

BoxFib%FirstAgent(:) = 0
BoxFib%LastAgent(:)  = 0
BoxFib%AgentBox(:)   = 0
BoxFib%VerletList(:) = 0




!-------- Initialisation des agents -------------------------------------------------------------------------------
if (Pmod%Ndirtot.ne.0) then
	open(unit = 2, file = trim(adjustl(Psim%nomdossier))//'/dirichlet_layer.dat',status='old',action='read')
	read(2,*) skip
	read(2,*) skip
	do i=Pmod%Nfib+1,Pmod%Nfib+Pmod%Ndirtot
		read(2,'(6(e23.17e2,1x))') Xfib(i,1), Xfib(i,2), Xfib(i,3), omega(i,1), omega(i,2), omega(i,3)
	end do
	close(2)
end if




!-------- Initialisation du fichier de sauvegarde -----------------------------------------------------------------
open(unit = 10, file = trim(adjustl(Psim%nomdossier))//'/analyse_alltime.dat')
write(10,*) "# time, Al (moy, et, min, max), Nlinks (moy, et, min, max), Nrep (moy, et, min, max),",&
	" Ncross (moy, et, min, max), Nneigh (moy, et, min, max), Nlinks_tot, Nrep_tot, Ncross_tot,",&
	" Nb of fibers with more than ",Nneigh_threshold," neighbours (absolute, percentage), chi_real"




!------------------------------------------------------------------------------------------------------------------
!--------------------------------- LOOP IN TIME -------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
do k = index0,indexmax,1
	
!---------------------------- Reinitialize all quantifiers --------------------------------------------------------	
	Nlinks(:) = 0
	Nrep(:)   = 0
	Ncross(:) = 0
	Nneigh(:) = 0
	Al(:)     = -1 ! Default value for fibers with less than the threshold number of neighbours
	
	Nlinks_tot  = 0
	Nneigh_mean = 0
	Nneigh_std  = 0
	Nfat = 0
	Al_mean = 0
	Al_std  = 0
	Al_min  = 1
	

	

!---------------------------- Find data files name ----------------------------------------------------------------
	write(ch_iter,*) k
	if (k==0) then
		name_fibers = trim(adjustl(Psim%nomdossier))//"/INPUTFIBERS.dat"
		name_links  = trim(adjustl(Psim%nomdossier))//"/INPUTLINKS.dat"
	else
		name_fibers = trim(adjustl(Psim%nomdossier))//'/fibers'//trim(adjustl(ch_iter))//'.dat'
		if (Pmod%LinkType.eq.0) then
			name_links = trim(adjustl(Psim%nomdossier))//'/links'//trim(adjustl(ch_iter))//'.dat'
		else
			name_links = trim(adjustl(Psim%nomdossier))//"/INPUTLINKS.dat"
		end if	
	end if
	
	
	
	
!---------------------------- Read fibers' data files and attribute boxes------------------------------------------
	open(2,file=name_fibers,status='old',action='read')
	do i=1,Pmod%Nfib
		read(2,'(6(e23.17e2,1x))') Xfib(i,:), omega(i,:)
	end do
	close(2)

	BoxFib%FirstAgent(:) = 0
	BoxFib%LastAgent(:)  = 0
	BoxFib%AgentBox(:)   = 0
	BoxFib%VerletList(:) = 0
	do i = 1,Pmod%Nfib
		call UpdateVerlet(i,Xfib(i,:),Pdom,BoxFib)
	end do
	call FillGhostBoxes(Pmod,Pdom,BoxFib)
	
	
	
	
!---------------------------- Read links' data files and compute linking statistics -------------------------------
	if ((Pmod%LinkType.eq.0).or.(Pmod%LinkType.eq.1)) then
		open(2,file=name_links,status='old',action='read')
		read(2,*,iostat=ios) i,j,lij,lji
		do while(ios==0)
			Nlinks_tot = Nlinks_tot + 1
			Nlinks(i) = Nlinks(i) + 1
			if(j<=Pmod%Nfib) then
				Nlinks(j) = Nlinks(j) + 1
			end if
			
			read(2,*,iostat=ios) i,j,lij,lji
		end do
		close(2)
	end if

    Nlinks_mean = sum(Nlinks)*1.0/Pmod%Nfib
    Nlinks_std  = sqrt( sum(Nlinks**2)*1.0/Pmod%Nfib - Nlinks_mean**2 )
    Nlinks_min  = minval(Nlinks)
    Nlinks_max  = maxval(Nlinks)
	
	
	
	
!---------------------------- Compute quantifiers for each fiber --------------------------------------
	!$OMP PARALLEL DO private(Mat_proj,Intermed,l,l1,l2,l3,ltemp,j,Xtemp,dij,lij,lji, &
	!$OMP& Pvec,Pval,it_num,rot_num,Pval_mean), num_threads(10)
	do i=1,Pmod%Nfib
		Intermed(1,:) = omega(i,:) 
		Mat_proj(:,:) = matmul(transpose(Intermed),Intermed)
		
		l = BoxFib%AgentBox(i)
		do l1 = -1,1,1 ! for each neighbour j of fiber i, compute fiber-fiber interaction
			do l2 = -1,1,1
				do l3 = -1,1,1
					ltemp = l + l1 + 2*l2*(Pdom%Nx+1) + 4*l3*(Pdom%Nx+1)*(Pdom%Ny+1) ! depends if we have ghost agents or not

					if ( (ltemp.lt.1).or.(ltemp.gt.8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)) ) then
						write(*,*) "BoxFib number is out of range at step ", k
						write(*,*) "Fiber ", i, " is in box ", l, " (with l1 =", l1, " l2 =", l2, " l3 =", l3, " and ltemp = ", ltemp, ")."
						write(*,*) "Its coordinates are ", Xfib(i,:)
						STOP
					end if
					
					j = BoxFib%FirstAgent(ltemp)
					do while (j.ne.0)
						if (j.ne.i) then
							! Take into account periodic boundary condition 
							if (BoxFib%BoxType(ltemp,1).eq.1) then
								Xtemp(:) = Xfib(j,:) + BoxFib%BoxType(ltemp,2:4)
							else
								Xtemp(:) = Xfib(j,:)
							end if
						
							! Compute distance between the two fibers
							call find_overlap_fibfib(Xfib(i,:),Xtemp,omega(i,:),omega(j,:),Pmod,dij,lij,lji)
							
							! Overlapping range
		        			if (dij.le.2*Pmod%Rf) then
								Nrep(i) = Nrep(i) + 1
		        			end if
							
							! Cross-linking range
		        			if (dij.le.2*Pmod%Rf+Pmod%eps) then
								Ncross(i) = Ncross(i) + 1
		        			end if
							
							! Correlation matrix
		        			if (dij.le.Rneigh) then
		        			    Nneigh(i) = Nneigh(i) + 1
								Intermed(1,:) = omega(j,:) 
		        			    Mat_proj = Mat_proj + matmul(transpose(Intermed),Intermed)
		        			end if	
						end if
						
						j = BoxFib%VerletList(j)
					end do
					
				end do
			end do
		end do
		
        if (Nneigh(i).ge.Nneigh_threshold) then
			!$OMP CRITICAL (update_Nfat)
			Nfat = Nfat + 1
			Nneigh_mean = Nneigh_mean + Nneigh(i)
			Nneigh_std  = Nneigh_std + Nneigh(i)**2
			!$OMP END CRITICAL (update_Nfat)
			
            Mat_proj = Mat_proj*1.0/(Nneigh(i)+1)
			call jacobi_eigenvalue(3,Mat_proj,20,Pvec,Pval,it_num,rot_num)
            Pval_mean = sum(Pval(:))/3.0
            Al(i) = sqrt(1.5)*sqrt(sum( (Pval(:) - Pval_mean)**2 ))/sqrt(sum( Pval(:)**2 ))
			
			!$OMP CRITICAL (update_Almean)
			Al_mean = Al_mean + Al(i)
			Al_std  = Al_std + Al(i)**2
			Al_min  = min(Al_min, Al(i))
			!$OMP END CRITICAL (update_Almean)
		end if
	end do
	!$OMP END PARALLEL DO
	
	
	
	
!---------------------------- Compute average quantifiers --------------------------------------
	Al_mean = Al_mean/Nfat
	Al_std  = sqrt( Al_std/Nfat - Al_mean**2 )
	Al_max  = maxval(Al)
    	
	Nrep_mean = sum(Nrep)*1.0/Pmod%Nfib
	Nrep_std  = sqrt( sum((Nrep - Nrep_mean)**2)*1.0/Pmod%Nfib )
	Nrep_min  = minval(Nrep)
	Nrep_max  = maxval(Nrep)
    
	Ncross_mean = sum(Ncross)*1.0/Pmod%Nfib
	Ncross_std  = sqrt( sum((Ncross - Ncross_mean)**2)*1.0/Pmod%Nfib )
	Ncross_min  = minval(Ncross)
	Ncross_max  = maxval(Ncross)
    
	Nneigh_mean = Nneigh_mean*1.0/Nfat ! sum(Nneigh)*1.0/Pmod%Nfib
	Nneigh_std  = sqrt( Nneigh_std*1.0/Nfat - Nneigh_mean**2 ) ! sqrt( sum((Nneigh - Nneigh_mean)**2)*1.0/Pmod%Nfib )
	Nneigh_min  = minval(Nneigh)
	Nneigh_max  = maxval(Nneigh)
	
	Nrep_tot  = sum(Nrep)/2.0
	Ncross_tot = sum(Ncross)/2.0
	chi_real = Nlinks_tot*1.0/Ncross_tot
	
	
	
	
!------------------------------------ Save quantifiers ------------------------------------
	write(10,1000) k*Psim%Tenr, Al_mean, Al_std, Al_min, Al_max,&
			Nlinks_mean, Nlinks_std, Nlinks_min, Nlinks_max,&
			Nrep_mean,   Nrep_std,   Nrep_min,   Nrep_max,&
			Ncross_mean, Ncross_std, Ncross_min, Ncross_max,&
			Nneigh_mean, Nneigh_std, Nneigh_min, Nneigh_max,&
			Nlinks_tot, Nrep_tot, Ncross_tot,&
			Nfat, Nfat*100.0/Pmod%Nfib, chi_real
	
	
	
			
!--------------------------------- Save state in csv file ---------------------------------
!!! Comment the if / end if lines if you want the detailed result file to be written for all time-step.
	if (k==indexmax) then
		open(unit = 2, file = trim(adjustl(Psim%nomdossier))//'/fibers'//trim(adjustl(ch_iter))//'_seg.csv',action='write')
		write(2,*)"X,Y,Z,wX,wY,wZ,Lfib,Alfib,Nlinks,Nrep,Ncross,Nneigh"
		do i=1,Pmod%Nfib
			write(2,2000) Xfib(i,:), omega(i,:), Pmod%Lf, Al(i), Nlinks(i), Nrep(i), Ncross(i), Nneigh(i)
		end do
		close(2)
	end if
		
		
			
end do
close(10)








1000 FORMAT(f10.2,1x,f10.8,1x,f12.10,1x,f10.8,1x,f10.8,1x,&
		4(f12.3,1x,f12.5,1x,i8,1x,i8,1x),&
		4(i10,1x),f8.4,1x,f6.4)
		
2000 FORMAT(8(f10.6,", "),4(i4,", "))

!------------------------------------------------------------------------------------------------------------------
!------------------------------ FINALISATION ----------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
deallocate(Xfib,omega)
deallocate(Nlinks,Nrep,Ncross,Nneigh,Al)
deallocate(BoxFib%BoxType)
deallocate(BoxFib%VerletList)
deallocate(BoxFib%AgentBox)
deallocate(BoxFib%FirstAgent)
deallocate(BoxFib%LastAgent)

end program main
