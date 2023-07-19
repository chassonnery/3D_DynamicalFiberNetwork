module fonctions

implicit none

save

contains

function gauss(n) 
! génère un nombre aléatoire choisi selon une loi normale 0,1 avec une "précision" n
	implicit none
	double precision    :: gauss
	integer, intent(in) :: n
	integer             :: i
	double precision    :: x
	
	gauss = 0
	do i=1,n
		call random_number(x)
		gauss = gauss + x
	enddo
	gauss = (gauss-n/2d0)/sqrt(n/12d0)
	
end function gauss


subroutine init_random_seed(yesorno)
! Initialize a pseudo-random number sequence, with an option for repeatability
	implicit none
	integer, intent(in)                :: yesorno
	integer                            :: i, n, clock
	integer, dimension(:), allocatable :: seed
	
	call Random_seed(size = n)
	allocate(seed(n))

	if (yesorno.eq.0) then
		call System_clock(COUNT=clock)
		seed = clock + 37 * (/ (i - 1, i = 1, n) /)
	else
		seed = 51748307 * (/ (i - 1, i = 1, n) /)
	end if
	
	call Random_seed(PUT = seed)
	deallocate(seed)

end subroutine


subroutine pick_and_reject(unif) 
! génère un vecteur aléatoire choisi selon une loi uniforme dans la sphère unité
	implicit none
	double precision, dimension(3), intent(out) :: unif
	integer :: IsInUnitSphere
	
	IsInUnitSphere = 0
	do while(IsInUnitSphere == 0)
		call random_number(unif(1))
		call random_number(unif(2))
		call random_number(unif(3))
		unif(:) = (unif(:) - 0.5)*2
	
		if (norm(unif)<=1.0 .and. norm(unif)>=0.5) then
			IsInUnitSphere = 1
		end if
		unif(:) = unif(:)*1./norm(unif)
	end do
	
end subroutine pick_and_reject


function norm(X)
	implicit none
	double precision, dimension(3), intent(in) :: X
	double precision                           :: norm

	norm = sqrt(X(1)**2+X(2)**2+X(3)**2)
end function


subroutine insemination(Xcell,Rcell,Pmod,Pdom)
use definitions
	implicit none
	double precision, dimension(:,:), intent(inout) :: Xcell
	double precision, dimension(:), intent(inout)   :: Rcell
	type(ParametersMod), intent(inout) :: Pmod
	type(ParametersDom), intent(in)    :: Pdom
	double precision                   :: ratio

	ratio = 1.0
	Pmod%Ncell = Pmod%Ncell + 1
	
	call random_number(Xcell(Pmod%Ncell,1))
	call random_number(Xcell(Pmod%Ncell,2))
	call random_number(Xcell(Pmod%Ncell,3))
	Xcell(Pmod%Ncell,1) = (2.*Pdom%xmaxR*Xcell(Pmod%Ncell,1) - Pdom%xmaxR)/ratio
	Xcell(Pmod%Ncell,2) = (2.*Pdom%ymaxR*Xcell(Pmod%Ncell,2) - Pdom%ymaxR)/ratio
	Xcell(Pmod%Ncell,3) = (2.*Pdom%zmaxR*Xcell(Pmod%Ncell,3) - Pdom%zmaxR)/ratio
	
	Rcell(Pmod%Ncell) = Pmod%Rmin
	
end subroutine


subroutine find_overlap_fibfib(Xi,Xk,wi,wk,Pmod,dik,lik,lki)
! Find the closest points segments and rays of two fibers
! http://geomalgorithms.com/a07-_distance.html
use definitions
	implicit none
	double precision, dimension(3), intent(in) :: Xi, Xk, wi, wk
	type(ParametersMod), intent(in)            :: Pmod
	double precision                           :: norme,prod,s,t
	double precision, dimension(3)             :: x
	double precision, intent(out)              :: dik,lik,lki

	x = Xi-Xk
	prod = (dot_product(wi,wi)*dot_product(wk,wk) - dot_product(wi,wk)**2)
	
	if (prod < 0.001) then ! lines almost parallel
		s = 0
		t = dot_product(wk,x)/dot_product(wk,wk)
	else
		s = (dot_product(wi,wk)*dot_product(wk,x) - dot_product(wk,wk)*dot_product(wi,x))/prod
		t = (dot_product(wi,wi)*dot_product(wk,x) - dot_product(wi,wk)*dot_product(wi,x))/prod
		if (s>=Pmod%Lf/2d0) then
			s = Pmod%Lf/2d0
			t = (dot_product(wk,x)+Pmod%Lf/2d0*dot_product(wk,wi))/dot_product(wk,wk)
		elseif (s<=-Pmod%Lf/2d0) then
			s = - Pmod%Lf/2d0
			t = (dot_product(wk,x)-Pmod%Lf/2d0*dot_product(wk,wi))/dot_product(wk,wk)
		endif
	endif
	
	if (t <= -Pmod%Lf/2d0) then
		t = - Pmod%Lf/2d0 ! recompute s for this edge
		s = - (dot_product(wi,x) + Pmod%Lf/2d0*dot_product(wi,wk))/dot_product(wi,wi)
		if (s >= Pmod%Lf/2d0) then
			s = Pmod%Lf/2d0
		elseif (s <= -Pmod%Lf/2d0) then
			s = -Pmod%Lf/2d0
		end if
	elseif (t >= Pmod%Lf/2d0) then
		t = Pmod%Lf/2d0 ! recompute s for this edge
		s = - (dot_product(wi,x) - Pmod%Lf/2d0*dot_product(wi,wk))/dot_product(wi,wi)
		if (s >= Pmod%Lf/2d0) then
			s = Pmod%Lf/2d0
		elseif (s <= -Pmod%Lf/2d0) then
			s = -Pmod%Lf/2d0
		end if
	endif
	
	dik = norm(Xi + s*wi - (Xk + t*wk))
	lik = s
	lki = t
end subroutine find_overlap_fibfib



subroutine find_overlap_cellfib(Xc,Xf,wf,Pmod,d,t)
! Find point on fiber i closest to the center of cell k
! http://geomalgorithms.com/a07-_distance.html
use definitions
    implicit none
	double precision, dimension(3), intent(in) :: Xc, Xf, wf
	type(ParametersMod), intent(in)            :: Pmod
    double precision, intent(out)              :: d,t

	t = - dot_product(wf,Xf-Xc)/dot_product(wf,wf)
	
	if (t <= -Pmod%Lf/2d0) then
	    t = -Pmod%Lf/2d0
	elseif (t >= Pmod%Lf/2d0) then
	    t = Pmod%Lf/2d0
	endif
	
	d = norm(Xf + t*wf - Xc)

end subroutine find_overlap_cellfib



subroutine grad_cellcell(grad,Xi,Xk,Ri,Rk,Pmod)
use definitions
	implicit none
	double precision, dimension(3), intent(inout) :: grad
	double precision, dimension(3), intent(in) :: Xi,Xk
	double precision, intent(in)               :: Ri,Rk
	type(ParametersMod), intent(in)            :: Pmod
	double precision                           :: dik
	double precision, dimension(3)             :: Vik

	
	dik = norm(Xi - Xk)
	if (dik<Ri+Rk) then ! cells i and k overlap
		if (dik<0.001) then
			call pick_and_reject(Vik)
		else
			Vik = (Xi - Xk)/dik
		end if
		
		grad(:) = grad(:) + Pmod%alpha_repCC*(Ri + Rk)**2 * (1.0 - dik/(Ri + Rk))**(1.5) * Vik(:) ! dimentioned force
	end if
		
end subroutine grad_cellcell



subroutine grad_cellfib(grad,Xi,Xk,Ri,omegak,Pmod)
! Contrary to the function grad_cellcell and grad_fibfib, which are antisymetric with respect to arguments i and k (they return the force/torque exerted by agent k on agent i), in this function the arguments i and k CAN NOT be reverted. 
! This function always return the force and torque exerted by the cell i on the fiber k. To obtain the force exerted by the fiber k on the cell i, call grad_cellfib(gradtemp,Xi,Xk,Ri,omegak,Pmod) then -gradtemp(1:3)
use definitions
	implicit none
	double precision, dimension(6), intent(inout) :: grad
	double precision, dimension(3), intent(in) :: Xi,Xk,omegak
	double precision, intent(in)               :: Ri
	type(ParametersMod), intent(in)            :: Pmod
	double precision                           :: dik,lki
	double precision, dimension(3)             :: Vki,torque_rep

	
	call find_overlap_cellfib(Xi,Xk,omegak,Pmod,dik,lki)

	if (dik < Ri+Pmod%Rf) then ! cell i and fiber k overlap
		if (dik<0.001) then
			call pick_and_reject(Vki)
		else
			Vki = (Xk + lki*omegak - Xi)/dik
		end if
		
		torque_rep(1) = lki*( Vki(1)*(omegak(2)**2 + omegak(3)**2) - Vki(2)*omegak(1)*omegak(2) - Vki(3)*omegak(1)*omegak(3) )  
		torque_rep(2) = lki*( Vki(2)*(omegak(1)**2 + omegak(3)**2) - Vki(1)*omegak(1)*omegak(2) - Vki(3)*omegak(2)*omegak(3) ) 
		torque_rep(3) = lki*( Vki(3)*(omegak(1)**2 + omegak(2)**2) - Vki(1)*omegak(1)*omegak(3) - Vki(2)*omegak(2)*omegak(3) ) 
		
		grad(1) = grad(1) + Pmod%alpha_repCF*(Ri + Pmod%Rf)**2 * (1.0 - dik/(Ri + Pmod%Rf))*(1.5) * Vki(1) ! dimentioned
		grad(2) = grad(2) + Pmod%alpha_repCF*(Ri + Pmod%Rf)**2 * (1.0 - dik/(Ri + Pmod%Rf))*(1.5) * Vki(2)
		grad(3) = grad(3) + Pmod%alpha_repCF*(Ri + Pmod%Rf)**2 * (1.0 - dik/(Ri + Pmod%Rf))*(1.5) * Vki(3)
		grad(4) = grad(4) + Pmod%alpha_repCF*(Ri + Pmod%Rf)**2 * (1.0 - dik/(Ri + Pmod%Rf))*(1.5) * torque_rep(1)
		grad(5) = grad(5) + Pmod%alpha_repCF*(Ri + Pmod%Rf)**2 * (1.0 - dik/(Ri + Pmod%Rf))*(1.5) * torque_rep(2)
		grad(6) = grad(6) + Pmod%alpha_repCF*(Ri + Pmod%Rf)**2 * (1.0 - dik/(Ri + Pmod%Rf))*(1.5) * torque_rep(3)
	end if
		
end subroutine grad_cellfib



subroutine grad_fibfib(grad,Xi,Xk,omegai,omegak,i,k,Pmod,link,distpoints)
use definitions
	implicit none
	double precision, dimension(6), intent(inout) :: grad
	double precision, dimension(:,:), intent(in)  :: link,distpoints
	double precision, dimension(3), intent(in)    :: Xi,Xk,omegai,omegak
	type(ParametersMod), intent(in)               :: Pmod
	integer, intent(in)                           :: i,k
	double precision                              :: norme, dik, lik, lki, c
	double precision, dimension(3)                :: Vik,vpoints,torque_rep,torque_rappel,torque_align,om,v
	double precision, dimension(3,3)              :: R,Rot
	double precision, parameter                   :: pi=3.141592653589793238

	
	call find_overlap_fibfib(Xi,Xk,omegai,omegak,Pmod,dik,lik,lki)
	! reference : Interparticle torques suppress motility-induced phase separation for rodlike particles, Damme, Rondenburg, 2019, arXiv

	if ( (link(i,k)==0).or.(Pmod%cumul_replink.eq.0) ) then
		if (dik<2*Pmod%Rf) then ! i and k overlap
			vpoints = Xi + lik*omegai - (Xk + lki*omegak)
			if (norm(Vpoints)<0.001) then
				call pick_and_reject(Vik)
!				write(*,*) "attention (i,k)= ", i, k
			else
				Vik = vpoints/norm(vpoints)
			end if
			
			torque_rep(1) = lik*( Vik(1)*(omegai(2)**2 + omegai(3)**2) - Vik(2)*omegai(1)*omegai(2) - Vik(3)*omegai(1)*omegai(3) )  
			torque_rep(2) = lik*( Vik(2)*(omegai(1)**2 + omegai(3)**2) - Vik(1)*omegai(1)*omegai(2) - Vik(3)*omegai(2)*omegai(3) ) 
			torque_rep(3) = lik*( Vik(3)*(omegai(1)**2 + omegai(2)**2) - Vik(1)*omegai(1)*omegai(3) - Vik(2)*omegai(2)*omegai(3) ) 
			
			grad(1) = grad(1) + Pmod%alpha_repFF*4.0*Pmod%Rf**2 * (1.0 - dik/(2.0*Pmod%Rf))**(1.5) * Vik(1) ! dimentioned
			grad(2) = grad(2) + Pmod%alpha_repFF*4.0*Pmod%Rf**2 * (1.0 - dik/(2.0*Pmod%Rf))**(1.5) * Vik(2)
			grad(3) = grad(3) + Pmod%alpha_repFF*4.0*Pmod%Rf**2 * (1.0 - dik/(2.0*Pmod%Rf))**(1.5) * Vik(3)
			grad(4) = grad(4) + Pmod%alpha_repFF*4.0*Pmod%Rf**2 * (1.0 - dik/(2.0*Pmod%Rf))**(1.5) * torque_rep(1)
			grad(5) = grad(5) + Pmod%alpha_repFF*4.0*Pmod%Rf**2 * (1.0 - dik/(2.0*Pmod%Rf))**(1.5) * torque_rep(2)
			grad(6) = grad(6) + Pmod%alpha_repFF*4.0*Pmod%Rf**2 * (1.0 - dik/(2.0*Pmod%Rf))**(1.5) * torque_rep(3)
		end if
	end if
		
	if (link(i,k)==1) then
		vpoints = Xi + distpoints(i,k)*omegai - (Xk + distpoints(k,i)*omegak)
		Vik = - (norm(vpoints) - Pmod%d0)*vpoints/norm(vpoints) ! /!\ dimentioned
		
		torque_rappel(1) = distpoints(i,k)*( Vik(1)*(omegai(2)**2 + omegai(3)**2) &
		                   - Vik(2)*omegai(1)*omegai(2) - Vik(3)*omegai(1)*omegai(3) )  
		torque_rappel(2) = distpoints(i,k)*( Vik(2)*(omegai(1)**2 + omegai(3)**2) &
		                   - Vik(1)*omegai(1)*omegai(2) - Vik(3)*omegai(2)*omegai(3) ) 
		torque_rappel(3) = distpoints(i,k)*( Vik(3)*(omegai(1)**2 + omegai(2)**2) &
		                   - Vik(1)*omegai(1)*omegai(3) - Vik(2)*omegai(2)*omegai(3) ) 
		
		! NEMATIC alignment force
		c = dot_product(omegai,omegak)
		if (acos(c) < pi/2d0) then
			v = omegak
		else
			v = -omegak
		endif
		c = dot_product(omegai,v)
		
		om(1) = omegai(2)*v(3) - omegai(3)*v(2)
		om(2) = omegai(3)*v(1) - omegai(1)*v(3)
		om(3) = omegai(1)*v(2) - omegai(2)*v(1)
		Rot=0
		if (abs(1.0 - c**2)>0.01) then
			R = 0
			R(1,2) = -om(3)
			R(1,3) =  om(2)
			R(2,1) =  om(3)
			R(2,3) = -om(1)
			R(3,1) = -om(2)
			R(3,2) =  om(1)
			Rot = R + (1-c)/(norm(om)**2) * matmul(R,R)
		endif
		
		torque_align(1) = Rot(1,1)*omegai(1) + Rot(1,2)*omegai(2) + Rot(1,3)*omegai(3)
		torque_align(2) = Rot(2,1)*omegai(1) + Rot(2,2)*omegai(2) + Rot(2,3)*omegai(3)
		torque_align(3) = Rot(3,1)*omegai(1) + Rot(3,2)*omegai(2) + Rot(3,3)*omegai(3)
		
		grad(1) = grad(1) + Pmod%alpha_rappel*Vik(1)
		grad(2) = grad(2) + Pmod%alpha_rappel*Vik(2)
		grad(3) = grad(3) + Pmod%alpha_rappel*Vik(3)
		grad(4) = grad(4) + Pmod%alpha_rappel*torque_rappel(1) + Pmod%alpha_align*torque_align(1)
		grad(5) = grad(5) + Pmod%alpha_rappel*torque_rappel(2) + Pmod%alpha_align*torque_align(2)
		grad(6) = grad(6) + Pmod%alpha_rappel*torque_rappel(3) + Pmod%alpha_align*torque_align(3)
	endif

end subroutine grad_fibfib



subroutine computeBoxType(Pdom,Pmod,Box)
use definitions
    implicit none
	type(ParametersDom),intent(in)    :: Pdom
	type(ParametersMod),intent(in)    :: Pmod
	type(ParametersBox),intent(inout) :: Box
	integer                           :: ix,iy,iz,l
		
	Box%BoxType(:,:) = 0.0
	! BoxType(i,1) = 0 --> real box
	!				 1 --> periodic ghost box
	!				 2 --> must-be-empty box (mixed border condition)
	!				 3 --> dirichlet box
	
	! Types are computed by increasing order of priority : if a box is on the edge between a dirichlet (3) or must-be-empty (2) border and a periodic border (1), then it is of type periodic (1). If a box is on the edge between a must-be-empty (2) and a dirichlet border (3), then it is of type must-be-empty (2).
	
	do ix = 0,2*Pdom%Nx+1
		do iy = 0,2*Pdom%Ny+1
			do iz = 0,2*Pdom%Nz+1
				l = 1 + ix + 2*iy*(Pdom%Nx+1) + 4*iz*(Pdom%Nx+1)*(Pdom%Ny+1)
				
				! Dirichlet boxes
				if ( (ix.eq.0).and.(Pmod%DirLim(1).eq.1) ) then
					Box%BoxType(l,1) = 3
				elseif ( (ix.eq.2*Pdom%Nx+1).and.(Pmod%DirLim(2).eq.1) ) then
					Box%BoxType(l,1) = 3
				end if
				
				if ( (iy.eq.0).and.(Pmod%DirLim(3).eq.1) ) then
					Box%BoxType(l,1) = 3
				elseif ( (iy.eq.2*Pdom%Ny+1).and.(Pmod%DirLim(4).eq.1) ) then
					Box%BoxType(l,1) = 3
				end if
				
				if ( (iz.eq.0).and.(Pmod%DirLim(5).eq.1) ) then
					Box%BoxType(l,1) = 3
				elseif ( (iz.eq.2*Pdom%Nz+1).and.(Pmod%DirLim(6).eq.1) ) then
					Box%BoxType(l,1) = 3
				end if

				! Must-be-empty boxes
				if ( (ix.eq.0).and.(Pmod%DirLim(1).eq.0) ) then
					Box%BoxType(l,1) = 2
				elseif ( (ix.eq.2*Pdom%Nx+1).and.(Pmod%DirLim(2).eq.0) ) then
					Box%BoxType(l,1) = 2
				end if
				if ( (iy.eq.0).and.(Pmod%DirLim(3).eq.0) ) then
					Box%BoxType(l,1) = 2
				elseif ( (iy.eq.2*Pdom%Ny+1).and.(Pmod%DirLim(4).eq.0) ) then
					Box%BoxType(l,1) = 2
				end if
				if ( (iz.eq.0).and.(Pmod%DirLim(5).eq.0) ) then
					Box%BoxType(l,1) = 2
				elseif ( (iz.eq.2*Pdom%Nz+1).and.(Pmod%DirLim(6).eq.0) ) then
					Box%BoxType(l,1) = 2
				end if
				
				! Periodic boxes
				if ( (ix.eq.0).and.( (Pmod%DirLim(1).eq.0).and.(Pmod%DirLim(2).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,2) = -2*Pdom%xmaxR
				else if ( (ix.eq.2*Pdom%Nx+1).and.( (Pmod%DirLim(1).eq.0).and.(Pmod%DirLim(2).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,2) = 2*Pdom%xmaxR
				end if
				
				if ( (iy.eq.0).and.( (Pmod%DirLim(3).eq.0).and.(Pmod%DirLim(4).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,3) = -2*Pdom%ymaxR
				else if ( (iy.eq.2*Pdom%Ny+1).and.( (Pmod%DirLim(3).eq.0).and.(Pmod%DirLim(4).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,3) = 2*Pdom%ymaxR
				end if
				
				if ( (iz.eq.0).and.( (Pmod%DirLim(5).eq.0).and.(Pmod%DirLim(6).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,4) = -2*Pdom%zmaxR
				else if ( (iz.eq.2*Pdom%Nz+1).and.( (Pmod%DirLim(5).eq.0).and.(Pmod%DirLim(6).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,4) = 2*Pdom%zmaxR
				end if
				
			end do
		end do
	end do
end subroutine computeBoxType



subroutine findBox(X,k,Pdom)
! Computes the index k of the box in which point X is located
! NB: this function INCLUDES THE CL (box 1 = left bottom with CL)
use definitions
	implicit none
	double precision, dimension(:), intent(in) :: X
	integer, intent(out)                       :: k
	type(ParametersDom), intent(in)            :: Pdom
	integer                                    :: ix,iy,iz
	
	ix = aint( (X(1)+Pdom%xmax)/Pdom%dx ) ! use fonction aint for the truncature
	iy = aint( (X(2)+Pdom%ymax)/Pdom%dy )
	iz = aint( (X(3)+Pdom%zmax)/Pdom%dz )
	
	k = 1 + ix + 2*iy*(Pdom%Nx+1) + 4*iz*(Pdom%Nx+1)*(Pdom%Ny+1)

end subroutine findBox



subroutine UpdateVerlet(i,X,Pdom,Box)
use definitions
	implicit none
	integer, intent(in)                       :: i
	double precision, dimension(3),intent(in) :: X
	type(ParametersDom),intent(in)            :: Pdom
	type(ParametersBox),intent(inout)         :: Box
	integer                                   :: k

	call findBox(X,k,Pdom)
	
	Box%AgentBox(i) = k
	
	if (Box%FirstAgent(k).eq.0) then ! The box is still empty, begin it with i
		Box%FirstAgent(k) = i 
		Box%LastAgent(k)  = i
	else ! add index i to the end of verletlist and change last agent
		Box%VerletList(Box%LastAgent(k)) = i
		Box%LastAgent(k) = i
	end if

end subroutine UpdateVerlet



subroutine cleanInsideDomain(Pmod,Pdom,Box)
! Function which clean the array FirstAgent and LastAgent for all boxes *inside* the domain 
use definitions
	implicit none
	type(ParametersMod), intent(in)    :: Pmod
	type(ParametersDom), intent(in)    :: Pdom
	type(ParametersBox), intent(inout) :: Box
	integer                            :: k,ix,iy,iz
	
	do ix = 1,2*Pdom%Nx
		do iy = 1,2*Pdom%Ny
			do iz = 1,2*Pdom%Nz
				k = 1 + ix + 2*iy*(Pdom%Nx+1) + 4*iz*(Pdom%Nx+1)*(Pdom%Ny+1)
				Box%FirstAgent(k) = 0
				Box%LastAgent(k) = 0
			end do
		end do
	end do
	
end subroutine cleanInsideDomain



subroutine FillGhostBoxes(Pmod,Pdom,Box)
! Function which fill the "ghost boxes" on the periodic borders of the domain
use definitions
	implicit none
	type(ParametersMod), intent(in)    :: Pmod
	type(ParametersDom), intent(in)    :: Pdom
	type(ParametersBox), intent(inout) :: Box
	integer                            :: k,kreal
	integer                            :: ix,iy,iz,ixreal,iyreal,izreal
	
	do ix = 0,2*Pdom%Nx+1
		do iy = 0,2*Pdom%Ny+1
			do iz = 0,2*Pdom%Nz+1
				
				k = 1 + ix + 2*iy*(Pdom%Nx+1) + 4*iz*(Pdom%Nx+1)*(Pdom%Ny+1)
				! If this box is a "ghost" one, then its first agent is equal to the first agent of its associated "real box" or "Dirichlet box"
				if (Box%BoxType(k,1).eq.1) then
					
					if ( (ix.eq.0).and.( (Pmod%DirLim(1).eq.0).and.(Pmod%DirLim(2).eq.0) )) then 
						ixreal = 2*Pdom%Nx
					elseif ( (ix.eq.2*Pdom%Nx+1).and.( (Pmod%DirLim(1).eq.0).and.(Pmod%DirLim(2).eq.0) )) then
						ixreal = 1
					else
						ixreal = ix
					end if
					
					if ( (iy.eq.0).and.( (Pmod%DirLim(3).eq.0).and.(Pmod%DirLim(4).eq.0) )) then 
						iyreal = 2*Pdom%Ny
					elseif ( (iy.eq.2*Pdom%Ny+1).and.( (Pmod%DirLim(3).eq.0).and.(Pmod%DirLim(4).eq.0) )) then
						iyreal = 1
					else
						iyreal = iy
					end if
					
					if ( (iz.eq.0).and.( (Pmod%DirLim(5).eq.0).and.(Pmod%DirLim(6).eq.0) )) then 
						izreal = 2*Pdom%Nz
					elseif ( (iz.eq.2*Pdom%Nz+1).and.( (Pmod%DirLim(5).eq.0).and.(Pmod%DirLim(6).eq.0) )) then
						izreal = 1
					else
						izreal = iz
					end if
				
					kreal = 1 + ixreal + 2*iyreal*(Pdom%Nx+1) + 4*izreal*(Pdom%Nx+1)*(Pdom%Ny+1)
					Box%FirstAgent(k) = Box%FirstAgent(kreal)
					Box%LastAgent(k) = Box%LastAgent(kreal)
					
					if (kreal.eq.k) then
						write(*,*) "Erreur : BoxType of box ", k, "is badly defined"
						STOP
					end if
				end if
				
			end do
		end do
	end do
	
end subroutine FillGhostBoxes



subroutine HandleCellBoundaryCondition(X,Pmod,Pdom,Box)
use definitions
	implicit none
	double precision, dimension(3), intent(inout) :: X
	type(ParametersMod), intent(in) :: Pmod
	type(ParametersDom), intent(in) :: Pdom
	type(ParametersBox), intent(in) :: Box
	integer                         :: ix, iy, iz, k

	call findBox(X,k,Pdom)
	
	! If a cell exit the computation domain by entering a periodic ghost box : make it come back from the opposite side.
	if (Box%BoxType(k,1).eq.1) then
		X(:) = X(:) - Box%BoxType(k,2:4)
		call findBox(X,k,Pdom)
		if (Box%BoxType(k,1).eq.1) then
			write(*,*) "Error : cell is still in a periodic box after applying boundary condition."
			STOP
		end if
	end if
	
	! For the cells, dirichlet boxes are the same as must-be-empty boxes.
	! If a cell exit the computation domain by entering a must-be-empty box : bring it back inside from the same side.
	if ((Box%BoxType(k,1).eq.2).or.(Box%BoxType(k,1).eq.3)) then
		ix = aint( (X(1)+Pdom%xmax)/Pdom%dx )
		iy = aint( (X(2)+Pdom%ymax)/Pdom%dy )
		iz = aint( (X(3)+Pdom%zmax)/Pdom%dz )
		
		if ( ((ix.eq.0).or.(ix.eq.2*Pdom%Nx+1)) .and. ((Pmod%DirLim(1).eq.1).or.(Pmod%DirLim(2).eq.1)) ) then
			X(1) = sign(1d0,X(1))*0.999*Pdom%xmaxR
		end if
		
		if ( ((iy.eq.0).or.(iy.eq.2*Pdom%Ny+1)) .and. ((Pmod%DirLim(3).eq.1).or.(Pmod%DirLim(4).eq.1)) ) then
			X(2) = sign(1d0,X(2))*0.999*Pdom%ymaxR
		end if
		
		if ( ((iz.eq.0).or.(iz.eq.2*Pdom%Nz+1)) .and. ((Pmod%DirLim(5).eq.1).or.(Pmod%DirLim(6).eq.1)) ) then
			X(3) = sign(1d0,X(3))*0.999*Pdom%zmaxR
		end if
		
		call findBox(X,k,Pdom)
		if ((Box%BoxType(k,1).eq.2).or.(Box%BoxType(k,1).eq.3)) then
			write(*,*) "Error : cell is still in a must-be-empty or dirichlet box after applying boundary condition."
			STOP
		end if
	end if
		
end subroutine



subroutine HandleFibBoundaryCondition(X,IsInDirichlet,Pmod,Pdom,Box)
use definitions
	implicit none
	double precision, dimension(3), intent(inout) :: X
	integer, intent(inout)          :: IsInDirichlet
	type(ParametersMod), intent(in) :: Pmod
	type(ParametersDom), intent(in) :: Pdom
	type(ParametersBox), intent(in) :: Box
	integer                         :: ix, iy, iz, k
	
	call findBox(X,k,Pdom)
	if ( (k.lt.1).or.(k.gt.8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)) ) then
		write(*,*) "Error : motion bring a fiber into the box ",k," which does not exist !" 
		STOP
	end if
	
	! If a fiber exit the computation domain by entering a periodic ghost box : make it come back from the opposite side.
	if (Box%BoxType(k,1).eq.1) then
		X(:) = X(:) - Box%BoxType(k,2:4)
		
		call findBox(X,k,Pdom)		
		if (Box%BoxType(k,1).eq.1) then
			write(*,*) "Error : fiber is still in a periodic box after applying boundary condition."
			STOP
		end if
	end if
	
	! If a fiber exit the computation domain by entering a must-be-empty box : bring it back inside from the same side.
	if (Box%BoxType(k,1).eq.2) then
		ix = aint( (X(1)+Pdom%xmax)/Pdom%dx )
		iy = aint( (X(2)+Pdom%ymax)/Pdom%dy )
		iz = aint( (X(3)+Pdom%zmax)/Pdom%dz )
		
		if ( (ix.eq.0).and.(Pmod%DirLim(1).eq.0).and.(Pmod%DirLim(2).eq.1) ) then
			X(1) = - 0.999*Pdom%xmaxR
		else if ( (ix.eq.2*Pdom%Nx+1).and.(Pmod%DirLim(1).eq.1).and.(Pmod%DirLim(2).eq.0) ) then
			X(1) = + 0.999*Pdom%xmaxR
		end if
		
		if ( (iy.eq.0).and.(Pmod%DirLim(3).eq.0).and.(Pmod%DirLim(4).eq.1) ) then
			X(2) = - 0.999*Pdom%ymaxR
		else if ( (iy.eq.2*Pdom%Ny+1).and.(Pmod%DirLim(3).eq.1).and.(Pmod%DirLim(4).eq.0) ) then
			X(2) = + 0.999*Pdom%ymaxR
		end if
		
		if ( (iz.eq.0).and.(Pmod%DirLim(5).eq.0).and.(Pmod%DirLim(6).eq.1) ) then
			X(3) = - 0.999*Pdom%zmaxR
		else if ( (iz.eq.2*Pdom%Nz+1).and.(Pmod%DirLim(5).eq.1).and.(Pmod%DirLim(6).eq.0) ) then
			X(3) = + 0.999*Pdom%zmaxR
		end if
		
		call findBox(X,k,Pdom)
		if (Box%BoxType(k,1).eq.2) then
			write(*,*) "Error : fiber is still in a must-be-empty box after applying boundary condition."
			STOP
		end if
	end if
	
	! If a fiber exit the computation domain by entering a dirichlet box : fix it here.
	if (Box%BoxType(k,1).eq.3) then
		if (IsInDirichlet.eq.0) then
			IsInDirichlet = 1
		else if (IsInDirichlet.eq.1) then
			IsInDirichlet = 2
		end if
	end if

end subroutine



subroutine RandomInitialization(Xcell,Rcell,Xfib,omega,link,distpoints,Pmod,Pdom,BoxCell,BoxFib)
use definitions
	implicit none
	double precision, dimension(:,:), intent(inout) :: Xcell, Xfib, omega, link, distpoints
	double precision, dimension(:), intent(inout)   :: Rcell
	type(ParametersMod), intent(in)                 :: Pmod
	type(ParametersDom), intent(in)                 :: Pdom
	type(ParametersBox), intent(inout)              :: BoxCell,BoxFib

	double precision, parameter :: pi=3.141592653589793238
	double precision :: alea1,alea2,alea3,chil,dik,lik,lki,Xtemp(3),ratio
	integer          :: i, k, l, l1, l2, l3, ltemp
	integer          :: ix, iy, iz, Ntemp, Ndirperbox

	! Cells
	ratio = 4.0
	do i=1,Pmod%Ncell
		! Positions homogeneously distributed inside the whole domain
		call random_number(alea1)
		call random_number(alea2)
		call random_number(alea3)
		Xcell(i,1) = (2d0*Pdom%xmaxR*alea1 - Pdom%xmaxR)/ratio
		Xcell(i,2) = (2d0*Pdom%ymaxR*alea2 - Pdom%ymaxR)/ratio
		Xcell(i,3) = (2d0*Pdom%zmaxR*alea3 - Pdom%zmaxR)/ratio
	
		! Security in case of a change in the initialization process
		if (abs(Xcell(i,1))>Pdom%xmaxR) Xcell(i,1) = Xcell(i,1) - sign(1d0,Xfib(i,1))*2d0*Pdom%xmaxR
		if (abs(Xcell(i,2))>Pdom%ymaxR) Xcell(i,2) = Xcell(i,2) - sign(1d0,Xfib(i,2))*2d0*Pdom%ymaxR
		if (abs(Xcell(i,3))>Pdom%zmaxR) Xcell(i,3) = Xcell(i,3) - sign(1d0,Xfib(i,3))*2d0*Pdom%zmaxR

		! Each and every cell begin with the minimal radius Rmin
		Rcell(i) = Pmod%Rmin	
	end do
	
	BoxCell%FirstAgent(:) = 0
	BoxCell%LastAgent(:)  = 0
	BoxCell%AgentBox(:)   = 0
	BoxCell%VerletList(:) = 0
	
	do i = 1,Pmod%Ncell
		call UpdateVerlet(i,Xcell(i,:),Pdom,BoxCell)
	end do
	call FillGhostBoxes(Pmod,Pdom,BoxCell)
	
	
	
	! System's fibers
	do i=1,Pmod%Nfib
		! Positions homogeneously distributed inside the whole domain
		call random_number(alea1)
		call random_number(alea2)
		call random_number(alea3)
		Xfib(i,1) = 2d0*Pdom%xmaxR*alea1 - Pdom%xmaxR
		Xfib(i,2) = 2d0*Pdom%ymaxR*alea2 - Pdom%ymaxR
		Xfib(i,3) = 2d0*Pdom%zmaxR*alea3 - Pdom%zmaxR
	
		! Security in case of a change in the initialization process
		if (abs(Xfib(i,1))>Pdom%xmaxR) Xfib(i,1) = Xfib(i,1) - sign(1d0,Xfib(i,1))*2d0*Pdom%xmaxR
		if (abs(Xfib(i,2))>Pdom%ymaxR) Xfib(i,2) = Xfib(i,2) - sign(1d0,Xfib(i,2))*2d0*Pdom%ymaxR
		if (abs(Xfib(i,3))>Pdom%zmaxR) Xfib(i,3) = Xfib(i,3) - sign(1d0,Xfib(i,3))*2d0*Pdom%zmaxR

		! Orientations homogeneously distributed in the sphere
		call pick_and_reject(omega(i,:))	
	end do
	
	
	! Dirichlet layer's fibers
	Ndirperbox = int(Pmod%DirFillingRate * Pdom%dx*Pdom%dy*Pdom%dz)
	Ntemp = Pmod%Nfib
	
	do ix = 1,2*Pdom%Nx+2
		do iy = 1,2*Pdom%Ny+2
			do iz = 1,2*Pdom%Nz+2
				l = ix + 2*(iy-1)*(Pdom%Nx+1) + 4*(iz-1)*(Pdom%Nx+1)*(Pdom%Ny+1)

				if (BoxFib%BoxType(l,1).eq.3) then
					k = 1
					do while(k.le.Ndirperbox)
						!  Positions homogeneously distributed inside box l = (ix,iy,iz)
						call random_number(alea1)
						call random_number(alea2)
						call random_number(alea3)
						Xfib(Ntemp+k,1) = alea1*Pdom%dx + (ix-1)*Pdom%dx - Pdom%xmax
						Xfib(Ntemp+k,2) = alea2*Pdom%dy + (iy-1)*Pdom%dy - Pdom%ymax
						Xfib(Ntemp+k,3) = alea3*Pdom%dz + (iz-1)*Pdom%dz - Pdom%zmax
        				
						! Orientation parallele to the border between box l and the inside of the domain
						call random_number(alea1)
						if ( (ix.eq.1).or.(ix.eq.2*Pdom%Nx+2) ) then
							omega(Ntemp+k,1) = 0.0
							omega(Ntemp+k,2) = cos(2*pi*alea1)
							omega(Ntemp+k,3) = sin(2*pi*alea1)
						elseif ( (iy.eq.1).or.(iy.eq.2*Pdom%Ny+2) ) then
							omega(Ntemp+k,1) = cos(2*pi*alea1)
							omega(Ntemp+k,2) = 0.0
							omega(Ntemp+k,3) = sin(2*pi*alea1)
						elseif ( (iz.eq.1).or.(iz.eq.2*Pdom%Nz+2) ) then
							omega(Ntemp+k,1) = cos(2*pi*alea1)
							omega(Ntemp+k,2) = sin(2*pi*alea1)
							omega(Ntemp+k,3) = 0.0
						end if

						k = k+1
					end do
					Ntemp = Ntemp + Ndirperbox
				end if
				
			end do
		end do
	end do	
	
	
	! Fill ALL fiber-boxes (including ghost and Dirichlet boxes)
	BoxFib%FirstAgent(:) = 0
	BoxFib%LastAgent(:)  = 0
	BoxFib%AgentBox(:)   = 0
	BoxFib%VerletList(:) = 0
	
	do i = 1,Pmod%Nfib+Pmod%Ndirtot
		call UpdateVerlet(i,Xfib(i,:),Pdom,BoxFib)
	end do
	call FillGhostBoxes(Pmod,Pdom,BoxFib)
	

	! Links
	link = 0
	distpoints = 0
	
	if (Pmod%LinkType.ne.2) then
	
		chil = Pmod%freq_link*1.0/(Pmod%freq_link + Pmod%freq_dlink)
		
		do i=1,Pmod%Nfib
			! Pour rechercher les voisins on utilise une méthode de localisation sur la grille numérique (appelée VerletMethod)
			l = BoxFib%AgentBox(i)
			do l1 = -1,1
				do l2 = -1,1
					do l3 = -1,1
						ltemp = l + l1 + 2*l2*(Pdom%Nx+1) + 4*l3*(Pdom%Nx+1)*(Pdom%Ny+1)
						k = BoxFib%FirstAgent(ltemp)
						
						do while (k /= 0)  ! i is real | k is ghost/real/dirichlet
							if ((k<=Pmod%Nfib).and.(k>i)) then ! k is a moving fiber (ghost or real) different from i
								
								if (BoxFib%BoxType(ltemp,1).eq.1) then
									Xtemp(:) = Xfib(k,:) + BoxFib%BoxType(ltemp,2:4)
								else
									Xtemp(:) = Xfib(k,:)
								end if
								
								call find_overlap_fibfib(Xfib(i,:),Xtemp,omega(i,:),omega(k,:),Pmod,dik,lik,lki)
								call random_number(alea1)
								
								if ( ((dik-2*Pmod%Rf-Pmod%eps)<0d0).and.(alea1<chil) ) then
									link(i,k) = 1
									link(k,i) = 1
									distpoints(i,k) = lik
									distpoints(k,i) = lki
								end if
								
							elseif (k>Pmod%Nfib) then ! k is a fiber from the Dirichlet layer
								
								call find_overlap_fibfib(Xfib(i,:),Xfib(k,:),omega(i,:),omega(k,:),Pmod,dik,lik,lki)
								call random_number(alea1)
								
								if ( ((dik-2*Pmod%Rf-Pmod%eps)<0d0).and.(alea1<chil) ) then
									link(i,k) = 1
									distpoints(i,k) = lik
									distpoints(k,i) = lki
								end if
								
							end if
							k = BoxFib%VerletList(k)
						end do
						
					end do
				end do
			end do
		end do
		
	end if

end subroutine RandomInitialization



end module fonctions