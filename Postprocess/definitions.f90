module definitions

type ParametersSim
	integer            :: InitType
	double precision   :: TF, dt, Tenr
	character(len=100) :: nomdossier
end type ParametersSim

type ParametersMod
	integer           :: Nfib, Ncell, Nmax, DirLim(6), Ndirtot, cumul_replink, LinkType
	double precision  :: Rmin, Rmax, Lf, Rf, eps, d0 
	double precision  :: alpha_repCC, alpha_repCF, alpha_repFF, alpha_rappel, alpha_align
	double precision  :: freq_ens, Kcroiss, freq_link, freq_dlink
	double precision  :: DirFillingRate
end type ParametersMod

type ParametersDom
	integer           :: Nx,Ny,Nz ! nombre de carres en x, en y et en z
	double precision  :: dx,dy,dz ! ca va de soi
	double precision  :: xmax,ymax,zmax,xmaxR,ymaxR,zmaxR ! demi longueur du domaine dans les directions x, y et z (avec et sans les boîtes fantômes à l'exterieur)
end type ParametersDom


type ParametersBox
	integer, dimension(:), allocatable :: AgentBox
	integer, dimension(:), allocatable :: VerletList
	integer, dimension(:), allocatable :: FirstAgent
	integer, dimension(:), allocatable :: LastAgent
	double precision, dimension(:,:), allocatable :: BoxType
end type ParametersBox


type SaveState
	double precision, dimension(:,:,:), allocatable :: CellState, FibState, LinkState
	double precision, dimension(:,:), allocatable   :: VariousReal ! time, tcpu
	integer, dimension(:,:), allocatable            :: VariousInt ! Nlinks_int, Nlinks_ext, Ncell, kbetween
end type SaveState

contains


end module definitions