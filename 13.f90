!in this version i've set Swi = 0 , Ect = 0 , Not3 = 0 , da_di_ds(4), in fact all of da_di_ds
!module esclec is to save and load parameter files and save and output tooth morphology files
!module coreop2d is the model itself. It includes the following subroutines
!      subroutine initialconditions  			specifies the default initial conditions
!      subroutine allocateinitialstate             	allocates the matrices for the initial conditions and sets them     
!      subroutine initializecellpositions               specifies the position of cells in the initial conditions
!      subroutine calculatemargins   			calculates the shape of each epithelial cell (this is the position of its margins)
!      subroutine applydiffusion    		        calculates diffusion of all the molecules between cells
!      subroutine applydifferentiation     		updates cells differentiation values  
!      subroutine EpGrowthBorderForce 			calculates pushing between cells resulting from epithelial growth and border growth
!      subroutine BoyForce            			calculates pushing between cells resulting from buoyancy
!      subroutine repulseneighbor       		calculates repulsion between neighboring cells
!      subroutine repelnonneigh   			checks if non-neigh cells get too close and applies repulsion between them (it never happens for seals)
!      subroutine applyborderbias     			applies BMP4 concentration in the buccal and lingual borders of the tooth
!      subroutine applynucleartraction  		calculates the nucleus traction by the cell borders
!      subroutine updatecellposition    		updates cell positions
!      subroutine addcell            			calculates where cell divisions occur and adds new cells accordingly
!      subroutine updatebordercells     		identifies which of the cells added in addcell are in the border of the tooth
!      subroutine iteration          			determines the other by which the subroutines are called (invariable)
!      subroutine setparams	     			loads parameter file

!***************************************************************************
!***************  MODULE ***************************************************
!***************************************************************************

module coreop2d   !WARNING: This module has been optimized; do not add extra things to it

implicit none
public :: iteration,calculatemargins,applydiffusion,allocateinitialstate,initialconditions

real*8, public, allocatable  :: positions(:,:)  ! the positions of the nodes x, y, z
real*8, public, allocatable  :: border(:,:,:)  ! internodes positions for each cell x,y,z; position of border btw cell i and neighbor j;  the distance between center of cell i and border shared with each neighboring cell
integer, public, allocatable :: neigh(:,:)
integer, public, allocatable :: knots(:)
integer, public, allocatable :: num_neigh(:)      ! How many neigh cell i has (just a count)
real*8, public, allocatable  :: q3d(:,:,:)    ! q3d(cell i,z=depth into tissue,species): 1=[Act]=[Activator], 2=[Inh]=[Inhibitor], 3=[Sec]=[FGF]. quantities that are 3d (activator,inh,fgf) ! species 4 and 5 unused (4 was ectodin, weakly used)
real*8, public, allocatable  :: forces(:,:),forcesnapshot(:,:) !forces = sum of all mechanical forces acting on each cell, forcesnapshot = preserved copy of this specific force contribution, 
real*8, public, allocatable  :: px(:),py(:),pz(:)
real*8, public, allocatable  :: DiffState(:) ! DiffState is a 1d array of the differentiation state for each cell (a list of diff. states).. Diffstate per-cell state variables that apply to whole cell and dont vary with depth z. Previously q2d(cell,1)

integer, public :: num_active_cells   ! number of active cells
integer, public :: num_all_cells   !  number of total cells (active and ghost) within the real radius
integer, public :: max_z_layers     ! number of depth z for the calculation of quantities
integer, parameter, public :: nvmax=30   !maximum number of neighbors
integer, public :: Rad
integer, public, parameter :: num_species_in_q3d=3	!this 3 is the indices used in q3d for concentration of Act, Inh, and FGF
integer, public :: temps,npas
real*8, public,  parameter :: la=1.      !original distance between node

! true parameters
real*8, public :: Egr                  ! Egr (epithelial proliferation rate, "tacre")
real*8, public :: Mgr                  ! Mgr (mesenchymal proliferation rate, "tahor")
real*8, public :: Rep                  ! Rep (Young's modulus, stiffnes, "elas")

real*8, public :: Adh                  ! Adh (traction between neigh, "Adh")
real*8, public :: Act                  ! Act (activator auto-activation, "acac")
real*8, public :: Inh                  ! Inh (inhibition of activator, "iMac")

real*8, public :: Sec                  ! Sec (growth factor secretion rate, "ih")

real*8, public :: Da  		       ! Activator (BMP4) diffusion rate
real*8, public :: Di                   ! Inhibitor (SHH) diffusion rate
real*8, public :: Ds                   ! FGF (Sec) diffusion rate

real*8, public :: Int                  ! Int (differentiation threshold, "us")
real*8, public :: Set                  ! Set (differentiation threshold, "ud")

real*8, public :: Boy                  ! mesenchymal buoyancy (Eq. 13). In subroutine BoyForce. Paper Eq. 13. Mech. resistance of mesench. to invagination of epithelium . 
real*8, public :: Dff                  ! differentiation rate (Eq. 6). Formerly tadif). Equation 6. ! Originally difq2d(3). 
real*8, public :: Bgr                  ! border growth to invagination coupling (Used In Subroutine "updatecellposition". Originally difq2d(4).)
real*8, public, parameter :: dmax=2    ! maximum distance before making a new node
real*8, public :: Abi,Pbi,Bbi,Lbi      ! bias posterior, anterior, lingual, buccal pabl (NOTE: Abi="bip", Pbi="bia", Bbi="bib", Lbi="bil")
real*8, public :: Deg                  ! Deg (protein degradation rate, "mu")
real*8, public :: Dgr                  ! Dgr (sharpness maxima, indicates how the epithelium pulls down if we do not have pressure of the mesenchyme, "tazmax")
                                       ! it doesn't affect much unless it's very low, 100 is a good value
real*8, public :: Ntr                  ! Ntr (mechanical traction from the borders to the nucleus. "radibi")
real*8, public :: Bwi                  ! radius of the center where we apply the bias ap
real*8, public :: ina                  ! initial activator concentration
real*8, public :: umgr                 ! basal mesenchymal prolif. rate (independent of Sec level)

!semi-constants
real*8, public :: umelas
integer, public :: maxcels

!implementation detail 
real*8, parameter ::  delta=0.005D1 , vmin=0.015D1 
integer,public :: nca,icentre,centre,first_border_cell,focus
real*8,public:: x,y,xx,yy
real*8, public :: csu,ssu,csd,ssd,cst,sst,csq,ssq,csc,ssc,css,sss
integer, public :: num_new_cells
integer, public :: num_border_cells,nmap
integer, public,allocatable :: border_cell_list(:),mmap(:)

!typical values
integer, public :: i,j,k,ii,jj,kk,iii,jjj,kkk,iiii,jjjj,kkkk,iiiii,jjjjj,kkkkk
real*8, public :: a,b,c,d,e,f,g,h,aa,bb,cc,dd,ee,ff,gg,hh,aaa,bbb,ccc,ddd,eee,fff,ggg,hhh,panic

!for visualization
integer, public :: vlinies,vrender,showborders,vvec,vvecx,vveck,vex,vn
integer, public :: nc
integer, public :: pin,pina !went to catch Act,inh,fgf,p
integer, public :: nivell,kko !nivell al que fem el tall
integer, public :: iti,iterationtotal
character*31, public ::  nfpro

!constants
real*8, public, parameter ::  pii = 31.41592653589793D-1

contains



subroutine initialconditions
  !core values
  max_z_layers=4
  temps=0
  npas=1
  !for visualization
  vlinies=0
  vrender=1
  vvec=0
  vn=0
  vveck=0
  vvecx=0
  showborders=0
  vex=0
  pin=1
  pina=1
  nivell=1
  panic=0
end subroutine initialconditions



subroutine allocateinitialstate
  integer, allocatable :: cv(:,:)
  real*8 , allocatable :: temp_positions(:,:)

  umelas=1-Rep

  !Variable initial values
  j=0
  do i=1,Rad   ; j=j+i ; end do ; num_all_cells=6*j+1             ! if j is zero then num_all_cells=1 to begin with 
  j=0
  do i=1,Rad-1 ; j=j+i ; end do ; num_active_cells=6*j+1
  a=pii*0.2D1/0.36D3
  csu=dsin(0*a)  ; ssu=dcos(0*a)
  csd=dsin(60*a)  ; ssd=dcos(60*a)
  cst=dsin(120*a) ; sst=dcos(120*a)
  csq=dsin(180*a) ; ssq=dcos(180*a)
  csc=dsin(240*a) ; ssc=dcos(240*a)
  css=dsin(300*a) ; sss=dcos(300*a)

  !allocations
  allocate(cv(num_all_cells,nvmax))
  allocate(temp_positions(num_all_cells,3))              ! temp_positions is a temporary copy of cell positions used during reindexing

  allocate(positions(num_all_cells,3))
  allocate(neigh(num_all_cells,nvmax))
  allocate(forces(num_all_cells,3))
  allocate(forcesnapshot(num_all_cells,3))
  allocate(border(num_all_cells,nvmax,8))
  allocate(knots(num_all_cells))
  allocate(num_neigh(num_all_cells))
  allocate(DiffState(num_all_cells))
  allocate(q3d(num_all_cells,max_z_layers,num_species_in_q3d))
  allocate(mmap(Rad))
  allocate(border_cell_list(Rad))

  !visualization matrices
  allocate(px(num_all_cells)) ; allocate(py(num_all_cells)) ; allocate(pz(num_all_cells))

  !values of zero
  neigh=0. ; num_neigh=0. ; positions=0. ; DiffState = 0. ; q3d=0. ; knots=0 ; forces=0. ; forcesnapshot=0.
  
  !initial values
  positions(1,1)=0. ; positions(1,2)=0. ; positions(1,3)=1.
  nca=1
  num_neigh=6

  !initial mesh values

  al: do icentre=1,num_active_cells
    x=positions(icentre,1) ; y=positions(icentre,2)

    xx=x+csu*la ; yy=y+ssu*la ; j=1 ; jj=4 ; call initializecellpositions
    xx=x+csd*la ; yy=y+ssd*la ; j=2 ; jj=5 ; call initializecellpositions
    xx=x+cst*la ; yy=y+sst*la ; j=3 ; jj=6 ; call initializecellpositions
    xx=x+csq*la ; yy=y+ssq*la ; j=4 ; jj=1 ; call initializecellpositions
    xx=x+csc*la ; yy=y+ssc*la ; j=5 ; jj=2 ; call initializecellpositions
    xx=x+css*la ; yy=y+sss*la ; j=6 ; jj=3 ; call initializecellpositions

  end do al
  do i=2,num_active_cells
    do j=1,nvmax
      if (neigh(i,j)>num_active_cells) then ; neigh(i,j)=num_all_cells ; end if 
    end do
  end do
  do k=1,3
    do i=2,num_active_cells
      do j=1,nvmax-1
        if (neigh(i,j)==num_all_cells.and.neigh(i,j+1)==num_all_cells) then 
          do jj=j,nvmax-1
            neigh(i,jj)=neigh(i,jj+1)
          end do
        end if 
      end do
    end do
  end do
  do i=2,num_active_cells
    k=0
    do j=1,nvmax
      if (neigh(i,j)==num_all_cells.and.k==0) then ; k=1 ; cycle ; end if
      if (neigh(i,j)==num_all_cells.and.k==1) then ; neigh(i,j)=0 ; exit ; end if
    end do
  end do
  positions=dnint(positions*1D14)*1D-14
  do i=1,num_active_cells
    do j=1,3
      if (abs(positions(i,j))<1D-14) positions(i,j)=0. ; end do ; end do

  !original distance calculation between nodes

  !investment so that the former are on the margins

  cv=neigh
  temp_positions=positions
  do i=num_active_cells,1,-1
    neigh(i,:)=cv(num_active_cells-i+1,:)
    positions(i,:)= temp_positions(num_active_cells-i+1,:)
  end do

  cv=neigh
  do i=num_active_cells,1,-1
    ii=num_active_cells-i+1
    do jj=1,num_active_cells
      do jjj=1,nvmax
        if (cv(jj,jjj)==i)  neigh(jj,jjj)=ii
      end do
    end do
  end do

  call calculatemargins
  num_neigh=3
  num_neigh(1)=6
  border(:,:,4:5)=la
  do i=1,num_active_cells
  end do
  centre=num_active_cells
  first_border_cell=(Rad-1)*6+1
  focus=first_border_cell/2-1
  if (Rad==2) focus=3
  border_cell_list=0 ; mmap=0
  !margin values ​​ATTENTION NOT SCALABLE WITH Rad, only for Rad=2
  do i=1,Rad 
    mmap(i)=i 
  end do
  ii=1
  do i=1,Rad 
    border_cell_list(ii)=first_border_cell/2+i 
    ii=ii+1 
  end do
  num_border_cells=Rad ; nmap=Rad
  call calculatemargins
  q3d=0
end subroutine allocateinitialstate



subroutine initializecellpositions
  al: do i=1,nca 
      if (i==icentre) cycle al
      if (dnint(1000000*positions(i,1))==dnint(1000000*xx).and.dnint(1000000*positions(i,2))==dnint(1000000*yy)) then ; 
        do ii=1,nvmax ; if (neigh(icentre,ii)==i) then ; return ; end if ;  end do
        neigh(icentre,j)=i ; neigh(i,jj)=icentre ; num_neigh(i)=num_neigh(i)+1 ; num_neigh(icentre)=num_neigh(icentre)+1 ; return
      end if
    end do al
    num_neigh(icentre)=num_neigh(icentre)+1 ; nca=nca+1 ; neigh(icentre,j)=nca ; neigh(nca,jj)=icentre ; positions(nca,1)=xx ; 
    positions(nca,2)=yy ; positions(nca,3)=1. ; num_neigh(nca)=num_neigh(nca)+1
end subroutine initializecellpositions



subroutine calculatemargins
  real*8 cont
  integer kl

  border(:,:,1:3)=0.
  do i=1,num_active_cells
    aa=0. ; bb=0. ; cc=0. ; kl=0
    do j=1,nvmax
      if (neigh(i,j)/=0.) then 
        a=0. ; b=0. ; c=0. ; cont=0
        iii=i 
        a=positions(i,1) ; b=positions(i,2) ; c=positions(i,3) ; cont=1
        ii=neigh(i,j)  
        if (ii>num_active_cells) then
          do jj=j-1,1,-1
            if (neigh(i,jj)/=0) then  
              if (neigh(i,jj)<num_active_cells+1) then 
                ii=neigh(i,jj)
                a=a+positions(ii,1) 
                b=b+positions(ii,2)
                c=c+positions(ii,3) 
                cont=cont+1
                goto 77
              else
                goto 77
              end if
            end if 
          end do 
          do jj=nvmax,j+1,-1
            if (neigh(i,jj)/=0) then  
              if (neigh(i,jj)<num_active_cells+1) then
                ii=neigh(i,jj)
                a=a+positions(ii,1) 
                b=b+positions(ii,2)
                c=c+positions(ii,3) 
                cont=cont+1
                goto 77
              else
                goto 77
              end if
            end if 
          end do 
          goto 77 
        end if
  66    if (ii==i) goto 77
        kl=kl+1 
        if (kl>100) then 
          do jj=j-1,1,-1 
            if (neigh(i,jj)/=0) then  
              if (neigh(i,jj)<num_active_cells+1) then
                ii=neigh(i,jj) 
                a=a+positions(ii,1) 
                b=b+positions(ii,2) 
                c=c+positions(ii,3) 
                cont=cont+1 
                goto 77 
              else 
                goto 77 
              end if
            end if 
          end do 
          do jj=nvmax,j+1,-1 
            if (neigh(i,jj)/=0) then  
              if (neigh(i,jj)<num_active_cells+1) then 
                ii=neigh(i,jj) 
                a=a+positions(ii,1) 
                b=b+positions(ii,2) 
                c=c+positions(ii,3) 
                cont=cont+1 
                goto 77 
              else 
                goto 77 
              end if
            end if 
          end do 
          goto 77 
        end if
        a=a+positions(ii,1) ; b=b+positions(ii,2) ; c=c+positions(ii,3) ; cont=cont+1
        do jj=1,nvmax
          if (neigh(ii,jj)==iii) then
            jjj=jj
            exit
          end if
        end do
        do jj=jjj+1,nvmax  !comencem la gira
          if (neigh(ii,jj)/=0) then 
            if (neigh(ii,jj)>num_active_cells) goto 77 
            iii=ii 
            ii=neigh(iii,jj) 
            goto 66 
          end if
        end do
        do jj=1,jjj-1  !comencem la gira
          if (neigh(ii,jj)/=0) then 
            if (neigh(ii,jj)>num_active_cells) goto 77 
            iii=ii 
            ii=neigh(iii,jj) 
            goto 66 
          end if
        end do
      end if
  77    border(i,j,1)=a/cont ; border(i,j,2)=b/cont ; border(i,j,3)=c/cont 
    end do
  end do
end subroutine calculatemargins



subroutine applydiffusion
  real*8 pes(num_all_cells,nvmax)   ! area of contact between i and neighbor (i, j)
  real*8 areap(num_all_cells,nvmax)
  real*8 suma,areabottom              ! fraction of a cell’s surface area that is in contact with the “bottom” (the z-direction neighbors / substrate) and is used to weight vertical diffusion.areabottom
  real*8 hq3d(num_all_cells,max_z_layers,num_species_in_q3d) ! num_all_cells = number of cells (within "real radius"), max_z_layers = z depth, num_species_in_q3d = 3 (constant)
  real*8 ux,uy,uz,dx,dy,dz,ua,ub,uc

  hq3d=0.

  do i=1,num_active_cells
    pes(i,:)=0. ; areap(i,:)=0.
  ui: do j=1,nvmax
      if (neigh(i,j)/=0.) then 
        ua=positions(i,1) ; ub=positions(i,2) ; uc=positions(i,3)
        do jj=j+1,nvmax
          if (neigh(i,jj)/=0.) then
  pes(i,j)=sqrt((border(i,j,1)-border(i,jj,1))**2+(border(i,j,2)-border(i,jj,2))**2+(border(i,j,3)-border(i,jj,3))**2)
            ux=border(i,j,1)-ua  ; uy=border(i,j,2)-ub  ;  uz=border(i,j,3)-uc 
            dx=border(i,jj,1)-ua ; dy=border(i,jj,2)-ub ; dz=border(i,jj,3)-uc
            areap(i,j)=0.05D1*sqrt((uy*dz-uz*dy)**2+(uz*dx-ux*dz)**2+(ux*dy-uy*dx)**2)
            cycle ui
          end if
        end do
  pes(i,j)=sqrt((border(i,j,1)-border(i,1,1))**2+(border(i,j,2)-border(i,1,2))**2+(border(i,j,3)-border(i,1,3))**2) 
        ux=border(i,j,1)-ua ; uy=border(i,j,2)-ub ; uz=border(i,j,3)-uc 
        dx=border(i,1,1)-ua ; dy=border(i,1,2)-ub ; dz=border(i,1,3)-uc
        areap(i,j)=0.05D1*sqrt((uy*dz-uz*dy)**2+(uz*dx-ux*dz)**2+(ux*dy-uy*dx)**2)
      end if
    end do ui
    areabottom=sum(areap(i,:))
    suma=sum(pes(i,:))+2*areabottom ; areabottom=areabottom/suma ; pes(i,:)=pes(i,:)/suma 
    do k=1,num_species_in_q3d !ATTENTION
      do kk=2,max_z_layers-1
        hq3d(i,kk,k)=hq3d(i,kk,k)+areabottom*(q3d(i,kk-1,k)-q3d(i,kk,k))
        hq3d(i,kk,k)=hq3d(i,kk,k)+areabottom*(q3d(i,kk+1,k)-q3d(i,kk,k))
        do j=1,nvmax
          if (neigh(i,j)/=0) then 
            ii=neigh(i,j)
            if (ii==num_all_cells) then
              hq3d(i,kk,k)=hq3d(i,kk,k)+pes(i,j)*(-q3d(i,kk,k)*0.044D1)      !sink     
            else
              hq3d(i,kk,k)=hq3d(i,kk,k)+pes(i,j)*(q3d(ii,kk,k)-q3d(i,kk,k)) 
            end if
          end if
        end do
      end do
      hq3d(i,max_z_layers,k)=areabottom*(-q3d(i,max_z_layers,k)*0.044D1) 
      hq3d(i,max_z_layers,k)= hq3d(i,max_z_layers,k)+areabottom*(q3d(i,max_z_layers-1,k)-q3d(i,max_z_layers,k))
      do j=1,nvmax
        if (neigh(i,j)/=0) then 
          ii=neigh(i,j)
          if (ii==num_all_cells) then
              hq3d(i,max_z_layers,k)=hq3d(i,max_z_layers,k)+pes(i,j)*(-q3d(i,max_z_layers,k)*0.044D1)   !sink     
          else
            hq3d(i,max_z_layers,k)=hq3d(i,max_z_layers,k)+pes(i,j)*(q3d(ii,max_z_layers,k)-q3d(i,max_z_layers,k)) 
          end if
        end if
      end do
    end do
    pes(i,:)=pes(i,:)*suma ; areabottom=areabottom*suma ; suma=suma-areabottom ; pes(i,:)=pes(i,:)/suma 
    areabottom=areabottom/suma
    do k=1,num_species_in_q3d !ATTENTION
      hq3d(i,1,k)=areabottom*(q3d(i,2,k)-q3d(i,1,k))
      do j=1,nvmax
        if (neigh(i,j)/=0) then 
          ii=neigh(i,j)
          if (ii==num_all_cells) then
              hq3d(i,1,k)=hq3d(i,1,k)+pes(i,j)*(-q3d(i,1,k)*0.044D1)     
          else
            hq3d(i,1,k)=hq3d(i,1,k)+pes(i,j)*(q3d(ii,1,k)-q3d(i,1,k)) 
          end if
        end if
      end do
    end do
  end do
  ! Activator diffusion
  q3d(:,:,1) = q3d(:,:,1) + delta * Da * hq3d(:,:,1)

  ! Inhibitor diffusion
  q3d(:,:,2) = q3d(:,:,2) + delta * Di * hq3d(:,:,2)

  ! FGF diffusion
  q3d(:,:,3) = q3d(:,:,3) + delta * Ds * hq3d(:,:,3)

  ! REACTION
  hq3d=0.
  do i=1,num_active_cells
    if (q3d(i,1,1)>1) then
      if (i>=first_border_cell) knots(i)=1 
    end if
    
    a = Act*q3d(i,1,1)  ! Act = Activator (activator auto-activation)
    if (a<0) a=0.
    hq3d(i,1,1) = a/(1+Inh*q3d(i,1,2)) - Deg*q3d(i,1,1) ! Eq. (14) sans diffusion
    if (DiffState(i)>Int) then     ! Int (initial inhibitor threshold)
      hq3d(i,1,2) = q3d(i,1,1)*DiffState(i) - Deg*q3d(i,1,2) ! Eq. (17). NOTE: DiffState <= 1.0
    else
      if (knots(i)==1) then
        hq3d(i,1,2) = q3d(i,1,1) - Deg*q3d(i,1,2)        ! Eq. (17)
      end if
    end if
    if (DiffState(i)>Set) then     ! Set (growth factor threshold)
      a = Sec*DiffState(i) - Deg*q3d(i,1,3)                   ! Eq. (18). NOTE: DiffState <= 1.0
      if(a<0.) a=0.
      hq3d(i,1,3) = a
    else
      if (knots(i)>Set) then
        a = Sec - Deg*q3d(i,1,3) ! Eq. (18). Sec = Sec (growth factor secretion rate)
        if(a<0.) a=0.
        hq3d(i,1,3) = a
      end if
    end if

  end do

  if (maxval(abs(hq3d(:,1,1:2)))>1D100) then ; 
    !print *,"PANIC OVERFLOW" ; 
    panic=1 
    return 
  end if
  do i=1,3 
    q3d(:,1,i) = q3d(:,1,i) + delta*hq3d(:,1,i)   ! explicit Euler method
  end do

  where(q3d<0.) q3d=0.  

end subroutine applydiffusion



subroutine applydifferentiation
  do i=1,num_active_cells
    DiffState(i) = DiffState(i) + Dff*(q3d(i,1,3))    ! paper Eq. 6
    if (DiffState(i)>1.) DiffState(i)=1.
  end do
end subroutine applydifferentiation



subroutine EpGrowthBorderForce
  real*8 uux,uuy,uuz,ua,ub,uc,uaa,ubb,uuux,uuuy,duux,duuy

  forcesnapshot=0.
  forces=0
  do i=first_border_cell,num_active_cells
    if (knots(i)==1) cycle
    ua=positions(i,1) ; ub=positions(i,2) ; uc=positions(i,3)
    aa=0 ; bb=0 ; cc=0
    do j=1,nvmax
      k=neigh(i,j)
      if (k==0.or.k>num_active_cells) cycle
      b=uc-positions(k,3)
      if (b<-1D-4) then
        uux=ua-positions(k,1)      ; uuy=ub-positions(k,2)      ; uuz=uc-positions(k,3) 
        d=sqrt(uux**2+uuy**2+uuz**2)
        d=1/d
        aa=aa-uux*d ; bb=bb-uuy*d ; cc=cc-uuz*d
      end if
    end do
    d = sqrt(aa**2+bb**2+cc**2)
    if (d>0) then
      d = Egr/d
      a = 1-DiffState(i) ; if (a<0) a=0.  ! DiffState(i) = differentiation state
      d = d*a
      forces(i,1) = aa*d      ! Eq. (5). aa,bb,cc are x,y,z components of u_ij 
      forces(i,2) = bb*d
      forces(i,3) = cc*d
    end if
  end do

  do i=1,first_border_cell-1
    aa=0. ; bb=0. ; a=-0.3 ; b=0. ; c=0. 
    ua=positions(i,1) ; ub=positions(i,2)
    do j=1,nvmax
      k=neigh(i,j)
      if (k<1.or.k>num_active_cells) cycle
      if (k>first_border_cell-1) then 
        uux=ua-positions(k,1)   ; uuy=ub-positions(k,2) 
        d=sqrt(uux**2+uuy**2)
        if (d>0) then
          c=acos(uux/d)
          if (uuy<0) c=2*pii-c !acos(uux/d)
        end if
      else
        uux=ua-positions(k,1)      ; uuy=ub-positions(k,2)
        d=sqrt(uux**2+uuy**2)
        if (d>0) then
        if (a==-0.3) then
          a=acos(uux/d)
          if (uuy<0) a=2*pii-a
          if (d>0) then
            dd=1/d
            uuux=-uuy*dd           ; uuuy=uux*dd
            uaa=acos(uuux)
            if (uuuy<0) uaa=2*pii-uaa
          end if
        else
          b=acos(uux/d)
          if (uuy<0) b=2*pii-b
          if (d>0) then
            dd=1/d
            duux=-uuy*dd           ; duuy=uux*dd
            ubb=acos(duux)
            if (duuy<0) ubb=2*pii-ubb!acos(duux) !ubb
          end if
        end if  
        end if
      end if
    end do

      if (a<b) then ; d=a ; a=b ; b=d ; end if
      if (c<a.and.c>b) then 
        if (uaa<a.and.uaa>b) then ; uuux=-uuux ; uuuy=-uuuy ; end if ! is on the inside side and then we have to invert it
        if (ubb<a.and.ubb>b) then ; duux=-duux ; duuy=-duuy ; end if
      else
        if (uaa>a.or.uaa<b) then ; uuux=-uuux ; uuuy=-uuuy ; end if ! is on the inside side and then we have to invert it
        if (ubb>a.or.ubb<b) then ; duux=-duux ; duuy=-duuy ; end if 
      end if    
      aa=-uuux-duux ; bb=-uuuy-duuy  

      !now let's see if it's outward from the tooth to the shabby
      a=ua+aa ; b=ub+bb
      c=ua-aa ; d=ub-bb
      dd=sqrt(a**2+b**2)
      ddd=sqrt(c**2+d**2)
      if (ddd>dd) then ; aa=-aa ; bb=-bb ; end if

    ! We also have the downward traction due to the adhesion to the mesenchyme.
    d=sqrt(aa**2+bb**2)
    if (d>0) then
      d = (d + Mgr*q3d(i,1,3) + umgr) / d   ! Eq. (12) + basal Mgr
      aa = aa*d				    ! Eq. 10
      bb = bb*d				    ! Eq. 11
    end if
    cc=Dgr
    d=sqrt(aa**2+bb**2+cc**2)
    if (d>0) then
      d=Egr/d
      a = 1-DiffState(i) ; if (a<0) a=0.  ! DiffState(i) = differentiation state
      d=d*a                           !this seems to be wrong and wrong
      forces(i,1)=aa*d
      forces(i,2)=bb*d
      forces(i,3)=cc*d
    end if
  end do

  forcesnapshot=forces

end subroutine EpGrowthBorderForce



subroutine BoyForce
  real*8 ax,ay

  !does it push perpendicularly? This is equation 13. 

  do i=1,num_active_cells
    ax=forces(i,1) ; ay=forces(i,2) 
    d=sqrt(ax**2+ay**2)
    if (d/=0) then
      c=forces(i,3)
      if (d>0) then
        a=sqrt(ax**2+ay**2+c**2) ; a=-c/a
        ax=ax*a ; ay=ay*a
        dd=sqrt(ax**2+ay**2+d**2) ; dd=Boy*q3d(i,1,3)/dd                         !!! I believe this is Eq. 13. Boy corresponds to k_Boy
        if (dd>0) then                                              !! q3d(i,1,3) = [Sec], DiffState(i) = d_i
        a=1-DiffState(i) ; if (a<0) a=0.                                !! this is the differentiation gate (DiffState(i) is DIFF state variable)
        ax=ax*dd*a ; ay=ay*dd*a ; d=d*dd*a
        forces(i,1)=forces(i,1)-ax ; forces(i,2)=forces(i,2)-ay ;forces(i,3)=forces(i,3)-d     
        end if
      end if
    end if
  end do

end subroutine BoyForce



subroutine repulseneighbor
  real*8 ux,uy,uz,ua,ub,uc,d,dr,rd
  real*8 persu(nvmax,3)
  !finite element roll
  do i=1,num_active_cells
    ua=positions(i,1) ; ub=positions(i,2) ; uc=positions(i,3)
    persu=0.
    do j=1,nvmax
      k=neigh(i,j)
      if (k>0.and.k<num_active_cells+1) then
        ux=positions(k,1)-ua ; uy=positions(k,2)-ub ; uz=positions(k,3)-uc 	!eq1: r_ij = p_j - p_i
        if (abs(ux)<1D-15) ux=0.
        if (abs(uy)<1D-15) uy=0.
        if (abs(uz)<1D-15) uz=0.
        dr=sqrt(ux**2+uy**2+uz**2)						!eq1: d_r = ||p_j - p_i||
        rd=border(i,j,5)							!eq1: rd = ||p_ij^0||
        if (dr<1D-8) dr=0.
        if (rd<1D-8) rd=0.
        if (knots(i)==1.and.knots(k)==1) then
          d=dr-rd 								!eq1: d= (||p_j - p_i|| -||p_ij^0||) (||p_j - p_i||)
          dr=d/dr 								!introduce normalization
          persu(j,1)=ux*dr ; persu(j,2)=uy*dr ; persu(j,3)=uz*dr		!persu = instantenaous "force" contribution vector. will be updated & stored as Force
        else
          if (dr<rd) then
            d=dr-rd 
            dr=d/dr 
            persu(j,1)=ux*dr ; persu(j,2)=uy*dr ; persu(j,3)=uz*dr		!normalized eq1: d= (||p_j - p_i|| -||p_ij^0||) (||p_j - p_i||)
          else
            if (i>first_border_cell-1) then					!the cells from "1" to "first_border_cell-1" are interior epithelial cells, so this is saying "if i IS a border cell!"
              persu(j,1)=ux*Adh ; persu(j,2)=uy*Adh ; persu(j,3)=uz*Adh 	!eq2 (without normalization) because ux = positions(k,1) - positions(i,1), uy = positions(k,2) - positions(i,2) , uz = positions(k,3) - positions(i,3)
            end if
          end if
        end if
      end if
    end do

    !fast version without sorting (possible biases for floats)
    c=Rep 
    if (c>1) c=1
    a=0. ; do j=1,nvmax ; a=a+persu(j,1) ; end do ;				! x direction forces
    forces(i,1)=forces(i,1)+a*c							! eq1 = kRep * (||p_j - p_i|| - ||p_o||)(p_j - p_i)
      a=0. ; do j=1,nvmax ; a=a+persu(j,2) ; end do ;				! y direction forces
    forces(i,2)=forces(i,2)+a*c
      a=0. ; do j=1,nvmax ; a=a+persu(j,3) ; end do ;				! z direction forces
    forces(i,3)=forces(i,3)+a*c

  end do

end subroutine repulseneighbor



subroutine repelnonneigh
  real*8 ux,uy,uz,ua,ub,uc,dd,d
  real*8, allocatable :: persu(:,:),cpersu(:,:)   !arbitrary CRITICAL OPTIMIZATION FACTOR
  integer conta,espai,espaia
  espai=20
  allocate(persu(espai,3))
  do i=1,num_active_cells
    ua=positions(i,1) ; ub=positions(i,2) ; uc=positions(i,3)
    persu=0. ; conta=0
  gg: do ii=1,num_active_cells
      if (ii==i) cycle 
      do j=1,nvmax ; if (neigh(i,j)==ii) cycle gg ; end do
      ux=positions(ii,1)-ua 
      if (ux>0.14D1) cycle
      uy=positions(ii,2)-ub 
      if (uy>0.14D1) cycle
      uz=positions(ii,3)-uc
      if (uz>0.14D1) cycle
      if (abs(ux)<1D-15) ux=0.
      if (abs(uy)<1D-15) uy=0.
      if (abs(uz)<1D-15) uz=0.
      d=sqrt(ux**2+uy**2+uz**2)
      if (d<0.14D1) then
        conta=conta+1
        if (conta>espai) then 
          espaia=espai
          espai=espai+20

          if (allocated(cpersu)) deallocate(cpersu)
          allocate(cpersu(espai,3))
          cpersu= 0.0d0
          cpersu(1:espaia,:)=persu
          deallocate(persu)
          allocate(persu(espai,3))
          persu=cpersu
          deallocate(cpersu)
        end if 
        dd=1/(d+10D-1)**8 ; d=dd/d ; d=aint(d*1D8)*1D-8
        persu(conta,1)=-ux*d ; persu(conta,2)=-uy*d  ; persu(conta,3)=-uz*d 
      end if
    end do gg

    !quick unsorted version (possible biases for floats)
    a=0. ; do j=1,espai ; a=a+persu(j,1) ; end do ;
  forces(i,1)=forces(i,1)+a*Rep
    a=0. ; do j=1,espai ; a=a+persu(j,2) ; end do ;
  forces(i,2)=forces(i,2)+a*Rep
    a=0. ; do j=1,espai ; a=a+persu(j,3) ; end do ;
  forces(i,3)=forces(i,3)+a*Rep
  end do

end subroutine repelnonneigh



subroutine applyborderbias
  do i=1,first_border_cell-1
    if (positions(i,2) < 0.0d0) then
        if (q3d(i,1,1)<Lbi) then  ! Won't reset activator concentration to value of bias; to make Ina work in borders.
            q3d(i,1,1)=Lbi
        end if
    else
      if (positions(i,2) > 0.0d0) then
        if (q3d(i,1,1)<Bbi) then  ! Same as above.
            q3d(i,1,1)=Bbi
        end if
      end if
    end if
  end do
end subroutine



subroutine initact  ! Adds initial activator concentration in each cell. Defined by Ina parameter.
  do i=1,num_active_cells
    q3d(i,1,1)=ina
  end do
end subroutine



subroutine applynucleartraction  !creates bias
  real*8 n
  real*8 positions_after_traction(num_all_cells,3)

  positions_after_traction=positions
  do i=first_border_cell,num_active_cells
    if (DiffState(i)==1) cycle
    a=0. ; b=0. ; c=0. ; n=0
    do j=1,nvmax
      k=neigh(i,j)
      if (k/=0.and.k<num_active_cells+1) then
        a=a+positions(k,1) ; b=b+positions(k,2) ; c=c+positions(k,3)
        n=n+1
      end if
    end do
    n=1/n
    a=a*n ; b=b*n ; c=c*n 
    a=a-positions(i,1)
    b=b-positions(i,2)
    c=c-positions(i,3)
    positions_after_traction(i,1)=positions(i,1)+delta*Ntr*a
    positions_after_traction(i,2)=positions(i,2)+delta*Ntr*b
    if (knots(i)==0) then
      a=1-DiffState(i) ; if (a<0) a=0.
      positions_after_traction(i,3)=positions(i,3)+delta*Ntr*c*a
    end if
  end do

  !for the margins
  do i=1,first_border_cell-1
    if (DiffState(i)==1) cycle
    a=0. ; b=0. ; c=0. ; n=0
    do j=1,nvmax
      k=neigh(i,j)
      if (k>0.and.k<first_border_cell.and.k<num_active_cells+1) then
        a=a+positions(k,1) ; b=b+positions(k,2) ; c=c+positions(k,3)
        n=n+1
      end if
    end do
    n=1/n
    a=a*n ; b=b*n ; c=c*n 
    a=a-positions(i,1)
    b=b-positions(i,2)
    c=c-positions(i,3)
    positions_after_traction(i,1)=positions(i,1)+delta*Ntr*a
    positions_after_traction(i,2)=positions(i,2)+delta*Ntr*b
    if (knots(i)==0) then
      a=1-DiffState(i) ; if (a<0) a=0.
      positions_after_traction(i,3)=positions(i,3)+delta*Ntr*c*a
    end if
  end do
  positions=positions_after_traction
end subroutine applynucleartraction



subroutine updatecellposition

  !we determine the extremes
  do i=1,first_border_cell-1
    if (abs(positions(i,2))<Bwi) then
      if (positions(i,1)>0) then ; forces(i,1)=forces(i,1)*Pbi ; 
        forces(i,3)=forces(i,3)*Bgr ; 
      end if
      if (positions(i,1)<0) then ; forces(i,1)=forces(i,1)*Abi ; 
        forces(i,3)=forces(i,3)*Bgr ; 
      end if
    end if    
  end do

  do i=1,num_active_cells
    if (forces(i,3)<0) forces(i,3)=0. !it is due to the pressure of the stelate
  end do

  do i=1,num_active_cells
    if (knots(i)==1) forces(i,3)=0.
  end do
  do i=1,num_active_cells
      positions(i,:)=positions(i,:)+delta*forces(i,:)  ! Eq3 kinda! This is the final position update (in the case where Swi = 0)
  end do
end subroutine updatecellposition



subroutine addcell

  real*8,  allocatable :: temp_positions(:,:)  ! tempposiitons stores positions while cell indices are reordered
  integer, allocatable :: temp_neigh(:,:),temp_new_neigh(:,:)
  integer, allocatable :: temp_num_neigh(:)
  integer, allocatable :: cknots(:)
  real*8, allocatable  :: cDiffState(:),temp_border(:,:,:)    ! temp_border is temporary border
  real*8,  allocatable :: cq3d(:,:,:)    ! quantities that are 3d: Activator,inh,fgf
  integer,  allocatable :: snapshot_neigh(:,:)
  integer,  allocatable :: snapshot_num_neigh(:)
  real*8 ,  allocatable :: postdivisiontemp_positions(:,:)! postdivisiontemp_positions is a temporary snapshot of cell positions used when swapping/reordering existing cells after new cells are added
  real*8 ,  allocatable :: scq3d(:,:,:)
  real*8, allocatable  :: scDiffState(:)
  real*8 ,  allocatable :: temp_borderofexternalcells(:,:,:) !temp_borderofexternalcells is temporary border of external cells
  integer ,  allocatable :: scknots(:)

  real*8 ua,ub,uc,ux,uy,uz

  integer new_cell_pairs(num_active_cells*nvmax,2),new_cell_is_external(num_active_cells*nvmax)
  integer pillats(nvmax),cpillats(nvmax)
  integer new_num_active_cells,prev_num_active_cells
  integer cj,ini,fi,sjj,ji,ij,jji

  num_new_cells=0
  new_cell_pairs=0
  new_cell_is_external=0

  !first we identify and name the new nodes and rescale the mesh matrix and see
  do i=1,num_active_cells
    kkk=0 ; ji=0
    do j=1,nvmax
      if (neigh(i,j)>num_active_cells) then ; ji=1 ; exit ; end if
    end do
    do j=1,nvmax
      k=neigh(i,j)
      ua=positions(i,1) ; ub=positions(i,2) ; uc=positions(i,3)
      if (k/=0.and.k>i.and.k<=num_active_cells) then
        ux=positions(k,1) ; uy=positions(k,2) ; uz=positions(k,3)
        ux=ux-ua ; uy=uy-ub ; uz=uz-uc
        a=sqrt(ux**2+uy**2+uz**2)
        a=dnint(a*1D9)*1D-9
        if (a>dmax) then  ! we add new node
          num_new_cells=num_new_cells+1 ; new_cell_pairs(num_new_cells,1)=i ; new_cell_pairs(num_new_cells,2)=k
          if (i<first_border_cell.and.k<first_border_cell) then ; new_cell_is_external(num_new_cells)=1 ; end if
        end if
      end if    
    end do 
  end do

  if (num_new_cells>0) then
    do i=1,num_all_cells
      do j=1,nvmax
        if (neigh(i,j)==0) then
          do jj=j,nvmax-1
            neigh(i,jj)=neigh(i,jj+1)
          end do
        end if
      end do
    end do

    new_num_active_cells=num_all_cells+num_new_cells

    allocate(temp_positions(new_num_active_cells,3))   ; allocate(temp_neigh(new_num_active_cells,nvmax))
    allocate(temp_num_neigh(new_num_active_cells))    ; allocate(cDiffState(new_num_active_cells))   
    allocate(cq3d(new_num_active_cells,max_z_layers,num_species_in_q3d))
    allocate(temp_new_neigh(num_new_cells,nvmax)) ; allocate(cknots(new_num_active_cells))     
    allocate(temp_border(new_num_active_cells,nvmax,8))

    temp_positions=0. ; temp_neigh=0. ; temp_num_neigh=0. ; cDiffState=0. ; cq3d=0. ; cknots=0 ; temp_border=0.

    do i=1,num_new_cells
      ii=new_cell_pairs(i,1) ; kk=new_cell_pairs(i,2)
    end do

    do i=new_num_active_cells,num_active_cells+num_new_cells+1,-1
      temp_positions(i,:)=positions(i-num_new_cells,:)     ; temp_neigh(i,:)=neigh(i-num_new_cells,:)
      temp_num_neigh(i)=num_neigh(i-num_new_cells)       ; cDiffState(i)=DiffState(i-num_new_cells)
      cq3d(i,:,:)=q3d(i-num_new_cells,:,:)     ; cknots(i)=knots(i-num_new_cells)
      temp_border(i,:,4:8)=border(i-num_new_cells,:,4:8)
    end do

    temp_positions(1:num_active_cells,:)=positions(1:num_active_cells,:)
    temp_neigh(1:num_active_cells,:)=neigh(1:num_active_cells,:) 
    temp_num_neigh(1:num_active_cells)=num_neigh(1:num_active_cells)
    cDiffState(1:num_active_cells)=DiffState(1:num_active_cells)
    cq3d(1:num_active_cells,:,:)=q3d(1:num_active_cells,:,:)     ; cknots(1:num_active_cells)=knots(1:num_active_cells)
    temp_border(1:num_active_cells,:,:)=border(1:num_active_cells,:,:)   

    do i=1,num_active_cells+num_new_cells ; do j=1,nvmax 
    if (temp_neigh(i,j)>num_active_cells) temp_neigh(i,j)=new_num_active_cells 
    end do ; end do

    temp_neigh(num_active_cells+num_new_cells+1:new_num_active_cells,:)=0

    do i=1,num_new_cells
      ii=new_cell_pairs(i,1) ; kk=new_cell_pairs(i,2) ; jj=num_active_cells+i 
      temp_neigh(jj,:)=0 ; temp_neigh(jj,1)=ii ; temp_neigh(jj,2)=kk

      a= temp_positions(ii,1)+ temp_positions(kk,1) ; b= temp_positions(ii,2)+ temp_positions(kk,2)
      d=sqrt(a**2+b**2)
      a=a/d ; b=b/d
      d=d/0.20D1
      d=dnint(d*1D10)*1D-10
      a=d*a ; b=d*b
      temp_positions(jj,1)=a ; temp_positions(jj,2)=b

      temp_positions(jj,3)=(positions(ii,3)+positions(kk,3))*0.05D1
      cq3d(jj,:,:)=(cq3d(ii,:,:)+cq3d(kk,:,:))*0.05D1 ; cDiffState(jj)=(cDiffState(ii)+cDiffState(kk))*0.05D1
      !The 0.333 instead of 0.25 is a trick because before we multiplied the parents by 3/4
    end do
 
    do i=1,num_new_cells
      ii=new_cell_pairs(i,1) ; kk=new_cell_pairs(i,2)
      jj=num_active_cells+i
      do j=1,nvmax ; if (temp_neigh(ii,j)==kk) then ; temp_neigh(ii,j)=jj ; exit ; end if ; end do
      do j=1,nvmax ; if (temp_neigh(kk,j)==ii) then ; temp_neigh(kk,j)=jj ; exit ; end if ; end do
    end do

    do i=1,num_new_cells
      pillats=0
      ii=new_cell_pairs(i,1) ; kk=new_cell_pairs(i,2)
      jj=num_active_cells+i
      !now I have to look at which parent is to follow it (it will be called ini) it is the one that has no external nodes towards j+1
      do j=1,nvmax ; if (neigh(ii,j)==kk) then ; jjj=j ; exit ; end if ; end do
      kkk=0
      do jjjj=jjj+1,nvmax
        if (neigh(ii,jjjj)>0) then
          if (neigh(ii,jjjj)<num_active_cells+1) then 
            ini=ii ; fi=kk ; kkk=1 ; exit
          else 
            ini=kk ; fi=ii ; kkk=1 ; exit 
          end if
        end if       
      end do
      if (kkk==0) then
        do jjjj=1,jjj-1
          if (neigh(ii,jjjj)>0) then
            if (neigh(ii,jjjj)<num_active_cells+1) then 
              ini=ii ; fi=kk ; exit
            else 
              ini=kk ; fi=ii ; exit 
            end if       
          end if
        end do
      end if
      iii=ini
      cj=1 ; pillats(cj)=iii
      !now we look for the j in which iii has jj (the cell whose neigh we are looking for)
      do j=1,nvmax
        if (temp_neigh(iii,j)==jj) then ; jjj=j ; exit ; end if
      end do
      !side j+1(right); we follow the neighbor towards the j+1 side of iii
      
      kkk=0
      do j=jjj+1,nvmax ; jji=temp_neigh(iii,j) 
        if (jji/=0.and.jji<num_active_cells+num_new_cells+1) then ; iiii=jji 
        kkk=1 ; exit  ; end if
      end do
      if (kkk==0) then !I couldn't find the neighbor and I need to go over it again.
        do j=1,jjj-1 ; jji=temp_neigh(iii,j) 
        if (jji/=0.and.jji<num_active_cells+num_new_cells+1) then ; iiii=jji 
        kkk=1 ; exit  ; end if ; end do
      end if

      cj=cj+1 ; if (cj>nvmax) then ; panic=1 ; return ; end if ; pillats(cj)=iiii
      do j=1,nvmax ; if (temp_neigh(iiii,j)==iii) then ; jjjj=j ; exit ; end if ; end do
  88  iii=iiii ; jjj=jjjj
      
      kkk=0
      do j=jjj+1,nvmax ; jji=temp_neigh(iii,j) 
        if (jji/=0.and.jji<num_active_cells+num_new_cells+1) then ; iiii=jji 
        kkk=1 ; exit  ; end if
      end do
      if (kkk==0) then !I couldn't find the neighbor and I need to go over it again.
        do j=1,jjj-1 ; jji=temp_neigh(iii,j) 
          if (jji/=0.and.jji<num_active_cells+num_new_cells+1) then ; iiii=jji 
          kkk=1 ; exit  ; end if
        end do
      end if

      cj=cj+1 ; if (cj>nvmax) then ; panic=1 ; return ; end if ; pillats(cj)=iiii
      do j=1,nvmax ; if (temp_neigh(iiii,j)==iii) then ; jjjj=j ; exit ; end if ; end do
      if (iiii==fi) then ! equinox
        kkk=0
        do kkkk=jjjj+1,nvmax
          if (temp_neigh(iiii,kkkk)/=0.and.kkk==1) then  
            if (temp_neigh(iiii,kkkk)>num_active_cells+num_new_cells) then
              iiii=ini ; kkk=2 ; cj=cj+1 ; exit
            else              
              sjj=kkkk ; kkk=2 ; exit 
            end if
          end if
          if (temp_neigh(iiii,kkkk)/=0.and.kkk==0) then; kkk=1 ; sjj=kkkk ;  end if
        end do
        if (kkk<2) then
          do kkkk=1,jjjj-1
            if (temp_neigh(iiii,kkkk)/=0.and.kkk==1) then
              if (temp_neigh(iiii,kkkk)/=0.and.kkk==1) then  
                if (temp_neigh(iiii,kkkk)>num_active_cells+num_new_cells) then
                  iiii=ini ; cj=cj+1 ; exit
                else              
                  sjj=kkkk ; exit 
                end if
              end if
            end if
            if (temp_neigh(iiii,kkkk)/=0.and.kkk==0) then; kkk=1 ; sjj=kkkk ; end if
          end do
        end if
        jjjj=sjj-1
      end if

      !we have gone all the way around
      if (iiii==ini) then
        cpillats=0
        if (cj>nvmax) then 
          panic=1 ; return 
        end if
        pillats(cj)=0 ; cj=cj-1
        do jjj=1,cj
          cpillats(cj-jjj+1)=pillats(jjj) 
        end do
        pillats=cpillats
        !now let's see what nodes can actually be
        jjj=0
        if (cj>nvmax) then 
          panic=1 ; return 
        end if ;
        do kkk=1,cj
          kkkk=pillats(kkk)
          if (kkkk>num_active_cells.and.kkkk<=num_active_cells+num_new_cells) jjj=jjj+1
        end do
        if (jjj==0) then  !we don't have new nodes on the sides then the crossing is impossible
          temp_new_neigh(i,:)=pillats  
          do j=1,nvmax
            k=temp_new_neigh(i,j)
            if (k/=0) then
              ii=new_cell_pairs(i,1) ; iiii=new_cell_pairs(i,2)
              if (k/=ii.and.k/=iiii.and.k<num_active_cells+num_new_cells+1) then ! It's one of those that I have to connect with.
                kkkk=0
  uu:           do kk=1,nvmax
  uuu:            do kkk=1,nvmax
                    if (temp_neigh(k,kk)/=0.and.temp_neigh(k,kk)==pillats(kkk).and.kkkk==1) then 
                      ji=kk ; exit uu
                    end if
                    if ( temp_neigh(k, kk) /= 0 .and. &
                        temp_neigh(k, kk) == pillats(kkk) .and. &
                        kkkk == 0 ) then
                      kkkk = 1
                      ij   = kk
                      exit uuu
                    end if
                  end do uuu
                end do uu
                if (ji-ij==1) then
                  do kk=nvmax,ji+1,-1 ; temp_neigh(k,kk)=temp_neigh(k,kk-1) 
                  temp_border(k,kk,4:8)=temp_border(k,kk-1,4:8) ; end do 
                  temp_neigh(k,ji)=jj
                else
                  temp_neigh(k,ji+1)=jj ; 
                end if
              end if
            end if
          end do     
        else                       !the txungu will be in sinomes we have a new one because then there will only be 3 neigh
          temp_new_neigh(i,:)=0
          temp_new_neigh(i,1)=ini
          kkkk=1
          if (cj>nvmax) then ; panic=1 ; return ; end if
rtt:      do kkk=1,cj
            jjj=pillats(kkk)
            if (jjjj==fi) then 
              kkkk=kkkk+1 ; temp_new_neigh(i,kkkk)=fi 
            else
              if (jjjj>num_active_cells) then
                kkkk=kkkk+1 ; temp_new_neigh(i,kkkk)=jjjj
              end if
            end if
          end do rtt
          goto 899
!           do j=1,nvmax
!             k=temp_new_neigh(i,j)
!             if (k/=0) then
!               ii=new_cell_pairs(i,1) ; iiii=new_cell_pairs(i,2)
!               if (k/=ii.and.k/=iiii.and.k<num_active_cells+1) then ! It's one of those that I have to connect with.
!                 kkkk=0
! uuuu:           do kk=1,nvmax
! uuuuu:            do kkk=1,nvmax
!                     if ( temp_neigh(k, kk) /= 0 .and. &
!                         temp_neigh(k, kk) == pillats(kkk) .and. &
!                         kkkk == 1 ) then
!                       ji = kk
!                       exit uuuu
!                     end if
!                     if ( temp_neigh(k, kk) /= 0 .and. &
!                         temp_neigh(k, kk) == pillats(kkk) .and. &
!                         kkkk == 0 ) then
!                       kkkk = 1
!                       ij   = kk
!                       exit uuuuu
!                     end if
!                   end do uuuuu
!                 end do uuuu
!                 !we have to connect between ij and ji
!                 if (ji-ij==1) then
!                   do kk=nvmax,ji+1,-1
!                     temp_neigh(k,kk)=temp_neigh(k,kk-1) 
!                     temp_border(k,kk,4:8)=temp_border(k,kk-1,4:8)
!                   end do
!                   temp_neigh(k,ji)=jj
!                 else
!                   temp_neigh(k,ji+1)=jj
!                 end if
!               end if
!             end if
!           end do
  899     continue     
        end if

        !now we need to add the connections to external nodes
        ii=new_cell_pairs(i,1) ; kk=new_cell_pairs(i,2)
        kkk=0 ; jjj=0
        do j=1,nvmax
          if (temp_neigh(ii,j)>num_active_cells+num_new_cells) then ; kkk=1 ; exit ; end if
        end do
        do j=1,nvmax
          if (temp_neigh(kk,j)>num_active_cells+num_new_cells) then ; kkk=kkk+1 ; exit ; end if
        end do
        if (kkk==2) then
          do j=1,nvmax
            if (temp_new_neigh(i,j)==ii) then ; ij=j ; exit ; end if
          end do
          do j=1,nvmax
            if (temp_new_neigh(i,j)==kk) then ; jjj=j ; exit ; end if
          end do 
          if (ij>jjj) then 
            ji=ij ; ij=jjj
          else
            ji=jjj
          end if
          !we have to connect between ij and ji
          if (ji-ij==1) then
            do kk=nvmax,ji+1,-1 ; temp_new_neigh(i,kk)=temp_new_neigh(i,kk-1) 
            temp_border(i,kk,4:8)=temp_border(i,kk-1,4:8) ; end do 
            temp_new_neigh(i,ji)=new_num_active_cells;
          else
            temp_new_neigh(i,ji+1)=new_num_active_cells ; !temp_border(i,ji+1,4:5)=0.
          end if
        end if 
        cycle 
      end if
      goto 88
    end do

    !now we need to add the external connections

    !now we replace
    temp_neigh(num_active_cells+1:num_active_cells+num_new_cells,:)=temp_new_neigh(1:num_new_cells,:)

    deallocate(positions)
    allocate(positions(new_num_active_cells,3))
    positions=temp_positions

    !we calculate the new basal distances of the new cells
    do i=num_active_cells+1,num_active_cells+num_new_cells
      ua=positions(i,1) ; ub=positions(i,2) ; uc=positions(i,3)
      do j=1,nvmax
        ii=temp_neigh(i,j)
        if (ii>0.and.ii<num_active_cells+num_new_cells+1) then
          ux=positions(ii,1) ; uy=positions(ii,2) ; uz=positions(ii,3) 
          ux=ux-ua ; uy=uy-ub ; uz=uz-uc
          if (abs(ux)<10D-14) ux=0.
          if (abs(uy)<10D-14) uy=0.
          if (abs(uz)<10D-14) uz=0.
          d=sqrt(ux**2+uy**2)
          temp_border(i,j,5)=d
          d=sqrt(ux**2+uy**2+uz**2)
          temp_border(i,j,4)=d
          temp_border(i,j,6)=ux ;  temp_border(i,j,7)=uy ;  temp_border(i,j,8)=uz 
        end if
      end do
    end do

    !we calculate the basal distances of the new connections between old and new cells
    do i=1,num_active_cells
      ua=positions(i,1) ; ub=positions(i,2) ; uc=positions(i,3)
      do j=1,nvmax
        ii=temp_neigh(i,j)
        if (ii>num_active_cells.and.ii<num_active_cells+num_new_cells+1) then
          ux=positions(ii,1) ; uy=positions(ii,2) ; uz=positions(ii,3) 
          ux=ux-ua ; uy=uy-ub ; uz=uz-uc
          if (abs(ux)<10D-14) ux=0.
          if (abs(uy)<10D-14) uy=0.
          if (abs(uz)<10D-14) uz=0.
          d=sqrt(ux**2+uy**2)
          temp_border(i,j,5)=d
          d=sqrt(ux**2+uy**2+uz**2)
          temp_border(i,j,4)=d
          temp_border(i,j,6)=ux ;  temp_border(i,j,7)=uy ;  temp_border(i,j,8)=uz 
        end if
      end do  
    end do  

    num_all_cells=new_num_active_cells
    prev_num_active_cells=num_active_cells
    num_active_cells=num_active_cells+num_new_cells

    deallocate(neigh)    ; deallocate(forces)  ; deallocate(forcesnapshot)
    deallocate(border) ;  deallocate(num_neigh) ; deallocate(DiffState)
    deallocate(q3d)   ; deallocate(px)    ;  deallocate(py)     ; deallocate(pz)      
    deallocate(knots)

    allocate(neigh(new_num_active_cells,nvmax)) ; allocate(forces(new_num_active_cells,3))  
    allocate(forcesnapshot(new_num_active_cells,3))
    allocate(border(new_num_active_cells,nvmax,8)) ;  allocate(num_neigh(new_num_active_cells))    
    allocate(DiffState(new_num_active_cells))  
    allocate(q3d(new_num_active_cells,max_z_layers,num_species_in_q3d))
    allocate(px(new_num_active_cells)) ; allocate(py(new_num_active_cells)) 
    allocate(pz(new_num_active_cells)) 
    allocate (knots(new_num_active_cells))

    neigh=temp_neigh ; num_neigh=temp_num_neigh ; DiffState=cDiffState ; q3d=cq3d ; knots=cknots 
    border=temp_border ; forces=0.
    deallocate(temp_positions)  ; deallocate(temp_neigh)    
    deallocate(temp_num_neigh) ; deallocate(cDiffState)  ; deallocate(cq3d)
    deallocate(temp_new_neigh)   ; deallocate(cknots); deallocate(temp_border) 

    !we remove the zeros
    do i=1,num_all_cells
      do j=1,nvmax
        if (neigh(i,j)==0) then
          do jj=j,nvmax-1
            neigh(i,jj)=neigh(i,jj+1)
          end do
        end if
      end do
    end do
    do i=1,num_all_cells
      ii=0
      do j=1,nvmax
        if (neigh(i,j)>0) ii=ii+1
      end do
      num_neigh(i)=ii
    end do
  end if   !END IF WE HAVE NEW CELLS

  do iii=1,num_new_cells
    if (new_cell_is_external(iii)==1) then  !we have a new external cell
      ii=prev_num_active_cells+iii
      if (first_border_cell==centre) centre=ii
      allocate(snapshot_neigh(num_active_cells,nvmax))
      allocate(snapshot_num_neigh(num_active_cells))
      allocate(postdivisiontemp_positions(num_active_cells,3))
      allocate(scq3d(num_active_cells,max_z_layers,num_species_in_q3d))
      allocate(scDiffState(num_active_cells))
      allocate(temp_borderofexternalcells(num_active_cells,nvmax,8))
      allocate(scknots(num_active_cells))
      postdivisiontemp_positions=positions
      snapshot_neigh=neigh
      scDiffState=DiffState
      scq3d=q3d
      snapshot_num_neigh=num_neigh
      temp_borderofexternalcells=border
      scknots=knots
      neigh(ii,:)=snapshot_neigh(first_border_cell,:)
      neigh(first_border_cell,:)=snapshot_neigh(ii,:)
      num_neigh(ii)=snapshot_num_neigh(first_border_cell)
      num_neigh(first_border_cell)=snapshot_num_neigh(ii)
      positions(ii,:)=postdivisiontemp_positions(first_border_cell,:)
      positions(first_border_cell,:)=postdivisiontemp_positions(ii,:)
      DiffState(ii)=scDiffState(first_border_cell)
      DiffState(first_border_cell)=scDiffState(ii)
      q3d(ii,:,:)=scq3d(first_border_cell,:,:)
      q3d(first_border_cell,:,:)=scq3d(ii,:,:)
      border(ii,:,:)=temp_borderofexternalcells(first_border_cell,:,:)
      border(first_border_cell,:,:)=temp_borderofexternalcells(ii,:,:)
      knots(ii)=scknots(first_border_cell)
      knots(first_border_cell)=scknots(ii)
      snapshot_neigh=neigh

      do i=1,num_active_cells
        do j=1,nvmax
          if (snapshot_neigh(i,j)==ii) neigh(i,j)=first_border_cell
        end do
      end do
      do i=1,num_active_cells
        do j=1,nvmax
         if (snapshot_neigh(i,j)==first_border_cell) neigh(i,j)=ii
        end do
      end do
      deallocate(snapshot_neigh)
      deallocate(snapshot_num_neigh)
      deallocate(postdivisiontemp_positions)
      deallocate(scq3d)
      deallocate(scDiffState)
      deallocate(temp_borderofexternalcells)
      deallocate(scknots)
      first_border_cell=first_border_cell+1                              ! wxpand the contiguous border-cell block by 1first_border_cell
    end if
  end do


  if (num_new_cells>0) then
    call updatebordercells
  end if

end subroutine addcell



subroutine updatebordercells
 !makes the new cells that are on the margin between two extremes for the bias extreme
 integer,allocatable :: old_border_cell_list(:)
 integer new_border_cells(first_border_cell)
 integer num_new_border_cells,old_num_border_cells
  !now let's correct the ends
  new_border_cells=0
  num_new_border_cells=0
  er:  do i=1,first_border_cell-1
       if (i==3.or.i==6) cycle !ATTENTION IS NOT SCALABLE WITH RAD
       kk=0
       do ii=1,num_border_cells
         iii=border_cell_list(ii)
         if (iii==i) cycle er
       end do
  err:   do j=1,nvmax
         k=neigh(i,j)
         if (k<first_border_cell) then
           do ii=1,num_border_cells
             iii=border_cell_list(ii)
             if (k==iii) then
               if (kk==1) then
                 num_new_border_cells=num_new_border_cells+1
                 new_border_cells(num_new_border_cells)=i    
                 cycle er
               else
                 kk=1
                 cycle err
               end if
             end if
           end do
         end if
       end do err
     end do er

     if (num_new_border_cells>0) then
       old_num_border_cells=num_border_cells
       allocate(old_border_cell_list(num_border_cells))
       old_border_cell_list=border_cell_list
       deallocate(border_cell_list)
       num_border_cells=num_border_cells+num_new_border_cells
       allocate(border_cell_list(num_border_cells))
       border_cell_list(1:num_border_cells-num_new_border_cells)=old_border_cell_list
       do i=1,num_new_border_cells
         border_cell_list(old_num_border_cells+i)=new_border_cells(i)
       end do 
       deallocate(old_border_cell_list)
     end if

     new_border_cells=0
     num_new_border_cells=0
  era:  do i=1,first_border_cell-1
       if (i==3.or.i==6) cycle !ATTENTION IS NOT SCALABLE WITH RAD
       kk=0
       do ii=1,nmap
         iii=mmap(ii)
         if (iii==i) cycle era
       end do
  erra:   do j=1,nvmax
         k=neigh(i,j)
         if (k<first_border_cell) then
           do ii=1,nmap
             iii=mmap(ii)
             if (k==iii) then
               if (kk==1) then
                 num_new_border_cells=num_new_border_cells+1
                 new_border_cells(num_new_border_cells)=i    
                 cycle era
               else
                 kk=1
                 cycle erra
               end if
             end if
           end do
         end if
       end do erra
     end do era

     if (num_new_border_cells>0) then
       old_num_border_cells=nmap
       allocate(old_border_cell_list(nmap))
       old_border_cell_list=mmap
       deallocate(mmap)
       nmap=nmap+num_new_border_cells
       allocate(mmap(nmap))
       mmap(1:nmap-num_new_border_cells)=old_border_cell_list
       do i=1,num_new_border_cells
         mmap(old_num_border_cells+i)=new_border_cells(i)
       end do 
       deallocate(old_border_cell_list)
     end if
end subroutine



subroutine iteration(tbu)
  integer tbu,ite
  do ite=1,tbu
    panic=0
    forces=0.
    call applydiffusion
    if (panic==1) return
    call applyborderbias
    call applydifferentiation
    call EpGrowthBorderForce
    call BoyForce
    call repelnonneigh
    call repulseneighbor
    call applynucleartraction
    call updatecellposition
    call addcell
    call calculatemargins
    temps=temps+1
  end do
end subroutine iteration



end module coreop2d

!***************************************************************************
!***************  MOQUL ***************************************************
!***************************************************************************

module esclec


use coreop2d
! public:: guardaforma
! character*50, public :: fifr,fivr,fipr,fifw,fivw,fipw 
integer, public :: file_opened_success,map,read_failed,pass,passs,maptotal,is,maxll
integer, parameter :: mamax=5000
real*8,  public :: positionsp(1000,3,mamax)  !attention if num_active_cells>1000 the system crashes when reading
real*8,  public :: parap(32,mamax)
character*50, public :: param_names(32)
integer, public :: knotsp(1000,mamax)
integer, public, allocatable :: neighp(:,:,:)
real*8, public,allocatable :: ma(:)
real*8, public :: vamax,vamin
character*30, public :: cac,cad,caq,cat,cas,cau,cass
contains



subroutine read_param_file
  character*50 param_name
  do i=3,32
    read (2,*,END=666,ERR=777)  a,param_name ;
    parap(i,map)=a ;
    param_names(i)=param_name 
  end do
  Return

  777 print *,"reading error for" ; read_failed=1 ; close(2) ; return
  666 read_failed=1 ; close(2) ; return
end subroutine read_param_file



subroutine setparams(imap)                 ! imap: index of parameter set (always 1 in current usage)
  integer imap                               ! Declare imap as a local integer variable, which is the subroutine argument passed in by the caller.
  temps=parap(1,imap)                      ! global counter of the number of simulation time steps that have elapsed, measured in number of iterations that have already been executed
  num_active_cells=parap(2,imap)           ! Current number of cells
  Egr=parap(3,imap) 
  Mgr=parap(4,imap) 
  Rep=parap(5,imap) 
  !Nothing
  Adh=parap(7,imap)
  Act=parap(8,imap)  
  Inh=parap(9,imap) 
  !Nothing
  Sec=parap(11,imap) 
  Da=parap(13,imap)    
  Di=parap(14,imap) 
  Ds=parap(15,imap) 
  !Nothing
  Int=parap(17,imap) 
  Set=parap(18,imap)
  !Nothing
  Boy=parap(19,imap)                      ! Boy (Eq. 13)
  Dff=parap(20,imap)                      ! Dff (differentiation rate). (Eq. 6) (formerly tadif)
  Bgr=parap(21,imap)                      ! Bgr (border growth factor) (softeq. "Growth Biases", subroutine updatecellposition). (formerly fac)
  Abi=parap(22,imap) 
  Pbi=parap(23,imap) 
  Bbi=parap(24,imap) 
  Lbi=parap(25,imap)
  Rad=parap(26,imap)
  Deg=parap(27,imap)  
  Dgr=parap(28,imap)
  Ntr=parap(29,imap) 
  Bwi=parap(30,imap)
  ina = parap(31,imap)
  umgr = parap(32,imap)
end subroutine setparams



subroutine initialize_from_parameter_file
  open(2,file=cac,status='old',iostat=i)
  map=1
  call read_param_file
  call setparams(map)
  close(2)
end subroutine initialize_from_parameter_file



subroutine guardaveinsoff_2(temp_neigh)
  integer temp_neigh(num_all_cells,nvmax)
  real*8 c(4),mic(4)
    !I have to go from neighboring to faces

  integer face(num_active_cells*20,5)
  integer nfaces
  integer bi,nop
  real*8 mamx

  allocate(ma(num_active_cells))
    call mat

  nfaces=0

    do i=1,num_active_cells
  ale: do j=1,nvmax
       bi=0
       ii=neigh(i,j) ; if (ii==0.or.ii>num_active_cells) cycle
  ele:   do k=1,nvmax
         iii=neigh(ii,k) ; if (iii==0.or.iii>num_active_cells.or.iii==i) cycle
         do kk=1,nvmax
           iiii=neigh(iii,kk) ; if (iiii==0.or.iiii>num_active_cells) cycle
           if (iiii==i) then !triangle trobat
             nfaces=nfaces+1          
             bi=bi+1
             nop=iii     
             if (bi==1) cycle ele
             cycle ale
           end if
         end do
       end do ele
      
       do k=1,nvmax
         iii=neigh(ii,k) ; if (iii==0.or.iii>num_active_cells.or.iii==i.or.iii==nop) cycle
         if (bi==0) cycle
         ! to by the squares
         do kk=1,nvmax
           iiii=neigh(iii,kk) ; 
           if (iiii==0.or.iiii>num_active_cells.or.iiii==ii.or.iiii==nop) cycle
           do kkk=1,nvmax
             jj=neigh(iiii,kkk)
             if (jj==i) then !triangle found
               nfaces=nfaces+1
               !write (2,*) nfaces,i,ii,iii,iiii
               cycle ale
             end if
           end do
         end do
       end do

    end do ale
  end do
  
  write (2,*) "COFF"  
  write (2,*) num_active_cells,nfaces,0
  do i=1,num_active_cells ; call get_rainbow_knot(ma(i),vamin,vamax,c,i) ; 
  write (2,*) positions(i,:),c ; end do

  nfaces=0

    do i=1,num_active_cells
  aale: do j=1,nvmax
       bi=0
       ii=neigh(i,j) ; if (ii==0.or.ii>num_active_cells) cycle
  aele:   do k=1,nvmax
         iii=neigh(ii,k) ; if (iii==0.or.iii>num_active_cells.or.iii==i) cycle
         do kk=1,nvmax
           iiii=neigh(iii,kk) ; if (iiii==0.or.iiii>num_active_cells) cycle
           if (iiii==i) then !triangle found
               nfaces=nfaces+1
               call get_rainbow(ma(i),vamin,vamax,c)
               mic=c ; mamx=ma(i)
               call get_rainbow(ma(ii),vamin,vamax,c)
               mic=mic+c
               call get_rainbow(ma(iii),vamin,vamax,c)
               mic=mic+c
               mic=mic/3.
             write (2,67) 3,i-1,ii-1,iii-1 
  67 format (4I10,4F10.6)
             bi=bi+1
             nop=iii          
             if (bi==1) cycle aele
             cycle aale
           end if
         end do
       end do aele
      
       do k=1,nvmax
         iii=neigh(ii,k) ; if (iii==0.or.iii>num_active_cells.or.iii==i.or.iii==nop) cycle
         if (bi==0) cycle
         ! to by the squares
         do kk=1,nvmax
           iiii=neigh(iii,kk) ; 
           if (iiii==0.or.iiii>num_active_cells.or.iiii==ii.or.iiii==nop) cycle
           do kkk=1,nvmax
             jj=neigh(iiii,kkk)
             if (jj==i) then !triangle found
               nfaces=nfaces+1
               call get_rainbow(ma(i),vamin,vamax,c)
               mic=c ; mamx=ma(i)
               call get_rainbow(ma(ii),vamin,vamax,c)
               mic=mic+c 
               call get_rainbow(ma(iii),vamin,vamax,c)
               mic=mic+c ; 
               call get_rainbow(ma(iiii),vamin,vamax,c)
               mic=mic+c ; 
               mic=c/4.
               write (2,68) 4,i-1,ii-1,iii-1,iiii-1
  68 format (5I10,4F10.6)
               cycle aale
             end if
           end do
         end do
       end do

    end do aale
  end do
  deallocate(ma)
end subroutine guardaveinsoff_2



subroutine mat
  ma=0
      do i=1,num_active_cells
        if (knots(i)==1) then ; ma(i)=1.0
        else
          if (DiffState(i)>Int) ma(i)=0.1
          if (DiffState(i)>Set) ma(i)=1.0
        end if
      end do
  vamax=maxval(ma)
  vamin=minval(ma)
end subroutine



subroutine get_rainbow(val,minval,maxval,c)
  real*8, intent(in) :: val,maxval,minval
  real*8, intent(out) :: c(4)

  real*8 :: f

  if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
  else ! probably maxval==minval
   f = 0.5
  endif

  if (f < .07) then
   c(1) = 0.6
   c(2) = 0.6
   c(3) = 0.6
   c(4) = 0.8
  elseif (f < .2) then
   c(1) = 1.0
   c(2) = f
   c(3) = 0.0
   c(4) = 0.5
  elseif (f < 1.0) then
   c(1) = 1.0
   c(2) = f*3
   c(3) = 0.0
   c(4) = 1.0
  else
   c(1) = 1.0
   c(2) = 1.0
   c(3) = 0.0
   c(4) = 1.0
  endif

end subroutine get_rainbow



subroutine get_rainbow_knot(val,minval,maxval,c,i)
  real*8, intent(in) :: val,maxval,minval
  real*8, intent(out) :: c(4)

  real*8 :: f

  if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
  else ! probably maxval==minval
   f = 0.5
  endif

  if (knots(i) == 1) then
   c(1) = 0
   c(2) = 1
   c(3) = 1
   c(4) = 0
  elseif (f < .07) then
   c(1) = 0.6
   c(2) = 0.6
   c(3) = 0.6
   c(4) = 0.8
  elseif (f < .2) then
   c(1) = 1.0
   c(2) = f
   c(3) = 0.0
   c(4) = 0.5
  elseif (f < 1.0) then
   c(1) = 1.0
   c(2) = f*3
   c(3) = 0.0
   c(4) = 1.0
  else
   c(1) = 1.0
   c(2) = 1.0
   c(3) = 0.0
   c(4) = 1.0
  endif

end subroutine get_rainbow_knot



end module esclec


!---------------------------------------------------------------------------

!***************************************************************************
!***************  PROGRAMA           ****************************************
!***************************************************************************

program tresdac

  use coreop2d
  use esclec
  implicit none

  integer   sstep,nt
  real*8    parapo(32)
  character*12 cu,cd,ct,cq
  character*12 acu,acd,act_str,acq
  character*2 dcu,dcd,dct
  character*6 dcq
  character*31 nff
  character*51 nfoff,nfes,accaufolder,caufolder
  character*2 dir_suffix
  integer paras(19)
  integer sis,siss
  integer prev_num_active_cells

  integer idi

  paras(1)=3 ; paras(2)=4 ;paras(3)=5 ;paras(4)=7 ;paras(5)=8 ;
  paras(6)=9 ; paras(7)=11 ;paras(8)=13 ;paras(9)=14 ;paras(10)=15 ;
  paras(11)=17 ; paras(12)=18 ;paras(13)=19 ;paras(14)=20 ;paras(15)=21 ;
  paras(16)=23 ; paras(17)=27 ;paras(18)=28 ;paras(19)=29 

  call getarg(1,cac) !input file
  call getarg(2,caufolder) !folder to save output
  call getarg(3,cau) !output file
  call getarg(4,cad) ! number of time steps per save (i.e for ./runt.e ./P2.txt . P2.txt 9000 1, it is the 9000). 
  call getarg(5,cass) !repetition of the save_blocks, it safe a file at this step, so an example input of ./runt.e ./P2.txt . P2.txt 9000 2, it is the "2". It will run from 9000x2 iterations, saving an output file every 9000 time steps

  close (4)

  READ(cad,*)iterationtotal
  READ(cass,*)sstep

  if (sstep<0) then ; siss=-1 ; else ;siss=1 ;end if

  max_z_layers=4
  call initialconditions
  call initialize_from_parameter_file
  call allocateinitialstate

  prev_num_active_cells=num_active_cells
  call setparams(1)
  num_active_cells=prev_num_active_cells
  
  call initact  ! Sets initial activation concentration (Ina).
  
  do iti=1,sstep 
    ii=0
    write (ct,*) (idi+ii)*sis
    write (cq,*) iti*iterationtotal
 
    acu=adjustl(cu)
    acd=adjustl(cd)
    act_str=adjustl(ct)
    acq=adjustl(cq)
    accaufolder=trim(adjustl(caufolder))

    dcu=acu ; dcd=acd ;dct=act_str ;dcq=acq

    nff = trim(adjustl(accaufolder)) // "/" // trim(adjustl(dcq)) // "_" // trim(adjustl(cau))

    do i=1,len(nff)
      if (nff(i:i)==" ") nff(i:i)="_"
    end do

    nfoff=nff(1:31)//"_"//".off"
    nfes=nff(1:26)//"_"//".txt"
    nfpro=nfpro(1:15)//"_progressbar"//".txt"

    temps=0

    nt=iterationtotal !*iti
    call iteration(nt)

    dir_suffix=trim(adjustl(dcu))

    open(2,file=nfoff,iostat=file_opened_success)
    call guardaveinsoff_2(neigh)
    close(2)
  end do

  parap(:,1)=parapo
  call setparams(1)

end program tresdac

!***************************************************************************
!***************  FI  PROGRAMA      ****************************************
!***************************************************************************