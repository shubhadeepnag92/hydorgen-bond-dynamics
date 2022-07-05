!	NON-EQUILIBRIUM MONTE CARLO PROGRAM TO SIMULATE THE MIXTURE OF HYDROCARBONS IN FLEXIBLE ZEOLITE NaY
!	Developer - Shubhadeep Nag
!	Latest Change Date - 23/7/2021

!	MODIFICATION FROM v2 
	!	neighbor list is implemented
	!	change com_p(:) to com_p1,com_p2,com_3 & com_new(:) to com_new1,com_new2,com_new3
	!	incorporate molecule transition to unit cell ( 1 x 1 x 1 ) during position calculation for MORE than ( 1 x 1 x 1 )
!	MODIFICATION FROM v3 
	!	inhomogeneity sampling is inplemented along [111]
	!	no PBC along x,y & z respectively
!	MODIFICATION FROM v5 
	!	write com as trajectory file
!	MODIFICATION FROM v6 
	!	put hot zone along [011] instead of [111]
	!	increase efficiency
!   MODIFICATION FROM v7 
    !   put hot zone along x-axis at every quarter of unit boxlength
!   MODIFICATION FROM v8 
    !   some minor change in variables etc.
!   MODIFICATION FROM v9
    !   put flexibility of host medium
!	MODIFICATION FROM v10
	!	put 50% acceptence criteria
	!	put long range correction of interaction energy
!	MODIFICATION FROM v11
	!	put one sided repulsive potential
!	MODIFICATION FROM v13
	!	little change in code to make faster (for 10x2x2 earlier it took 0.48sec, now it takes 0.42sec)
!	MODIFICATION FROM v20
	!	Intra-angular motion included
!	MODIFICATION FROM v21
	!	put angle rotation to atomic position in a subroutine and put a bondlength array to calculate upto 16th decimal
!   MODIFICATION FROM v22
    !   change the transition probability : new transition from one to another thermal bath will divided into three steps ; (a) old
    !   position to boundary - delta (b) boundary - delta to boundary + delta and (c) boundary + delta to new position
!	MODIFICATION IN v25
	!	make the transition of each step independent
	!	accept or reject three independent steps independently
!	MODIFICATION IN v26
	!	can make different number of sites of two hydrocarbons
!	MODIFICATION IN v27
	!	gaussian temperature profile

	module parameter

	real*8 :: boxlength,boxlengthx,boxlengthy,boxlengthz,halfboxlengthx,halfboxlengthy,halfboxlengthz,density
	real*8 :: pe1(-3000:3000),pe2(-3000:3000),pe1_gg(-3000:3000),pe1_gh(-3000:3000),pe2_gg(-3000:3000),pe2_gh(-3000:3000)    
        real*8 :: E_g11,E_g12,E_g22,E_g1h,E_g2h,guestguest_E,guesthost_E
	real*8 :: dx,dy,dz,d,rc,s12,s6,E,E_totalinter,dxaa,dyaa,dzaa,ratio,rsq,E_old,dr_max,E_new,bond_length
	real*8 :: xrand,yrand,zrand,RT,prob,xt,R,T,dE,accp,rej,phi,psi,rand,sin_psi,sin_phi,two_pie,sin_theta,cos_phi,cos_psi
	real*8 :: rot(3,3),onemcospsi,nz,nzsq,nx,nxsq,ny,nysq,cos_theta,nxy,nyz,nxz,dpsi_max,E_present,ls,rs,Th,Tl,RTh,RTl
	real*8 :: xn1,yn1,zn1,xn2,yn2,zn2,rn1,rn2,dotneu,dotden,dotratio,ndhangle,cos_dh,sin_dh,dh_max,detrot,onemcosdh,pt_distance
	real*8 :: com_new1,com_new2,com_new3,com_p1,com_p2,com_p3,pt_distance_old,pt_distance_new,E_vdw
	real*8 :: dx1,dx2,dy1,dy2,dz1,dz2,E_totalintra,E_bond,E_angle,d1,d2,halfboxlength,angle,equ_distSiO,force_constSiO
	real*8 :: force_constOSiO,equ_angleOSiO,force_constSiOSi,equ_angleSiOSi
	real*8 :: accphh,accptr,accprot,accpdih,rejhh,E_ub,force_constub,equ_ub,start,finish
	real*8 :: dr_maxhh,hostaccpratio,guest1accpratiotrans,guest1accpratiorot,guest1accpratiodih
	real*8 :: guest2accpratiotrans,guest2accpratiorot,guest2accpratiodih,accpgg1trans,accpgg1rot,accpgg1dih
	real*8 :: accpgg2trans,accpgg2rot,accpgg2dih,rejgg1trans,rejgg1rot,rejgg1dih,rejgg2trans,rejgg2rot,rejgg2dih
        real*8 :: accpgg1trans1,accpgg1trans2,accpgg1trans3,accpgg2trans1,accpgg2trans2,accpgg2trans3
        real*8 :: rejgg1trans1,rejgg1trans2,rejgg1trans3,rejgg2trans1,rejgg2trans2,rejgg2trans3
        real*8 :: dr_maxgg1trans,dr_maxgg1rot,dr_maxgg1dih,dr_maxgg2trans,dr_maxgg2rot,dr_maxgg2dih,separation_factor,y
	real*8 :: c0(1:3),c1(1:3),c2(1:3),c3(1:3),initial_rg_atomic(1:3,1:10,1:1000)
	real*8 :: equi_angle(2),angle_parm(2),nxab,nyab,nzab
	real*8 :: uxab,uyab,uzab,uxbc,uybc,uzbc,theta_angle,Eintra_angle,Eold_angle,Enew_angle,ntheta_angle,dangle_max
	real*8 :: nxbc,nybc,nzbc,nxybc,nxzbc,nyzbc,nxsqbc,nysqbc,nzsqbc,dtheta_angle,bond_angle_length,ux,uy,uz
	real*8 :: dE_x,dE_y,dE_z,prob_x,prob_x1,prob_x2,prob_x3,prob_y,prob_z,Ex_new,Ex_old,Ey_new,Ey_old,Ez_new,Ez_old,RThomo,df_pos
	real*8 :: xtemp(1:10000),boundary,leftneo,rightneo,leftnp,rightnp,alpha
        real*8 :: equ_distAlO,force_constAlO,equ_angleOAlO,force_constOAlO,equ_angleAlOSi,force_constAlOSi,equ_ubAl,force_constubAl
	real*8 :: com_inter_x

	real*8,allocatable :: com(:,:),rg_atomic(:,:,:),sigmag(:),epsg(:),sigmasqgg(:,:),fourepsgg(:,:),sigmasqhh(:,:),fourepshh(:,:)
	real*8,allocatable :: rg_natomic(:,:),atomic_p(:,:),rh_atomic(:,:),sigmah(:),epsh(:),sigmasqgh(:,:),fourepsgh(:,:)
	real*8,allocatable :: dhangle(:),vx(:,:),vy(:,:),vz(:,:),ratom(:,:),rcage(:,:),cagetocomdist(:),initial_rh_atomic(:,:)
	real*8,allocatable :: aw(:,:),bw(:,:),cw(:,:),dw(:,:),com_wind_dist(:),rhatomic_p(:),rh_natomic(:),temp(:),bondlength(:,:,:)

	integer,allocatable :: nlistgg(:),listgg(:,:),nlistgh(:),listgh(:,:),gmoltype(:),gatomtype(:,:),hatomtype(:),dhtype(:)
	integer,allocatable :: listcagewindow(:,:),bond_nlist(:),bond_list(:,:),nangle_list(:),vdw_nlist(:),vdw_list(:,:)

	integer :: m,j,igatom,jgatom,igmol,jgmol,ngmol,igmoltype,ngmoltype,natom,igatomtype,jgatomtype,ngatomtype,cage_com_pos
	integer :: lx,ly,lz,translation,imc,rotation,ngatom,ihatom,nhatom,nhatomtype,ihatomtype,jhatomtype,idh,dh_rotation,jhatom
	integer :: i,nmc,ibin1,ibin2,count1(-5000:5000),count2(-5000:5000),neq,ncagecenter,icagecenter,nlistcagewind,iwindow,nwindow
	integer :: jcagecenter,ilistcagewind,jlistcagewind,k,icagecenter_old,icagecenter_new,ilistcagewind_old,ilistcagewind_new
	integer :: khatom,movement,iwinfo,nwinfo,count1left,count1right,count2left,count2right,nsteady,simulation,restart,ix,iangle
	integer :: ngatom1,ngatom2,ndh(1:2),idhatom(1:4,1:10,1:2),igatom1,igatom2,igatom3,igatom4,nangle(2),index_angle(1:3,1:10,1:2)
	integer :: iangle1,iangle2,iangle3,jangle,part,cartesian,itemp
	integer :: nigatom,njgatom,temp_profile_type

	character(10) :: a

	end module parameter

	program mcmixhydro

!	MONTECARLO PROGRAM

	use parameter
	implicit none

        start = 0.0d0 ; finish = 0.0d0  

        call cpu_time(start)

!	OPEN FILES TO READ AND  WRITE

	open(unit=1,file='inputmcal.dat')
	open(unit=2,file='host.xyz')
	open(unit=3,file='guestposition.xyz')
	open(unit=4,file='atomic_hydrocarbon.xyz')
	open(unit=5,file='neighborlist.dat')
	open(unit=16,file='cagecentreszeoY.xyz')
	open(unit=7,file='window_equation.dat')
	open(unit=8,file='trajectory.xyz')
	open(unit=9,file='com_trajectory.xyz')
	open(unit=10,file='output.dat')
	open(unit=11,file='energy_output.dat')
        open(unit=12,file='final_guestcom.xyz')
	open(unit=13,file='final_guestatomic.xyz')
	open(unit=14,file='acceptance_ratio.dat')
        open(unit=15,file='final_hostatom.xyz')
	open(unit=17,file='proper_traj.dat')
	open(unit=18,file='temperature_data.dat')
        open(unit=21,file='separation_factor.dat')	

!	READ & CALCUALTE INITIAL VARIABLES

	read(1,*)
	read(1,*) simulation,restart,temp_profile_type
	read(1,*)
	read(1,*) lx,ly,lz
	read(1,*)
	read(1,*) dr_maxhh,nwinfo
	read(1,*)
	read(1,*) boxlength,rc,nhatomtype,ngmoltype,ngatomtype,dr_max,R
	read(1,*)
	read(1,*) nmc,nsteady,dpsi_max,rs,ls,T
	read(1,*)
	read(1,*) ndh(1),ndh(2)
	read(1,*)
	do i=1,ndh(1)
	read(1,*) idhatom(:,i,1)
	enddo !i
	read(1,*)
	do i=1,ndh(2)
	read(1,*) idhatom(:,i,2)
	write(*,*) 'i'
	enddo !i
	read(1,*)
	read(1,*) dh_max,bond_length
	read(1,*)
	read(1,*)
	read(1,*) c0(1),c1(1),c2(1),c3(1)
	read(1,*) c0(2),c1(2),c2(2),c3(2)

	dr_maxgg1trans = dr_max   ; dr_maxgg2trans = dr_max
	dr_maxgg1rot   = dpsi_max ; dr_maxgg2rot   = dpsi_max
	dr_maxgg1dih   = dh_max   ; dr_maxgg2dih   = dh_max
	boxlengthx = boxlength * lx; boxlengthy = boxlength * ly; boxlengthz = boxlength * lz 
	halfboxlengthx = boxlengthx * 0.5d0;  halfboxlengthy = boxlengthy * 0.5d0;  halfboxlengthz = boxlengthz * 0.5d0
	RT = R * T; RThomo = R * T
        count1left = 0; count1right = 0 ;count2left = 0; count2right =0
	leftneo = 0.0d0; rightneo = 0.0d0; leftnp = 0.0d0; rightnp = 0.0d0
!        density = ( hostmass * lx * ly * lz + guestmass * ngmol ) / ( boxlengthx * boxlengthy * boxlengthz * 10E-30 )

!	READ COORDINATES OF HOST ZEOLITE ATOMS

	read(2,*) nhatom
	read(2,*)

!	READ THE CENTER OF MASS COORDINATES OF MIXTURE OF MOLECULES

	read(3,*) ngmol
	read(3,*)

	read(4,*) ngatom1,ngatom2
	read(4,*)

	read(16,*) ncagecenter
	read(16,*)

	if(ngatom1.ge.ngatom2) ngatom = ngatom1
	if(ngatom1.lt.ngatom2) ngatom = ngatom2

	allocate ( hatomtype(1:nhatom),rh_atomic(1:3,1:nhatom),ratom(1:3,1:ngatom),rcage(1:3,1:ncagecenter) )
	allocate ( gmoltype(1:ngmol),gatomtype(1:ngatom,1:ngmol),rg_natomic(1:3,1:ngatom) )
	allocate ( com(1:3,1:ngmol),rg_atomic(1:3,1:ngatom,1:ngmol),nlistgh(1:ngmol),listgh(1:nhatom,1:ngmol) )
	allocate ( sigmag(1:ngatomtype),epsg(1:ngatomtype),sigmasqgg(1:ngatomtype,1:ngatomtype),fourepsgg(1:ngatomtype,1:ngatomtype) )
	allocate ( nlistgg(1:ngmol),listgg(1:ngmol,1:ngmol),atomic_p(1:3,1:ngatom),dhangle(1:ndh(2)) )
	allocate ( sigmah(1:nhatomtype),epsh(1:nhatomtype),sigmasqgh(1:ngatomtype,1:nhatomtype),fourepsgh(1:ngatomtype,1:nhatomtype) )
	allocate ( vx(1:ngatom,1:ngatom),vy(1:ngatom,1:ngatom),vz(1:ngatom,1:ngatom),temp(1:25000) )
	allocate ( cagetocomdist(1:ncagecenter),vdw_nlist(1:nhatom),vdw_list(1:nhatom,1:nhatom) )
	allocate ( bond_nlist(1:nhatom),bond_list(1:nhatom,1:nhatom),nangle_list(1:nhatom),rhatomic_p(1:3),rh_natomic(1:3) )
	allocate ( sigmasqhh(1:nhatomtype,1:nhatomtype),fourepshh(1:nhatomtype,1:nhatomtype),dhtype(1:ngmol) )
	allocate ( bondlength(1:ngatom,1:ngatom,1:ngmol),initial_rh_atomic(1:3,1:nhatom) )

!	READ THE HOST ATOMIC COORDINATES

	do ihatom=1,nhatom
		read(2,*) hatomtype(ihatom),rh_atomic(:,ihatom)
                initial_rh_atomic(:,ihatom) = rh_atomic(:,ihatom)
	enddo !ihatom

	do igmol=1,ngmol
	   read(3,*) gmoltype(igmol),dhtype(igmol),com(:,igmol)
	enddo !imol

!	READ THE ATOMIC COORDINATES OF ATOMS OF MOLECULE WRT COM #1 TYPE

	do igmol=1,int(ngmol/2)

	   do igatom=1,ngatom1
		read(4,*) gatomtype(igatom,igmol),rg_atomic(:,igatom,igmol)
		initial_rg_atomic(:,igatom,igmol) = rg_atomic(:,igatom,igmol)
	   enddo !igatom

!	CALCULATTION OF BOND-LENGTH UPTO 16th DECIMAL

	do igatom=1,ngatom1
	    do jgatom = 1,ngatom1
		if(igatom.ne.jgatom) then
			d = sqrt((rg_atomic(1,igatom,igmol)-rg_atomic(1,jgatom,igmol))**2+(rg_atomic(2,igatom,igmol) &
				- rg_atomic(2,jgatom,igmol))**2 + (rg_atomic(3,igatom,igmol)-rg_atomic(3,jgatom,igmol))**2 )
			bondlength(igatom,jgatom,igmol) = d
			!write(*,*) igmol,igatom,jgatom,bondlength(igatom,jgatom,igmol)
		endif
	    enddo !jgatom
	enddo !igatom

	enddo !igmol

!	READ THE ATOMIC COORDINATES OF ATOMS OF MOLECULE WRT COM #2 TYPE

	do igmol=int(ngmol/2)+1,ngmol

	   do igatom=1,ngatom2
		read(4,*) gatomtype(igatom,igmol),rg_atomic(:,igatom,igmol)
		initial_rg_atomic(:,igatom,igmol) = rg_atomic(:,igatom,igmol)
	   enddo !igatom

!	CALCULATTION OF BOND-LENGTH UPTO 16th DECIMAL

	do igatom=1,ngatom2
	    do jgatom = 1,ngatom2
		if(igatom.ne.jgatom) then
			d = sqrt((rg_atomic(1,igatom,igmol)-rg_atomic(1,jgatom,igmol))**2+(rg_atomic(2,igatom,igmol) &
				- rg_atomic(2,jgatom,igmol))**2 + (rg_atomic(3,igatom,igmol)-rg_atomic(3,jgatom,igmol))**2 )
			bondlength(igatom,jgatom,igmol) = d
			!write(*,*) igmol,igatom,jgatom,bondlength(igatom,jgatom,igmol)
		endif
	    enddo
	enddo

	enddo !igmol

	write(10,*)


!	READ THE CAGE CENTER POSITION

	do icagecenter = 1,ncagecenter
	   read(16,*) rcage(:,icagecenter)
	enddo !icagecenter

	read(16,*)
	read(16,*) nlistcagewind

	allocate ( listcagewindow(1:nlistcagewind,1:ncagecenter),com_wind_dist(1:nlistcagewind),dw(1:nlistcagewind,1:icagecenter) )
	allocate ( aw(1:nlistcagewind,1:icagecenter),bw(1:nlistcagewind,1:icagecenter),cw(1:nlistcagewind,1:icagecenter) )

	do icagecenter = 1,ncagecenter

		read(16,*) i,(listcagewindow(j,icagecenter),j=1,nlistcagewind)

	enddo !icagecenter

!	write(20,*) listcagewindow(:,:)


!	READ WINDOW CENTER POSITION AND EQUATION

	read(7,*)

	do icagecenter = 1,ncagecenter

	   do ilistcagewind = 1,nlistcagewind
		k=listcagewindow(ilistcagewind,icagecenter)
		
      	read(7,*) j,i,aw(ilistcagewind,icagecenter),bw(ilistcagewind,icagecenter),cw(ilistcagewind,icagecenter)&
		,dw(ilistcagewind,icagecenter)
	      !write(20,*) j,i,aw(k,icagecenter),bw(k,icagecenter),cw(k,icagecenter),dw(k,icagecenter)

	   enddo !ilistcagewindow

		!write(20,*) i,listcagewindow(:,icagecenter)

	enddo !icagecenter

	!write(20,*) '2',(listcagewindow(j,2),j=1,nlistcagewindow)



!	WRITE IN OUTPUT

	write(10,1) '|********************MONTE CARLO SIMULATION OF MIXTURE OF HYDROCARBONS**********************|'
1	format(20x,A)
	write(10,*) 
	if(simulation.eq.0) write(10,2) '$********************EQUILIBRIUM   SIMULATION********************$'
2	format(40x,A)
	if(simulation.eq.1) write(10,3) '$*****NON-EQUILIBRIUM SIMULATION (with hot-zone perpendicular to cage-window distance)****$'
3	format(20x,A)
	if(simulation.eq.2) write(10,4) '$*****NON-EQUILIBRIUM SIMULATION (with hot-zone perpendicular to X-axis)*******$'
4	format(10x,A)
	write(10,*)
	write(10,5) '$*****INITIAL CONFIGURATION DETAILS*******$  '
5	format(45x,A)
	write(10,*)
	write(10,6) 'No of types of molecules',ngmoltype
6	format(2x,A,40x,i5)
	write(10,7) 'No of types of atoms',ngatomtype
7	format(2x,A,44x,i5)
	write(10,8) 'No of molecules',ngmol
8	format(2x,A,49x,i5)
	write(10,9) 'No of atoms in first molecule of mixture of hydrocarbon',ngatom
9	format(2x,A,9x,i5)
	write(10,10) 'No of atoms in second molecule of mixture of hydrocarbon',ngatom
10	format(2x,A,8x,i5)
	write(10,11) 'No of unitcell along x,y and z respectively',lx,ly,lz
11	format(2x,A,22x,i4,10x,i4,10x,i4)
	write(10,12) 'Boxlength along x,y and z respectively',boxlengthx,real(boxlengthy),real(boxlengthz)
12	format(2x,A,30x,f13.8,3x,f11.8,3x,f11.8)
	write(10,*)
	write(10,*)
	write(10,21) '$**********************SIMULATION DETAILS**************************$'
21	format(35x,A)
	write(10,*)
	write(10,22) 'No of MC steps',nmc
22	format(2x,A,40x,i10)
	write(10,23) 'No of production steps',(nmc-nsteady)
23	format(2x,A,32x,i10)
	write(10,34) 'Interval of printing data',nwinfo 
34	format(2x,A,34x,i5)
	write(10,24) 'Cut-off radius',real(rc)
24	format(2x,A,46x,f14.8)
	write(10,25) 'Gas constant(kJ/mol)',real(R)
25	format(2x,A,39x,f14.8)
	write(10,26) 'Ambient Temperature(K)',real(Tl)
26	format(2x,A,39x,f14.8)
	write(10,27) 'Hot-zone Temperature(K)',real(Th)
27	format(2x,A,38x,f14.8)
	write(10,28) 'Right side of hot-zone from window',real(rs)
28	format(2x,A,25x,f14.8)
 	write(10,29) 'Left side of hot-zone from window',real(ls)
29	format(2x,A,26x,f14.8)
	write(10,*)
        write(11,39) 'No of MC  steps','E_g1g1','E_g1g2','E_g2g2','E_g1h','E_g2h'
39      format(2x,A,14x,A,14x,A,14x,A,14x,A,14x,A)
	write(14,42) 'NMC','host_accp','guest1_trans_accp','guest1_rot_accp','guest2_trans_accp','guest2_rot_accp',&
	'guest2_dih_accp'
42	format(A,6x,A,6x,A,6x,A,6x,A,6x,A,6x,A)	

!	READ THE FORCE-FILED PARAMETER

	write(10,13) '$********************FORCE-FIELD PARAMETERS**********************$'
13	format(35x,A)
	write(10,*)
	write(10,*)

!	READ BOND LENGTH PARAMETER

	write(10,14) ' ~~~BOND PARAMETER~~~ '
14	format(10x,A)
	write(10,*)
	write(10,15) 'equilibrium distance(Angs)','force constant(kJ/mol/A^2)'
15	format(20x,A,10x,A)

	read(1,*)
	read(1,*)
	read(1,*)	
	read(1,*) a,equ_distSiO,force_constSiO
	write(10,16) a,real(equ_distSiO),real(force_constSiO)
16	format(10x,A,5x,f14.8,10x,f20.8)
        read(1,*) a,equ_distAlO,force_constAlO
        write(10,16) a,real(equ_distAlO),real(force_constAlO)

	write(10,*)

!	READ ANGLE LENGTH PARAMETER

	write(10,14) ' ~~~ANGLE PARAMETER~~~ '
	write(10,*)
	write(10,15) 'equilibrium angle(deg)','force constant(kJ/mol/rad^2)'
	write(10,*)

	read(1,*)
	read(1,*)
	read(1,*)
	read(1,*) a,equ_angleOSiO,force_constOSiO
	write(10,16) a,real(equ_angleOSiO),real(force_constOSiO)
	read(1,*) a,equ_angleOAlO,force_constOAlO
        write(10,16) a,real(equ_angleOAlO),real(force_constOAlO)
        read(1,*) a,equ_angleSiOSi,force_constSiOSi
	write(10,16) a,real(equ_angleSiOSi),real(force_constSiOSi)
        read(1,*) a,equ_angleAlOSi,force_constAlOSi
        write(10,16) a,real(equ_angleAlOSi),real(force_constAlOSi)

	write(10,*)
	

!	READ UREY-BRADLEY PARAMETER

	write(10,14) ' ~~~UREY-BRADLEY PARAMETER~~~ '
	write(10,*)
	write(10,15) 'equilibrium distance(Angs)','force constant(kJ/mol/A^2)'
	write(10,*)

	read(1,*)
	read(1,*)
	read(1,*)
	read(1,*) a,equ_ub,force_constub
	write(10,16) a,real(equ_ub),real(force_constub)
        read(1,*) a,equ_ubAl,force_constubAl
        write(10,16) a,real(equ_ubAl),real(force_constubAl)
	write(10,*)

!	READ THE LJ PARAMETER

	write(10,14) ' ~~~LENNARD - JONES PARAMETER~~~ '
	write(10,*)
	write(10,17) 'atomtype','index','sigma(A)','epsilon(KJ/mol)'
17	format(10x,A,10x,A,10x,A,10x,A)
	write(10,*)
	write(10,*) '		GUEST			'
	write(10,*)

	read(1,*)
	read(1,*)
	read(1,*)
	read(1,*) 

!	GUEST PARAMETERS
	
	do igatomtype=1,ngatomtype
	   read(1,*) a,sigmag(igatomtype),epsg(igatomtype)
	   write(10,18) a,igatomtype,real(sigmag(igatomtype)),real(epsg(igatomtype))
18	   format(10x,A,10x,i1,10x,f14.8,10x,f14.8)
	enddo !iatomtype

!	HOST PARAMETERS

	read(1,*)
	read(1,*)
	read(1,*)

	write(10,*)
	write(10,*) '		HOST		'
	write(10,*)

	do ihatomtype=1,nhatomtype
	   read(1,*) a,sigmah(ihatomtype),epsh(ihatomtype)
	   write(10,18) a,ihatomtype,real(sigmah(ihatomtype)),real(epsh(ihatomtype))
	enddo !ihatomtype

!	CALCULATION OF INTERACTION LJ PARAMETER

!	HOST-HOST

	do ihatomtype = 1,nhatomtype
	   do jhatomtype = 1,nhatomtype
	sigmasqhh(ihatomtype,jhatomtype)=0.25d0*(sigmah(ihatomtype)+sigmah(jhatomtype))*(sigmah(ihatomtype)+sigmah(jhatomtype))
		fourepshh(ihatomtype,jhatomtype)=4.0d0*sqrt(epsh(ihatomtype)*epsh(jhatomtype))
		write(10,*) 'host-host',ihatomtype,jhatomtype,sigmasqhh(ihatomtype,jhatomtype),fourepshh(ihatomtype,jhatomtype)
	   enddo !jatomtype
	enddo !iatomtype

!	GUEST-GUEST

	do igatomtype = 1,ngatomtype
	   do jgatomtype = 1,ngatomtype
	sigmasqgg(igatomtype,jgatomtype)=0.25d0*(sigmag(igatomtype)+sigmag(jgatomtype))*(sigmag(igatomtype)+sigmag(jgatomtype))
		fourepsgg(igatomtype,jgatomtype)=4.0d0*sqrt(epsg(igatomtype)*epsg(jgatomtype))
		write(10,*) 'guest-guest',igatomtype,jgatomtype,sigmasqgg(igatomtype,jgatomtype),fourepsgg(igatomtype,jgatomtype)
	   enddo !jatomtype
	enddo !iatomtype

!	GUEST-HOST

	do igatomtype = 1,ngatomtype
	   do jhatomtype = 1,nhatomtype
	sigmasqgh(igatomtype,jhatomtype)=0.25d0*(sigmag(igatomtype)+sigmah(jhatomtype))*(sigmag(igatomtype)+sigmah(jhatomtype))
		fourepsgh(igatomtype,jhatomtype)=4.0d0*sqrt(epsg(igatomtype)*epsh(jhatomtype))
		write(10,*) 'guest-host',igatomtype,jhatomtype,sigmasqgh(igatomtype,jhatomtype),fourepsgh(igatomtype,jhatomtype)
	   enddo !jatomtype
	enddo !iatomtype

	write(10,*)
	write(10,14) ' ~~~DIHEDRAL PARAMETER~~~ '
	write(10,*)
	write(10,19) 'C0','C1','C2','C3'
19	format(20x,A,20x,A,20x,A,20x,A)
	write(10,*)
	write(10,20) c0(1),c1(1),c2(1),c3(1)
20	format(10x,f14.8,10x,f14.8,10x,f14.8,10x,f14.8)
	write(10,*)

!	read angle value

	read(1,*)
	read(1,*)
	read(1,*)
	read(1,*) nangle(1),nangle(2),dangle_max

!	read angle parameters

	read(1,*)
	read(1,*)
	read(1,*)

	read(1,*) equi_angle(1),angle_parm(1)
	read(1,*) equi_angle(2),angle_parm(2)

	write(10,*) equi_angle(1),angle_parm(1)
	write(10,*) equi_angle(2),angle_parm(2)

	read(1,*)

!	write the angle of molecules

	write(10,*)
	write(10,*) 'No of angles of two types of hydrocarbon'
	write(10,*) nangle(1),nangle(2)
	write(10,*)

!	read the angle parameter
	write(10,*)
	write(10,*) 'Atom index of #1 hydrocarbon'
	read(1,*)
	do iangle = 1,nangle(1)
	read(1,*) index_angle(1:3,iangle,1)!,indexrot_angle1(iangle)
	write(10,*) '1',iangle,index_angle(1:3,iangle,1)
	enddo !iangle
	write(10,*)

	write(10,*)
	write(10,*) 'Atom index of #2 hydrocarbon'
	read(1,*)
	do iangle = 1,nangle(2)
	read(1,*) index_angle(1:3,iangle,2)!,indexrot_angle2(iangle)
	write(10,*) '2',iangle,index_angle(1:3,iangle,2)
	enddo !iangle
	write(10,*)

!	store temperature data

!	read(16,*)

	do ix = 1,25000!int(boxlengthx/0.01d0)

		read(18,*) y,temp(ix)
                
	enddo !ix


!	CREATE THE NEIGHBOR LIST

	call neighborlist

!	CREATE THE NEIGHBOR LIST FOR FLEXIBILITY

	call neighborlistflex

!	CALCULATE THE TOTAL ENERGY

	call totalenergy

	call intramolenergy

	do igmol = 1, ngmol
	   call energyith
	   E_present = E_present + E
	enddo !igmol

	!write(*,*) E_present*0.5d0


	E_present = E_totalinter

	write(10,*)
	write(10,*) 'Total Initial Intramolecular Energy in KJ/mol',E_totalintra
	write(10,*)
	write(10,*) 'Total Initial Intermolecular Energy in KJ/mol',E_totalinter

        write(10,44) 'No of MC steps','Total Host-Host Energy','Total Guest-Host energy','Separation Factor'
44      format(A,10x,A,10x,A,10x,A)

						!	MONTE - CARLO SIMULATION	!

	accphh = 0.0d0 ; rejhh = 0.0d0 
	accpgg1trans = 0.0d0 ; rejgg1trans = 0.0 ; accpgg1rot = 0.0d0 ; rejgg1rot = 0.0 ;accpgg1dih = 0.0d0 ; rejgg1dih = 0.0
	accpgg2trans = 0.0d0 ; rejgg2trans = 0.0 ; accpgg2rot = 0.0d0 ; rejgg2rot = 0.0 ;accpgg2dih = 0.0d0 ; rejgg2dih = 0.0
        accpgg1trans1 = 0.0d0; accpgg1trans2 = 0.0d0; accpgg1trans3 = 0.0d0
        accpgg2trans1 = 0.0d0; accpgg2trans2 = 0.0d0; accpgg2trans3 = 0.0d0


        write(17,*) 'file to store the character of the molecule, &
                        x,y,z-coordinates, guest-guest energy, guest-host energy at every 1000 MC steps'
        write(17,*)

	do imc=1,nmc

	write(*,*) imc

        if(mod(imc,1000).eq.0) then

	write(17,*) 'No of MC Steps',imc
	write(17,*)

        endif


!	Only If you can not keep the coordinates correct after rotation, angular and dihedral motion

	if(mod(imc,5000).eq.0) then
        do ihatom = 1,nhatom
                rh_atomic(:,ihatom) = initial_rh_atomic(:,ihatom)
        enddo!ihatom
	do igmol=1,int(ngmol/2)
	   do igatom=1,ngatom1
	!	rg_atomic(:,igatom,igmol) = initial_rg_atomic(:,igatom,igmol)
	   enddo
	enddo
	do igmol=int(ngmol/2)+1,ngmol
	   do igatom=1,ngatom2
	!	rg_atomic(:,igatom,igmol) = initial_rg_atomic(:,igatom,igmol)
	   enddo
	enddo
	endif

		!-----------------------------------FLEXIBILITY TEST OF HOST SYSTEM---------------------!

!	CREATE THE NEIGHBOR LIST FOR FLEXIBILITY TEST

	if(mod(imc,10).eq.0)        call neighborlistflex

	do ihatom = 1,0!nhatom

       ! if(bond_nlist(ihatom).ne.0) then

!	CALCULATION OF OLD ENERGY

	   movement = 0
	   call intramolenergyith
	   E_old = E_bond+E_angle+E_vdw+E_ub
!	   if(ihatom.eq.577) write(30,*) E_old,E_bond,E_angle,E_vdw,E_ub

!	RANDOM MOVEMENT OF ATOMS

	   !rand=0.0d0
	   call random_number(xrand)
	   rh_natomic(1)= rh_atomic(1,ihatom) + (2*xrand - 1)*dr_maxhh
	   if(rh_natomic(1).gt.boxlengthx)  rh_natomic(1) = rh_natomic(1) - boxlengthx
	   if(rh_natomic(1).lt.0.0d0) 	    rh_natomic(1) = rh_natomic(1) + boxlengthx
	
	   !rand=0.0d0
	   call random_number(yrand)
	   rh_natomic(2)= rh_atomic(2,ihatom) + (2*yrand - 1)*dr_maxhh
	   if(rh_natomic(2).gt.boxlengthy)  rh_natomic(2) = rh_natomic(2) - boxlengthy
	   if(rh_natomic(2).lt.0.0d0) 	    rh_natomic(2) = rh_natomic(2) + boxlengthy

	   !rand=0.0d0
	   call random_number(zrand)
	   rh_natomic(3)= rh_atomic(3,ihatom) + (2*zrand - 1)*dr_maxhh
	   if(rh_natomic(3).gt.boxlengthz)  rh_natomic(3) = rh_natomic(3) - boxlengthz
	   if(rh_natomic(3).lt.0.0d0) 	    rh_natomic(3) = rh_natomic(3) + boxlengthz

!        if(ihatom.eq.577) write(30,*) ihatom,real(rh_atomic(:,ihatom)),xrand,yrand,zrand
!        if(ihatom.eq.577)  write(30,*) real(rh_natomic(:))

!	CALCULATION OF NEW ENERGY

	   movement = 1
	   call intramolenergyith
	   E_new = E_bond+E_angle+E_vdw+E_ub
!	    if(ihatom.eq.577) write(30,*) E_new,E_bond,E_angle,E_vdw,E_ub

!	IMPORTANCE SAMPLING

!        write(30,*) ihatom,E_old,E_new

	   dE=E_new-E_old
           RTl = R * temp(int(rh_atomic(1,ihatom)/0.01d0) )
	   if(dE.gt.0.0d0) then
	      	prob=exp(-dE/RTl)
      	call random_number(rand)
	   endif

	   if(dE.lt.0.0d0.or.rand.lt.prob) then
		rh_atomic(:,ihatom)=rh_natomic(:)
		accphh=accphh+1.0d0
!                if(ihatom.eq.577)  write(30,*) 'accp'
	   else
		rejhh=rejhh+1.0d0
!                if(ihatom.eq.577)  write(30,*) 'rej'
	   endif

!          dE = 0.0d0; prob = 0.0d0; xt = 0.0d0
!          dE = E_new - E_old

!	  T1 = temp(int(rh_atomic(1,ihatom)/0.01d0) ) ; T2 = temp(int(rh_natomic(1)/0.01d0))
!	  RT1 = R * T1 ; RT2 = R * T2
	
!          if(T1.eq.T2) then

!               HOMOGENEOUS TRANSITION
!
!                        if (dE.gt.0.0d0) then
!                        prob = exp ( - dE/RT1 )
!                        call random_number(xt)
!                        endif

!                        if (dE.lt.0.0d0.or.xt.lt.prob) then
!                                rh_atomic(:,ihatom)=rh_natomic(:)
!		                accphh=accphh+1.0d0
!	                else
!		                rejhh=rejhh+1.0d0
!	                endif

!          endif

!          if(T1.gt.T2) then

!               HOT TO COLD TRANSITION

!                        prob=exp(+E_new/RT2-E_old/RT1+dE*(1/RT2-1/RT1))
!                        call random_number(xt)

!                        if (xt.lt.prob) then
!                                rh_atomic(:,ihatom)=rh_natomic(:)
!		                accphh=accphh+1.0d0
!	                else
!		            rejhh=rejhh+1.0d0
!	                endif

!          endif

!          if(T1.lt.T2) then

!               COLD TO HOT TRANSITION

!                        prob=exp(-E_new/RT2+E_old/RT1-dE*(1/RT2-1/RT1))
!                        call random_number(xt)

!                        if (xt.lt.prob) then
!                                rh_atomic(:,ihatom)=rh_natomic(:)
!                      		accphh=accphh+1.0d0
!	                else
!		                rejhh=rejhh+1.0d0
!	                endif

!	   endif

	   !write(30,*) ihatom,E_old,E_new,dE,prob,rand!,accp,rej

!	endif

	enddo !ihatom

			!-------------------------------GUEST-GUEST & GUEST-HOST SIMULATION---------------------------!


!	CREATE THE NEIGHBOR LIST

	if(mod(imc,5).eq.0) call neighborlist

	do igmol=1,ngmol

	  E_old = 0.0d0;E_new=0.0d0;dE=0.0d0;prob=0.0d0;xt=0.0d0

!	TRANSLATION OF COM
	
	   translation = 0
	   rotation	= 0
	   call energyith
	   E_old = E

	   call translationalmotion
	
	   !write(20,*) com(1,igmol)
	   !write(20,*) com_new1!,com_new2,com_new3

	   translation = 1
	   rotation	= 0
	   call energyith
	   E_new = E

		!write(*,*) imc,com_new1,com_new2,com_new3
	  
	  if(temp_profile_type.eq.0) then		!	square temperature profile

	   	if(simulation.eq.0) call transhomosampling
	   	if(simulation.eq.1) call perp_cage_wind_inhomosampling
           	if(simulation.eq.2) call perp_X_inhomosampling

	   endif

	   if(temp_profile_type.eq.1) then		!	gaussian temperature profile
		
		call gaussian_sampling

	   endif
	
	   !write(*,*) imc,igmol,com(:,igmol)  

!	ROTATION OF ATOMS WRT COM

	   E_old = 0.0d0;E_new=0.0d0;dE=0.0d0;prob=0.0d0;xt=0.0d0

	   translation = 0
	   rotation	= 0
	   call energyith
	   E_old = E

		!write(20,*) 'rot',imc,E_old

	   call rotationalmotion

	   translation = 0
	   rotation	= 1
	   call energyith
	   E_new = E

		!write(20,*) 'rot',imc,E_new

	   call rothomosampling

	! ------------------------------------------- INTRA-ANGULAR MOTION --------------------------------------	!

	   do iangle = 1,nangle(gmoltype(igmol))

		translation=0
		rotation=0
		dh_rotation=0

		iangle1 = index_angle(1,iangle,gmoltype(igmol))
		iangle2 = index_angle(2,iangle,gmoltype(igmol)) 
		iangle3 = index_angle(3,iangle,gmoltype(igmol))

		call intramol_angle
		!write(*,*) theta_angle

!	energy associate with the angle of interest

		Eold_angle = angle_parm(gmoltype(igmol)) * 0.5d0 * (( theta_angle - equi_angle(gmoltype(igmol))) ) &
		* (( theta_angle - equi_angle(gmoltype(igmol))) )
		!write(30,*) theta_angle,Eold_angle

!	energy of gg + gh

		!E_old=0.0d0
		!call energyith
		!E_old = E

!	total old energy

		if(gmoltype(igmol).eq.1) then
			call intraenergy_angle
			E_old =  Eold_angle 
			do jangle = 1,nangle(gmoltype(igmol))
				if(iangle.ne.jangle) then
					iangle1 = index_angle(1,jangle,gmoltype(igmol))
					iangle2 = index_angle(2,jangle,gmoltype(igmol)) 
					iangle3 = index_angle(3,jangle,gmoltype(igmol))
					call intramol_angle
					E_old = E_old + angle_parm(gmoltype(igmol)) * 0.5d0 * (( theta_angle - equi_angle(gmoltype(igmol))) ) &
							* (( theta_angle - equi_angle(gmoltype(igmol))) )
				endif
			enddo !jangle
			!write(30,*) E_old,Eintra_angle,Eold_angle,E_old
		endif
	
		if(gmoltype(igmol).eq.2) then
			!call intraenergy_angle
			E_old = Eold_angle !+ E_old
		endif

!	random movement of intra angle
		!write(*,*) theta_angle
		call random_number(rand)
		ntheta_angle = theta_angle + (2*rand-1)*dangle_max
		dtheta_angle = - ( ntheta_angle - theta_angle )

		!write(30,*) 'dtheta_angle',dtheta_angle

		!write(*,*) 'change in angle',dtheta_angle,theta_angle,ntheta_angle

!	rotation of atoms of hydrocarbon


		rotation = 1

!	for type #1	BRANCHED HYDROCARBON

		if(gmoltype(igmol).eq.1) then

		  iangle1 = index_angle(1,iangle,gmoltype(igmol))
		  iangle2 = index_angle(2,iangle,gmoltype(igmol)) 
		  iangle3 = index_angle(3,iangle,gmoltype(igmol))

			do igatom = 1,ngatom1

			if (igatom.ne.iangle3) rg_natomic(:,igatom) = rg_atomic(:,igatom,igmol)

			enddo !igatom

			call intraangle_rotation

			ratom(1,iangle3) = rg_atomic(1,iangle3,igmol) - rg_atomic(1,iangle2,igmol)
			ratom(2,iangle3) = rg_atomic(2,iangle3,igmol) - rg_atomic(2,iangle2,igmol)
			ratom(3,iangle3) = rg_atomic(3,iangle3,igmol) - rg_atomic(3,iangle2,igmol)

			!write(*,*) ratom(:,3)

	

			rg_natomic(1,iangle3) = rot(1,1)*ratom(1,iangle3) + rot(1,2)*ratom(2,iangle3) + rot(1,3)*ratom(3,iangle3)
 			rg_natomic(1,iangle3) = rg_natomic(1,iangle3) + rg_atomic(1,iangle2,igmol)
				
			rg_natomic(2,iangle3) = rot(2,1)*ratom(1,iangle3) + rot(2,2)*ratom(2,iangle3) + rot(2,3)*ratom(3,iangle3)
			rg_natomic(2,iangle3) = rg_natomic(2,iangle3) + rg_atomic(2,iangle2,igmol)
	
			rg_natomic(3,iangle3) = rot(3,1)*ratom(1,iangle3) + rot(3,2)*ratom(2,iangle3) + rot(3,3)*ratom(3,iangle3)
			rg_natomic(3,iangle3) = rg_natomic(3,iangle3) + rg_atomic(3,iangle2,igmol)

			!write(40,*) '1',rg_natomic(:,iangle3)

			!d = sqrt((rg_natomic(1,iangle3)-rg_natomic(1,iangle2))**2+(rg_natomic(2,iangle3)-rg_natomic(2,iangle2))**2 &
			!+	(rg_natomic(3,iangle3)-rg_natomic(3,iangle2))**2 )
			!if(mod(imc,1000).eq.0) write(20,*) d

		endif



!	for type #2	LINEAR HYDROCARBON

		if(gmoltype(igmol).eq.2) then

		  iangle1 = index_angle(1,iangle,gmoltype(igmol))
		  iangle2 = index_angle(2,iangle,gmoltype(igmol)) 
		  iangle3 = index_angle(3,iangle,gmoltype(igmol))

		  do igatom = 1,ngatom2

		     if (igatom.le.iangle2) then
			rg_natomic(:,igatom) = rg_atomic(:,igatom,igmol)
		     endif

		     if (igatom.ge.iangle3) then

			call intraangle_rotation

			ratom(1,igatom) = rg_atomic(1,igatom,igmol) - rg_atomic(1,iangle2,igmol)
			ratom(2,igatom) = rg_atomic(2,igatom,igmol) - rg_atomic(2,iangle2,igmol)
			ratom(3,igatom) = rg_atomic(3,igatom,igmol) - rg_atomic(3,iangle2,igmol)

			!write(*,*) ratom(:,3)

	

	     		rg_natomic(1,igatom) = rot(1,1)*ratom(1,igatom) + rot(1,2)*ratom(2,igatom) + rot(1,3)*ratom(3,igatom)
 	     		rg_natomic(1,igatom) = rg_natomic(1,igatom) + rg_atomic(1,iangle2,igmol)
				
	     		rg_natomic(2,igatom) = rot(2,1)*ratom(1,igatom) + rot(2,2)*ratom(2,igatom) + rot(2,3)*ratom(3,igatom)
	     		rg_natomic(2,igatom) = rg_natomic(2,igatom) + rg_atomic(2,iangle2,igmol)

	     		rg_natomic(3,igatom) = rot(3,1)*ratom(1,igatom) + rot(3,2)*ratom(2,igatom) + rot(3,3)*ratom(3,igatom)
	     		rg_natomic(3,igatom) = rg_natomic(3,igatom) + rg_atomic(3,iangle2,igmol)

			!write(*,*) rg_natomic(:,igatom)

			!d = sqrt((rg_natomic(1,igatom)-rg_natomic(1,igatom-1))**2+(rg_natomic(2,igatom)-rg_natomic(2,igatom-1))**2 &
			!+	(rg_natomic(3,igatom)-rg_natomic(3,igatom-1))**2 )
			!if(mod(imc,1000).eq.0) write(20,*) d
	     

		    endif

		  enddo !igatom

		endif

		!call intramol_angle

		Enew_angle = 0.0d0
		Enew_angle = angle_parm(gmoltype(igmol)) * 0.5d0 * (( ntheta_angle - equi_angle(gmoltype(igmol))) ) &
		* ( (ntheta_angle - equi_angle(gmoltype(igmol))) )

!	energy of gg + gh

		!E_new=0.0d0
		!call energyith
		!E_new = E

		!write(*,*) Enew_angle

!	total new energy

		if(gmoltype(igmol).eq.1) then
			call intraenergy_angle
			E_new = Enew_angle 
			do jangle = 1,nangle(gmoltype(igmol))
				if(iangle.ne.jangle) then
					iangle1 = index_angle(1,jangle,gmoltype(igmol))
		  			iangle2 = index_angle(2,jangle,gmoltype(igmol)) 
		  			iangle3 = index_angle(3,jangle,gmoltype(igmol))
					call intramol_angle
					E_new = E_new + angle_parm(gmoltype(igmol)) * 0.5d0 * (( ntheta_angle - equi_angle(gmoltype(igmol))) ) &
							* ( (ntheta_angle - equi_angle(gmoltype(igmol))) )
				endif
			enddo !jangle
		endif
	
		if(gmoltype(igmol).eq.2) then
			!call intraenergy_angle
			E_new =  Enew_angle !+ E_new
		endif

		!write(30,*) igmol,E_old,theta_angle,E_new,ntheta_angle,equi_angle(gmoltype(igmol))

		dE = 0.0d0; prob = 0.0d0; xt = 0.0d0
		dE = E_new - E_old

        	RTl = R*temp(int((com(1,igmol))/0.01d0))

		if (dE.gt.0.0d0) then
			prob = exp ( - dE/RTl )
			call random_number(xt)
		endif

		if (dE.lt.0.0d0.or.xt.lt.prob) then
			if(gmoltype(igmol).eq.1) then
				do igatom = 1,ngatom1
					rg_atomic(:,igatom,igmol) = rg_natomic(:,igatom)
				enddo
			endif
			if(gmoltype(igmol).eq.2) then
				do igatom = 1,ngatom2
					rg_atomic(:,igatom,igmol) = rg_natomic(:,igatom)
				enddo
			endif
			!write(30,*) imc,'accp',dE,xt,prob
		else
			!if(gmoltype(igmol).eq.1) rejgg1rot = rejgg1rot   + 1
			!if(gmoltype(igmol).eq.2) rejgg2rot = rejgg2rot   + 1
		endif


	enddo !iangle

!  ---------- DIHEDRAL MOTION ----------	!

	    do idh=1,ndh(gmoltype(igmol))
		if(gmoltype(igmol).eq.2) then
		translation=0
		rotation=0
		dh_rotation=0
		E_old=0.0d0
		call energyith
		E_old = E
		call caldihedralangle
		call dihedenergy
		E_old = E_old + E
		call dihedralmotion
		rotation=1
		dh_rotation=1
		E_new = 0.0d0
		call energyith
		E_new = E
		call caldihedralangle
		call dihedenergy
		E_new = E_new + E
		call dihhomosampling
		endif
	   enddo !idh


	enddo !igmol

	translation = 0
	rotation = 0
	dh_rotation = 0


	do igmol = 1,ngmol

		translation = 0
		rotation = 0

             if(gmoltype(igmol).eq.1) then
                call energyith                  ! calcualte energy of igmol

                call cagetocomdistance          ! calculate perp distance from nearest window to com of igmol

                !ibin1 = nint(pt_distance/0.1d0)                !
                !if you want to bin perp distance from nearest window
                !plane to com of a molecule
                !if(pt_distance.ge.-0.05.and.pt_distance.le.0.05)
                !write(60,*) igmol,pt_distance,ibin1,E

                ibin1 = nint(com(1,igmol)/0.1d0)                ! if you want to bin com position of a molecule

                if(ibin1.gt.0.0d0.and.ibin1.lt.boxlength*10.0d0) count1left = count1left + 1
                if(ibin1.gt.boxlength*(lx-1)*10.0d0.and.ibin1.lt.boxlength*lx*10.0) count1right = count1right + 1

                if(E.le.0) count1(ibin1)=count1(ibin1)+1
                if(E.le.0) pe1(ibin1) = pe1(ibin1) + E
                if(E.le.0) pe1_gg(ibin1) = pe1_gg(ibin1) + guestguest_E
                if(E.le.0) pe1_gh(ibin1) = pe1_gh(ibin1) + guesthost_E

           endif

           if(gmoltype(igmol).eq.2) then
                call energyith                ! calcualte energy of igmol

		call cagetocomdistance
                !ibin2 = nint(pt_distance/0.1d0)        !       if you
                !want to bin perp distance from nearest window plane to
                !com of a molecule
                ibin2 = nint(com(1,igmol)/0.1d0)                ! if you want to bin com position of a molecule

                !if(pt_distance.ge.-0.05.and.pt_distance.le.0.05)
                !write(70,*) igmol,pt_distance,ibin2,E


                if(ibin2.gt.0.0d0.and.ibin2.lt.boxlength*10.0d0) count2left = count2left + 1
                if(ibin2.gt.boxlength*(lx-1)*10.0d0.and.ibin2.lt.boxlength*lx*10.0d0) count2right = count2right + 1

                if(E.le.0)count2(ibin2)=count2(ibin2)+1
                if(E.le.0) pe2(ibin2) = pe2(ibin2) + E
                if(E.le.0) pe2_gg(ibin2) = pe2_gg(ibin2) + guestguest_E
                if(E.le.0) pe2_gh(ibin2) = pe2_gh(ibin2) + guesthost_E


           endif

        enddo !igmol

						!*****************OUTPUT*********************!

        do igmol = 1,ngmol
                if(igmol.le.int(ngmol*0.5d0)) then
                        if(com(1,igmol).le.24.8536d0) leftneo = leftneo + 1
                        if(com(1,igmol).ge.74.5608d0) rightneo = rightneo + 1
                endif
                if(igmol.gt.int(ngmol*0.5d0)) then
                        if(com(1,igmol).le.24.8536d0) leftnp = leftnp + 1
                        if(com(1,igmol).ge.74.5608d0) rightnp = rightnp + 1
                endif
        enddo !igmol

        !       WRITE THE CENTER OF MASS COORDINATES OF GUEST MOLECULES ALONG TRAJECTORY

        if(iwinfo.eq.nwinfo-1) then

        write(9,*) ngmol
        write(9,51) 'id','x','y','z'!,'tot_E','E_gg','E_gh'
51      format(3x,A,3x,A,5x,A,6x,A)!,5x,A,3x,A,3x,A)

        do igmol=1,ngmol

          write(9,50) gmoltype(igmol),com(:,igmol)!,E,guestguest_E,guesthost_E
50      format(1x,i3,1x,f6.2,1x,f6.2,1x,f6.2)!,1x,f8.2,1x,f8.2,1x,f8.2)
        enddo !igmol
        endif

	iwinfo = iwinfo + 1
	if(iwinfo.ge.nwinfo) then
	iwinfo = 0

        alpha = (leftneo / leftnp) / (rightneo / rightnp)
        write(21,*) imc,alpha,log(alpha)/log(10.0d0),leftneo,rightneo,leftnp,rightnp
        leftneo = 0.0d0; rightneo = 0.0d0; leftnp = 0.0d0; rightnp = 0.0d0

!	CALCULATE THE TOTAL ENERGY

	call totalenergy
        call intramolenergy

        !write(*,*) count1left,count1right,count2left,count2right

        !separation_factor = (count1left/count2left) / (count1right/count2right)

 	
	write(10,40) imc,real(E_totalintra),real(E_totalinter)!,count1left,count1right,count2left,count2right
40      format(i7,4x,f20.8,4x,f16.8)!,4x,i8,4x,i8,4x,i8,4x,i8)
	write(11,41) imc,real(E_g11),real(E_g12),real(E_g22),real(E_g1h),real(E_g2h),real(E_totalinter)
41      format(2x,i7,4x,f16.8,4x,f16.8,4x,f16.8,4x,f16.8,4x,f16.8,4x,f16.8)
	write(14,43)imc,real(accphh/(accphh+rejhh)),real(accpgg1trans1/(accpgg1trans1+rejgg1trans1)), &
	real(accpgg1trans2/(accpgg1trans2+rejgg1trans2)),real(accpgg1trans3/(accpgg1trans3+rejgg1trans3)), &
        real(accpgg2trans1/(accpgg2trans1+rejgg2trans1)),real(accpgg2trans2/(accpgg2trans2+rejgg2trans2)), &
        real(accpgg2trans3/(accpgg2trans3+rejgg2trans3)), &
	real(accpgg1rot/(accpgg1rot+rejgg1rot)),real(accpgg2rot/(accpgg2rot+rejgg2rot)),real(accpgg2dih/(accpgg2dih+rejgg2dih))
43      format(2x,i7,2x,f16.8,2x,f16.8,2x,f16.8,2x,f16.8,2x,f16.8,2x,f16.8,2x,f16.8,2x,f16.8,2x,f16.8,2x,f16.8 )

!	CHECK THE MAXIMUM DISPLACEMENT OF HOST ATOMS THROUGH ACCEPTANCE CRITERIA
	
	hostaccpratio = 0.0d0
	hostaccpratio = accphh/(accphh+rejhh)
	!if(hostaccpratio.gt.0.6) dr_maxhh = dr_maxhh * 1.05d0
	!if(hostaccpratio.lt.0.6) dr_maxhh = dr_maxhh * 0.95d0
	accphh = 0.0d0 ; rejhh = 0.0d0
	!write(*,*) dr_maxhh

!	CHECK THE MAXIMUM DISPLACEMENT OF GUEST ATOMS THROUGH ACCEPTANCE CRITERIA

						!	GUEST MOLECULE #1	!
!	TRANSLATION
	
	guest1accpratiotrans = 0.0d0
	guest1accpratiotrans = accpgg1trans/(accpgg1trans+rejgg1trans)
	if(guest1accpratiotrans.gt.0.5) dr_maxgg1trans = dr_maxgg1trans * 1.05d0
	if(guest1accpratiotrans.lt.0.5) dr_maxgg1trans = dr_maxgg1trans * 0.95d0
	accpgg1trans = 0.0d0 ; rejgg1trans = 0.0d0
	!write(*,*) dr_maxgg1trans

!	ROTATION

	guest1accpratiorot = 0.0d0
	guest1accpratiorot = accpgg1rot/(accpgg1rot+rejgg1rot)
	if(guest1accpratiorot.gt.0.5) dr_maxgg1rot = dr_maxgg1rot * 1.05d0
	if(guest1accpratiorot.lt.0.5) dr_maxgg1rot = dr_maxgg1rot * 0.95d0
	accpgg1rot = 0.0d0 ; rejgg1rot = 0.0d0
	!write(*,*) dr_maxgg1rot

!	DIHEDRAL

	guest1accpratiodih = 0.0d0
	guest1accpratiodih = accpgg1dih/(accpgg1dih+rejgg1dih)
	if(guest1accpratiodih.gt.0.5) dr_maxgg1dih = dr_maxgg1dih * 1.05d0
	if(guest1accpratiodih.lt.0.5) dr_maxgg1dih = dr_maxgg1dih * 0.95d0
	accpgg1dih = 0.0d0 ; rejgg1dih = 0.0d0
	!write(*,*) dr_maxgg1dih

						!	GUEST MOLECULE #2	!

!	TRANSLATION

	guest2accpratiotrans = 0.0d0
	guest2accpratiotrans = accpgg2trans/(accpgg2trans+rejgg2trans)
	if(guest2accpratiotrans.gt.0.5) dr_maxgg2trans = dr_maxgg2trans * 1.05d0
	if(guest2accpratiotrans.lt.0.5) dr_maxgg2trans = dr_maxgg2trans * 0.95d0
	accpgg2trans = 0.0d0 ; rejgg2trans = 0.0d0
	!write(*,*) dr_maxgg2trans

!	ROTATION

	guest2accpratiorot = 0.0d0
	guest2accpratiorot = accpgg2rot/(accpgg2rot+rejgg2rot)
	if(guest2accpratiorot.gt.0.5) dr_maxgg2rot = dr_maxgg2rot * 1.05d0
	if(guest2accpratiorot.lt.0.5) dr_maxgg2rot = dr_maxgg2rot * 0.95d0
	accpgg2rot = 0.0d0 ; rejgg2rot = 0.0d0
	!write(*,*) dr_maxgg2rot

!	DIHEDRAL

	guest2accpratiodih = 0.0d0
	guest2accpratiodih = accpgg2dih/(accpgg2dih+rejgg2dih)
	if(guest2accpratiodih.gt.0.8) dr_maxgg2dih = dr_maxgg2dih * 1.05d0
        if(guest2accpratiodih.lt.0.8) dr_maxgg2dih = dr_maxgg2dih * 0.95d0
	accpgg2dih = 0.0d0 ; rejgg2dih = 0.0d0
	!write(*,*) dr_maxgg2dih



!	WRITE THE TRAJECTORY OF HOST & GUEST MOLECULES

					!	HOST		!

	write(8,*) int(ngmol/2)*(ngatom1+ngatom2) + nhatom
	write(8,*)

	do ihatom=1,nhatom
	   if(hatomtype(ihatom).eq.1) write(8,*) 'Si',rh_atomic(:,ihatom)
           if(hatomtype(ihatom).eq.2) write(8,*) 'Al',rh_atomic(:,ihatom)
	   if(hatomtype(ihatom).eq.3) write(8,*) 'O',rh_atomic(:,ihatom)
	   if(hatomtype(ihatom).eq.4) write(8,*) 'Na',rh_atomic(:,ihatom)
	enddo !ihatom

					!	GUEST		!

	do igmol=1,ngmol

	   do igatom=1,ngatom1
	   	if(gmoltype(igmol).eq.1) write(8,*) 'c',com(:,igmol)+rg_atomic(:,igatom,igmol)
	   enddo !iatom

	   do igatom=1,ngatom2
	   	if(gmoltype(igmol).eq.2) write(8,*) 'p',com(:,igmol)+rg_atomic(:,igatom,igmol)
	   enddo 

	enddo !igmol

!	WRITE THE CENTER OF MASS COORDINATES OF GUEST MOLECULES ALONG TRAJECTORY

!	write(9,*) ngmol
!        write(9,*)

!        do igmol=1,ngmol
!           write(9,50) gmoltype(igmol),com(:,igmol)
!50      format(1x,i3,1x,f10.4,1x,f10.4,1x,f10.4)
!        enddo !igmol

        endif

        if(mod(imc,100000).eq.0) then

!	WRITE THE DENSITY PROFILE

	!do i = -3000,3000
	!  write(40+imc,*) i*0.1d0,pe1(i)/count1(i),count1(i),count1(i)*1.0d0/(100000)
	!  write(50+imc,*) i*0.1d0,pe2(i)/count2(i),count2(i),count2(i)*1.0d0/(100000)
	!enddo !i

	endif
	

	enddo !imc

	 do i = -300,2000
          write(40,*) i*0.1d0,pe1(i)/count1(i),pe1_gg(i)/count1(i),pe1_gh(i)/count1(i),count1(i)
          write(50,*) i*0.1d0,pe2(i)/count2(i),pe2_gg(i)/count2(i),pe2_gh(i)/count2(i),count2(i)
        enddo !i

				!----------------------	WRITE THE OUTPUT DATA	--------------------------!

!	WRITE THE DISTRIBUTION AND ENERGY LANDSCAPE 0F GUEST MOLECULES

        !write(10,*) 'distribution of population'

        !write(10,*) '1',count1left,count1right
        !write(10,*) '2',count2left,count2right 

!       WRITE THE FINAL CONFIGURATION OF HOST ATOMS

        write(15,*) nhatom
        write(15,*)

        do ihatom=1,nhatom
           if(hatomtype(ihatom).eq.1) write(15,*) '1',rh_atomic(:,ihatom)
           if(hatomtype(ihatom).eq.2) write(15,*) '2',rh_atomic(:,ihatom)
           if(hatomtype(ihatom).eq.3) write(15,*) '3',rh_atomic(:,ihatom)
        enddo !ihatom

!	WRITE THE FINAL CONFIGURATION OF GUEST CENTER OF MASS COORDINATES

        write(12,*) ngmol
        write(12,*)

        do igmol=1,ngmol
           write(12,*) gmoltype(igmol),gmoltype(igmol),com(:,igmol)
        enddo !igmol

!	WRITE THE FINAL CONFIGURATION OF GUEST ATOMIC COORDINATES WRT CENTER OF MASS

	write(13,*) ngmol*5
	write(13,*)
	
	do igmol=1,ngmol
	   do igatom=1,ngatom
	   if(gmoltype(igmol).eq.1) write(13,*) gatomtype(igatom,igmol),rg_atomic(:,igatom,igmol)
	   enddo !iatom

	   do igatom=1,ngatom
	   if(gmoltype(igmol).eq.2) write(13,*) gatomtype(igatom,igmol),rg_atomic(:,igatom,igmol)
	   enddo 
	enddo !igmol

        call cpu_time(finish)

        write(10,*) 'Total Time Taken in Seconds',finish-start

	end

			!.............. SUBROUTINE TO CREATE NEIGHBOR LIST (GUEST-GUEST + GUEST-HOST)...................!

	subroutine neighborlist
	
	use parameter
	implicit none

	do igmol=1,ngmol
	   nlistgg(igmol)=0
	   nlistgh(igmol) = 0
	   
	   do jgmol=1,ngmol
	
!	avoinding self-interaction
	       if (igmol.ne.jgmol) then

		d=0.0d0;dx=0.0d0;dy=0.0d0;dz=0.0d0

!	calculation of distance between center of mass of molecules
		dx = com(1,jgmol) - com(1,igmol);dy = com(2,jgmol) - com(2,igmol);dz = com(3,jgmol) - com(3,igmol)

!	applying minimum boundary condition

		call minimumdisplacementcondition

		d = sqrt (dx*dx+dy*dy+dz*dz)

		!write(20,*) igmol,jgmol,d

!	applying cut-off sphere radius
	        if (d.le.rc) then
		   nlistgg(igmol) = nlistgg(igmol) + 1
		   listgg(nlistgg(igmol),igmol)= jgmol
	        endif
	      endif
	     enddo !jgmol

!	write(5,*) igmol,nlistgg(igmol),(listgg(j,igmol),j=1,nlistgg(igmol))

	 do ihatom=1,nhatom

			d=0.0d0;dx=0.0d0;dy=0.0d0;dz=0.0d0

			dx = rh_atomic(1,ihatom) - com(1,igmol) !- rg_atomic(1,igatom,igmol)
			dy = rh_atomic(2,ihatom) - com(2,igmol) !- rg_atomic(2,igatom,igmol)
			dz = rh_atomic(3,ihatom) - com(3,igmol) !- rg_atomic(3,igatom,igmol)

!	applying minimum boundary condition

			call minimumdisplacementcondition

			d = sqrt (dx*dx+dy*dy+dz*dz)

!	applying cut-off sphere radius
	        	if (d.le.rc) then
		   		nlistgh(igmol) = nlistgh(igmol) + 1
		   		listgh(nlistgh(igmol),igmol)= ihatom
	        	endif
	      		
	     enddo !ihatom

!	write(5,*) igmol,nlistgh(igmol),(listgh(j,igmol),j=1,nlistgh(igmol))

	enddo !igmol

	end subroutine
			!................. SUBROUTINE TO CALCULATE ENERGY ........................!

	subroutine totalenergy

	use parameter
	implicit none

!	real :: start,finish

!	call cpu_time(start)

	E_totalinter = 0.0d0 ; E_g11 = 0.0d0 ; E_g12 = 0.0d0 ; E_g22 = 0.0d0 ; E_g1h = 0.0d0 ; E_g2h = 0.0d0

!	GUEST-GUEST

	do igmol=1,ngmol

!	   do jgmol=1,ngmol
!		if(igmol.ne.jgmol) then

	  do j=1,nlistgg(igmol)
	     jgmol = listgg(j,igmol)

		E=0.0d0;dx=0.0d0;dy=0.0d0;dz=0.0d0;dxaa=0.0d0;dyaa=0.0d0;dzaa=0.0d0

		dx = com(1,jgmol) - com(1,igmol);dy = com(2,jgmol) - com(2,igmol);dz = com(3,jgmol) - com(3,igmol)

!	applying minimum boundary condition

		call minimumdisplacementcondition

!	calculation of inter & intra molecular energy

		if(gmoltype(igmol).eq.1) nigatom = ngatom1
		if(gmoltype(igmol).eq.2) nigatom = ngatom1

		if(gmoltype(jgmol).eq.1) njgatom = ngatom1
		if(gmoltype(jgmol).eq.2) njgatom = ngatom2

		   do igatom=1,nigatom
		      igatomtype = gatomtype(igatom,igmol)

		      do jgatom=1,njgatom
		      	jgatomtype = gatomtype(jgatom,jgmol)

		      	dxaa =  rg_atomic(1,jgatom,jgmol) + dx - rg_atomic(1,igatom,igmol)
		      	dyaa =  rg_atomic(2,jgatom,jgmol) + dy - rg_atomic(2,igatom,igmol)
		      	dzaa =  rg_atomic(3,jgatom,jgmol) + dz - rg_atomic(3,igatom,igmol)

		      	rsq=dxaa*dxaa+dyaa*dyaa+dzaa*dzaa

		      	ratio = sigmasqgg(igatomtype,jgatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

		      	E = fourepsgg(igatomtype,jgatomtype) * ( s12 - s6 )
		      	E_totalinter = E_totalinter + E

			if(gmoltype(igmol).eq.1.and.gmoltype(jgmol).eq.1) E_g11 = E + E_g11
			if(gmoltype(igmol).eq.1.and.gmoltype(jgmol).eq.2) E_g12 = E + E_g12
			if(gmoltype(igmol).eq.2.and.gmoltype(jgmol).eq.2) E_g12 = E + E_g12
			if(gmoltype(igmol).eq.2.and.gmoltype(jgmol).eq.2) E_g22 = E + E_g22

			!write(20,*) imol,jmol,iatom,iatomtype,jatom,jatomtype
			!write(20,*) rsq,sigmasq(iatomtype,jatomtype),ratio,E_totalinter
			!write(20,*)

		      enddo !jatom

		   enddo !iatom 

!		endif

	     enddo !j

	enddo !imol

	E_totalinter = E_totalinter*0.5d0

!	GUEST-HOST

	do igmol=1,ngmol

	do j=1,nlistgh(igmol)
	     ihatom = listgh(j,igmol)
 
!	   do ihatom=1,nhatom
	   ihatomtype = hatomtype(ihatom)
	
		if(gmoltype(igmol).eq.1) nigatom = ngatom1
		if(gmoltype(igmol).eq.2) nigatom = ngatom2

		do igatom=1,nigatom
		igatomtype = gatomtype(igatom,igmol)

			E=0.0d0;dx=0.0d0;dy=0.0d0;dz=0.0d0

			dx = rh_atomic(1,ihatom) - com(1,igmol) - rg_atomic(1,igatom,igmol)
			dy = rh_atomic(2,ihatom) - com(2,igmol) - rg_atomic(2,igatom,igmol)
			dz = rh_atomic(3,ihatom) - com(3,igmol) - rg_atomic(3,igatom,igmol)

!	applying minimum boundary condition

			call minimumdisplacementcondition

			rsq=dx*dx+dy*dy+dz*dz

		      	ratio = sigmasqgh(igatomtype,ihatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

		      	E = fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 )
		      	E_totalinter = E_totalinter + E

			if(gmoltype(igmol).eq.1) E_g1h = E + E_g1h
			if(gmoltype(igmol).eq.2) E_g2h = E + E_g2h

		enddo !igatom

	    enddo !ihatom

	enddo !igmol

	dh_rotation = 0

	do igmol=1,ngmol

			do idh=1,ndh(gmoltype(igmol))

			call caldihedralangle
	
			call dihedenergy
	
			E_totalinter = E_totalinter + E
			!write(*,*) E

			!E_g22 = E + E_g22

			enddo !idh

	enddo !igmol

!	call cpu_time(finish)

!	write(30,*) 'all',finish - start
		

	end subroutine

				!................. SUBROUTINE TO ENERGY OF iTH MOLECULES........................!

	subroutine energyith

	use parameter
	implicit none

!	real :: start,finish

!	call cpu_time(start)

!	if(translation.eq.0) com_p(:)=com(:,igmol)		! p stands for present
!	if(translation.eq.1) com_p(:)=com_new(:)

	if(translation.eq.0) then	
		com_p1 = com(1,igmol) ; com_p2 = com(2,igmol) ; com_p3 = com(3,igmol)		! p stands for present
	endif
	if(translation.eq.1) then
		com_p1 = com_new1 ; com_p2 = com_new2 ; com_p3 = com_new3
	endif

	if(rotation.eq.0) atomic_p(:,:)=rg_atomic(:,:,igmol)
	if(rotation.eq.1) atomic_p(:,:)=rg_natomic(:,:)

	E = 0.0d0; guestguest_E = 0.0d0; guesthost_E = 0.0d0

	   do j=1,nlistgg(igmol)
		jgmol=listgg(j,igmol)

!	    do jgmol=1,ngmol
!		if(igmol.ne.jgmol) then

		dx=0.0d0;dy=0.0d0;dz=0.0d0;dxaa=0.0d0;dyaa=0.0d0;dzaa=0.0d0

		dx = com(1,jgmol) - com_p1;dy = com(2,jgmol) - com_p2;dz = com(3,jgmol) - com_p3

!	applying minimum boundary condition

		call minimumdisplacementcondition

!	calculation of inter & intra molecular energy

		if(gmoltype(igmol).eq.1) nigatom = ngatom1
		if(gmoltype(igmol).eq.2) nigatom = ngatom1

		if(gmoltype(jgmol).eq.1) njgatom = ngatom1
		if(gmoltype(jgmol).eq.2) njgatom = ngatom2
		
		   do igatom=1,nigatom
		      igatomtype = gatomtype(igatom,igmol)

		      do jgatom=1,njgatom
		      	jgatomtype = gatomtype(jgatom,jgmol)

		      	dxaa =  rg_atomic(1,jgatom,jgmol) + dx - atomic_p(1,igatom)
		      	dyaa =  rg_atomic(2,jgatom,jgmol) + dy - atomic_p(2,igatom)
		      	dzaa =  rg_atomic(3,jgatom,jgmol) + dz - atomic_p(3,igatom)

		      	rsq=dxaa*dxaa+dyaa*dyaa+dzaa*dzaa

		      	ratio = sigmasqgg(igatomtype,jgatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

			guestguest_E = guestguest_E + fourepsgg(igatomtype,jgatomtype) * ( s12 - s6 )

		      	E = E + fourepsgg(igatomtype,jgatomtype) * ( s12 - s6 )

		      enddo !jatom

		   enddo !iatom

!		endif 

	     enddo !j

!	do ihatom=1,nhatom
!	   ihatomtype = hatomtype(ihatom)

	do j=1,nlistgh(igmol)
	     ihatom = listgh(j,igmol)
	   ihatomtype = hatomtype(ihatom)

		if(gmoltype(igmol).eq.1) nigatom = ngatom1
		if(gmoltype(igmol).eq.2) nigatom = ngatom2

		do igatom=1,nigatom
		igatomtype=gatomtype(igatom,igmol)

			dx=0.0d0;dy=0.0d0;dz=0.0d0

			dx = rh_atomic(1,ihatom) - com_p1 - atomic_p(1,igatom)
			dy = rh_atomic(2,ihatom) - com_p2 - atomic_p(2,igatom)
			dz = rh_atomic(3,ihatom) - com_p3 - atomic_p(3,igatom)

!	applying minimum boundary condition

			call minimumdisplacementcondition

			rsq=dx*dx+dy*dy+dz*dz

		      	ratio = sigmasqgh(igatomtype,ihatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

			guesthost_E = guesthost_E + fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 )

		      	E = E + fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 )

		enddo !igatom

	    enddo !ihatom

!	call cpu_time(finish)

!	write(30,*) 'ith',finish-start

	end subroutine

				!................. SUBROUTINE TO CALCULATE BOUNDARY ENERGY ........................!

	subroutine boundaryenergy

	use parameter
	implicit none

!	real :: start,finish

!	call cpu_time(start)

!	if(translation.eq.0) com_p(:)=com(:,igmol)		! p stands for present
!	if(translation.eq.1) com_p(:)=com_new(:)

	if(translation.eq.0) then	
		com_p1 = com(1,igmol) ; com_p2 = com(2,igmol) ; com_p3 = com(3,igmol)		! p stands for present
	endif
	if(translation.eq.1) then
		com_p1 = com_new1 ; com_p2 = com_new2 ; com_p3 = com_new3
	endif

	if(rotation.eq.0) atomic_p(:,:)=rg_atomic(:,:,igmol)
	if(rotation.eq.1) atomic_p(:,:)=rg_natomic(:,:)

	E = 0.0d0; guestguest_E = 0.0d0; guesthost_E = 0.0d0

	   do j=1,nlistgg(igmol)
		jgmol=listgg(j,igmol)

!	    do jgmol=1,ngmol
!		if(igmol.ne.jgmol) then

		dx=0.0d0;dy=0.0d0;dz=0.0d0;dxaa=0.0d0;dyaa=0.0d0;dzaa=0.0d0

		dx = com(1,jgmol) - com_p1;dy = com(2,jgmol) - com_p2;dz = com(3,jgmol) - com_p3

!	applying minimum boundary condition

		call minimumdisplacementcondition

!	calculation of inter & intra molecular energy
		if(gmoltype(igmol).eq.1) nigatom = ngatom1
		if(gmoltype(igmol).eq.2) nigatom = ngatom1

		if(gmoltype(jgmol).eq.1) njgatom = ngatom1
		if(gmoltype(jgmol).eq.2) njgatom = ngatom2

		   do igatom=1,nigatom
		      igatomtype = gatomtype(igatom,igmol)

		      do jgatom=1,njgatom
		      	jgatomtype = gatomtype(jgatom,jgmol)

		      	dxaa =  rg_atomic(1,jgatom,jgmol) + dx - atomic_p(1,igatom)
		      	dyaa =  rg_atomic(2,jgatom,jgmol) + dy - atomic_p(2,igatom)
		      	dzaa =  rg_atomic(3,jgatom,jgmol) + dz - atomic_p(3,igatom)

		      	rsq=dxaa*dxaa+dyaa*dyaa+dzaa*dzaa

		      	ratio = sigmasqgg(igatomtype,jgatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

			guestguest_E = guestguest_E + fourepsgg(igatomtype,jgatomtype) * ( s12 - s6 )

		      	E = E + fourepsgg(igatomtype,jgatomtype) * ( s12 - s6 )

		      enddo !jatom

		   enddo !iatom

!		endif 

	     enddo !j

!	do ihatom=1,nhatom
!	   ihatomtype = hatomtype(ihatom)

	do j=1,nlistgh(igmol)
	     ihatom = listgh(j,igmol)
	   ihatomtype = hatomtype(ihatom)

		if(gmoltype(igmol).eq.1) nigatom = ngatom1
		if(gmoltype(igmol).eq.2) nigatom = ngatom2

		do igatom=1,nigatom
		igatomtype=gatomtype(igatom,igmol)

			dx=0.0d0;dy=0.0d0;dz=0.0d0

			dx = rh_atomic(1,ihatom) - com_p1 - atomic_p(1,igatom)
			dy = rh_atomic(2,ihatom) - com_p2 - atomic_p(2,igatom)
			dz = rh_atomic(3,ihatom) - com_p3 - atomic_p(3,igatom)

!	applying minimum boundary condition

			call minimumdisplacementcondition

			rsq=dx*dx+dy*dy+dz*dz

		      	ratio = sigmasqgh(igatomtype,ihatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

			guesthost_E = guesthost_E + fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 )

		      	E = E + fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 )

		enddo !igatom

	    enddo !ihatom

!	call cpu_time(finish)

!	write(30,*) 'ith',finish-start

	end subroutine







		!.............	TRANSLATE THE CENTER OF MASS OF MOLECULES ..............!

	subroutine translationalmotion
	use parameter
	implicit none

	if(gmoltype(igmol).eq.1) dr_max = dr_maxgg1trans
	if(gmoltype(igmol).eq.2) dr_max = dr_maxgg2trans

	call random_number(xrand)
	com_new1 = com(1,igmol) + (2*xrand-1)*dr_max	! random displacement of centre of mass coordinates

	if(simulation.eq.0) then
        	if(com_new1.gt.boxlengthx)     com_new1 = com_new1 - boxlengthx
		if(com_new1.lt.0.0d0)          com_new1 = com_new1 + boxlengthx
	endif

	if(simulation.eq.1.or.simulation.eq.2) then

		do while (com_new1.gt.boxlengthx)
            	call random_number(xrand)
        	com_new1=com(1,igmol)+(2*xrand-1)*dr_max
        	end do ! while

        	do while (com_new1.lt.0.0d0)
        	call random_number(xrand)
        	com_new1=com(1,igmol)+(2*xrand-1)*dr_max
        	end do ! while

	endif
	      
	call random_number(yrand) 
	com_new2 = com(2,igmol) + (2*yrand-1)*dr_max	! along x,y and z respectively
	if(com_new2.gt.boxlengthy)     com_new2 = com_new2 - boxlengthy
	if(com_new2.lt.0.0d0)         com_new2 = com_new2 + boxlengthy

	!do while (com_new2.gt.boxlengthy)
        !call random_number(yrand)
        !com_new2=com(2,igmol)+(2*yrand-1)*dr_max
        !end do ! while

        !do while (com_new2.lt.0.0d0)
        !call random_number(yrand)
        !com_new2=com(2,igmol)+(2*yrand-1)*dr_max
        !end do ! while
	     
	call random_number(zrand)
	com_new3 = com(3,igmol) + (2*zrand-1)*dr_max
	if(com_new3.gt.boxlengthz)     com_new3 = com_new3 - boxlengthz
	if(com_new3.lt.0.0d0)         com_new3 = com_new3 + boxlengthz

	!do while (com_new3.gt.boxlengthz)
        !call random_number(zrand)
        !com_new3=com(3,igmol)+(2*zrand-1)*dr_max
        !end do ! while

        !do while (com_new3.lt.0.0d0)
        !call random_number(zrand)
        !com_new3=com(3,igmol)+(2*zrand-1)*dr_max
        !end do ! while

	end subroutine

		!---------------------------- SUBROUTINE TO ROTATE THE ATOMS -----------------------------!

	subroutine rotationalmotion
	use parameter
	implicit none

	if(gmoltype(igmol).eq.1) dpsi_max = dr_maxgg1rot
	if(gmoltype(igmol).eq.2) dpsi_max = dr_maxgg2rot

	if(gmoltype(igmol).eq.1) ngatom = ngatom1
	if(gmoltype(igmol).eq.2) ngatom = ngatom2
 
	two_pie = 2.0d0 * 3.141592654d0

!	rotational
	
	rand=0.0d0
	call random_number(rand)
	phi 	   = rand * two_pie
	cos_phi = dcos( phi )
	sin_phi = dsin( phi )

	rand=0.0d0
	call random_number(rand)
	cos_theta  = 2.0d0 * rand - 1.0d0
	sin_theta  = sqrt(1.0d0 - cos_theta * cos_theta)

	rand=0.0d0
	call random_number(rand)
	psi 	   = (2.0d0 * rand - 1.0d0) * dpsi_max
	cos_psi    = dcos( psi )
	sin_psi    = dsin( psi ) 
	onemcospsi = 1 - cos_psi

	nx = sin_theta * cos_phi
	ny = sin_theta * sin_phi
	nz = cos_theta

	nxy  = nx * ny; nyz  = ny * nz; nxz  = nx * nz   
	nxsq = nx * nx; nysq = ny * ny; nzsq = nz * nz

!	Rotational Matrix
! 	rot(i,j) = rotational matrix element at ith row & jth column
 
	rot(1,1) = cos_psi + nxsq * onemcospsi
	rot(1,2) = nxy * onemcospsi + nz * sin_psi
	rot(1,3) = nxz * onemcospsi - ny * sin_psi

	rot(2,1) = nxy * onemcospsi - nz * sin_psi
	rot(2,2) = cos_psi + nysq*onemcospsi
	rot(2,3) = nyz * onemcospsi + nx * sin_psi

	rot(3,1) = nxz * onemcospsi + ny * sin_psi
	rot(3,2) = nyz * onemcospsi - nx * sin_psi
	rot(3,3) = cos_psi + nzsq*onemcospsi

!	calculation of determinant

!	detrot=rot(1,1)*(rot(2,2)*rot(3,3)-rot(2,3)*rot(3,2))
!	detrot=detrot+rot(1,2)*(rot(3,1)*rot(2,3)-rot(2,1)*rot(3,3))
!	detrot=detrot+rot(1,3)*(rot(2,1)*rot(3,2)-rot(3,1)*rot(2,2))

	!write(*,*) detrot

!	if(detrot.lt.0.9999.or.detrot.gt.1.0001) then
!		write(10,*) 'deh','determinant is not 1',detrot
!		!stop
!	endif

	do igatom = 1,ngatom

	  rg_natomic(1,igatom)=rot(1,1)*rg_atomic(1,igatom,igmol)+rot(1,2)*rg_atomic(2,igatom,igmol)&
				+rot(1,3)*rg_atomic(3,igatom,igmol)
 						
	  rg_natomic(2,igatom)=rot(2,1)*rg_atomic(1,igatom,igmol) + rot(2,2)*rg_atomic(2,igatom,igmol) &
				+ rot(2,3)*rg_atomic(3,igatom,igmol)
	 	
	  rg_natomic(3,igatom)=rot(3,1)*rg_atomic(1,igatom,igmol) + rot(3,2)*rg_atomic(2,igatom,igmol) &
				+ rot(3,3)*rg_atomic(3,igatom,igmol)
					
	enddo !igatom

	end subroutine

		!------------ SUBROUTINE TO SAMPLE ISOTHERMAL MOTION (IMPORTANCE SAMPLING) ----------------!

	subroutine transhomosampling
	use parameter
	implicit none

	dE = 0.0d0; prob = 0.0d0; xt = 0.0d0
	dE = E_new - E_old

	if (dE.gt.0.0d0) then
		prob = exp ( - dE/RTl )
		call random_number(xt)
	endif

	if (dE.lt.0.0d0.or.xt.lt.prob) then
		com(1,igmol) = com_new1 ;  com(2,igmol) = com_new2 ;  com(3,igmol) = com_new3
		if(gmoltype(igmol).eq.1) accpgg1trans = accpgg1trans + 1
		if(gmoltype(igmol).eq.2) accpgg2trans = accpgg2trans + 1
	else
		if(gmoltype(igmol).eq.1) rejgg1trans = rejgg1trans   + 1
		if(gmoltype(igmol).eq.2) rejgg2trans = rejgg2trans   + 1
	endif

	end subroutine

	subroutine rothomosampling
	use parameter
	implicit none

!        if(E_new.gt.0.0d0) goto 80

	dE = 0.0d0; prob = 0.0d0; xt = 0.0d0
	dE = E_new - E_old

        RTl = R*temp(int(com(1,igmol)/0.01d0))


	if (dE.gt.0.0d0) then
		prob = exp ( - dE/RTl )
		call random_number(xt)
	endif

	if (dE.lt.0.0d0.or.xt.lt.prob) then
		if(gmoltype(igmol).eq.1) then
				do igatom = 1,ngatom1
					rg_atomic(:,igatom,igmol) = rg_natomic(:,igatom)
				enddo
			endif
			if(gmoltype(igmol).eq.2) then
				do igatom = 1,ngatom2
					rg_atomic(:,igatom,igmol) = rg_natomic(:,igatom)
				enddo
			endif
		if(gmoltype(igmol).eq.1) accpgg1rot = accpgg1rot + 1
		if(gmoltype(igmol).eq.2) accpgg2rot = accpgg2rot + 1
	else
		if(gmoltype(igmol).eq.1) rejgg1rot = rejgg1rot   + 1
		if(gmoltype(igmol).eq.2) rejgg2rot = rejgg2rot   + 1
	endif

	end subroutine

	subroutine dihhomosampling
	use parameter
	implicit none

	dE = 0.0d0; prob = 0.0d0; xt = 0.0d0
	dE = E_new - E_old

        RTl = R * temp(int(com(1,igmol)/0.01d0))

	if (dE.gt.0.0d0) then
		prob = exp ( - dE/RTl )
		call random_number(xt)
	endif

	if (dE.lt.0.0d0.or.xt.lt.prob) then
		if(gmoltype(igmol).eq.1) then
				do igatom = 1,ngatom1
					rg_atomic(:,igatom,igmol) = rg_natomic(:,igatom)
				enddo
			endif
			if(gmoltype(igmol).eq.2) then
				do igatom = 1,ngatom2
					rg_atomic(:,igatom,igmol) = rg_natomic(:,igatom)
				enddo
			endif
		if(gmoltype(igmol).eq.1) accpgg1dih = accpgg1dih + 1
		if(gmoltype(igmol).eq.2) accpgg2dih = accpgg2dih + 1
	else
		if(gmoltype(igmol).eq.1) rejgg1dih = rejgg1dih   + 1
		if(gmoltype(igmol).eq.2) rejgg2dih = rejgg2dih   + 1
	endif

	end subroutine

	

!	if (dE.lt.0.0d0.or.xt.lt.prob) then
!		if (translation.eq.1) then
!		    com(1,igmol) = com_new1 ;  com(2,igmol) = com_new2 ;  com(3,igmol) = com_new3
!		    goto 50
!		endif
!		if (rotation.eq.1) then
!		    rg_atomic(:,:,igmol) = rg_natomic(:,:)
!		    goto 50
!		endif
!50		E_present = E_present + dE
!		accp = accp + 1
!	else
!		rej = rej + 1
!	endif



		!................. SUBROUTINE TO PERFORM MINIMUM DISPLACEMENT CONDITION .............!

	subroutine minimumdisplacementcondition

	use parameter
	implicit none


	if (dx.gt.halfboxlengthx)   dx = dx - boxlengthx
	if ( dx.lt.-halfboxlengthx) dx = dx + boxlengthx

	if (dy.gt.halfboxlengthy)   dy = dy - boxlengthy
	if ( dy.lt.-halfboxlengthy) dy = dy + boxlengthy

	if (dz.gt.halfboxlengthz)   dz = dz - boxlengthz
	if ( dz.lt.-halfboxlengthz) dz = dz + boxlengthz


	end subroutine

		!--------------------------- SUBROUTINE TO CALCULATE DIHEDRAL ENERGY ---------------------------!

	subroutine dihedenergy
	use parameter
	implicit none

	E=0.0d0

	!write(*,*) dhangle(idh),E
	
	E = c0(dhtype(igmol))
	E = E + c1(dhtype(igmol)) * ( 1 + dcos(dhangle(idh)))
	E = E + c2(dhtype(igmol)) * ( 1 - dcos(2.0d0*dhangle(idh)))
	E = E + c3(dhtype(igmol)) * ( 1 + dcos(3.0d0*dhangle(idh)))

	!write(13,*) E,c0,c1,c2,c3

	end subroutine

		!------------------------------ SUBROUTINE TO CALCULATE DIHEDRAL ANGLE--------------------!

	subroutine caldihedralangle
	use parameter
	implicit none
	

	if(dh_rotation.eq.0) atomic_p(:,:)=rg_atomic(:,:,igmol)
	if(dh_rotation.eq.1) atomic_p(:,:)=rg_natomic(:,:)

	igatom1 = idhatom(1,idh,gmoltype(igmol)); igatom2 = idhatom(2,idh,gmoltype(igmol))
	igatom3 = idhatom(3,idh,gmoltype(igmol)); igatom4 = idhatom(4,idh,gmoltype(igmol))

	vx(1,2)=atomic_p(1,igatom1)-atomic_p(1,igatom2)
	vy(1,2)=atomic_p(2,igatom1)-atomic_p(2,igatom2)
	vz(1,2)=atomic_p(3,igatom1)-atomic_p(3,igatom2)

	vx(2,3)=atomic_p(1,igatom2)-atomic_p(1,igatom3)
	vy(2,3)=atomic_p(2,igatom2)-atomic_p(2,igatom3)
	vz(2,3)=atomic_p(3,igatom2)-atomic_p(3,igatom3)

	vx(3,4)=atomic_p(1,igatom3)-atomic_p(1,igatom4)
	vy(3,4)=atomic_p(2,igatom3)-atomic_p(2,igatom4)
	vz(3,4)=atomic_p(3,igatom3)-atomic_p(3,igatom4)	

	 !do idh=1,ndh
	
!	cross-product of unit vectors n1 and n2
	
	 xn1 = (vy(1,2)*vz(2,3)-vy(2,3)*vz(1,2))
	 yn1 = (vx(1,2)*vz(2,3)-vx(2,3)*vz(1,2))
	 zn1 = (vx(1,2)*vy(2,3)-vx(2,3)*vy(1,2))
	 rn1 = sqrt(xn1*xn1+yn1*yn1+zn1*zn1)

	 xn2 = (vy(2,3)*vz(3,4)-vy(3,4)*vz(2,3))
	 yn2 = (vx(2,3)*vz(3,4)-vx(3,4)*vz(2,3))
	 zn2 = (vx(2,3)*vy(3,4)-vx(3,4)*vy(2,3))
	 rn2 = sqrt(xn2*xn2+yn2*yn2+zn2*zn2)

!	dot-product of the unit vector to calculate the dihedral angle

	 dotneu = xn1*xn2+yn1*yn2+zn1*zn2
	 dotden = rn1*rn2

	 dotratio = real(dotneu/dotden)

	 dhangle(idh) = dacos(dotratio)

	!write(*,*) dhangle(idh)

	end subroutine

	!---------------------------- SUBROUTINE TO ROTATE THE DIHEDRAL ANGLES -----------------------------!

	subroutine dihedralmotion
	use parameter
	implicit none

40	if(gmoltype(igmol).eq.1) dh_max = dr_maxgg1dih
	if(gmoltype(igmol).eq.2) dh_max = dr_maxgg2dih

	igatom1 = idhatom(1,idh,gmoltype(igmol)); igatom2 = idhatom(2,idh,gmoltype(igmol))
	igatom3 = idhatom(3,idh,gmoltype(igmol)); igatom4 = idhatom(4,idh,gmoltype(igmol))


!	calculation of unit vector

	nx = (rg_atomic(1,igatom3,igmol)-rg_atomic(1,igatom2,igmol))/bondlength(igatom3,igatom2,igmol)
	ny = (rg_atomic(2,igatom3,igmol)-rg_atomic(2,igatom2,igmol))/bondlength(igatom3,igatom2,igmol)
 	nz = (rg_atomic(3,igatom3,igmol)-rg_atomic(3,igatom2,igmol))/bondlength(igatom3,igatom2,igmol)


	nxsq=nx*nx;nysq=ny*ny;nzsq=nz*nz;nxy=nx*ny;nxz=nx*nz;nyz=ny*nz
	

!	random movement of dihedral angle

	call random_number(rand)
	ndhangle = (2*rand-1)*dh_max

!	calculation of matrix 

	cos_dh = dcos(ndhangle)
	sin_dh = dsin(ndhangle)
	onemcosdh = 1 - cos_dh

	rot(1,1) = cos_dh + nxsq * onemcosdh
	rot(1,2) = nxy * onemcosdh + nz * sin_dh
	rot(1,3) = nxz * onemcosdh - ny * sin_dh

	rot(2,1) = nxy * onemcosdh - nz * sin_dh
	rot(2,2) = cos_dh + nysq*onemcosdh
	rot(2,3) = nyz * onemcosdh + nx * sin_dh

	rot(3,1) = nxz * onemcosdh + ny * sin_dh
	rot(3,2) = nyz * onemcosdh - nx * sin_dh
	rot(3,3) = cos_dh + nzsq*onemcosdh

	!write(20,*) cos_dh
	!write(20,*) sin_dh
	!write(20,*) onemcosdh
	!write(20,*) idh,nx,ny,nz,ndhangle,sqrt(nxsq+nysq+nzsq)
	!write(20,*) idh,rot(1,:)
	!write(20,*) idh,rot(2,:)
	!write(20,*) idh,rot(3,:)

!	calculation of determinant

!	detrot=rot(1,1)*(rot(2,2)*rot(3,3)-rot(2,3)*rot(3,2))
!	detrot=detrot+rot(1,2)*(rot(3,1)*rot(2,3)-rot(2,1)*rot(3,3))
!	detrot=detrot+rot(1,3)*(rot(2,1)*rot(3,2)-rot(3,1)*rot(2,2))

	!write(*,*) detrot

!	if(detrot.lt.0.999999.or.detrot.gt.1.000001) then
!		write(10,*) 'deh','determinant is not 1',detrot
		!stop
!	endif

!	apply rotation matrix to rotate the atomic coordinates

!			------------------------------- LINEAR HYDROCARBON --------------------------		!

	if(dhtype(igmol).eq.2) then

	if(gmoltype(igmol).eq.1) ngatom = ngatom1
	if(gmoltype(igmol).eq.2) ngatom = ngatom2

	do igatom=1,idh+2
	     
	     rg_natomic(:,igatom) = rg_atomic(:,igatom,igmol)
	
	enddo !iatom

	do igatom=idh+3,ngatom

	ratom(1,igatom) = rg_atomic(1,igatom,igmol) - rg_atomic(1,idh+2,igmol)
	ratom(2,igatom) = rg_atomic(2,igatom,igmol) - rg_atomic(2,idh+2,igmol)
	ratom(3,igatom) = rg_atomic(3,igatom,igmol) - rg_atomic(3,idh+2,igmol)

	enddo

	do igatom = idh+3,ngatom

	     rg_natomic(1,igatom) = rot(1,1)*ratom(1,igatom) + rot(1,2)*ratom(2,igatom) + rot(1,3)*ratom(3,igatom)
 	     rg_natomic(1,igatom) = rg_natomic(1,igatom) + rg_atomic(1,idh+2,igmol)
				
	     rg_natomic(2,igatom) = rot(2,1)*ratom(1,igatom) + rot(2,2)*ratom(2,igatom) + rot(2,3)*ratom(3,igatom)
	     rg_natomic(2,igatom) = rg_natomic(2,igatom) + rg_atomic(2,idh+2,igmol)

	     rg_natomic(3,igatom) = rot(3,1)*ratom(1,igatom) + rot(3,2)*ratom(2,igatom) + rot(3,3)*ratom(3,igatom)
	     rg_natomic(3,igatom) = rg_natomic(3,igatom) + rg_atomic(3,idh+2,igmol)
	     
	enddo !iatom

	endif

!			----------------------------- BRANCHED HYDROCARBON ------------------------		!

	if(dhtype(igmol).eq.1) then

	if(gmoltype(igmol).eq.1) ngatom = ngatom1
	if(gmoltype(igmol).eq.2) ngatom = ngatom2

	do igatom=1,ngatom
	     
	     if (igatom.ne.igatom4) rg_natomic(:,igatom) = rg_atomic(:,igatom,igmol)
		!write(30,*) rg_natomic(:,igatom),rg_atomic(:,igatom,igmol)
	
	enddo !iatom

	ratom(1,igatom4) = rg_atomic(1,igatom4,igmol) - rg_atomic(1,igatom3,igmol)
	ratom(2,igatom4) = rg_atomic(2,igatom4,igmol) - rg_atomic(2,igatom3,igmol)
	ratom(3,igatom4) = rg_atomic(3,igatom4,igmol) - rg_atomic(3,igatom3,igmol)


	rg_natomic(1,igatom4) = rot(1,1)*ratom(1,igatom4) + rot(1,2)*ratom(2,igatom4) + rot(1,3)*ratom(3,igatom4)
 	rg_natomic(1,igatom4) = rg_natomic(1,igatom4) + rg_atomic(1,igatom3,igmol)
				
	rg_natomic(2,igatom4) = rot(2,1)*ratom(1,igatom4) + rot(2,2)*ratom(2,igatom4) + rot(2,3)*ratom(3,igatom4)
	rg_natomic(2,igatom4) = rg_natomic(2,igatom4) + rg_atomic(2,igatom3,igmol)

	rg_natomic(3,igatom4) = rot(3,1)*ratom(1,igatom4) + rot(3,2)*ratom(2,igatom4) + rot(3,3)*ratom(3,igatom4)
	rg_natomic(3,igatom4) = rg_natomic(3,igatom4) + rg_atomic(3,igatom3,igmol)

	do igatom = 1,ngatom

	if (igatom.ne.igatom1.and.igatom.ne.igatom2.and.igatom.ne.igatom3.and.igatom.ne.igatom4) then

	dx = rg_natomic(1,igatom4) - rg_natomic(1,igatom)
	dy = rg_natomic(2,igatom4) - rg_natomic(2,igatom)
	dz = rg_natomic(3,igatom4) - rg_natomic(3,igatom)
	d = sqrt ( dx*dx + dy*dy + dz*dz )
	!write(*,*) igatom4,igatom,rg_natomic(:,igatom4),rg_natomic(:,igatom),d
	if (d.lt.2.4999745d0.or.d.gt.2.5200255d0) goto 40 

	endif

	enddo !igatom

	endif

!	do igatom=1,ngatom-1
!		dx=rg_natomic(1,igatom+1)-rg_natomic(1,igatom)
!		dy=rg_natomic(2,igatom+1)-rg_natomic(2,igatom)
!		dz=rg_natomic(3,igatom+1)-rg_natomic(3,igatom)
!		d=dx*dx+dy*dy+dz*dz
!		!write(*,*) iatom,iatom+1,sqrt(d)
!		if(sqrt(d).lt.1.535.or.sqrt(d).gt.1.545) then
!
!			write(10,*) 'bond',d,detrot
!			!stop
!		endif
!	enddo !iatom

	end subroutine	


		!----------------------------SUBROUTINE TO CALCULATE COM TO WINDOW DISTANCE-----------------------!


	subroutine cagetocomdistance
	use parameter
	implicit none

	real*8 :: com1,com2,com3

	integer :: nmolx,nmoly,nmolz		! nmolx,nmoly,nmolz = cell position of x,y,z-coordinate of com of molecule

	halfboxlength = boxlength * 0.5d0

	if(translation.eq.0) then	
		com_p1 = com(1,igmol) ; com_p2 = com(2,igmol) ; com_p3 = com(3,igmol)		! p stands for present
	endif
	if(translation.eq.1) then
		com_p1 = com_new1 ; com_p2 = com_new2 ; com_p3 = com_new3
	endif

!	write(20,*) com_p1,com_p2,com_p3

	nmolx = int(com_p1 / boxlength); nmoly = int(com_p2 / boxlength) ; nmolz = int(com_p3 / boxlength )

!	write(20,*) nmolx,nmoly,nmolz,boxlength

	com_p1 = com_p1 - boxlength * nmolx ; com_p2 = com_p2 - boxlength * nmoly ; com_p3 = com_p3 - boxlength * nmolz

!	write(20,*) com_p1,com_p2,com_p3

	if(com_p1.gt.boxlength.or.com_p2.gt.boxlength.or.com_p3.gt.boxlength) then

		write(30,*) com_p1,com_p2,com_p3
		stop

	endif

!	write(20,*) com(:,igmol)

!	CALCULATE THE DISTANCE OF igmol FROM EVERY CAGE CENTER

	do icagecenter = 1,ncagecenter

	   dx=com_p1-rcage(1,icagecenter); dy=com_p2-rcage(2,icagecenter); dz=com_p3-rcage(3,icagecenter)

	   if (dx.gt.halfboxlength)   dx = dx - boxlength
	   if ( dx.lt.-halfboxlength) dx = dx + boxlength

	   if (dy.gt.halfboxlength)   dy = dy - boxlength
	   if ( dy.lt.-halfboxlength) dy = dy + boxlength

	   if (dz.gt.halfboxlength)   dz = dz - boxlength
	   if ( dz.lt.-halfboxlength) dz = dz + boxlength


	   cagetocomdist(icagecenter) = sqrt ( dx*dx + dy*dy + dz*dz )

	enddo !icagecenter



!	FIND OUT WHICH CAGE CENTER IS CLOSEST

	do icagecenter = 1,ncagecenter

	   do jcagecenter = 1,ncagecenter

		if(icagecenter.ne.jcagecenter) then

			if( cagetocomdist(icagecenter).gt.cagetocomdist(jcagecenter) ) goto 30
			
		endif
	   enddo !jcagecenter

	   cage_com_pos = icagecenter
	   goto 50

30	enddo !icagecenter

!50	write(20,*) cage_com_pos

50	dx=com_p1-rcage(1,icagecenter); dy=com_p2-rcage(2,icagecenter); dz=com_p3-rcage(3,icagecenter)

	   if (dx.gt.halfboxlength)   dx = dx - boxlength
	   if ( dx.lt.-halfboxlength) dx = dx + boxlength

	   if (dy.gt.halfboxlength)   dy = dy - boxlength
	   if ( dy.lt.-halfboxlength) dy = dy + boxlength

	   if (dz.gt.halfboxlength)   dz = dz - boxlength
	   if ( dz.lt.-halfboxlength) dz = dz + boxlength

	com1 = rcage(1,icagecenter) + dx ; com2 = rcage(2,icagecenter) + dy ; com3 = rcage(3,icagecenter) + dz

!	write(20,*) (cagetocomdist(i),i=1,8)
!	write(20,*) com1,com2,com3,icagecenter

!	write(20,*) (listcagewindow(j,icagecenter),j=1,nlistcagewind)

	do ilistcagewind = 1,nlistcagewind

	d = com1*aw(ilistcagewind,icagecenter)+com2*bw(ilistcagewind,icagecenter) & 
		+ com3*cw(ilistcagewind,icagecenter) - dw(ilistcagewind,icagecenter)

	d = d / sqrt( aw(ilistcagewind,icagecenter)**2+bw(ilistcagewind,icagecenter)**2+cw(ilistcagewind,icagecenter)**2 )

	com_wind_dist(ilistcagewind) = d 

!	write(20,*) aw(ilistcagewind,icagecenter),bw(ilistcagewind,icagecenter),cw(ilistcagewind,icagecenter) &
 !		,dw(ilistcagewind,icagecenter),com_wind_dist(ilistcagewind)
  
	enddo

	do ilistcagewind = 1,nlistcagewind

	   do jlistcagewind = 1,nlistcagewind

		if(ilistcagewind.ne.jlistcagewind) then

			if( abs(com_wind_dist(ilistcagewind)).gt.abs(com_wind_dist(jlistcagewind)) ) goto 60
			
		endif
	   enddo !jcagecenter

	   pt_distance = com_wind_dist(ilistcagewind)
       !    write(20,*) pt_distance
	   goto 70

60	enddo !icagecenter

!	write(20,*) com_p1,com_p2,com_p3,pt_distance - 21.8d0

!70	write(20,*) pt_distance 

!70	pt_distance = pt_distance + 21.8d0

!	write(*,*) pt_distance

70	end subroutine

			!................. SUBROUTINE TO ENERGY OF iTH MOLECULES FOR SAMPLING........................!

	subroutine samplingenergy

	use parameter
	implicit none

!	real :: start,finish

!	call cpu_time(start)

	if(translation.eq.0) then

	if(cartesian.eq.1) then
		if(part.eq.0) then
			com_p1 = com(1,igmol); com_p2 = com(2,igmol); com_p3 = com(3,igmol)	! p stands for present
		endif
		if(part.eq.1) then
			com_p1 = com(1,igmol); com_p2 = com(2,igmol); com_p3 = com(3,igmol)	! p stands for present
		endif
		if(part.eq.3) then
			com_p1 = boundary + 0.0001; com_p2 = com(2,igmol); com_p3 = com(3,igmol)	! p stands for present
		endif
		
		
	endif
	if(cartesian.eq.2) then
		com_p1 = com(1,igmol); com_p2 = com(2,igmol); com_p3 = com(3,igmol)
	endif
	if(cartesian.eq.3) then
		com_p1 = com(1,igmol); com_p2 = com(2,igmol); com_p3 = com(3,igmol)
	endif

	!write(*,*) translation,cartesian,com_p1,com_p2,com_p3

	endif

	if(translation.eq.1) then

	if(cartesian.eq.1) then
		if(part.eq.0) then
			com_p1 = com_new1; com_p2 = com(2,igmol); com_p3 = com(3,igmol)	! p stands for present
		endif
		if(part.eq.1) then
			com_p1 = boundary - 0.0001; com_p2 = com(2,igmol); com_p3 = com(3,igmol)	! p stands for present
		endif
		if(part.eq.3) then
			com_p1 = com_new1; com_p2 = com(2,igmol); com_p3 = com(3,igmol)	! p stands for present
		endif
		
		
	endif
	if(cartesian.eq.2) then
		com_p1 = com(1,igmol); com_p2 = com_new2; com_p3 = com(3,igmol) 
	endif
	if(cartesian.eq.3) then
		com_p1 = com(1,igmol); com_p2 = com(2,igmol); com_p3 = com_new3 
	endif

	!write(*,*) translation,cartesian,com_p1,com_p2,com_p3

	endif


	if(rotation.eq.0) atomic_p(:,:)=rg_atomic(:,:,igmol)
	if(rotation.eq.1) atomic_p(:,:)=rg_natomic(:,:)

	E = 0.0d0; guestguest_E = 0.0d0; guesthost_E = 0.0d0

	   do j=1,nlistgg(igmol)
		jgmol=listgg(j,igmol)

!	    do jgmol=1,ngmol
!		if(igmol.ne.jgmol) then

		dx=0.0d0;dy=0.0d0;dz=0.0d0;dxaa=0.0d0;dyaa=0.0d0;dzaa=0.0d0

		dx = com(1,jgmol) - com_p1;dy = com(2,jgmol) - com_p2;dz = com(3,jgmol) - com_p3

!	applying minimum boundary condition

		call minimumdisplacementcondition

!	calculation of inter & intra molecular energy

		if(gmoltype(igmol).eq.1) nigatom = ngatom1
		if(gmoltype(igmol).eq.2) nigatom = ngatom1

		if(gmoltype(jgmol).eq.1) njgatom = ngatom1
		if(gmoltype(jgmol).eq.2) njgatom = ngatom2

		   do igatom=1,nigatom
		      igatomtype = gatomtype(igatom,igmol)

		      do jgatom=1,njgatom
		      	jgatomtype = gatomtype(jgatom,jgmol)

		      	dxaa =  rg_atomic(1,jgatom,jgmol) + dx - atomic_p(1,igatom)
		      	dyaa =  rg_atomic(2,jgatom,jgmol) + dy - atomic_p(2,igatom)
		      	dzaa =  rg_atomic(3,jgatom,jgmol) + dz - atomic_p(3,igatom)

		      	rsq=dxaa*dxaa+dyaa*dyaa+dzaa*dzaa

		      	ratio = sigmasqgg(igatomtype,jgatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

			guestguest_E = guestguest_E + fourepsgg(igatomtype,jgatomtype) * ( s12 - s6 )

		      	E = E + fourepsgg(igatomtype,jgatomtype) * ( s12 - s6 )

		      enddo !jatom

		   enddo !iatom

!		endif 

	     enddo !j

!	do ihatom=1,nhatom
!	   ihatomtype = hatomtype(ihatom)

	do j=1,nlistgh(igmol)
	     ihatom = listgh(j,igmol)
	   ihatomtype = hatomtype(ihatom)

		if(gmoltype(igmol).eq.1) nigatom = ngatom1
		if(gmoltype(igmol).eq.2) nigatom = ngatom2

		do igatom=1,nigatom
		igatomtype=gatomtype(igatom,igmol)

			dx=0.0d0;dy=0.0d0;dz=0.0d0

			dx = rh_atomic(1,ihatom) - com_p1 - atomic_p(1,igatom)
			dy = rh_atomic(2,ihatom) - com_p2 - atomic_p(2,igatom)
			dz = rh_atomic(3,ihatom) - com_p3 - atomic_p(3,igatom)

!	applying minimum boundary condition

			call minimumdisplacementcondition

			rsq=dx*dx+dy*dy+dz*dz

		      	ratio = sigmasqgh(igatomtype,ihatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

			guesthost_E = guesthost_E + fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 )

		      	E = E + fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 )

		enddo !igatom

	    enddo !ihatom

!	call cpu_time(finish)

!	write(30,*) 'ith',finish-start

	end subroutine

        			!------------------ INHOMOGENEOUS SAMPLING (hot zones are placed along x-axis) -----------!

        subroutine perp_X_inhomosampling

        use parameter
        implicit none

        real :: T1,T2,RT1,RT2,pos_hot_zone1,pos_hot_zone2
        
        dE = 0.0d0; prob = 0.0d0; xt = 0.0d0

!	transition from (x,y,z) to (x',y',z')
!	First transition probability (x,y,z) to (x,y,z')

	dE_z = 0.0d0

	translation = 0
	cartesian = 3
	call samplingenergy
	Ez_old = E

	translation = 1
	call samplingenergy
	Ez_new = E

	dE_z = Ez_new - Ez_old
	prob_z = exp ( - dE_z / RThomo ) 
	!write(*,*) RThomo
	call random_number(xt)

	if (dE_z.lt.0.0d0.or.xt.lt.prob_z) then
		com(3,igmol) = com_new3
                if(gmoltype(igmol).eq.1) accpgg1trans = accpgg1trans + 1
                if(gmoltype(igmol).eq.2) accpgg2trans = accpgg2trans + 1
        else
                if(gmoltype(igmol).eq.1) rejgg1trans = rejgg1trans   + 1
                if(gmoltype(igmol).eq.2) rejgg2trans = rejgg2trans   + 1
        endif

	!write(*,*) com(:,igmol),dE_z,Ez_new,Ez_old
	!write(*,*) com_new1,com_new2,com_new3

!	Second transition probability (x,y,z') to (x,y',z')

	dE_y = 0.0d0

	translation = 0
	cartesian = 2
	call samplingenergy
	Ey_old = E

	translation = 1
	call samplingenergy
	Ey_new = E

	RT = R * temp(int(com(2,igmol)/0.01d0) )
	dE_y = Ey_new - Ey_old
	prob_y = exp ( - dE_y / RThomo ) 
	call random_number(xt)

	if (dE_y.lt.0.0d0.or.xt.lt.prob_y) then
		com(2,igmol) = com_new2 
                if(gmoltype(igmol).eq.1) accpgg1trans = accpgg1trans + 1
                if(gmoltype(igmol).eq.2) accpgg2trans = accpgg2trans + 1
        else
                if(gmoltype(igmol).eq.1) rejgg1trans = rejgg1trans   + 1
                if(gmoltype(igmol).eq.2) rejgg2trans = rejgg2trans   + 1
        endif

	!write(*,*) com(:,igmol),dE_y,Ey_new,Ey_old

!	Third transition probability (x,y',z') to (x',y',z')

	dE_x = 0.0d0

	T1 = temp(int(com(1,igmol)/0.01d0) ) ; T2 = temp(int(com_new1/0.01d0))
	if(T1.le.0.01) T1 = 300
	if(T2.le.0.01) T2 = 300
	RT1 = R * T1 ; RT2 = R * T2
	!write(20,*) T1,T2

	if(T1.eq.T2) then

		translation = 0
		part = 0
		cartesian = 1
		call samplingenergy
		Ex_old = E

		translation = 1
		call samplingenergy
		Ex_new = E

		dE_x = Ex_new - Ex_old
		prob_x = exp ( - dE_x / RT1 ) 

		if (dE_x.lt.0.0d0.or.xt.lt.prob_x) then
			com(1,igmol) = com_new1 
                        if(gmoltype(igmol).eq.1) accpgg1trans = accpgg1trans + 1
                        if(gmoltype(igmol).eq.2) accpgg2trans = accpgg2trans + 1
              else
                        if(gmoltype(igmol).eq.1) rejgg1trans = rejgg1trans   + 1
                        if(gmoltype(igmol).eq.2) rejgg2trans = rejgg2trans   + 1
              endif

		!write(*,*) com(:,igmol),dE_x,Ex_new,Ex_old,prob_x,xt


	endif

	if(T1.gt.T2) then

		!write(20,*) 'hot to cold'

		!df_pos = com_new1 - com(1,igmol)
		!write(*,*) df_pos,com_new1,com(1,igmol)
		!do itemp = 0,int(df_pos/0.01d0)
		!   	xtemp(itemp)= temp(int((itemp*0.01d0+com(1,igmol))/0.01d0))
	        !	xtemp(itemp+1)= temp(int(((itemp+1)*0.01d0+com(1,igmol))/0.01d0))
		!       	if(xtemp(itemp+1).lt.xtemp(itemp)) then
		!		boundary = (itemp+1)*0.01d0 + com(1,igmol)
		!		!write(*,*) boundary
		!		goto 50
		!	endif
		!write(20,*) xtemp(itemp),xtemp(itemp+1),itemp,int((itemp*0.01d0+com(1,igmol))/0.01d0),int(((itemp+1)*0.01d0+com(1,igmol))/0.01d0)
		!enddo

                pos_hot_zone1 = int(com(1,igmol)/6.2134d0) + 1
                pos_hot_zone2 = int(com_new1/6.2134d0) + 1
        
                if(pos_hot_zone1.eq.pos_hot_zone2) then
                        rs = 6.2134*pos_hot_zone1 - 0.1d0; ls = 6.2134*pos_hot_zone1 - 0.9d0
                        if(abs(com_new1-rs).le.abs(com_new1-ls)) then
                                boundary = rs
                        else
                                boundary = ls
                        endif
                else
                        boundary = 6.2134*pos_hot_zone1-0.1d0
                endif

                !write(20,*) com(1,igmol),com_new1,boundary,rs,ls,pos_hot_zone1,pos_hot_zone2

                 !write(20,*) pos_hot_zone1,pos_hot_zone2,rs,ls,boundary,com(1,igmol)-rs,com_new1-ls

		translation = 0
		cartesian = 1
		part = 1
		call samplingenergy
		Ex_old = E

		translation = 1
		cartesian = 1
		part = 1
		call samplingenergy
		Ex_new = E

                dE_x = 0.0d0        
		dE_x = Ex_new - Ex_old
                xt = 0.0d0
		prob_x1 = exp ( - dE_x / RT1 )
                call random_number(xt)

                if(dE_x.lt.0.0d0.or.xt.lt.prob_x1) then

                        if(gmoltype(igmol).eq.1) accpgg1trans1= accpgg1trans1+ 1
                        if(gmoltype(igmol).eq.2) accpgg2trans1= accpgg2trans1+ 1
                        
                        xt = 0.0d0
		        part = 2
		        prob_x2 = (T1/T2)
                        call random_number(xt)

                        if(xt.lt.prob_x2) then

                                if(gmoltype(igmol).eq.1) accpgg1trans2= accpgg1trans2+ 1
                                if(gmoltype(igmol).eq.2) accpgg2trans2= accpgg2trans2+ 1

		                translation = 0
		                cartesian = 1
		                part = 3
		                call samplingenergy
		                Ex_old = E

		                translation = 1
		                cartesian = 1
		                part = 3
		                call samplingenergy
		                Ex_new = E

                                xt = 0.0d0
                                dE_x = 0.0d0
		                dE_x = Ex_new - Ex_old
		                prob_x3 = exp ( - dE_x / RT2 )
                                call random_number(xt)

		!prob_x = prob_x1 * prob_x2 * prob_x3
		!call random_number(xt)

                                if (dE_x.lt.0.0d0.or.xt.lt.prob_x3) then
			                com(1,igmol) = com_new1 ;  com(2,igmol) = com_new2 ;  com(3,igmol) = com_new3
                                        if(gmoltype(igmol).eq.1) accpgg1trans3 = accpgg1trans3 + 1
                                        if(gmoltype(igmol).eq.2) accpgg2trans3 = accpgg2trans3 + 1
                                else
                                        if(gmoltype(igmol).eq.1) rejgg1trans3 = rejgg1trans3 + 1
                                        if(gmoltype(igmol).eq.2) rejgg2trans3 = rejgg2trans3 + 1
                                endif

                        else
                                if(gmoltype(igmol).eq.1) rejgg1trans2 = rejgg1trans2 + 1
                                if(gmoltype(igmol).eq.2) rejgg2trans2 = rejgg2trans2 + 1
                        endif

                else

                        if(gmoltype(igmol).eq.1) rejgg1trans1 = rejgg1trans1 + 1
                        if(gmoltype(igmol).eq.2) rejgg2trans1 = rejgg2trans1 + 1

                endif

	endif

	if(T1.lt.T2) then

		!write(20,*) 'cold to hot'

		!df_pos = com_new1 - com(1,igmol)
		!write(20,*) df_pos,com_new1,com(1,igmol)
		!do itemp = 0,int(df_pos/0.01d0)
		!   	xtemp(itemp)= temp(int((itemp*0.01d0+com(1,igmol))/0.01d0))
	        !		xtemp(itemp+1)= temp(int(((itemp+1)*0.01d0+com(1,igmol))/0.01d0))
		!	if(xtemp(itemp+1).gt.xtemp(itemp)) then
		!		boundary = (itemp+1)*0.01d0 + com(1,igmol)
		!		!write(*,*) boundary
		!		goto 60
		!	endif
		!write(20,*) xtemp(itemp),xtemp(itemp+1),itemp,int((itemp*0.01d0+com(1,igmol))/0.01d0),int(((itemp+1)*0.01d0+com(1,igmol))/0.01d0)
		!enddo


                pos_hot_zone1 = int(com(1,igmol)/6.2134d0) + 1
                pos_hot_zone2 = int(com_new1/6.2134d0) + 1

                if(pos_hot_zone1.eq.pos_hot_zone2) then
                        rs = 6.2134*pos_hot_zone1 - 0.1d0; ls = 6.2134*pos_hot_zone1 - 0.9d0
                        if(abs(com(1,igmol)-rs).le.abs(com(1,igmol)-ls)) then
                                boundary = rs
                        else
                                boundary = ls
                        endif
                else
                        boundary = 6.2134*pos_hot_zone2 - 0.1d0
                endif

                 !write(20,*) com(1,igmol),com_new1,boundary,rs,ls,pos_hot_zone1,pos_hot_zone2

                 !write(20,*) pos_hot_zone1,pos_hot_zone2,rs,ls,boundary,com(1,igmol)-rs,com_new1-ls

		translation = 0
		cartesian = 1
		part = 1
		call samplingenergy
		Ex_old = E

		translation = 1
		cartesian = 1
		part = 1
		call samplingenergy
		Ex_new = E

                dE_x = 0.0d0
		dE_x = Ex_new - Ex_old
                xt = 0.0d0
		prob_x1 = exp ( - dE_x / RT1 )
	        call random_number(xt)

                if (dE_x.lt.0.0d0.or.xt.lt.prob_x1) then

			if(gmoltype(igmol).eq.1) accpgg1trans1= accpgg1trans1+ 1
                        if(gmoltype(igmol).eq.2) accpgg2trans1= accpgg2trans1+ 1

		        part = 2
                        xt = 0.0d0
                        call random_number(xt)
		        prob_x2 = (T1/T2)
		        !write(*,*) prob_x2
                        if (xt.lt.prob_x2) then

				if(gmoltype(igmol).eq.1) accpgg1trans2= accpgg1trans2+ 1
                                if(gmoltype(igmol).eq.2) accpgg2trans2= accpgg2trans2+ 1

		                translation = 0
		                cartesian = 1
		                part = 3
		                call samplingenergy
		                Ex_old = E

		                translation = 1
		                cartesian = 1
		                part = 3
		                call samplingenergy
		                Ex_new = E
                                
                                dE_x = 0.0d0
		                dE_x = Ex_new - Ex_old
                                xt = 0.0d0
		                prob_x3 = exp ( - dE_x / RT2 )
				call random_number(xt)
		                !write(*,*) prob_x3,dE_x,Ex_new,Ex_old

		!prob_x = prob_x1 * prob_x2 * prob_x3
		!call random_number(xt)
		!write(20,*) prob_x , prob_x1 , prob_x2 , prob_x3,xt


                                if (dE_x.lt.0.0d0.or.xt.lt.prob_x3) then
			                com(1,igmol) = com_new1 ;  com(2,igmol) = com_new2 ;  com(3,igmol) = com_new3
                                        if(gmoltype(igmol).eq.1) accpgg1trans = accpgg1trans + 1
                                        if(gmoltype(igmol).eq.2) accpgg2trans = accpgg2trans + 1
                                else
                                        if(gmoltype(igmol).eq.1) rejgg1trans = rejgg1trans   + 1
                                        if(gmoltype(igmol).eq.2) rejgg2trans = rejgg2trans   + 1
                                endif
		!write(20,*) com(1,igmol),com_new1

                        else
                                if(gmoltype(igmol).eq.1) rejgg1trans2 = rejgg1trans2 + 1
                                if(gmoltype(igmol).eq.2) rejgg2trans2 = rejgg2trans2 + 1
                        endif

                else

                        if(gmoltype(igmol).eq.1) rejgg1trans1 = rejgg1trans1 + 1
                        if(gmoltype(igmol).eq.2) rejgg2trans1 = rejgg2trans1 + 1

                endif
                

	endif
                

        end subroutine

           !---------------- SAMPLING OF COM IN INHOMOGENEOUs MEDIUM ( PERPENDICULAR TO CAGE TO WINDOW )-----------------------!

        subroutine inhomosampling

        use parameter
        implicit none

        dE = 0.0d0; prob = 0.0d0; xt = 0.0d0
        dE = E_new - E_old

!       pt_distance_old =   ( com(1,igmol) + com(2,igmol) + com(3,igmol)!    ) * 0.577350269 - 21.52384897 
!       pt_distance_new =   ( com_new1 + com_new2 + com_new3 ) *       0.577350269 - 21.52384897

        !write(20,*) imc,igmol,icagecenter,ilistcagewind,pt_distance_old,pt_distance_new

        if(pt_distance_old.gt.-ls.and.pt_distance_old.lt.-rs) then
                if(pt_distance_new.gt.-ls.and.pt_distance_new.lt.-rs) then

                 !write(20,*)
                 !j,imc,igmol,icagecenter,ilistcagewind,pt_distance_old,pt_distance_new,'hh'
!
!               HOT  TO HOT TRANSITION
!
                        if (dE.gt.0.0d0) then
                        prob = exp ( - dE/RTh )
                        call random_number(xt)
                        endif

                        if (dE.lt.0.0d0.or.xt.lt.prob) then
                                com(1,igmol) = com_new1 ;  com(2,igmol) = com_new2 ;  com(3,igmol) = com_new3
                                if(gmoltype(igmol).eq.1) accpgg1trans = accpgg1trans + 1
                                if(gmoltype(igmol).eq.2) accpgg2trans = accpgg2trans + 1
                        else
                                if(gmoltype(igmol).eq.1) rejgg1trans = rejgg1trans   + 1
                                if(gmoltype(igmol).eq.2) rejgg2trans = rejgg2trans   + 1
                        endif

                else

                 !write(20,*)
                 !j,imc,igmol,icagecenter,ilistcagewind,pt_distance_old,pt_distance_new,'hc'

!
!               HOT TO COLD TRANSITION
!
                        !write(*,*) 'HC'
                        prob=exp(-dE*(1/RTl+1/RTh))
                        call random_number(xt)

                        if (xt.lt.prob) then
                                com(1,igmol) = com_new1 ;  com(2,igmol) = com_new2 ;  com(3,igmol) = com_new3
                                if(gmoltype(igmol).eq.1) accpgg1trans = accpgg1trans + 1
                                if(gmoltype(igmol).eq.2) accpgg2trans = accpgg2trans + 1
                        else
                                if(gmoltype(igmol).eq.1) rejgg1trans = rejgg1trans   + 1
                                if(gmoltype(igmol).eq.2) rejgg2trans = rejgg2trans   + 1
                        endif

                endif
        
        else

        
                if(pt_distance_new.gt.-ls.and.pt_distance_new.lt.-rs) then
!               if(pt_distance_new.gt.-1.6d0.and.pt_distance_new.lt.-0.6d0)
!               then

                ! write(20,*)
                ! j,imc,igmol,icagecenter,ilistcagewind,pt_distance_old,pt_distance_new,'ch'

!               COLD TO HOT TRANSITION

                        prob=exp(-dE*(1/RTh+1/RTl))
                        call random_number(xt)

                        if (xt.lt.prob) then
                                com(1,igmol) = com_new1 ;  com(2,igmol) = com_new2 ;  com(3,igmol) = com_new3
                                if(gmoltype(igmol).eq.1) accpgg1trans = accpgg1trans + 1
                                if(gmoltype(igmol).eq.2) accpgg2trans = accpgg2trans + 1
                        else
                                if(gmoltype(igmol).eq.1) rejgg1trans = rejgg1trans   + 1
                                if(gmoltype(igmol).eq.2) rejgg2trans = rejgg2trans   + 1
                        endif

                else

                 !write(20,*)
                 !j,imc,igmol,icagecenter,ilistcagewind,pt_distance_old,pt_distance_new,'cc'

!               COLD TO COLD TRANSITION


                        if (dE.gt.0.0d0) then
                        prob = exp ( - dE/RTl )
                        call random_number(xt)
                        endif

                        if (dE.lt.0.0d0.or.xt.lt.prob) then
                                com(1,igmol) = com_new1 ;  com(2,igmol) = com_new2 ;  com(3,igmol) = com_new3
                                if(gmoltype(igmol).eq.1) accpgg1trans = accpgg1trans + 1
                                if(gmoltype(igmol).eq.2) accpgg2trans = accpgg2trans + 1
                        else
                                if(gmoltype(igmol).eq.1) rejgg1trans = rejgg1trans   + 1
                                if(gmoltype(igmol).eq.2) rejgg2trans = rejgg2trans   + 1
                        endif

                endif

        endif

        end subroutine


		!----------------CALCULATION OF TRANSLATION MOTION TO CHOOSE BETWEEN HOMO OR INHOMO SAMPLING------!

	subroutine perp_cage_wind_inhomosampling

	use parameter
	implicit none

	!write(*,*) icagecenter,ilistcagewind

	if(icagecenter.eq.1) then
		if(ilistcagewind.eq.1.or.ilistcagewind.eq.2) then
			call inhomosampling
			!write(*,*) 'inhomo'
		else
			call transhomosampling
			!write(*,*) 'homo'
		endif
	endif
	if(icagecenter.eq.2) then
		if(ilistcagewind.eq.1.or.ilistcagewind.eq.4) then
			call inhomosampling
			!write(*,*) 'inhomo'
		else
			call transhomosampling
			!write(*,*) 'homo'
		endif
	endif
	if(icagecenter.eq.3) then
		if(ilistcagewind.eq.2.or.ilistcagewind.eq.4) then
			call inhomosampling
			!write(*,*) 'inhomo'
		else
			call transhomosampling
			!write(*,*) 'homo'
		endif
	endif
	if(icagecenter.eq.4) then
		if(ilistcagewind.eq.2.or.ilistcagewind.eq.4) then
			call inhomosampling
			!write(*,*) 'inhomo'
		else
			call transhomosampling
			!write(*,*) 'homo'
		endif
	endif
	if(icagecenter.eq.5) then
		if(ilistcagewind.eq.3.or.ilistcagewind.eq.4) then
			call inhomosampling
			!write(*,*) 'inhomo'
		else
			call transhomosampling
			!write(*,*) 'homo'
		endif
	endif
	if(icagecenter.eq.6) then
		if(ilistcagewind.eq.3.or.ilistcagewind.eq.4) then
			call inhomosampling
			!write(*,*) 'inhomo'
		else
			call transhomosampling
			!write(*,*) 'homo'
		endif
	endif
	if(icagecenter.eq.7) then
		if(ilistcagewind.eq.1.or.ilistcagewind.eq.2) then
			call inhomosampling
			!write(*,*) 'inhomo'
		else
			call transhomosampling
			!write(*,*) 'homo'
		endif
	endif
	if(icagecenter.eq.8) then
		if(ilistcagewind.eq.1.or.ilistcagewind.eq.2) then
			call inhomosampling
			!write(*,*) 'inhomo'
		else
			call transhomosampling
			!write(*,*) 'homo'
		endif
	endif
!	if(icagecenter.eq.5.or.icagecenter.eq.2.or.icagecenter.eq.7.or.icagecenter.eq.4) call homosampling

!	if(icagecenter.eq.1.or.icagecenter.eq.2.or.icagecenter.eq.3.or.icagecenter.eq.4) then
!		if(ilistcagewind.eq.1.or.icagecenter.eq.8.or.icagecenter.eq.12.or.icagecenter.eq.16) then
!			write(20,*) icagecenter,ilistcagewind
!			call inhomosampling
			!write(*,*) 'inhomosampling'
!		else
!			call homosampling
!		endif
!	else
!		call homosampling
!	endif

	end subroutine


		!**************************************** FLEXIBILITY STUDY ******************************************!

			!--------------- SUBROUTINE TO CALCULATE THE NEIGHBOR LIST ----------------!

	subroutine neighborlistflex
	use parameter
	implicit none

	real*8 :: rch

!	VDW NEIGHBORLIST

	do ihatom=1,nhatom
	   if(hatomtype(ihatom).eq.3.or.hatomtype(ihatom).eq.4) then
	   vdw_nlist(ihatom)=0
	   do jhatom=1,nhatom
	   if(hatomtype(jhatom).eq.4.or.hatomtype(jhatom).eq.3) then
	   if(ihatom.ne.jhatom) then
		dx = rh_atomic(1,ihatom) - rh_atomic(1,jhatom)
		dy = rh_atomic(2,ihatom) - rh_atomic(2,jhatom)
		dz = rh_atomic(3,ihatom) - rh_atomic(3,jhatom)

!	APPLYING MINIMUM BONDARY CONDITION

		call minimumdisplacementcondition
			d = sqrt ( dx*dx + dy*dy + dz*dz )

		rch = 2.0d0

!	applying cut-off sphere radius
	        if (d.le.rch+2.0d0) then
		   vdw_nlist(ihatom) = vdw_nlist(ihatom) + 1
		   vdw_list(vdw_nlist(ihatom),ihatom)= jhatom
	!		write(*,*) ihatom,jhatom,d
	        endif
	      endif
	      endif
	     enddo !j
	endif
	enddo !i

	do ihatom=1,nhatom
	!   write(5,*) ihatom,vdw_nlist(ihatom),(vdw_list(jhatom,ihatom),jhatom=1,vdw_nlist(ihatom))
	enddo !ihatom
	   	
!	LINEAR BOND NEIGHBORLIST

        if(imc.eq.1) then	
	do ihatom = 1,nhatom
	   bond_nlist(ihatom)=0
	   do jhatom = 1,nhatom
		if(ihatom.ne.jhatom) then

		dx = rh_atomic(1,ihatom) - rh_atomic(1,jhatom)
		dy = rh_atomic(2,ihatom) - rh_atomic(2,jhatom)
		dz = rh_atomic(3,ihatom) - rh_atomic(3,jhatom)

!	APPLYING MINIMUM BONDARY CONDITION

		call minimumdisplacementcondition
		d = sqrt ( dx*dx + dy*dy + dz*dz )
			
		if (d.lt.1.8d0) then
			bond_nlist(ihatom) = bond_nlist(ihatom) + 1
		   	bond_list(bond_nlist(ihatom),ihatom)= jhatom
		endif
		        
	   endif
	   enddo ! jhatom
	enddo ! ihatom
        

!	WRITE NEIGHBOR LIST OF LINEAR BONDED ATOMS

	!write(4,*) ' List of Linear Bonded atoms and Bond lengths'
	!write(4,*)

	do ihatom=1,nhatom
	  ! write(4,*) ihatom,bond_nlist(ihatom),(bond_list(jhatom,ihatom),jhatom=1,bond_nlist(ihatom))
	enddo !ihatom

!	CREATION OF NEIGHBOR LIST OF ANGULAR BONDED ATOMS

	do ihatom = 1,nhatom
	   nangle_list(ihatom) = 0
	   if(bond_nlist(ihatom).ne.0) then
	   do j = 1,bond_nlist(ihatom)
	    
	     do k = j+1,bond_nlist(ihatom)

		nangle_list(ihatom) = nangle_list(ihatom) + 1

!		write(4,*) jatom,iatom,katom

	       enddo !k

	     enddo !j

	     endif
	  ! write(4,*) ihatom,nangle_list(ihatom)

	enddo !ihatom
        endif

	end subroutine

				!--------------- SUBROUTINE TO CALCULATE INTRAMOLECULE ENERGY -------------------!

	subroutine intramolenergy
	use parameter
	implicit none
	
	E_totalintra = 0.0d0; E_vdw = 0.0d0; E_bond = 0.0d0; E_angle = 0.0d0; E_ub = 0.0d0

!	VAN-DER-WAALS ENERGY

	E=0.0d0
	do ihatom = 1,nhatom

        if(vdw_nlist(ihatom).ne.0) then

	   do j = 1,vdw_nlist(ihatom)
		jhatom = vdw_list(j,ihatom)

		dx = rh_atomic(1,jhatom) - rh_atomic(1,ihatom)
		dy = rh_atomic(2,jhatom) - rh_atomic(2,ihatom)
		dz = rh_atomic(3,jhatom) - rh_atomic(3,ihatom)

!	APPLYING MINIMUM BONDARY CONDITION

		call minimumdisplacementcondition
		d = sqrt ( dx*dx + dy*dy + dz*dz )

		ratio = sigmasqhh(ihatomtype,jhatomtype)/(d*d)
		s6 = ratio * ratio * ratio
		s12 = s6 * s6

		E = fourepshh(ihatomtype,jhatomtype) * ( s12 - s6 )
		E_vdw = E_vdw + E

	  enddo !j

        endif

!	LINEAR BONDED ENERGY


	   if(bond_nlist(ihatom).ne.0) then
	   do j = 1,bond_nlist(ihatom)
		jhatom = bond_list(j,ihatom)

		dx = rh_atomic(1,jhatom) - rh_atomic(1,ihatom)
		dy = rh_atomic(2,jhatom) - rh_atomic(2,ihatom)
		dz = rh_atomic(3,jhatom) - rh_atomic(3,ihatom)

!	APPLYING MINIMUM BONDARY CONDITION

		call minimumdisplacementcondition
		d = sqrt ( dx*dx + dy*dy + dz*dz )
		call linearbondenergy
		        
	   enddo ! j

	   endif


!	ANGULAR BONDED ENERGY

	   if(bond_nlist(ihatom).ne.0) then
	   do j = 1,bond_nlist(ihatom)
	     jhatom = bond_list(j,ihatom) 

	     dx = rh_atomic(1,ihatom) - rh_atomic(1,jhatom)
	     dy = rh_atomic(2,ihatom) - rh_atomic(2,jhatom)
	     dz = rh_atomic(3,ihatom) - rh_atomic(3,jhatom)

!	APPLYING MINIMUM IMAGE CONDITION

	     call minimumdisplacementcondition
	     d = sqrt ( dx*dx + dy*dy + dz*dz )
	     dx1=dx;dy1=dy;dz1=dz;d1=d

	     do k = j+1,bond_nlist(ihatom)

		khatom = bond_list(k,ihatom) 

		dx = rh_atomic(1,ihatom) - rh_atomic(1,khatom)
		dy = rh_atomic(2,ihatom) - rh_atomic(2,khatom)
		dz = rh_atomic(3,ihatom) - rh_atomic(3,khatom)

!	APPLYING MINIMUM IMAGE CONDITION

	     	call minimumdisplacementcondition
		d = sqrt ( dx*dx + dy*dy + dz*dz )
		dx2=dx;dy2=dy;dz2=dz;d2=d
	
		angle = dacos( (dx1*dx2 + dy1*dy2 + dz1*dz2) / ( d1*d2 ) )!*57.29d0

		call angularbondedenergy

	       enddo !k

	     enddo !j

	  endif

	  if(bond_nlist(ihatom).ne.0) then

		do j = 1,bond_nlist(ihatom)

	      	   jhatom = bond_list(j,ihatom)

		   dx = rh_atomic(1,jhatom) - rh_atomic(1,ihatom)
		   dy = rh_atomic(2,jhatom) - rh_atomic(2,ihatom)
		   dz = rh_atomic(3,jhatom) - rh_atomic(3,ihatom)

		   call minimumdisplacementcondition
		d = sqrt ( dx*dx + dy*dy + dz*dz )
	     	   dx1=dx;dy1=dy;dz1=dz;d1=d

	      	   do k = 1,bond_nlist(jhatom)

			khatom = bond_list(k,jhatom)
			if(khatom.ne.ihatom) then

			dx = rh_atomic(1,jhatom) - rh_atomic(1,khatom)
			dy = rh_atomic(2,jhatom) - rh_atomic(2,khatom)
			dz = rh_atomic(3,jhatom) - rh_atomic(3,khatom)
		   

	     	   	call minimumdisplacementcondition
			d = sqrt ( dx*dx + dy*dy + dz*dz )
		   	dx2=dx;dy2=dy;dz2=dz;d2=d
	
		   	angle = dacos( (dx1*dx2 + dy1*dy2 + dz1*dz2) / ( d1*d2 ) ) !*57.29d0
		   	call angularbondedenergy
			endif
	       enddo !k

	enddo !j

	endif

	enddo !iatom

	E_totalintra = E_bond*0.5d0 + E_angle + E_vdw*0.5d0 + E_ub

	end subroutine

			!------------- SUBROUTINE TO CALCULATE ENERGY BETWEEN iTH ATOM AND ALL OTHER ATOMS ---------------!

	subroutine intramolenergyith
	use parameter
	implicit none

	E_vdw = 0.0d0; E_bond = 0.0d0;	E_angle = 0.0d0; E_ub = 0.0d0

	if(movement.eq.0) rhatomic_p(:)=rh_atomic(:,ihatom)
	if(movement.eq.1) rhatomic_p(:)=rh_natomic(:)

!	VDW ENERGY

        if(vdw_nlist(ihatom).ne.0) then
	do j=1,vdw_nlist(ihatom)
	jhatom=vdw_list(j,ihatom)

        ihatomtype = hatomtype(ihatom)
        jhatomtype = hatomtype(jhatom)

        if(ihatomtype.eq.3.and.jhatomtype.eq.3) goto 40
        if(ihatomtype.eq.4.and.jhatomtype.eq.4) goto 40

        !write(40,*) j,ihatom,jhatom,vdw_nlist(ihatom)
	
	dx = rh_atomic(1,jhatom) - rhatomic_p(1)
	dy = rh_atomic(2,jhatom) - rhatomic_p(2)
	dz = rh_atomic(3,jhatom) - rhatomic_p(3)

!	APPLYING MINIMUM BONDARY CONDITION

	call minimumdisplacementcondition
	d = sqrt ( dx*dx + dy*dy + dz*dz )
	  
	ratio = sigmasqhh(ihatomtype,jhatomtype)/(d*d)
	s6 = ratio * ratio * ratio
	s12 = s6 * s6

	E = fourepshh(ihatomtype,jhatomtype) * ( s12 - s6 )
	E_vdw = E_vdw + E

        !write(40,*) j,E,E_vdw,d,(s12-s6),fourepshh(ihatomtype,jhatomtype),ihatomtype,jhatomtype

40	enddo !j
        endif

!	LINEAR BONDED ENERGY

	if(bond_nlist(ihatom).ne.0) then
	   do j = 1,bond_nlist(ihatom)
		jhatom = bond_list(j,ihatom)

		dx = rhatomic_p(1) - rh_atomic(1,jhatom)
		dy = rhatomic_p(2) - rh_atomic(2,jhatom)
		dz = rhatomic_p(3) - rh_atomic(3,jhatom)

!	APPLYING MINIMUM BONDARY CONDITION

		call minimumdisplacementcondition
		d = sqrt ( dx*dx + dy*dy + dz*dz )
		call linearbondenergy
		        
	   enddo ! jatom

	endif

!	ANGULAR BONDED ENERGY
	
	   if(bond_nlist(ihatom).ne.0) then
	   do j = 1,bond_nlist(ihatom)
	     jhatom = bond_list(j,ihatom) 

	     dx = rhatomic_p(1) - rh_atomic(1,jhatom)
	     dy = rhatomic_p(2) - rh_atomic(2,jhatom)
	     dz = rhatomic_p(3) - rh_atomic(3,jhatom)

!	APPLYING MINIMUM IMAGE CONDITION

	     call minimumdisplacementcondition
	     d = sqrt ( dx*dx + dy*dy + dz*dz )
	     dx1=dx;dy1=dy;dz1=dz;d1=d

	     do k = j+1,bond_nlist(ihatom)

		khatom = bond_list(k,ihatom) 

		dx = rhatomic_p(1) - rh_atomic(1,khatom)
		dy = rhatomic_p(2) - rh_atomic(2,khatom)
		dz = rhatomic_p(3) - rh_atomic(3,khatom)

!	APPLYING MINIMUM IMAGE CONDITION

	     	call minimumdisplacementcondition
		d = sqrt ( dx*dx + dy*dy + dz*dz )
		dx2=dx;dy2=dy;dz2=dz;d2=d
	
		angle = dacos( (dx1*dx2 + dy1*dy2 + dz1*dz2) / ( d1*d2 ) )!*57.29d0

		call angularbondedenergy

	       enddo !k

	     enddo !j

	  endif

 	if(bond_nlist(ihatom).ne.0) then

		do j = 1,bond_nlist(ihatom)

	      	   jhatom = bond_list(j,ihatom)

		   dx = rh_atomic(1,jhatom) - rhatomic_p(1)
		   dy = rh_atomic(2,jhatom) - rhatomic_p(2)
		   dz = rh_atomic(3,jhatom) - rhatomic_p(3)

		   call minimumdisplacementcondition
		   d = sqrt ( dx*dx + dy*dy + dz*dz )
	     	   dx1=dx;dy1=dy;dz1=dz;d1=d

	      	   do k = 1,bond_nlist(jhatom)

			khatom = bond_list(k,jhatom)
			if(khatom.ne.ihatom) then

			dx = rh_atomic(1,jhatom) - rh_atomic(1,khatom)
			dy = rh_atomic(2,jhatom) - rh_atomic(2,khatom)
			dz = rh_atomic(3,jhatom) - rh_atomic(3,khatom)
		   

	     	   	call minimumdisplacementcondition
			d = sqrt ( dx*dx + dy*dy + dz*dz )
		   	dx2=dx;dy2=dy;dz2=dz;d2=d
	
		   	angle = dacos( (dx1*dx2 + dy1*dy2 + dz1*dz2) / ( d1*d2 ) )!*57.29d0
		   	call angularbondedenergy

			dx = rhatomic_p(1) - rh_atomic(1,khatom)
			dy = rhatomic_p(2) - rh_atomic(2,khatom)
			dz = rhatomic_p(3) - rh_atomic(3,khatom)
		   

	     	   	call minimumdisplacementcondition
			d = sqrt ( dx*dx + dy*dy + dz*dz )
	
			call ubenergy
			
			endif
	       enddo !k

	enddo !j

	endif

	!write(30,*) iatom,E_bond,E_angle

	end subroutine
		
		

			!----------------- SUBROUTINE TO CALCULATE LINEAR BONDED ENERGY ---------------!

	subroutine linearbondenergy

	use parameter
	implicit none

	!E_bond=0.0d0

	if(hatomtype(ihatom).eq.1.and.hatomtype(jhatom).eq.3) then
	    E_bond = E_bond + force_constSiO * (equ_distSiO - d) * (equ_distSiO - d)
	endif

        if(hatomtype(ihatom).eq.2.and.hatomtype(jhatom).eq.3) then
            E_bond = E_bond + force_constSiO * (equ_distSiO - d) * (equ_distSiO - d)
        endif

	if(hatomtype(ihatom).eq.3.and.hatomtype(jhatom).eq.1) then
	    E_bond = E_bond + force_constSiO * (equ_distSiO - d) * (equ_distSiO - d)
	endif

        if(hatomtype(ihatom).eq.3.and.hatomtype(jhatom).eq.2) then
            E_bond = E_bond + force_constSiO * (equ_distSiO - d) * (equ_distSiO - d)
        endif

	!write(40,*) hatomtype(ihatom),hatomtype(jhatom)
	!write(40,*) ihatom,jhatom,E_bond,force_constSiO,equ_distSiO,d

	end subroutine

		!---------------------- SUBROUTINE TO CALCULATE ANGULAR BONDED ENERGY ---------------------!

	subroutine angularbondedenergy

	use parameter
	implicit none

	!E_angle = 0.0d0

	if(hatomtype(jhatom).eq.3.and.hatomtype(ihatom).eq.1.and.hatomtype(khatom).eq.3) then
	    E_angle = E_angle + force_constOSiO * (equ_angleOSiO - angle) * (equ_angleOSiO - angle)
	endif 

        if(hatomtype(jhatom).eq.3.and.hatomtype(ihatom).eq.2.and.hatomtype(khatom).eq.3) then
            E_angle = E_angle + force_constOSiO * (equ_angleOSiO - angle) * (equ_angleOSiO - angle)
        endif

	if(hatomtype(jhatom).eq.3.and.hatomtype(ihatom).eq.1.and.hatomtype(khatom).eq.1) then
	    E_angle = E_angle + force_constSiOSi * (equ_angleSiOSi - angle) * (equ_angleSiOSi - angle) 
	endif 

        if(hatomtype(jhatom).eq.3.and.hatomtype(ihatom).eq.2.and.hatomtype(khatom).eq.1) then
            E_angle = E_angle + force_constSiOSi * (equ_angleSiOSi - angle) * (equ_angleSiOSi - angle) 
        endif

         if(hatomtype(jhatom).eq.3.and.hatomtype(ihatom).eq.1.and.hatomtype(khatom).eq.2) then
            E_angle = E_angle + force_constSiOSi * (equ_angleSiOSi - angle) * (equ_angleSiOSi - angle) 
        endif

	!write(50,*) jhatom,ihatom,khatom,E_angle,angle!,d1,d2

	end subroutine

		!---------------------- SUBROUTINE TO CALCULATE UREY-BRADLEY ENERGY ---------------------!

	subroutine ubenergy

	use parameter
	implicit none

	!E_ub = 0.0d0

	if(hatomtype(jhatom).eq.3.and.hatomtype(ihatom).eq.1.and.hatomtype(khatom).eq.1) then
	    E_ub = E_ub + force_constub * ( equ_ub - d ) * ( equ_ub - d )
	endif 

        if(hatomtype(jhatom).eq.3.and.hatomtype(ihatom).eq.2.and.hatomtype(khatom).eq.1) then
            E_ub = E_ub + force_constub * ( equ_ub - d ) * ( equ_ub - d )
        endif

        if(hatomtype(jhatom).eq.3.and.hatomtype(ihatom).eq.1.and.hatomtype(khatom).eq.2) then
            E_ub = E_ub + force_constub * ( equ_ub - d ) * ( equ_ub - d )
        endif

	!write(60,*) jhatom,ihatom,khatom,E_ub,d, equ_ub

	end subroutine

		!----------------------SUBROUTINE TO CALCULATE THE INTRA ANGULAR ENERGY OF HYDROCARBONS-----------!

	subroutine intraenergy_angle

	use parameter
	implicit none

	if(gmoltype(igmol).eq.1) ngatom = ngatom1
	if(gmoltype(igmol).eq.2) ngatom = ngatom2
	
	if(rotation.eq.0) atomic_p(:,:)=rg_atomic(:,:,igmol)
	if(rotation.eq.1) atomic_p(:,:)=rg_natomic(:,:)

	Eintra_angle = 0.0d0

	do jangle = 1,nangle(gmoltype(igmol))

	  if(iangle.ne.jangle) then

	  E = 0.0d0
	   
	  iangle1 = index_angle(1,jangle,gmoltype(igmol))
	  iangle2 = index_angle(2,jangle,gmoltype(igmol))
	  iangle3 = index_angle(3,jangle,gmoltype(igmol))

!	calculate the angle

	  uxab = ( atomic_p(1,iangle1) - atomic_p(1,iangle2) ) / bondlength(iangle1,iangle2,igmol)
	  uyab = ( atomic_p(2,iangle1) - atomic_p(2,iangle2) ) / bondlength(iangle1,iangle2,igmol)
	  uzab = ( atomic_p(3,iangle1) - atomic_p(3,iangle2) ) / bondlength(iangle1,iangle2,igmol)

	  uxbc = -( atomic_p(1,iangle2) - atomic_p(1,iangle3) ) / bondlength(iangle2,iangle3,igmol)
	  uybc = -( atomic_p(2,iangle2) - atomic_p(2,iangle3) ) / bondlength(iangle2,iangle3,igmol)
	  uzbc = -( atomic_p(3,iangle2) - atomic_p(3,iangle3) ) / bondlength(iangle2,iangle3,igmol)

	  theta_angle = dacos ( uxab*uxbc + uyab*uybc + uzab*uzbc )

	  E=angle_parm(gmoltype(igmol))*0.5d0*(theta_angle-equi_angle(gmoltype(igmol)))*(theta_angle-equi_angle(gmoltype(igmol))) 

!	calculate the energy of angle

	  Eintra_angle = Eintra_angle + E

	 ! write(*,*) igmol,jangle,iangle1,iangle2,iangle3,Eintra_angle,theta_angle

	  endif

	 enddo !jangle


	end subroutine

!			---------INTRA MOLECULAR ANGLE CALCULATION ---------		!

	subroutine intramol_angle

	use parameter
	implicit none

	if(gmoltype(igmol).eq.1) ngatom = ngatom1
	if(gmoltype(igmol).eq.2) ngatom = ngatom2
	
	if(rotation.eq.0) atomic_p(:,:)=rg_atomic(:,:,igmol)
	if(rotation.eq.1) atomic_p(:,:)=rg_natomic(:,:)

	theta_angle = 0.0d0
	uxab=0.0d0;uyab=0.0d0;uzab=0.0d0;uxbc=0.0d0;uybc=0.0d0;uzbc=0.0d0;theta_angle=0.0d0

	iangle1 = index_angle(1,iangle,gmoltype(igmol))
	iangle2 = index_angle(2,iangle,gmoltype(igmol))
	iangle3 = index_angle(3,iangle,gmoltype(igmol))

!	calculate the angle

	uxab = ( atomic_p(1,iangle1) - atomic_p(1,iangle2) ) / bondlength(iangle1,iangle2,igmol)
	uyab = ( atomic_p(2,iangle1) - atomic_p(2,iangle2) ) / bondlength(iangle1,iangle2,igmol)
	uzab = ( atomic_p(3,iangle1) - atomic_p(3,iangle2) ) / bondlength(iangle1,iangle2,igmol)

	uxbc = -( atomic_p(1,iangle2) - atomic_p(1,iangle3) ) / bondlength(iangle2,iangle3,igmol)
	uybc = -( atomic_p(2,iangle2) - atomic_p(2,iangle3) ) / bondlength(iangle2,iangle3,igmol)
	uzbc = -( atomic_p(3,iangle2) - atomic_p(3,iangle3) ) / bondlength(iangle2,iangle3,igmol)

	theta_angle = dacos ( uxab*uxbc + uyab*uybc + uzab*uzbc )

	!write(20,*) uxab,uyab,uzab,uxab**2+uyab**2+uzab**2,bondlength(iangle1,iangle2,igmol)
	!write(20,*) uxbc,uybc,uzbc,uxbc**2+uybc**2+uzbc**2,bondlength(iangle2,iangle3,igmol)
	!write(20,*) theta_angle,uxab*uxbc + uyab*uybc + uzab*uzbc

	end subroutine

!	----------------   ROTATION OF INTRA ATOMIC POSITION   --------------	!

	subroutine intraangle_rotation

	use parameter
	implicit none

	nx = 0.0d0; ny = 0.0d0; nz = 0.0d0
	nxab = 0.0d0; nyab = 0.0d0; nzab = 0.0d0
	nxbc = 0.0d0; nybc = 0.0d0; nzbc = 0.0d0	
	uxab = 0.0d0; uyab = 0.0d0; uzab = 0.0d0
	uxbc = 0.0d0; uybc = 0.0d0; uzbc = 0.0d0

	nxab = ( rg_atomic(1,iangle1,igmol) - rg_atomic(1,iangle2,igmol) ) / bondlength(iangle1,iangle2,igmol)
	nyab = ( rg_atomic(2,iangle1,igmol) - rg_atomic(2,iangle2,igmol) ) / bondlength(iangle1,iangle2,igmol)
	nzab = ( rg_atomic(3,iangle1,igmol) - rg_atomic(3,iangle2,igmol) ) / bondlength(iangle1,iangle2,igmol)


	!write(20,*) 'unit vector',nxab,nyab,nzab,sqrt(nxab**2+nyab**2+nzab**2)

	nxbc = -( rg_atomic(1,iangle2,igmol) - rg_atomic(1,iangle3,igmol) ) / bondlength(iangle2,iangle3,igmol)
	nybc = -( rg_atomic(2,iangle2,igmol) - rg_atomic(2,iangle3,igmol) ) / bondlength(iangle2,iangle3,igmol)
	nzbc = -( rg_atomic(3,iangle2,igmol) - rg_atomic(3,iangle3,igmol) ) / bondlength(iangle2,iangle3,igmol)

	!write(20,*) 'unit vector',nxbc,nybc,nzbc,sqrt(nxbc**2+nybc**2+nzbc**2)

	ux = nyab*nzbc - nzab*nybc
	uy = nxbc*nzab - nxab*nzbc
	uz = nxab*nybc - nxbc*nyab

	nx = ux / sqrt(ux**2+uy**2+uz**2)
	ny = uy / sqrt(ux**2+uy**2+uz**2)
	nz = uz / sqrt(ux**2+uy**2+uz**2)

	nxsq=nx*nx;nysq=ny*ny;nzsq=nz*nz;nxy=nx*ny;nxz=nx*nz;nyz=ny*nz

	!write(20,*) 'cross vector',nx,ny,nz,sqrt(nxsq+nysq+nzsq)

!	calculation of matrix 
	
	!dtheta_angle = 1.5708d0
	cos_dh = dcos(dtheta_angle)
	sin_dh = dsin(dtheta_angle)
	onemcosdh = 1 - cos_dh

	!write(20,*) dtheta_angle,cos_dh,sin_dh,onemcosdh

	rot(1,1) = cos_dh + nxsq * onemcosdh
	rot(1,2) = nxy * onemcosdh + nz * sin_dh
	rot(1,3) = nxz * onemcosdh - ny * sin_dh
	
	rot(2,1) = nxy * onemcosdh - nz * sin_dh
	rot(2,2) = cos_dh + nysq*onemcosdh
	rot(2,3) = nyz * onemcosdh + nx * sin_dh
	
	rot(3,1) = nxz * onemcosdh + ny * sin_dh
	rot(3,2) = nyz * onemcosdh - nx * sin_dh
	rot(3,3) = cos_dh + nzsq*onemcosdh

	!write(20,*) rot(:,1)
	!write(20,*) rot(:,2)
	!write(20,*) rot(:,3)

!	calculation of determinant

	detrot=0.0d0
	detrot=rot(1,1)*(rot(2,2)*rot(3,3)-rot(2,3)*rot(3,2))
	detrot=detrot+rot(1,2)*(rot(3,1)*rot(2,3)-rot(2,1)*rot(3,3))
	detrot=detrot+rot(1,3)*(rot(2,1)*rot(3,2)-rot(3,1)*rot(2,2))

	!if(mod(imc,1).eq.0) write(20,*) imc,igmol,detrot

	end subroutine


				!------------------ INHOMOGENEOUS SAMPLING (hot zones are placed perpendicular to x-axis) -----------!


	subroutine gaussian_sampling

        use parameter
        implicit none

        real :: T1,T2,RT1,RT2,pos_hot_zone1,pos_hot_zone2
        
        dE = 0.0d0; prob = 0.0d0; xt = 0.0d0

				!$	transition from (x,y,z) to (x',y',z') by random displacments	$!

!	First transition probability (x,y,z) to (x,y,z')

	dE_z = 0.0d0

	cartesian = 3
	translation = 0
	call samplingenergy_gauss
	Ez_old = E

	translation = 1
	call samplingenergy_gauss
	Ez_new = E

	!write(20,*) '******z***'
	!write(20,*) com_new3,com(3,igmol),dE_z
	!write(20,*) Ez_new,Ez_old

	dE_z = Ez_new - Ez_old
	prob_z = exp ( - dE_z / RThomo ) 
	!write(*,*) RThomo
	call random_number(xt)

	if (dE_z.lt.0.0d0.or.xt.lt.prob_z) com(3,igmol) = com_new3

!	Second transition probability (x,y,z') to (x,y',z')

	dE_y = 0.0d0

	cartesian = 2
	translation = 0
	call samplingenergy_gauss
	Ey_old = E

	translation = 1
	call samplingenergy_gauss
	Ey_new = E

	!write(20,*) '******y***'
	!write(20,*) imc,igmol,com_new2,com(2,igmol),dE_y
	!write(20,*) Ey_new,Ey_old

	dE_y = Ey_new - Ey_old
	prob_y = exp ( - dE_y / RThomo ) 
	call random_number(xt)

	if (dE_y.lt.0.0d0.or.xt.lt.prob_y) com(2,igmol) = com_new2 

!	Third transition probability (x,y',z') to (x',y',z')

	prob_x = 0.0d0

	
	!write(20,*) prob_x
	T1 = temp(int(com(1,igmol)/0.01d0) ); T2 = temp(int(com_new1/0.01d0) )
	if(abs(T1-T2).ge.1) then
	RT1 = R*T1; RT2 = R * T2
	!if(RT2.ne.RT1) then
	!write(20,*) '******x***'
	!write(20,*) imc,igmol,com_new1,com(1,igmol),T2,T1
	call gaussian_prob
	!write(20,*) prob_x,xt,RT2,RT1
	prob_x = RT1/RT2*(exp(prob_x))
	call random_number(xt)

	!write(20,*) '******x***'
	!write(20,*) imc,igmol,com_new1,com(1,igmol),RT2,RT1
	!write(20,*) prob_x,xt,RT2,RT1
	
        if (xt.lt.prob_x) then
        if(com_new1.gt.1.0) then        
         com(1,igmol) = com_new1
        endif
        endif

	else
	
	!write(20,*) '******x***'
	!write(20,*) imc,com_new1,com(1,igmol),RT2,RT1

	dE_x = 0.0d0

	translation = 0
        part = 0
        cartesian = 1
        call samplingenergy
        Ex_old = E

        translation = 1
        call samplingenergy
        Ex_new = E

        dE_x = Ex_new - Ex_old
        prob_x = exp ( - dE_x / RThomo )

	call random_number(xt)

	if (dE_x.lt.0.0d0.or.xt.lt.prob_x) then
        if(com_new1.gt.1.0d0) then
         com(1,igmol) = com_new1
        endif
        endif

	endif
	 

        end subroutine

					!................. SUBROUTINE TO CALCULATE GAUSSIAN PROBABILITY..................!

	subroutine gaussian_prob

	use parameter
	implicit none

	integer*8 :: a0,b0,istep,nstep
	real*8 	  :: gauss_prob(0:1000),delx,Temp_gauss,energy_xh,energy_x

	a0 = com(1,igmol)*100
	b0 = com_new1*100

!	when the random displacement is less than 1E-3 Angstrom, the transition will take place automatically

	if(abs(b0-a0).lt.1E-3) then
	goto 20	
	prob_x = 1
	endif

	delx = 0.01!*1000

	nstep = abs(com(1,igmol)-com_new1)*100

	!write(20,*) nstep,a0,b0

	do istep = 0,nstep
	gauss_prob(istep) = 0
	enddo
	

	do istep = 0,nstep

	T=0.0d0; E=0.0d0; gauss_prob(istep) = 0.0

	com_inter_x = com(1,igmol) + delx*istep
	Temp_gauss = R*temp(int(com_inter_x/0.01d0) )
	call gaussian_prob_potential_energy
	energy_xh = E
	!write(20,*) com_inter_x,E
	com_inter_x = com_inter_x + 0.0001
	call gaussian_prob_potential_energy
	energy_x = E
	!write(20,*) com_inter_x,E
	E = (energy_xh - energy_x)/0.0001
	!write(20,*) E,temp_gauss/R
	gauss_prob(istep) = E/Temp_gauss
	!write(20,*) istep,com_inter_x,delx*istep,temp_gauss,gauss_prob(istep),prob_x

!	using trapezoidal rule to calculate the value

	if(istep.eq.0) then
		prob_x = prob_x + gauss_prob(istep)
	elseif(istep.eq.nstep) then
		prob_x = prob_x + gauss_prob(istep)
	else
		prob_x = prob_x + 2*gauss_prob(istep)
	endif

	!write(20,*) istep,com_inter_x,prob_x,gauss_prob(istep),E,temp_gauss

	enddo

	prob_x = prob_x * 0.01 * 0.5d0	 

20	end subroutine


			!................ SUBROUTINE TO CALCULATE POTETNIAL ENERGY AT INTERMEDIATE STEPS DURING NUMERICAL INTEGRATION ..........!

	subroutine gaussian_prob_potential_energy

	use parameter
	implicit none

	com_p1 = com_inter_x; com_p2 = com(2,igmol); com_p3 = com(3,igmol)
	!write(30,*) com_p1,com_p2,com_p3
	atomic_p(:,:)=rg_atomic(:,:,igmol)

	E = 0.0d0; guestguest_E = 0.0d0; guesthost_E = 0.0d0

	   do j=1,nlistgg(igmol)
		jgmol=listgg(j,igmol)

		dx=0.0d0;dy=0.0d0;dz=0.0d0;dxaa=0.0d0;dyaa=0.0d0;dzaa=0.0d0

		dx = com(1,jgmol) - com_p1;dy = com(2,jgmol) - com_p2;dz = com(3,jgmol) - com_p3

!	applying minimum boundary condition

		call minimumdisplacementcondition

!	calculation of inter & intra molecular energy

		if(gmoltype(igmol).eq.1) nigatom = ngatom1
                if(gmoltype(igmol).eq.2) nigatom = ngatom2

                if(gmoltype(jgmol).eq.1) njgatom = ngatom1
                if(gmoltype(jgmol).eq.2) njgatom = ngatom2

		   do igatom=1,nigatom
		      igatomtype = gatomtype(igatom,igmol)

		      do jgatom=1,njgatom
		      	jgatomtype = gatomtype(jgatom,jgmol)

		      	dxaa =  rg_atomic(1,jgatom,jgmol) + dx - atomic_p(1,igatom)
		      	dyaa =  rg_atomic(2,jgatom,jgmol) + dy - atomic_p(2,igatom)
		      	dzaa =  rg_atomic(3,jgatom,jgmol) + dz - atomic_p(3,igatom)

		      	rsq=dxaa*dxaa+dyaa*dyaa+dzaa*dzaa

		      	ratio = sigmasqgg(igatomtype,jgatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

			guestguest_E = guestguest_E + fourepsgg(igatomtype,jgatomtype) * ( s12 - s6 )

		      	E = E + fourepsgg(igatomtype,jgatomtype) * ( s12 - s6 )
			!write(30,*) igmol,jgmol,igatom,jgatom,E

		      enddo !jatom

		   enddo !iatom

!		endif 

	     enddo !j

	do j=1,nlistgh(igmol)
	     ihatom = listgh(j,igmol)
	     ihatomtype = hatomtype(ihatom)

		if(gmoltype(igmol).eq.1) nigatom = ngatom1
                if(gmoltype(igmol).eq.2) nigatom = ngatom2

		do igatom=1,nigatom
		igatomtype=gatomtype(igatom,igmol)

			dx=0.0d0;dy=0.0d0;dz=0.0d0

			dx = rh_atomic(1,ihatom) - com_p1 - atomic_p(1,igatom)
			dy = rh_atomic(2,ihatom) - com_p2 - atomic_p(2,igatom)
			dz = rh_atomic(3,ihatom) - com_p3 - atomic_p(3,igatom)

!	applying minimum boundary condition

			call minimumdisplacementcondition

			rsq=dx*dx+dy*dy+dz*dz

		      	ratio = sigmasqgh(igatomtype,ihatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

			guesthost_E = guesthost_E + fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 )

		      	E = E + fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 )
	if(ihatomtype.eq.5.and.E.gt.100) then
	!write(30,*) igmol,igatom,ihatom,igatomtype,ihatomtype,E,rsq
	!write(30,*) fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 ),fourepsgh(igatomtype,ihatomtype),sigmasqgh(igatomtype,ihatomtype)
	endif
			!write(30,*) igmol,igatom,ihatom,igatomtype,ihatomtype,E,rsq
			!write(30,*) rh_atomic(1,ihatom),com_p1+atomic_p(1,igatom)
			!write(30,*) rh_atomic(2,ihatom),com_p2+atomic_p(2,igatom)
			!write(30,*) rh_atomic(3,ihatom),com_p3+atomic_p(3,igatom)

		enddo !igatom

	    enddo !ihatom

	

	end subroutine

				!................. SUBROUTINE TO ENERGY OF iTH MOLECULES FOR SAMPLING........................!

	subroutine samplingenergy_gauss

	use parameter
	implicit none

!	real :: start,finish

!	call cpu_time(start)

	if(translation.eq.0) then

	if(cartesian.eq.1) then
		com_p1 = com(1,igmol); com_p2 = com(2,igmol); com_p3 = com(3,igmol)	! p stands for present
	endif
	if(cartesian.eq.2) then
		com_p1 = com(1,igmol); com_p2 = com(2,igmol); com_p3 = com(3,igmol)
	endif
	if(cartesian.eq.3) then
		com_p1 = com(1,igmol); com_p2 = com(2,igmol); com_p3 = com(3,igmol)
	endif

	!write(*,*) translation,cartesian,com_p1,com_p2,com_p3

	endif

	if(translation.eq.1) then

	if(cartesian.eq.1) then
		com_p1 = com_new1; com_p2 = com(2,igmol); com_p3 = com(3,igmol)		! p stands for present
	endif
	if(cartesian.eq.2) then
		com_p1 = com(1,igmol); com_p2 = com_new2; com_p3 = com(3,igmol) 
	endif
	if(cartesian.eq.3) then
		com_p1 = com(1,igmol); com_p2 = com(2,igmol); com_p3 = com_new3 
	endif

	!write(*,*) translation,cartesian,com_p1,com_p2,com_p3

	endif

!	atomic positions

	atomic_p(:,:)=rg_atomic(:,:,igmol)

	E = 0.0d0; guestguest_E = 0.0d0; guesthost_E = 0.0d0

	   do j=1,nlistgg(igmol)
		jgmol=listgg(j,igmol)

		dx=0.0d0;dy=0.0d0;dz=0.0d0;dxaa=0.0d0;dyaa=0.0d0;dzaa=0.0d0

		dx = com(1,jgmol) - com_p1;dy = com(2,jgmol) - com_p2;dz = com(3,jgmol) - com_p3

!	applying minimum boundary condition

		call minimumdisplacementcondition

!	calculation of inter & intra molecular energy

		if(gmoltype(igmol).eq.1) nigatom = ngatom1
		if(gmoltype(igmol).eq.2) nigatom = ngatom1

		if(gmoltype(jgmol).eq.1) njgatom = ngatom1
		if(gmoltype(jgmol).eq.2) njgatom = ngatom2

		   do igatom=1,nigatom
		      igatomtype = gatomtype(igatom,igmol)

		      do jgatom=1,njgatom
		      	jgatomtype = gatomtype(jgatom,jgmol)

		      	dxaa =  rg_atomic(1,jgatom,jgmol) + dx - atomic_p(1,igatom)
		      	dyaa =  rg_atomic(2,jgatom,jgmol) + dy - atomic_p(2,igatom)
		      	dzaa =  rg_atomic(3,jgatom,jgmol) + dz - atomic_p(3,igatom)

		      	rsq=dxaa*dxaa+dyaa*dyaa+dzaa*dzaa

		      	ratio = sigmasqgg(igatomtype,jgatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

			guestguest_E = guestguest_E + fourepsgg(igatomtype,jgatomtype) * ( s12 - s6 )

		      	E = E + fourepsgg(igatomtype,jgatomtype) * ( s12 - s6 )

		      enddo !jatom

		   enddo !iatom

!		endif 

	     enddo !j

	do j=1,nlistgh(igmol)
	     ihatom = listgh(j,igmol)
	     ihatomtype = hatomtype(ihatom)

		if(gmoltype(igmol).eq.1) nigatom = ngatom1
		if(gmoltype(igmol).eq.2) nigatom = ngatom2

		do igatom=1,nigatom
		igatomtype=gatomtype(igatom,igmol)

			dx=0.0d0;dy=0.0d0;dz=0.0d0

			dx = rh_atomic(1,ihatom) - com_p1 - atomic_p(1,igatom)
			dy = rh_atomic(2,ihatom) - com_p2 - atomic_p(2,igatom)
			dz = rh_atomic(3,ihatom) - com_p3 - atomic_p(3,igatom)

!	applying minimum boundary condition

			call minimumdisplacementcondition

			rsq=dx*dx+dy*dy+dz*dz

		      	ratio = sigmasqgh(igatomtype,ihatomtype)/rsq
		      	s6 = ratio * ratio * ratio
		      	s12 = s6 * s6

			guesthost_E = guesthost_E + fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 )

		      	E = E + fourepsgh(igatomtype,ihatomtype) * ( s12 - s6 )

		enddo !igatom

	    enddo !ihatom

!	call cpu_time(finish)

	end subroutine
