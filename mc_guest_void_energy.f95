!	Fortran Program to extract the data of tetra-tetra, octa-tetra, tetra-octa and octa-octa hopping of solute diffusion in bcc solvent
!	Developer - Shubhadeep Nag

	module parameter

	implicit none
	
	real*8,allocatable :: rg(:,:,:),scale_rg(:,:),rh(:,:,:),rv(:,:),energy_guest(:,:)
	integer,allocatable :: gtraj_void(:,:),void_index(:,:)

	real*8 :: count_tt,count_to,count_oo,count_ot
	real*8 :: energy_tt_lammps(0:1000),energy_tt_cal(0:1000),tt_count(0:1000),energy_ot_cal(0:1000)
	real*8 :: energy_to_cal(0:1000),to_count(0:1000),ot_count(0:1000),energy_oo_cal(0:1000),oo_count(0:1000)
	real*8 :: coord_tt(0:1000),coord_to(0:1000),coord_ot(0:1000),coord_oo(0:1000)
	real*8 :: dx,dy,dz,sqd,d,boxlength,halfboxlength,E,sigma,d_gv,d_gh,d_v,d_vv
	real*8 :: dx1,dx2,dy1,dy2,dz1,dz2,reac,reac_num,const

	integer*8 :: igatom,ngatom,ihatom,nhatom,ivoid,nvoid,itraj,ntraj,jvoid,jump_count(1:588)
	integer*8 :: k,count,ibingv,ibinreac,jtraj,ktraj,ireac,i,j,htraj
	integer*8 :: count_tt_1st,count_to_1st,count_oo_1st,count_ot_1st
	integer*8 :: count_tt_2nd,count_to_2nd,count_oo_2nd,total_count,count_hatom,kgatom,ibin1,ibin2,count_dist(1:3000,1:3000)
	integer*8 :: coordination,n


	character(4),allocatable :: c(:)

	end module parameter

	program bcc

	use parameter
	implicit none

!	open files to read and write
	
	open(unit=2,file='guest_1p1.xyz')
	open(unit=3,file='dump_pe.guest_1p1')
	open(unit=4,file='void_bcc.xyz')
	open(unit=8,file='host_1p1.xyz')
	open(unit=11,file='count_jump.dat',access='append')
	open(unit=17,file='tt_to_ot_oo_energy_1p1_rmin.dat')
	open(unit=18,file='tt_to_ot_oo_coord_1p1_rmin.dat')

!	read the trajectory points

	ntraj = 20000; ngatom = 588; sigma = 1.8; nhatom = 5488
	boxlength = 88.9d0; halfboxlength = 0.5 * boxlength
	count_tt = 0; count_to = 0; count_ot = 0; count_oo = 0

    do ibinreac = 0,1000
        tt_count(ibinreac) =0; to_count(ibinreac) = 0;ot_count(ibinreac) = 0
        oo_count(ibinreac) =0; energy_tt_lammps(ibinreac) = 0; energy_tt_cal(ibinreac) = 0; tt_count(ibinreac) = 0
		coord_tt(ibinreac) = 0; coord_to(ibinreac) = 0; coord_ot(ibinreac) = 0; coord_oo(ibinreac) = 0
        energy_to_cal(ibinreac) = 0; to_count(ibinreac) = 0; energy_oo_cal(ibinreac) = 0; oo_count(ibinreac) = 0
	enddo !ibingv

	do ibin1 = 1,3000
	do ibin2 = 1,3000
	count_dist(ibin1,ibin2) = 0
	enddo
	enddo

!	read guest positions

	allocate ( rg(1:3,1:ngatom,1:ntraj),scale_rg(1:3,1:ngatom) )

	count = 0

	do itraj = 1,ntraj
       ! write(*,*) itraj
		read(2,*) 
		read(2,*)

		count = count + 1

		do igatom = 1,ngatom
             
			read(2,*) k,rg(:,igatom,itraj)

		enddo !iatom

	enddo !itraj

!	read host positions

	allocate ( rh(1:3,1:nhatom,1:ntraj) )

	do itraj = 1,ntraj
       ! write(*,*) itraj
	read(8,*) 
	read(8,*)

	do ihatom = 1,nhatom
		
		read(8,*) k,rh(:,ihatom,itraj)

	enddo !ihatom

	enddo !itaj

!	read energy from lammps output

	allocate ( energy_guest(1:ngatom+nhatom,1:ntraj) )

	do itraj = 1,ntraj
      !  write(*,*) itraj 
	do i=1,9
	read(3,*)
	enddo

	do igatom = 1,ngatom
       
	read(3,*) kgatom,energy_guest(kgatom,itraj)
	enddo !itraj

	enddo !iatom

!	read void center positions

	read(4,*) nvoid
	read(4,*)
        
	allocate ( rv(1:3,1:nvoid),c(1:nvoid) )

	do ivoid = 1,nvoid
	!	write(*,*) ivoid
		read(4,*) c(ivoid),rv(:,ivoid)

	enddo !ivoid

!	locate the guest inside or outside the void

	allocate ( gtraj_void(1:ngatom,1:ntraj), void_index(1:ngatom,1:nvoid) )

        do igatom = 1,ngatom
        do itraj = 1,ntraj
        gtraj_void(igatom,itraj) = 0
        enddo
        do ivoid = 1,nvoid
        void_index(igatom,ivoid) = 0
        enddo 
        enddo


	c(0) = 'A'	!	assigninng a dummy character

	do i = 1,ngatom		!	initialize the variables			
     	!write(*,*) i
	jump_count(i) = 0
	enddo

	do igatom = 1,ngatom
		write(*,*) igatom
		do itraj = 1,ntraj
			!write(*,*) itraj
			do ivoid = 1,nvoid
				
				dx = 0.0d0; dy = 0.0d0; dz = 0.0d0; sqd = 0.0d0; d_gv = 0.0d0
				dx = rg(1,igatom,itraj) - rv(1,ivoid)
				dy = rg(2,igatom,itraj) - rv(2,ivoid)
				dz = rg(3,igatom,itraj) - rv(3,ivoid)

				call minimumimageconvention

				sqd = ( dx*dx + dy*dy + dz*dz )
				d_gv = sqrt(sqd)
				
				if(d_gv.lt.0.8d0) then						

					if(c(ivoid)=='T') then
						jump_count(igatom) = jump_count(igatom) + 1
						gtraj_void(igatom,jump_count(igatom)) = itraj
						void_index(igatom,jump_count(igatom)) = ivoid
                           !                     write(30,*) igatom,jump_count(igatom), &
                            !    gtraj_void(igatom,jump_count(igatom)),void_index(igatom,jump_count(igatom)),itraj,ivoid
					endif

					if(d_gv.lt.0.42d0) then
						if(c(ivoid)=='O') then
							jump_count(igatom) = jump_count(igatom) + 1
							gtraj_void(igatom,jump_count(igatom)) = itraj
							void_index(igatom,jump_count(igatom)) = ivoid
						endif
					endif

				endif

				

			enddo !ivoid

		enddo !itraj

	enddo !igatom

!	CALCULATE THE NUMBER OF JUMPS


	total_count=0;count_tt=0;count_to=0;count_ot=0;count_oo=0
        
        do igatom = 1,ngatom
        write(*,*) igatom
        do i=1,jump_count(igatom)-1
        j = i+ 1
        ivoid = void_index(igatom,i)
        jvoid = void_index(igatom,j)
        itraj = gtraj_void(igatom,i)
        jtraj = gtraj_void(igatom,j)


	if(ivoid.ne.jvoid) then
        if(itraj.ne.jtraj) then

					!!!!!!!!!!!	calculate the TETRA-TETRA jump		!!!!!!!!!!!!!

	if(c(ivoid)=='T'.and.c(jvoid)=='T') then

		call voidvoiddistance
		d_vv = d
		dx1 = dx; dy1 = dy; dz1 = dz	
		!write(*,*) d,dx1,dy1,dz1
		
		if(d_vv.ge.0.1) then 	! to avoid same void jump computation
			
		if (d_vv.le.2.25d0) then

			count_tt = count_tt + 1	
		        ktraj = 0	
			itraj = gtraj_void(igatom,i)
			jtraj = gtraj_void(igatom,j)
			
			do ktraj = itraj,jtraj
				
				E = energy_guest(igatom,ktraj)
				

!	computation of reaction coordinate
				call guestvoiddistance
				dx2 = dx; dy2 = dy; dz2 = dz
				
				reac_num = dx1*dx2 / abs(d_vv) + dy1*dy2 / abs(d_vv) + dz1*dz2 / abs(d_vv)
				reac = abs(reac_num) / abs(d_vv)
				

				ibinreac = nint(reac*100.0d0)
                               
				tt_count(ibinreac) = tt_count(ibinreac) + 1
				energy_tt_lammps(ibinreac) = energy_tt_lammps(ibinreac) + E
				
!	compute the energy
				call guestenergy
				!write(*,*) E
				energy_tt_cal(ibinreac) = energy_tt_cal(ibinreac) + E
				coord_tt(ibinreac) = coord_tt(ibinreac) + ihatom
				
			write(60,*) igatom,ibinreac,ihatom	
			enddo !ktraj

		endif

		endif	! same void
			 
	endif
		
					!!!!!!!!!!!	calculate the TETRA-OCTA jump		!!!!!!!!!!!!!

	if(c(ivoid)=='T'.and.c(jvoid)=='O') then
		
		call voidvoiddistance
		d_vv = d
		dx1 = dx; dy1 = dy; dz1 = dz
		
		if(d_vv.ge.0.1) then 	! to avoid same void jump computation
			
		if (d_vv.lt.1.59d0) then

			count_to = count_to + 1	
			
			itraj = gtraj_void(igatom,i)
			jtraj = gtraj_void(igatom,j)
			
                       do ktraj = itraj,jtraj

                                E = energy_guest(igatom,ktraj)

!       computation of reaction coordinate
                                call guestvoiddistance
                               dx2 = dx; dy2 = dy; dz2 = dz
                                reac_num = dx1*dx2 / abs(d_vv) + dy1*dy2 / abs(d_vv) + dz1*dz2 / abs(d_vv)
                                reac = abs(reac_num) / abs(d_vv)

                                ibinreac = nint(reac*100.0d0)
                                to_count(ibinreac) = to_count(ibinreac) + 1
                           
!       compute the energy
                                call guestenergy
                                energy_to_cal(ibinreac) = energy_to_cal(ibinreac) + E
                                coord_to(ibinreac) = coord_to(ibinreac) + ihatom
                                write(70,*) igatom,ibinreac,ihatom               
                        enddo !ktraj

		endif

	        endif	! same void
	endif

					!!!!!!!!!!!	calculate the OCTA-TETRA jump		!!!!!!!!!!!!!

	if(c(ivoid)=='O'.and.c(jvoid)=='T') then
		
		call voidvoiddistance
		d_vv = d
		dx1 = dx; dy1 = dy; dz1 = dz
		
		if(d_vv.ge.0.1) then 	! to avoid same void jump computation
			
		if (d_vv.lt.1.59d0) then

			count_ot = count_ot + 1	
			
			itraj = gtraj_void(igatom,i)
			jtraj = gtraj_void(igatom,j)

                       do ktraj = itraj,jtraj

                                E = energy_guest(igatom,ktraj)

!       computation of reaction coordinate
                                call guestvoiddistance
                               dx2 = dx; dy2 = dy; dz2 = dz
                                reac_num = dx1*dx2 / abs(d_vv) + dy1*dy2 / abs(d_vv) + dz1*dz2 / abs(d_vv)
                                reac = abs(reac_num) / abs(d_vv)

                                ibinreac = nint(reac*100.0d0)
                                ot_count(ibinreac) = ot_count(ibinreac) + 1
                           
!       compute the energy
                                call guestenergy
                                energy_ot_cal(ibinreac) = energy_ot_cal(ibinreac) + E
                                coord_ot(ibinreac) = coord_ot(ibinreac) + ihatom
                                
                        enddo !ktraj

		endif

        	endif	! same void
	endif

					!!!!!!!!!!!	calculate the OCTA-OCTA jump		!!!!!!!!!!!!!

	if(c(ivoid)=='O'.and.c(jvoid)=='O') then
		
		call voidvoiddistance
		d_vv = d
		dx1 = dx; dy1 = dy; dz1 = dz		
		

		if(d_vv.ge.0.1) then 	! to avoid same void jump computation
			
		if (d_vv.lt.3.17d0) then
						
			count_oo = count_oo + 1

			itraj = gtraj_void(igatom,i)
			jtraj = gtraj_void(igatom,j)

			do ktraj = itraj,jtraj

                              E = energy_guest(igatom,ktraj)

!       computation of reaction coordinate
                                call guestvoiddistance
                               dx2 = dx; dy2 = dy; dz2 = dz
                                reac_num = dx1*dx2 / abs(d_vv) + dy1*dy2 / abs(d_vv) + dz1*dz2 / abs(d_vv)
                                reac = abs(reac_num) / abs(d_vv)

                                ibinreac = nint(reac*100.0d0)
                                oo_count(ibinreac) = oo_count(ibinreac) + 1
                           
!       compute the energy
                                call guestenergy
                                energy_oo_cal(ibinreac) = energy_oo_cal(ibinreac) + E
                                coord_oo(ibinreac) = coord_oo(ibinreac) + ihatom
                                
                        enddo !ktraj

		endif

		endif	! same void

	endif

	endif

        endif

	enddo !ivoid

	enddo !iatom


	write(11,*) sigma, count_tt, count_to, count_ot, count_oo,count_tt+count_to+count_ot+count_oo 

	do ibinreac = 0,100
!	multiply with 4.2 to make it Kjoule from Kcal
		write(17,*) ibinreac*0.01,energy_tt_cal(ibinreac)*4.2/tt_count(ibinreac),energy_to_cal(ibinreac)*4.2/to_count(ibinreac) &
		,energy_ot_cal(ibinreac)*4.2/ot_count(ibinreac),energy_oo_cal(ibinreac)*4.2/oo_count(ibinreac)

		write(18,*) ibinreac*0.01,coord_tt(ibinreac)/tt_count(ibinreac),coord_to(ibinreac)/to_count(ibinreac) &
			,coord_ot(ibinreac)/ot_count(ibinreac),coord_oo(ibinreac)/oo_count(ibinreac)
	enddo! ibinreac

	end
	

!********************	compute the distance between two void center ***************!


	subroutine voidvoiddistance

	use parameter
	implicit none 

	ivoid = void_index(igatom,i)
	jvoid = void_index(igatom,j)

	dx = 0; dy = 0; dz= 0
	dx = rv(1,ivoid) - rv(1,jvoid); dy = rv(2,ivoid) - rv(2,jvoid); dz = rv(3,ivoid) - rv(3,jvoid)
	
	call minimumimageconvention

	sqd = dx*dx + dy*dy + dz*dz
	d = sqrt(sqd)

	!write(*,*) rv(:,ivoid),rv(:,jvoid),d

	end subroutine
	
!********************	compute the distance between guest from void center ***************!

	subroutine guestvoiddistance

	use parameter
	implicit none 
		
	dx = 0; dy = 0; dz= 0
	dx = rv(1,ivoid) - rg(1,igatom,ktraj); dy = rv(2,ivoid) - rg(2,igatom,ktraj); dz = rv(3,ivoid) - rg(3,igatom,ktraj)

	call minimumimageconvention

	end subroutine
	
!********************	compute the distance using Minimum Image Convection ***************!

	subroutine minimumimageconvention

	use parameter
	implicit none

	if(dx.gt.halfboxlength) dx = dx - boxlength
	if(dx.lt.-halfboxlength) dx = dx + boxlength

	if(dy.gt.halfboxlength) dy = dy - boxlength
	if(dy.lt.-halfboxlength) dy = dy + boxlength

	if(dz.gt.halfboxlength) dz = dz - boxlength
	if(dz.lt.-halfboxlength) dz = dz + boxlength

	end subroutine

!********************	compute the interaction energy of the solute moving from one void to another void center ***************!

	subroutine guestenergy

	use parameter
	implicit none

	real*8,allocatable :: guest_host_ene(:)
	real*8 :: temp,sum
	integer*8 :: jhatom

	allocate ( guest_host_ene(1:nhatom+ngatom) )

	do ihatom = 1,nhatom
	guest_host_ene(ihatom) = 0	
	enddo

	E = 0; coordination = 0
	do ihatom = 1,nhatom
		
		guest_host_ene(ihatom) = 0
		dx = 0.0d0; dy = 0.0d0; dz = 0.0d0; sqd = 0.0d0; d_gh = 0.0d0					
		dx = rg(1,igatom,ktraj) - rh(1,ihatom,ktraj)
		dy = rg(2,igatom,ktraj) - rh(2,ihatom,ktraj)
		dz = rg(3,igatom,ktraj) - rh(3,ihatom,ktraj)

		call minimumimageconvention

		sqd = ( dx*dx + dy*dy + dz*dz )
		d_gh = sqrt(sqd)
		if(d_gh.le.2.85d0) then
		coordination = coordination + 1
		guest_host_ene(coordination) = 4*0.717*( (sigma/d_gh)**12 - (sigma/d_gh)**6 )
		E = E + guest_host_ene(coordination)
		endif

	enddo !ihatom

!	arrange the interaction energy in descending order

	do 10 ihatom = 1,coordination-1
	   do 30 jhatom = ihatom + 1, coordination
		if(guest_host_ene(ihatom).gt.guest_host_ene(jhatom)) then
	 	temp=guest_host_ene(ihatom)
       		guest_host_ene(ihatom)=guest_host_ene(jhatom)
       		guest_host_ene(jhatom)=temp
		endif
30	continue
10	continue

	sum = 0.0
	do ihatom = 1,coordination
		sum = guest_host_ene(ihatom) + sum
	enddo 

	ihatom = ihatom - 1

20	end subroutine
