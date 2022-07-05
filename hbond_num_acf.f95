!	Program to compute the number of hydrogen bonds, continuous HB and intermittent HB correlation function
!	Date - 05/07/22

!	Declare Variables

	module parameter
	
	real,allocatable :: x(:,:),y(:,:),z(:,:),O(:,:),H2(:,:), hb_intermittent(:), hb_continuous(:), Chb(:), Shb(:), count_hb(:)
	real :: dOOx,dOOy,dOOz,dOHx,dOHy,dOHz, dOO, dOH, dOO_OH, theta, start, finish
	real :: boxlengthx,boxlengthy,boxlengthz,halflengthx,halflengthy,halflengthz
	integer :: itraj, ntraj, imol, jmol, nmol, iatom, natom, icorr, ncorr
	integer,allocatable :: h(:,:),  count_corr(:),  HB(:,:)
	character :: a
	
	end module parameter

	program HBcorreltaion

	use parameter
	implicit none
	
!	Open Files to Read & Write
	
	open (1,file='input_HB.dat')
	open (2,file='npt_500fs.xyz')
	open (3,file='npt_sort.xyz')
	open (4,file='hblist_OO.dat')
	open (6,file='hblist_OH.dat')
	open (7,file='hblist_OOH.dat')
	open (9,file='hblist_mol.dat')
	open (10,file='hbnum.dat')
	open (11,file='hbac.dat')
	
	call cpu_time(start)	! start calculating the time
	
!	Read the input variables

	read(1,*)
	read(1,*)
	read(1,*) ntraj, nmol, natom, ncorr
	read(1,*)
	read(1,*)
	read(1,*) boxlengthx, boxlengthy, boxlengthz	

!	Write the headings of each output files
!	writing data to file no. 3, 4, 6, 7 and 9 are optional, although necessary to check the program	

	call writeheadingofoutputfiles
	
	

!	calculate the half of the boxlength	
	halflengthx = 0.5*boxlengthx; halflengthy = 0.5*boxlengthy; halflengthz = 0.5*boxlengthz
	
!	Allocate the variables
	
	allocate ( x(1:natom,1:nmol),y(1:natom,1:nmol),z(1:natom,1:nmol),O(1:3,1:nmol),H2(1:3,1:nmol) )
	allocate ( h(1:nmol,0:ntraj), HB(1:nmol,0:ntraj), hb_intermittent(1:ntraj), count_corr(1:ntraj),  count_hb(1:ntraj) )
	allocate ( hb_continuous(1:ntraj), Chb(1:ntraj), Shb(1:ntraj) )
	
!	Initialise the variables
		
	do itraj = 1,ntraj
		hb_intermittent(itraj) = 0; hb_continuous(itraj) = 0; count_corr(itraj) = 0; Chb(itraj) = 0; Shb(itraj)  = 0
	enddo !itraj

!	Read the trajectories
	
	do itraj = 1,ntraj	!	ntraj = no. of trajectories
	
	read(2,*)
	read(2,*)
	
	do imol = 1,nmol	!	nmol = no. of molecules
		
		do iatom = 1,natom	!	natom = no. of atoms
		
			read(2,*) a,x(iatom,imol),y(iatom,imol),z(iatom,imol)

!	Store the oxygen and hydrogen atoms of OH group of MeOH
			
			if (iatom.eq.2) then
				O(1,imol) = x(iatom,imol); O(2,imol) = y(iatom,imol); O(3,imol) = z(iatom,imol)  
			endif
			
			if (iatom.eq.6) then
				H2(1,imol) = x(iatom,imol); H2(2,imol) = y(iatom,imol); H2(3,imol) = z(iatom,imol)  
			endif
			
		enddo !iatom
		
!	Write the hydrogen and oxygen atoms of all molecules at every trajectory
		
		write(3,*) itraj,imol,'	O',O(:,imol)
		write(3,*) itraj,imol,'	H2',H2(:,imol)
		
	enddo !imol
	
!	Start Calculating the Hydrogen bond
	
	do imol = 1,nmol-1
	
		do jmol = imol+1,nmol
		
			dOOx = O(1,jmol) - O(1,imol); dOOy = O(1,jmol) - O(1,imol); dOOz = O(1,jmol) - O(1,imol)
			
			if (dOOx.gt.halflengthx)  dOOx = dOOx - boxlengthx
			if (dOOx.lt.-halflengthx)  dOOx = dOOx + boxlengthx
			
			if (dOOy.gt.halflengthy)  dOOy = dOOy - boxlengthy
			if (dOOy.lt.-halflengthy)  dOOy = dOOy + boxlengthy
			
			if (dOOz.gt.halflengthz)  dOOz = dOOz - boxlengthz
			if (dOOz.lt.-halflengthz)  dOOz = dOOz + boxlengthz
	
			dOO = sqrt (dOOx*dOOx+dOOy*dOOy+dOOz*dOOz)
			
			if(dOO.le.3.5) then	!	First minima of RDF of Oxygen-Oxygen atoms of MeOH = 3.5 Angs

!	Write the molecules (A, B) at every time frame where the distance O_{A}O_{B} is less than 3.5 Angs 
			
				write(4,*) itraj, imol,jmol, dOO
				
				dOHx = H2(1,imol) - O(1,jmol); dOHy = H2(2,imol) - O(2,jmol); dOHz = H2(3,imol) - O(3,jmol)
				
				if (dOHx.gt.halflengthx)  dOHx = dOHx - boxlengthx
				if (dOHx.lt.-halflengthx)  dOHx = dOHx + boxlengthx
			
				if (dOHy.gt.halflengthy)  dOHy = dOHy - boxlengthy
				if (dOHy.lt.-halflengthy)  dOHy = dOHy + boxlengthy
			
				if (dOHz.gt.halflengthz)  dOHz = dOHz - boxlengthz
				if (dOHz.lt.-halflengthz)  dOHz = dOHz + boxlengthz
	
				dOH = sqrt (dOHx*dOHx+dOHy*dOHy+dOHz*dOHz)
				
				if(dOH.le.2.45) then	!	First minima of RDF of Oxygen-Hydrogen atoms of MeOH = 2.45 Angs

!	Write the molecules (A, B) at every time frame where the distance H_{A}O_{B} is less than 2.45 Angs
				
					write(6,*) itraj, imol,jmol, dOO, dOH

					dOHx = 0; dOHy = 0; dOHz = 0

					dOHx = H2(1,imol) - O(1,imol); dOHy = H2(2,imol) - O(2,imol); dOHz = H2(3,imol) - O(3,imol)
				
					if (dOHx.gt.halflengthx)  dOHx = dOHx - boxlengthx
					if (dOHx.lt.-halflengthx)  dOHx = dOHx + boxlengthx
			
					if (dOHy.gt.halflengthy)  dOHy = dOHy - boxlengthy
					if (dOHy.lt.-halflengthy)  dOHy = dOHy + boxlengthy
			
					if (dOHz.gt.halflengthz)  dOHz = dOHz - boxlengthz
					if (dOHz.lt.-halflengthz)  dOHz = dOHz + boxlengthz
				
					dOH = sqrt (dOHx*dOHx+dOHy*dOHy+dOHz*dOHz)

					dOO_OH = (dOHx*dOOx + dOHy*dOOy + dOHz*dOOz)
					theta = acos (dOO_OH / (dOO * dOH)) * 57.29	! calculate degree using dot product and convert theta to degree
					
					if(theta.le.30) then	! critical angle value = 30 degree

!	Write the molecules at every time frame where the angle HO_{A}O_{B} is than 30 degree					
					
						write(7,*) itraj, imol, jmol, dOO, dOH, theta
						count_hb(itraj) = count_hb(itraj) + 1	! count hydrogen bond
						h(imol,itraj) = 1						! define h(t) function of every molecule at every timeframe
					else
						h(imol,itraj) = 0						! define h(t) function of every molecule at every timeframe

					endif					
					
				endif
				
			endif
			
		enddo !jmol

!	Calculate the function H(t)

!	at the first trajectory frame, if there is A hydrogen bond between A and B molecule, assign H(t) = 1
!	at the first trajectory frame, if there is NO hydrogen bond between A and B molecule, assign H(t) = 0		

		if(itraj.eq.1) then
			if(h(imol,itraj).eq.1) HB(imol,itraj) = 1
			if(h(imol,itraj).eq.0) HB(imol,itraj) = 0
		endif
		
!	for trajectory frames, t >= 2, first check the value of H(t-1)
!	if the value of H(t-1) = 1; check further the value of h(t); if h(t) = 1, then put H(t) = 1; but if h(t) = 0, put H(t) = 0
!	if the value of H(t-1) = 0; check no further and put H(t) = 0	
		
		if(itraj.ge.2) then
			if(HB(imol,itraj-1).eq.1) then
				if(h(imol,itraj).eq.1) then
					HB(imol,itraj) = 1
				else
					HB(imol,itraj) = 0
				endif
			else
				HB(imol,itraj) = 0
			endif
		endif
				
!	Store the value of h(t) and H(t) of every molecule at every timeframe
		
		write(9,*) itraj,imol,h(imol,itraj), HB(imol,itraj)
		
	enddo !imol

!	Write the average number of hydrogen bonds at every timeframe
	
	count_hb(itraj) = count_hb(itraj)*1.0/ nmol
	write(10,*) itraj, count_hb(itraj)
	
	enddo !itraj	
	
!	End reading the trajectory	

!	Start Calculating the correlation function

	do icorr = 1,ncorr	! loop between 1 to correlation time

	   do itraj = 1,ntraj		! loop between 1 and number of trajectory
	   		if (itraj+icorr.le.ntraj) then

				do imol = 1,nmol

!	Calculation of the intermittent HB correlation function C_{HB}(t)    

					hb_intermittent(icorr) = hb_intermittent(icorr) + ( h(imol,itraj+icorr)*h(imol,itraj) * 1.0 )
		
!	Calculation of the continuous HB correlation function S_{HB}(t)

					hb_continuous(icorr) = hb_continuous(icorr) + ( HB(imol,itraj+icorr)*h(imol,itraj) * 1.0 )

!	Count the event

					count_corr(icorr) = count_corr(icorr) + 1
	
				enddo !imol

			endif

!	Compute C_{HB}(t) and S_{HB}(t)
		
			Chb(icorr) = hb_intermittent(icorr) * 1.0 / ( count_hb(icorr) * count_corr(icorr) )
			Shb(icorr) = hb_continuous(icorr) * 1.0 / ( count_hb(icorr) * count_corr(icorr) )

		enddo !itraj
		
		Chb(icorr) = Chb(icorr) / Chb(1)
		Shb(icorr) = Shb(icorr) / Shb(1)

!	Writing the functions C_{HB}(t) and S_{HB}(t) after normalization

		write(11,*) icorr,Chb(icorr),Shb(icorr)

	enddo !icorr

    call cpu_time(finish)	! calculate the final time
    write(11,*)
	write(11,*) 'Time required to run this program',finish-start,'	(sec)'
	
	end program HBcorreltaion
	
	subroutine	writeheadingofoutputfiles
	
	use parameter
	implicit none
	
	write(3,*) 'The hydrogen and oxygen atoms of all molecules at every trajectory obtained from the MD simulation'
	write(3,*)
	write(3,*) 'no. of trajectory at every 10 fs','	molecule ID',' X, Y and Z Coordinates of Oxygen and Hydrogen of OH group'
	write(3,*)
	
	write(4,*) 'Molecules IDs (A, B) at every time frame where the distance O_{B}O_{A} is less than 3.5 Angs '
	write(4,*)
	write(4,*) 'no. of trajectory at every 10 fs','	molecule ID A','	molecule ID B','	distance O_{B}O_{A} (A)'
	write(4,*)

	write(6,*) 'Molecules IDs (A, B) at every time frame where the distance H_{A}O_{B} is less than 2.45 Angs '
	write(6,*)
	write(6,*) 'no. of trajectory at every 10 fs','	molecule ID A','	molecule ID B','	distance O_{B}O_{A} (A)', &
	'	distance H_{A}O_{B} (A)'
	write(6,*)

	write(7,*) 'Molecules IDs (A, B) at every time frame where the angle H_{A}O_{A}O_{B} is less than 30 degree'
	write(7,*)
	write(7,*) 'no. of trajectory at every 10 fs','	molecule ID A','	molecule ID B','	distance O_{B}O_{A} (A)', &
	'	distance H_{A}O_{B} (A)','	Theta of H_{A}O_{A}O_{B} (degree)'
	write(7,*)
	
	write(9,*) 'The value of h(t) and H(t) of every molecule at every timeframe'
	write(9,*)
	write(9,*) 'no. of trajectory at every 10 fs','	molecule ID','	h(t)','	H(t)',' Reference: A. Chandra PRL. 85, 768, 2000'
	write(9,*)
	
	write(10,*) 'The average number of hydrogen bonds at every timeframe'
	write(10,*)
	write(10,*) 'no. of trajectory at every 10 fs','	Avg. Hydrogen Bond'
	write(10,*)	
	
	write(11,*)	'Input variables to calcualte HB bond'
	write(11,*)	'No. of trajectory','	No. of Molecules','	No. of Atoms of Each Molecule','	No. of Correlation'
	write(11,*)	ntraj, nmol, natom, ncorr
	write(11,*)
	write(11,*)	'	Boxlength along x		y		z in Angs'
	write(11,*)	boxlengthx, boxlengthy, boxlengthz
	write(11,*)
	write(11,*) 'The functions C_{HB}(t) and S_{HB}(t) at every correlation time'
	write(11,*)
	write(11,*) 'correlation time (x 10 fs)','	C_{HB}(t)','	S_{HB}(t)'
	write(11,*)	
	
	end subroutine
