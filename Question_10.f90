	!-------------------------------------------------------------------------------------------!
	!									Question 10												!
	!-------------------------------------------------------------------------------------------!
	!###########################################################################################!
	!-------------------------------------------------------------------------------------------!
	!									Main Function											!
	!-------------------------------------------------------------------------------------------!
	real x(100),f(100),fp,vs
	real t_id(100),p_id(100),vd(100)		 ! vd = diode voltage (in volt)
	real ab_err,pr_err(100)
	real is,kb,t,q,eta,id_out
	!-------------------------------------------------------------------------------------------!
	!									Given data												!
	!-------------------------------------------------------------------------------------------!
	is = 1.5e-3								 ! is = reversed biased saturation current (in mA)
	kb = 1.38e-23							 ! kb = boltzmann constant (in JK-1))
	t = 300.0								 ! t = temperature (in K)
	eta = 2.0								 ! n = ideality factor
	q = 1.602e-19							 ! q = electric charge (in C)
	!-------------------------------------------------------------------------------------------!
	open(unit=1,file="theoin.dat")
	write(*,*)"values of vd are :"
	do i=1,3
		read(1,*) vd(i)
		write(*,*) i," value is: ",vd(i)
	enddo
	close(1)
	write(*,*)"----------------------------------------------------------------"
	write(*,*)"----------------------------------------------------------------"
	!-------------------------------------------------------------------------------------------!
	!							calculate theoritical value										!
	!-------------------------------------------------------------------------------------------!
	do i=1,3
		vs=vd(i)
		call id(is,kb,t,q,eta,vs,id_out)
		t_id(i)=id_out
		write(*,*) "At ",vd(i)," theoritical value is: ",t_id(i)
	enddo
	write(*,*)"----------------------------------------------------------------"
	write(*,*)"----------------------------------------------------------------"
	!-------------------------------------------------------------------------------------------!
	!							calculate practical value										!
	!-------------------------------------------------------------------------------------------!
	open(unit=3,file="input.dat")
	read(3,*) n
	do i=1,n
	read(3,*) x(i), f(i)
	enddo
	close(3)
	open(unit=4,file="output.dat")
	do i=1,3
	  vs=vd(i)
	  call lagrangian(n,x,f,vs,fp)
	  p_id(i)=fp
	  write(*,*) "At ",vd(i)," practical value is: ",p_id(i)
	enddo
	write(*,*)"----------------------------------------------------------------"
	!-------------------------------------------------------------------------------------------!
	!							Calculating error in both case									!
	!-------------------------------------------------------------------------------------------!	
	write(*,*)"----------------------------------------------------------------"
	write(*,*)"    "
	open(unit=5,file='final_output.dat')
	write(*,*)"	vd	   theoritical	   practical	Percentage error	"
	write(5,*)"	vd	   theoritical	   practical	Percentage error	"
	do i=1,3
		ab_err=t_id(i)-p_id(i)					! absolute error
		pr_err(i)=abs(ab_err/t_id(i))*100.		! percentage error
		write(*,*) vd(i),t_id(i),p_id(i),pr_err(i)
		write(5,*) vd(i),t_id(i),p_id(i),pr_err(i)
	enddo	
	stop
	end
	!-------------------------------------------------------------------------------------------!
	!									End Main Function										!
	!-------------------------------------------------------------------------------------------!
	!###########################################################################################!
	!-------------------------------------------------------------------------------------------!
	!							subroutine for lagrangian										!
	!-------------------------------------------------------------------------------------------!
	subroutine lagrangian(n,x,f,xp,fp)
	integer n,i,j
	real x(100),f(100),xp,fp,sum,lf
	sum=0.0
	do i=1,n
		lf=1.0
		do j=1,n
			if (i.ne.j)	then
				lf=lf*(xp-x(j))/(x(i)-x(j))

			endif
		enddo
		sum=sum+lf*f(i)
	enddo
	fp=sum
	return
	end
	!-------------------------------------------------------------------------------------------!
	!					Defining function to calculate id (diode current)						!
	!-------------------------------------------------------------------------------------------!
	subroutine id(is,kb,t,q,n,vd,id_out)
	real id_out,vd,is,kb,t,q,n				 ! vd = diode voltage (in volt)
	real sum1, sum2
	sum1=q*vd
	sum2=n*kb*t
	id_out=is*(exp(sum1/sum2)-1.)
	return 
	end