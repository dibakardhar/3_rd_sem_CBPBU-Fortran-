	!-------------------------------------------------------------------------------------------!
	!										Question 8											!
	!-------------------------------------------------------------------------------------------!	
	!-------------------------------------------------------------------------------------------!
	!###############################	   Main function		################################!
	!-------------------------------------------------------------------------------------------!	
	!-------------------------------------------------------------------------------------------!
	!									Declare all variables									!
	!-------------------------------------------------------------------------------------------!	
	real f,w,pi,tp,pd,h,p
	real mass
	real x(100),y(100)
	real xx(6),yy(6)
	real dely(6,6),sum1(6),sum
	real dy(100), Tdy(100), erv
	real ke,Tke, erke, dh
	integer n,i,j,np,nt,ninte
	!-------------------------------------------------------------------------------------------!
	!							All given data from Question 8									!
	!-------------------------------------------------------------------------------------------!	
	mass=15						 ! mass of the particle
	nt=4						 ! number of time period
	ninte=20					 ! number of time interval in 1 tp
	pi=4.*atan(1.0)				 ! pi value
	w=15.						 !omega frequency
	pd=(2.*pi)/w				 ! time period
	tp=nt*pd					 ! total time
	n=nt*ninte					 ! total no. of interval
	h=tp/float(n)				 ! time step
	np=4						 ! number of points used to obtain derivative
	!-------------------------------------------------------------------------------------------!
	!								Create input Data files										!
	!-------------------------------------------------------------------------------------------!	
	open(unit=1,file="input.dat")
	do i=1,n+np
	x(i)=float(i-1)*h
	y(i)=f(x(i))
	write(1,*) x(i), y(i)
	enddo
	close(1)
	do i=1,n
	dh=h/2.5
	sum=x(i)-dh
	do j=1,np
	xx(j)=sum+(j-1)*dh
	yy(j)=f(xx(j))
	enddo
	call tab(np,yy,dely)
	p=(x(i)-xx(1))/h  
	call diff(p,np,sum1)
	call cal_dy(np,sum1,dely,dh,dy(i))
	Tdy(i)=15*6.5*cos(15*x(i))				! velocity array
	enddo
	!-------------------------------------------------------------------------------------------!
	!								Calculating Output file:									!
	!-------------------------------------------------------------------------------------------!	
	!-------------------------------------------------------------------------------------------!
	!										velocity											!
	!-------------------------------------------------------------------------------------------!	
	open(unit=2,file="velocity.dat")
	open(unit=3,file="KE.dat")
	open(unit=4,file='all_data.dat')
	do i=1,n
		erv=(Tdy(i)-dy(i))/Tdy(i)*100
		write(2,11) x(i),y(i),dy(i),Tdy(i),erv
		ke=mass*dy(i)*dy(i)/2.
		Tke=mass*Tdy(i)*Tdy(i)/2.
		erke=(Tke-ke)/Tke*100.
		write(3,12) x(i),y(i),ke,Tke,erke
		write(4,12) x(i),y(i),dy(i),Tdy(i),erv,ke,Tke,erke
	enddo
11	format(7(f14.7,x))
12	format(21(e21.7,x))
	write(*,*) "All done plot Data.... :)"
	stop
	end
	!-------------------------------------------------------------------------------------------!
	!##########################		 End Main Function		####################################!
	!-------------------------------------------------------------------------------------------!	
	!-------------------------------------------------------------------------------------------!
	!							Subroutine to calculate derivative of y							!
	!-------------------------------------------------------------------------------------------!	
	subroutine cal_dy(np,sum1,dely,h,dy)
	integer np
	real sum1(np), dely(np,np),h,dy , sum2
	integer n , i
	sum1=0.
	n=np-1
	sum2=0.
	do i=1,n
	!sum2=sum2+sum1(i)*dely(i,1)
	sum2=sum2+float((-1)**(n-1))*dely(i,1)/float(i)

	enddo
	dy=sum2/h
	return
	end
	!-------------------------------------------------------------------------------------------!
	!							Subroutine to calculate Differentiate							!
	!-------------------------------------------------------------------------------------------!	
	subroutine diff(p,np,sum1)
	real p,sum1(6),product
	integer i,j,np,n
	n=np-1
	product =1

	do k=1,n
	sum1(k)=0
	if(k.ne.1) then
		do i=0,k-1
			do j=0,k-1
			  if (i.ne.j) then
				product=product*(p-j)
				endif
			enddo
	sum1(k)=sum1(k)+product
	enddo
	else
	sum1(k)=1
	endif
	sum1(k)=sum1(k)/fact(k)

	enddo
	return
	end
	!-------------------------------------------------------------------------------------------!
	!					Subroutine to calculate forward difference table						!
	!-------------------------------------------------------------------------------------------!	
	subroutine tab(n,y,dy)
	integer n, i,j
	real dy(n,n),y(100)

	do i=1,n-1
	do j=1,n-i
	  if(i.eq.1) dy(i,j)=y(j+1)-y(j)
	if(i.ne.1) dy(i,j)=dy(i-1,j+1)-dy(i-1,j)
	enddo
	enddo
	return
	end


	!-------------------------------------------------------------------------------------------!
	!								Function for Acceleration									!
	!-------------------------------------------------------------------------------------------!	
	function f(t)
	real f,t
	f=6.5*sin(15.*t)
	return
	end
	!-------------------------------------------------------------------------------------------!
	!							Function to calculate Factorial									!
	!-------------------------------------------------------------------------------------------!	
	function fact(n)
	real fact, product
	integer n,i
	product=1
	do i=1,n
	product=product*float(i)
	enddo
	fact=product
	return
	end
	!-------------------------------------------------------------------------------------------!
	!##################################		End Program		####################################!
	!-------------------------------------------------------------------------------------------!	
