	!-------------------------------------------------------------------------------------------!
	!										Question 7											!
	!-------------------------------------------------------------------------------------------!
	real	m,k,pi,w,t,tp,dt,t1
	real	time(0:1000),pos,t_pos(0:1000),t_pe(0:1000),v(0:1000)
	real	p_pos(0:1000),p_pe(0:1000)!,p_time(0:1000)
	real	err_pos(0:1000),err_pe(0:1000),ae_pos,ae_pe
	real	ke(0:1000)
	integer i,j
	!-------------------------------------------------------------------------------------------!
	m=2.					! mass in kg
	k=2.5					! force constant N/m
	pi=4*atan(1.)			! value of pi
	w=sqrt(k/m)				! angular frequency
	t=2*pi/w				! single time period
	tp=3*t					! total time period
	dt=t/20					! time interval 
	!-------------------------------------------------------------------------------------------!
	!						Calculation of theoretical value									!
	!-------------------------------------------------------------------------------------------!
	i=0
	t1=0.
10	time(i)=t1
	
	t_pos(i)=pos(t1)
	ke(i)=(vel(t1)**2)*m*0.5
	t_pe(i)=.5*m*w*w*(t_pos(i)-1.5)**2
	t1=t1+dt
	i=i+1
	if(t1.lt.tp) goto 10
	!-------------------------------------------------------------------------------------------!
	!						Calculation of practical value										!
	!-------------------------------------------------------------------------------------------!
	!----------------------			Position points				--------------------------------!
	!-------------------------------------------------------------------------------------------!
	t1=0.							! iniial time
	nx=int((tp-t1)/dt)				! no. of positons to be calculated
	open(unit=1,file="output.dat")
	do i=1,nx
		n=2*i
		t2=time(i)
		hv=(t2-t1)/float(n)	
		do j=0,n
			t2=t1+float(j)*hv
			v(j)=vel(t2)
		enddo
		call simpson(n,hv,v,sum)
		p_pos(i)=sum-3.5/w+1.5
	enddo
	call meanpos(p_pos,nx,amean)
	write(*,*) "mean position is: ",amean
	do i=1,nx
		p_pe(i)=.5*k*(p_pos(i)-amean)**2
		write(1,11) time(i),t_pos(i),p_pos(i),t_pe(i),p_pe(i),ke(i),p_pe(i)+ke(i)
	enddo
	close(1)		
	!-------------------------------------------------------------------------------------------!
	!			Calculating percentage error in -- position -- kinetic energy					!
	!-------------------------------------------------------------------------------------------!
	open(unit=2,file="error.dat")
	write(*,*) "	time		t_pos		p_pos		t_pe		p_pe		err_pos		err_pe"
	do i=1,nx
		ae_pos=abs(t_pos(i)-p_pos(i))				! absolute error in position				
		err_pos(i)=abs(ae_pos/t_pos(i))*100.				! percentage error in position
		ae_pe=abs(t_pe(i)-p_pe(i))					! absolute error in position
		err_pe(i)=abs(ae_pe/t_pe(i))*100.				! percentage error in position
		write(2,*) time(i),err_pos(i),err_pe(i)
		write(*,11) time(i),t_pos(i),p_pos(i),t_pe(i),p_pe(i),err_pos(i),err_pe(i)
	enddo
11	format(21(e14.6,3x))
	!-------------------------------------------------------------------------------------------!
	stop
	end

	subroutine meanpos(x,n,a)
	integer n,i
	real x(n),a,minx,maxx
	minx=100000
	maxx=-100000
	do i=1,n
		    if(maxx.lt.x(i)) then
			maxx=x(i)
			endif
			if(minx.gt.x(i)) then
			minx=x(i)
			endif
	enddo
	write(*,*) "maxx is: ",maxx
	write(*,*) "minx is: ",minx
	a=(maxx+minx)/2.
	return
	end
	!-------------------------------------------------------------------------------------------!
	!						Subroutine for Simpsons 1/3 rd rule									!
	!-------------------------------------------------------------------------------------------!
	subroutine simpson(n,h,y,fx)
	real h,fx
	integer n, i
	real y(0:n)
	fx=y(0)+y(n)
	do i=1,n-1
	if(mod(i,2).ne.0.) then
	fx=fx+4.*y(i)
	else
	fx=fx+2.*y(i)
	endif
	enddo
	fx=fx*h/3.
	return
	end
	!-------------------------------------------------------------------------------------------!

	!-------------------------------------------------------------------------------------------!
	!								 velocity function											!
	!-------------------------------------------------------------------------------------------!
	function vel(t)
	real vel,t,w,k,m
	m=2.
	k=2.5
	w=sqrt(k/m)
	vel=3.5*sin(w*t)
	return
	end
	!-------------------------------------------------------------------------------------------!
	!-------------------------------------------------------------------------------------------!
	!							position function (theorical value)								!
	!-------------------------------------------------------------------------------------------!
	function pos(t)
	real pos,t,w,k,m
	m=2.
	k=2.5
	w=sqrt(k/m)
	x=3.5*cos(w*t)
	pos=1.5-x/w
	return
	end
	!-------------------------------------------------------------------------------------------!





