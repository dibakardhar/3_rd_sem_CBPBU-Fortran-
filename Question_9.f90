		!-------------------------------------------------------------------------!
		!								Question 9								  !
		!-------------------------------------------------------------------------!
		!							Least Square Fitting						  !
		!-------------------------------------------------------------------------!
		!						   decalre all variables						  !
		!-------------------------------------------------------------------------!
		real x,y,sx,sy,sxy,sx2,sx3,sx4,sx2y,sum
		integer	n,i,k,j,m
		real sm1, sm2, h, fxh, fx 
		dimension x(10),y(10),a(100,100),p(100)
		!-------------------------------------------------------------------------!
		n=3
		open(unit=1,file="input.txt")
		read(1,*)m
		do i=1,m
			read(1,*)x(i),y(i)
		enddo
		sx=0
		sy=0
		sxy=0
		sx2=0
		sx3=0
		sx4=0
		sx2y=0
		do i=1,m
			sx=sx+x(i)
			sy=sy+y(i)
			sxy=sxy+x(i)*y(i)
			sx2=sx2+x(i)**2
			sx3=sx3+x(i)**3
			sx4=sx4+x(i)**4
			sx2y=sx2y+((x(i)**2)*y(i))	
		enddo
		open(unit=1,file="output95.txt")
		write(1,*)m,sx,sx2,sy
		write(1,*)sx,sx2,sx3,sxy							
		write(1,*)sx2,sx3,sx4,sx2y
		endfile 1
		rewind 1
		do k=1,n
			read(1,*)(a(k,j),j=1,n+1)
		enddo
		do k=1,n-1
			do i=k+1,n
				do j=k+1,n+1
					a(i,j)=a(i,j)-a(i,k)/(a(k,k)*a(k,j))
				enddo
			enddo
		enddo
		p(n)=a(n,n+1)/a(n,n)
		do k=n-1,1,-1
			sum=0
			do j=k+1,n
				sum=sum+a(k,j)*p(j)
			enddo
			p(k)=1/a(k,k)*(a(k,n+1)-sum)
		enddo
		open(unit=8,file="output93.txt")
		write(8,*)(p(i),i=1,n)
		open(unit=20,file="opt.dat")
		do i=1,m
			sm1=(p(1)+p(2)*x(i)+p(3)*x(i)*x(i))
			sm2=(1.5+0.5*x(i)+5*x(i)*x(i))
			write(20,*) x(i),y(i),sm1, sm2	
			write(20,*) (sm2-sm1)/sm2*100
		enddo
		open(unit=30,file="velo.dat")
		h=0.001
		do i=1,m
			fx=(p(1)+p(2)*x(i)+p(3)*x(i)*x(i))
			fxh=(p(1)+p(2)*(x(i)+h)+p(3)*(x(i)+h)*(x(i)+h))
			sm1=(fxh-fx)/h
			sm2=(0.5+10*x(i))
			write(30,*) x(i),sm1, sm2	
			write(30,*) (sm2-sm1)/sm2*100
		enddo
		stop
		end