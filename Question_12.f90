!question 12
			  real x,y,dy,sum,num,err,o,val,b,m,sx,sy,sxy,sx2,A0,A1,pre,tim
				dimension x(10),p(10)
				dimension y(10)
				dimension dy(10,10)
				integer n,i,j
				 write(*,*)"Enter the temperature value"
			 read(*,*) o
			  write(*,*)"Enter the pressure value"
			 read(*,*) m

				 open(unit=1,file="input.txt")
			 read(1,*)n
				 do i=1,n
				read(1,*) x(i),y(i)
				enddo
			endfile 1
!Temperature Calculator Calculation

!Newton divide difference formula

            do i=1,n-1
            dy(1,i)=(y(i+1)-y(i))/(x(i+1)-x(i))
            enddo

            do j=2,n-1
             do i=1,n-j
            dy(j,i)=(dy(j-1,i+1)-dy(j-1,i))/(x(i+j)-x(i))
            enddo
            enddo
      sum=y(1)
      do i=1,n-1
      num=1.
      do j=1,i
                num=num*(o-x(j))
               enddo
                sum=sum+num*dy(i,1)
        enddo
!Lagrange method
         val=0.
            do i=1,n
                b=1.
                do j=1,n
                    if((i).ne.(j))then
            b=b*((o-x(j))/(x(i)-x(j)))
            endif
            enddo
           val=val+(b*y(i))
            enddo
              open(9,file="opt.txt")
        err=(sum-val)/sum*100
          

         write(9,*)"the value is",sum,val
         write(9,*) err

        
!pressure Calculation
           rewind 1
                    read(1,*)n
                       do i=1,n
                      read(1,*) y(i),x(i)
                       enddo
             endfile 1
 !Newton divide difference formula

          do i=1,n-1
            dy(1,i)=(y(i+1)-y(i))/(x(i+1)-x(i))
            enddo

            do j=2,n-1
             do i=1,n-j
            dy(j,i)=(dy(j-1,i+1)-dy(j-1,i))/(x(i+j)-x(i))
            enddo
            enddo
      sum=y(1)
      do i=1,n-1
      num=1.
      do j=1,i
                num=num*(m-x(j))
               enddo
                sum=sum+num*dy(i,1)
        enddo


        !Lagrange method
         val=0.
            do i=1,n
                b=1.
                do j=1,n
                    if((i).ne.(j))then
            b=b*((m-x(j))/(x(i)-x(j)))
            endif
            enddo
           val=val+(b*y(i))
            enddo
              open(9,file="opt.txt")
        err=(sum-val)/sum*100
          

         write(9,*)"the value is",sum,val
         write(9,*) err


    !Fitting of Pressure Vs Temperature curve
              rewind 1
                    read(1,*)n
                       do i=1,n
                      read(1,*) x(i),y(i)
                       enddo
             endfile 1
             sx=0
             sy=0
             sxy=0
             sx2=0
             do i=1,n
             sx=sx+x(i)
             sy=sy+y(i)
             sxy=sxy+x(i)*y(i)
             sx2=sx2+x(i)**2
             enddo
             A1=(sy*sx-n*sxy)/((sx)**2-n*sx2)
       
                       A0=(sxy*sx-sy*sx2)/(sx**2-n*sx2)         
              write(9,*) "the fitted equation is y=",A1,"*x+",A0
            open(2,file="fit.dat")
            do i=1,n
            p(i)=A1*x(i)+A0
            write(2,*) x(i),y(i),p(i)

              enddo
              endfile 2
              pre=A1*(52.5)+A0
             tim=(7.8-A0)/A1

             write(9,*)"The pressure from fitted curve",pre
              write(9,*)"The temperature from fitted curve",tim
      
             stop
             end