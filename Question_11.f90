!      Question 11
  real mass, vel, bmag
  real qch, rl, nt
  real nint, pi, w
  real t_per, tot_time
  real t0, n_t, h, np

  
  real x(110),y(110)
  real t(110)
  real vx(110),vy(110)
  real mx(110),my(110)
  real tvx(110),tvy(110)
  real tmx(110),tmy(110)
  real errvx(110),errvy(110)
  real errmx(110),errmy(110)

  mass=1.67e-27      ! mass of the particle
  vel=1.5e6        !velocity
  bmag=5000.        ! magnetic feild
  qch=1.6e-19        ! ion charge
  rl=mass*vel/bmag/qch  !radius
  nt=2.          ! number of time period
  nint=50.        ! number of time interval in 1 tp
  pi=4.*atan(1.0)      ! pi value
  w=(vel/rl)        !omega frequency
  t_per=(2.*pi)/w      ! time period
  tot_time=nt*t_per    ! total time
  n_t=nt*nint        ! total no. of interval
  h=tot_time/n_t      ! time step
  np=4.          ! number of points used to obtain derivative
  t0=0          ! innitial time

  call calculate_input(t0,rl,w,n_t,h,x,y)
  call save_data(n_t,x,y)
  call calculate_all(n_t,rl,w,np,h,mass,t0,t,vx,vy,mx,my)
  call calculate_theory(n_t,rl,w,mass,t,tvx,tvy,tmx,tmy)
  call calculate_error(n_t,vx,tvx,errvx)
  call calculate_error(n_t,vy,tvy,errvy)
  call calculate_error(n_t,mx,tmx,errmx)
  call calculate_error(n_t,my,tmy,errmy)

  call write_output(1,n_t,t,x,y,vx,vy,tvx,tvy,errvx,errvy)
  call write_output(2,n_t,t,x,y,mx,my,tmx,tmy,errmx,errmy)

  stop
  end



  subroutine write_output(m,n,t,x,y,vx,vy,tvx,tvy,errvx,errvy)

  integer m
  real n
  real x(110),y(110)
  real t(110)
  real vx(110),vy(110)  
  real tvx(110),tvy(110)
  real errvx(110),errvy(110)

   if(m.eq.1) then
   open(unit=1,file="vx-vy.dat")
   open(unit=2,file="v-err.dat")
   else
   open(unit=1,file="mx-my.dat")
   open(unit=2,file="m-err.dat")
   endif

   do i=1,int(n)
   write(1,11)t(i),x(i),y(i),vx(i),vy(i),tvx(i),tvy(i)
   write(2,11)t(i),vx(i),vy(i),tvx(i),tvy(i),errvx(i),errvy(i)

   enddo
   close(1)
   close(2)
11  format(19(e14.7,3x))
  return
  end

  subroutine calculate_error(n,x,tx,ex)
  real n,x(110),tx(110),ex(110) 

  do i=2, int(n)
  ex(i)=(tx(i)-x(i))/tx(i)*100
  enddo
  return
  end

  subroutine calculate_theory(n,r,w,mass,t,tvx,tvy,tmx,tmy)
  real n,r,w,mass
  real t(110)
  real tvx(110),tvy(110)
  real tmx(110),tmy(110)

  do i=1,int(n)
  tvx(i)=-r*w*cos(w*t(i))
  tvy(i)=-r*w*sin(w*t(i))
  tmx(i)=mass*tvx(i)
  tmy(i)=mass*tvy(i)
  enddo
  return
  end







  subroutine calculate_all(n_t,rl,w,np,h,mass,t0,t,vx,vy,mx,my)
  real t0,rl,w,np
  real hp,h,mass,n_t
  real t(110)
  real vx(110),vy(110)
  real mx(110),my(110)
 
  do i=1,int(n_t)
  t(i)=t0+(i-1)*h
  hp=h
  call calculate_diff(rl,w,t(i),np,hp,vx(i),vy(i))
  mx(i)=mass*vx(i)
  my(i)=mass*vy(i)
  enddo
  return
  end



  subroutine calculate_diff(r,w,t0,n,h,vx,vy)
  real r,w,t0,n,h,vx,vy
  real x(4),y(4)
  real del_x(3), del_y(3)


  call calculate_input(t0,r,w,n,h,x,y)
  call New_table(n,x,del_x)
  call New_table(n,y,del_y)
  vx=0
  vy=0
  do i=1,3
  vx=vx+((-1)*(i-1))*del_x(i)/float(i)
  vy=vy+((-1)*(i-1))*del_y(i)/float(i)
  enddo
  vx=vx/h
  vy=vy/h
  return
  end

  subroutine New_table(n,x,dx)
  real n,x(4),dx(3),dy(3,3)

  do i=1,int(n)-1
  do j=1,int(n)-i
  if(i.eq.1) dy(i,j)=x(j+1)-x(j)
  if(i.ne.1) dy(i,j)=dy(i-1,j+1)-dy(i-1,j)

  enddo
  enddo

  do i=1,3
   dx(i)=dy(i,1)
  enddo

  return
  end


  subroutine save_data(n,x,y)
  real n,x(110),y(110)
  open(unit=1,file="known.dat")
  do i=1,int(n)
  write(1,*) x(i),y(i)
  enddo
  close(1)
  return
  end



  subroutine calculate_input(t,r,w,n,h,x,y)
  real r,w,n,h,x(110),y(110),t

  do i=1,int(n)
  x(i)=0.014-r*sin(w*(t+(i-1)*h))
  y(i)=0.025+r*cos(w*(t+(i-1)*h))
  enddo
  return
  end