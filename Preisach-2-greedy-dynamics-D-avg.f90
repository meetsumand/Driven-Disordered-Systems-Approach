
module strings
implicit none
contains
  

function w(x)
double precision,intent(in)::x
double precision::w

w=1.d0/(1.d0+x**2)

end function w

function findxiplus(r,m)
double precision,intent(in)::r,m
double precision::findxiplus
findxiplus=m*sqrt((1.d0-r)/(r*(m**2)-1.d0))
end function findxiplus

 
function fitness(r,m,x)
double precision, intent(in)::r,m,x
double precision::fitness

fitness=r*w(x/m)

end function fitness


function intxy(r1,m1,r2,m2)
double precision,intent(in)::r1,m1,r2,m2
double precision::intxy
intxy=sqrt((r1-r2)/(r2/m1**2-r1/m2**2))
end function intxy


function bitstr(i,N)
  integer,intent(in)::N,i
  integer::k,j,bitstr(1:N)
    
  do j=1,N
     k=j
     bitstr(j)=(mod(i,2**(k))-mod(i,2**(k-1)))
     bitstr(j)=bitstr(j)/(2**(k-1))
  enddo


end function bitstr



function insertbit(str,i,ii,N)
  integer,intent(in)::N,i,ii
  integer::str(1:N),insertbit(1:N),k,j

  insertbit=str
  
  do j=2,i
     k=str(j)
     insertbit(j-1)=k
  enddo

  insertbit(i)=ii

end function insertbit


function dec(str,N)
  integer,intent(in)::N,str(1:N)
  integer::i,dec
!common N

  dec=0
  do i=1,n
     dec=dec+str(i)*(2**(i-1))
  enddo
  
end function dec

  function dec2(str,N)
  integer,intent(in)::N,str(1:N)
  integer::i,dec2
!common N

  dec2=0
  do i=1,n
     dec2=dec2+str(n-i+1)*(2**(i-1))
  enddo


end function dec2



function neighbor(str,i,N)
    integer,intent(in)::N,i,str(1:N)
    integer::neighbor(1:N)

    neighbor=str
    neighbor(i)=abs(1-str(i))


  end function neighbor

function bitflip(str,i,N)
    integer,intent(in)::N,i,str(1:N)
    integer::bitflip(1:N)

    bitflip=str
    bitflip(i)=abs(1-str(i))

  end function bitflip


    function strsize(str,N)
    integer,intent(in)::N
    integer::str(1:N),i,strsize
strsize=0
    do i=1,N
       strsize=strsize+1*str(i)
    enddo

  end function strsize

  
    
       
FUNCTION ran2(idum)
!    gerador de números aleatórios baseado na combinacao
!    de dois geradores lineares congruenciais tendo um periodo
!    maior do que 2x10^18. A saida e' "baralhada" 
!     -> valor inicial de idum=-1  (Numerical recipes 2a. edicao)
IMPLICIT NONE
INTEGER, PARAMETER :: im1=2147483563,im2=2147483399,imm1=im1-1,&
     ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,&
     ir2=3791,ntab=32,ndiv=1+imm1/ntab
DOUBLE PRECISION , PARAMETER ::   am=1.d0/im1,eps=1.d-14,rnmx=1.d0-eps
DOUBLE PRECISION :: ran2
INTEGER, DIMENSION(ntab) :: iv
INTEGER :: idum,idum2,j,k,iy

save iv,iy,idum2
data idum2/123456789/,iv /ntab*0/,iy /0/
      
if(idum.le.0) then
  idum=max(-idum,1)
  idum2=idum
  do j=ntab+8,1,-1
    k=idum/iq1
    idum=ia1*(idum-k*iq1)-ir1*k
    if(idum.lt.0) idum=idum+im1
    if(j.le.ntab) iv(j)=idum
   end do
   iy=iv(1)
endif

k=idum/iq1
idum=ia1*(idum-k*iq1)-ir1*k
if(idum.lt.0) idum=idum+im1
k=idum2/iq2
idum2=ia2*(idum2-k*iq2)-ir2*k

if(idum2.lt.0) idum2=idum2+im2

j=1+iy/ndiv
iy=iv(j)-idum2
iv(j)=idum
if (iy.lt.1)iy=iy+imm1
ran2=min(am*iy,rnmx)

END FUNCTION ran2


end module










!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM











program main
  
  use strings
  
  implicit none

  integer,parameter::N=20,NN=20,kexp=1,ns=200,npop=1000,xnm=1600
  double precision::reacount(1:xnm)=0.d0,nmut(1:xnm)=0.d0,ftns(1:xnm)=0.d0,xpts(1:xnm)=0.d0  
  integer::i,j,ii,m,a,ik,jj,rcntr=0,rea=10000,fg=0,kp=100,ini=0,fnbr,jmin,jmax,ssi,flg=0,comp=0,ini2=0
  integer::stability(0:2**n-1)=0,sstate(1:2**N)=0,sslabel(0:2**N-1)=0,ili=0
  double precision::xiplus(1:N)=0.d0,ximinus(1:N)=0.d0,xplus(1:2**N)=0.d0,xminus(1:2**N)=0.d0,x0plus=0.d0,x0minus=0.d0
  integer::str(1:N)=0,seed=-9971046,npk=0,nstr(1:N)=0,pkindicator(0:2**n-1)=0,str1(1:N)=0,str2(1:N)=0
 double precision::rx1,mx1,rx2,mx2,xint(1:n*2**(n-1))=0.d0,epsmin=0.d0,xiplusmax=0.d0,xi2,fcurrent=0.d0,fcompare=0.d0
integer::nbrr=0,ilis=0,start=0,inistep=0,flt=0,flpdni=0,lsi(1:2**n,1:2)=0
double precision::upcti=0.d0,upctvar=0.d0,upctser=0.d0
 
  double precision::rsingle(1:N)=1.d0,micw=1.d0,micsingle(1:N)=1.d0,r0(0:2**N-1)=1.d0,r(0:2**N-1)=1.d0,rh=0.d0,y(1:2)=0.d0,dx=0.d0
 double precision::rch(1:NN,1:2**n)=0.d0,ptot(1:NN)=0.d0,rch2(1:NN,1:2**n)=0.d0,peaklist(1:2**N,1:5)=0.d0,x1,x2,mmax=1.d0,upct=0.d0

  double precision::mic0(0:2**N-1)=1.d0,x=0.d0,rr(0:1)=0.d0,eps=.0000001d0,rnum,pk(1:NN)=0.d0,cdav=0.d0,cs(1:NN)=0.d0,pk2=0.d0
double precision::pi=acos(-1.d0)


30 format(f48.9,3x,f16.7)
character(len=12)::flnm

do ik=1,xnm
xpts(ik)=1.d0*N*0.001d0*(1.d0*ik-1.d0)
enddo

upct=0.d0
reacount=0.d0
upctvar=0.d0


do rcntr=1,rea,1
!print*,"rcntr",rcntr
!************ Assigning values to r_i 
    do i=1,N
     x1=ran2(seed)
     x2=ran2(seed)
     rsingle(i)=exp(-abs(sqrt(-2.d0*log(x1))*sin(2.d0*pi*x2)))
  end do
!************



!*********** Assigning values to m_i
mmax=1.d0
    do i=1,N
flg=0
      do while(flg==0)
     x1=ran2(seed)
     x2=ran2(seed)
     micsingle(i)=exp(1.d0*abs(sqrt(-2.d0*log(x1))*cos(2.d0*pi*x2)))
     !micsingle(i)=micsingle(i)*1.d0/sqrt(rsingle(i))
     if (micsingle(i)>1.d0/sqrt(rsingle(i))) flg=1
     enddo
enddo
!***********





!************ Fining xiplus and ximinus for mutations
do i=1,n
xiplus(i)=findxiplus(rsingle(i),micsingle(i))
ximinus(i)=xiplus(i)/micsingle(i)
!print*,"xi p m",xiplus(i),ximinus(i)
enddo
!*****************************************

!**************** rearranging indices to get x_1<x_2<...<x_n
do i=1,n
do j=1,n-1
if (xiplus(j)>xiplus(j+1)) then 
xi2=xiplus(j)
xiplus(j)=xiplus(j+1)
xiplus(j+1)=xi2

xi2=ximinus(j)
ximinus(j)=ximinus(j+1)
ximinus(j+1)=xi2


xi2=rsingle(j)
rsingle(j)=rsingle(j+1)
rsingle(j+1)=xi2

xi2=micsingle(j)
micsingle(j)=micsingle(j+1)
micsingle(j+1)=xi2

endif
enddo
enddo


xiplusmax=xiplus(n)


!**************** Assigning r and m to genotypes
r0=1.d0
mic0=1.d0
r0(0)=1.d0
mic0(0)=1.d0

  do i=1,2**n
  
     str=bitstr(i,N)
     do m=1,n
       if (str(m)==1)  r0(i)=r0(i)*rsingle(m)
     enddo


     do m=1,n
         if (str(m)==1) mic0(i)=mic0(i)*micsingle(m)
     enddo


  enddo
!***************** end

!************************** calculating minm intersection dist epsmin

!ili=0
!do i=0,2**n-1

!do j=1,n
!nstr=neighbor(bitstr(i,N),j,N)
!nbrr=dec(nstr,N)
!if (strsize(nstr,N)>strsize(bitstr(i,N),N)) then
!rx1=r0(i)
!mx1=mic0(i)

!rx2=r0(nbrr)
!mx2=mic0(nbrr)
!ili=ili+1
!xint(ili)=intxy(rx1,mx1,rx2,mx2)
!endif
!enddo

!enddo

!epsmin=1.d0
!do i=1,n*(2**(n-1))-1
!do j=i+1,n*(2**(n-1))
!if (abs(xint(i)-xint(j))<epsmin) epsmin=abs(xint(i)-xint(j))
!enddo
!enddo


!*********************************** Making stable states list

!epsmin=0.00001d0

stability=0
ssi=0
do i=0,2**n-1
  str=bitstr(i,n)
 x0minus=0.d0
 x0plus=xiplusmax+1.d0
!********** Finding minimun of x_i^+ and mximum of x_i^-
   do j=1,N
	if (str(j)==0) then 
!print*,"xiplus",xiplus(j)
         if (xiplus(j)<x0plus) then 
		jmin=j
		x0plus=xiplus(j)
!if (i==0) print*,"wt",xiplus,j
	 endif 
! print*,"x0plus",x0plus
	else
    	    if (ximinus(j)>x0minus) then 
		jmax=j
		x0minus=ximinus(j)
		endif
	endif 
   enddo
if (i==0) x0minus=0.d0
if (i==2**N-1) x0plus=(x0minus)*2.d0

if (x0plus>x0minus) then 
ssi=ssi+1
lsi(ssi,1)=jmin
lsi(ssi,2)=jmax
sstate(ssi)=i
stability(i)=1
sslabel(i)=ssi
xplus(ssi)=x0plus
xminus(ssi)=x0minus
else 
stability(i)=0
endif

!**********************************************************


enddo
!*************************************


!************************************ dynamics under changing concentration 

ini=2**n-1
x=0.000d0
upcti=0.d0
!.001227d0
!print*,"epsmin",epsmin

x=xplus(sslabel(ini))
x2=xminus(sslabel(ini))

!print*,"vals",x,x2,xpts(xnm)
do ik=1,xnm
flt=0
if ((x2.le.xpts(ik)).and.(x.ge.xpts(ik))) flt=1
if (ini==2**n-1.and.xpts(ik)>x) flt=1
if (flt==1) then 
nmut(ik)=nmut(ik)+strsize(bitstr(ini,N),N)
ftns(ik)=ftns(ik)+fitness(r0(ini),mic0(ini),xpts(ik))
reacount(ik)=reacount(ik)+1.d0
!print*,"ik x xl xu mut",ik,xpts(ik),x,x2,ini
endif
enddo 


do while(ini.ne.0)!************************* start of dynamics

flg=0
!dx=0.000005d0*(1.d0+xplus(sslabel(ini)))


x=xminus(sslabel(ini))!+dx
flpdni=lsi(sslabel(ini),2)
ini=dec(bitflip(bitstr(ini,N),flpdni,N),N)
x2=xminus(sslabel(ini))
!print*,ini,x
do ik=1,xnm
flt=0
if ((x2.le.xpts(ik)).and.(x.ge.xpts(ik))) flt=1
if (ini==2**n-1.and.xpts(ik)>x) flt=1
if (flt==1) then 
nmut(ik)=nmut(ik)+strsize(bitstr(ini,N),N)
ftns(ik)=ftns(ik)+fitness(r0(ini),mic0(ini),xpts(ik))
reacount(ik)=reacount(ik)+1.d0
!print*,"ik x xl xu mut",ik,xpts(ik),x,x2,ini
!print*,"if",nmut(ik),strsize(bitstr(ini,N),N)
endif
enddo 



enddo !************************************ end of dynamics


upctvar=upctvar+(upcti)**2
upct=upct+upcti
enddo
upct=upct/(1.d0*rea)
upctvar=sqrt(upctvar/(1.d0*rea)-upct**2)
upctser=upctvar/(sqrt(1.d0*rea))
do ik=1,xnm 
print*,xpts(ik),nmut(ik)/reacount(ik),ftns(ik)/reacount(ik)
!print*,xpts(ik),ftns(ik),nmut(ik)
enddo

!print*,rcntr
!print*,"no of trans",N,upct,upctser,"stdev",upctvar,"rea",rea

end program main




  
  

  
