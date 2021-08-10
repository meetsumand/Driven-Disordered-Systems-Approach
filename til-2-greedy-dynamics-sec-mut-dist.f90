
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
    integer::str(1:N),i
    double precision::strsize
    strsize=0.d0
    do i=1,N
       strsize=strsize+1.d0*str(i)
    enddo

  end function strsize


    function shdist(str1,str2,N)
    integer,intent(in)::N,str1(1:N),str2(1:N)
    integer::i
    double precision::shdist
    shdist=0.d0
    do i=1,N
       shdist=shdist+1.d0*str2(i)-1.d0*str1(i)
    enddo
 end function shdist

  function hdist(str1,str2,N)
    integer,intent(in)::N,str1(1:N),str2(1:N)
    integer::i
    double precision::hdist
    hdist=0.d0
    do i=1,N
       hdist=hdist+abs(1.d0*str2(i)-1.d0*str1(i))
    enddo
  end function hdist


       
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

  integer,parameter::N=20,NN=20,kexp=1,ns=200,npop=1000,xnm=400
  double precision::reacount(1:xnm,1:2)=0.d0,nmut(1:xnm,1:2)=0.d0,nmut2(1:xnm,1:2)=0.d0,ftns(1:xnm)=0.d0,xpts(1:xnm)=0.d0  
  integer::i,j,ii,m,a,ik,jj,rcntr=0,rea=100000,fg=0,kp=100,ini=0,fnbr,jmin,jmax,ssi,flg=0,comp=0,ini2=0
  integer::stability(0:2**n-1)=0,sstate(1:2**N)=0,sslabel(0:2**N-1)=0,transstate(0:2**N-1,1:2)=0.d0,ili=0
  double precision::xiplus(1:N)=0.d0,ximinus(1:N)=0.d0,xplus(1:2**N)=0.d0,xminus(1:2**N)=0.d0,x0plus=0.d0,x0minus=0.d0
  integer::str(1:N)=0,seed=-9227205,npk=0,nstr(1:N)=0,pkindicator(0:2**n-1)=0,str1(1:N)=0,str2(1:N)=0
  double precision::rx1,mx1,rx2,mx2,xint(1:n*2**(n-1))=0.d0,epsmin=0.d0,xiplusmax=0.d0,xi2,fcurrent=0.d0,fcompare=0.d0
  integer::nbrr=0,ilis=0,start=0,inistep=0,flt=0,st=0,initsize=0
  double precision::upcti=0.d0,upctvar=0.d0,upctser=0.d0
 
  double precision::rsingle(1:N)=1.d0,micw=1.d0,micsingle(1:N)=1.d0,r0(0:2**N-1)=1.d0,r(0:2**N-1)=1.d0,rh=0.d0,y(1:2)=0.d0,dx=0.d0
  double precision::rch(1:NN,1:2**n)=0.d0,ptot(1:NN)=0.d0,rch2(1:NN,1:2**n)=0.d0,peaklist(1:2**N,1:5)=0.d0,x1,x2,mmax=1.d0,upct=0.d0
  double precision::msig(1:xnm)=0.d0,msig2(1:xnm)=0.d0,rsig(1:xnm)=0.d0,sx=0.d0,sx2=0.d0,sreacount(1:xnm)=0.d0,sxpts(1:xnm)=0.d0
  double precision::snmut(1:xnm)=0.d0,snmut2(1:xnm)=0.d0,smcfone(0:n,1:2)=0.d0,smcftot(0:n,1:2)=0.d0
  double precision::smcfone2(0:n,1:2)=0.d0,smcftot2(0:n,1:2)=0.d0
  double precision::kdist(1:xnm,1:2)=0.d0,kdist2(1:xnm,1:2)=0.d0,distmin=0.d0
  integer::stc=0,statelist(0:2**n-1,1:2)=0,stepother=0,stcmax(1:2)=0,smcn(0:n,0:n,1:2)=0

  integer::step=0,ugen(1:xnm,1:2)=5,sig1=0,sig2=0,targ=0,ini0=0,smc=0
  double precision::distav(1:xnm)=0.d0,dist2(1:xnm)=0.d0,dist=0.d0,xl=0.d0,xu=0.d0,mutchange=0.d0,&
mutc(1:xnm,1:2)=0.d0,mutc2(1:xnm,1:2)=0.d0,smcnrea(0:n,1:2)=0,cent=0.d0,dent=0.d0
  double precision::smcdtot(0:n,1:2)=0.d0,smcdtotrea(1:2)=0.d0

  double precision::mic0(0:2**N-1)=1.d0,x=0.d0,rr(0:1)=0.d0,eps=.0000001d0,rnum,pk(1:NN)=0.d0,cdav=0.d0,cs(1:NN)=0.d0,pk2=0.d0
  double precision::pi=acos(-1.d0)


30 format(f48.9,3x,f16.7)
character(len=12)::flnm

open(unit=10,file='N20udmain.d',status='UNKNOWN')
open(unit=20,file='N20udsmdist.d',status='UNKNOWN')
open(unit=30,file='N20nmutetc.d',status='UNKNOWN')

do ik=1,xnm
xpts(ik)=(exp(1.d0))**(2.d0*N*ik/(1.d0*xnm)-3.5d0)
sxpts(ik)=20.d0*(1.d0*ik-0.999d0)/(1.d0*xnm)
enddo


reacount=0.d0



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




!*********************************** Making stable states list

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
targ=2**n-1
		jmax=j
		x0minus=ximinus(j)
		endif
	endif 
   enddo
if (i==0) x0minus=0.d0
if (i==2**n-1) x0plus=x0minus*5.d0

if (x0plus>x0minus) then 
ssi=ssi+1
sstate(ssi)=i
stability(i)=1
sslabel(i)=ssi
xplus(ssi)=mic0(i)*x0plus
xminus(ssi)=mic0(i)*x0minus
else 
stability(i)=0
endif

!**********************************************************


enddo
!*************************************

!print*,"allmut",xminus(sslabel(2**n-1)),xplus(sslabel(2**n-1)),sslabel(2**n-1)
!************************************ dynamics under changing concentration 

x=0.d0
do step=1,2

stc=0
if (step==1) then 
ini=0
x=xminus(sslabel(ini))
targ=2**n-1
!x2=xplus(sslabel(ini))
else
 ini=2**n-1
 targ=0
 x=xplus(sslabel(ini))
 !x2=xminus(sslabel(ini))
endif
stc=stc+1
statelist(stc,step)=ini

st=0
do !************************* start of dynamics


flg=0
ini2=ini
fnbr=ini

inistep=ini
smc=0

if (st.ne.0) then

   do while (flg==0)
  !start=start+1 
  ini=fnbr
!print*,"ini",ini
      do ik=1,N
!print*,ik,fnbr
           nstr=neighbor(bitstr(ini,N),ik,N)
           fcurrent=fitness(r0(fnbr),mic0(fnbr),x)
           fcompare=fitness(r0(dec(nstr,N)),mic0(dec(nstr,N)),x)
	   if (fnbr==inistep) fcurrent=fcurrent*.9d0
	   if (dec(nstr,N)==inistep) fcompare=fcompare*.9d0
		!print*,ik,dec(nstr,N),fitness(r0(dec(nstr,N)),mic0(dec(nstr,N)),x),fnbr,fcurrent
              if (fcompare.ge.fcurrent) then
              fnbr=dec(nstr,N)
              endif
!print*,"fnbr",fnbr
        enddo

   smc=smc+1
 !  print*,ini,smc
   if (fnbr==ini) flg=1

   enddo 



if (fnbr==inistep) then 
print*,"attention",ini2,rcntr,x,step,xminus(sslabel(2**n-1)),xplus(sslabel(2**n-1))
!print*,"xplus",xplus
exit
endif


ini2=fnbr
endif



if (step==1) then 
x2=xplus(sslabel(ini))!+dx
else
x2=xminus(sslabel(ini))
endif

if (st.ne.0) then
!print*,mutchange
stc=stc+1
statelist(stc,step)=ini
endif


if (st.ne.0) then 

initsize=int(strsize(bitstr(inistep,N),N))
smcn(initsize,smc-2,step)=smcn(initsize,smc-2,step)+1
smcnrea(initsize,step)=smcnrea(initsize,step)+1
smcdtot(smc-2,step)=smcdtot(smc-2,step)+1
smcdtotrea(step)=smcdtotrea(step)+1
!mutchange=1.d0*strsize(bitstr(ini,N),N)-1.d0*strsize(bitstr(inistep,N),N)
!mutc(ik,step)=mutc(ik,step)+mutchange
!mutc2(ik,step)=mutc2(ik,step)+1.d0*mutchange**2

smcftot(initsize,step)=smcftot(initsize,step)+1.d0*(1.d0*smc-2.d0)
 smcftot2(initsize,step)=smcftot2(initsize,step)+1.d0*((1.d0*smc-2.d0)**2)

if (smc>2) then 
smcfone(initsize,step)=smcfone(initsize,step)+1.d0
smcfone2(initsize,step)=smcfone2(initsize,step)+1.d0**2
reacount(ik,step)=reacount(ik,step)+1.d0
!print*,"step",ini,step,mutchange,mutc(ik,step)
!reacount(ik,step)=reacount(ik,step)+1.d0
endif

endif

!if (st.ne.0) then
do ik=1,xnm

if (step==1) then 
xu=x2
xl=x
else
xu=x
xl=x2
endif
flt=0
if ((xu.ge.xpts(ik)).and.(xl.le.xpts(ik))) flt=1
if (ini==2**n-1.and.xpts(ik)>xl) flt=1
if (flt==1) then 
nmut(ik,step)=nmut(ik,step)+1.d0*strsize(bitstr(ini,N),N)
!print*,"mid",ini,x,x2,"step",step,1.d0*strsize(bitstr(ini,N),N),"nmut",nmut(ik,step)
nmut2(ik,step)=nmut2(ik,step)+(1.d0*strsize(bitstr(ini,N),N)**2)
Ugen(ik,step)=ini
endif
enddo


 
!endif



!endif

st=1
x=x2
!print*,ini,x,1.d0*strsize(bitstr(ini,N),N)
if (ini==targ) exit
enddo !************************************ end of dynamics

stcmax(step)=stc

enddo !************************************* end of step




do ik=1,xnm

sig1=Ugen(ik,1)
sig2=Ugen(ik,2)
str1=bitstr(sig1,n)
str2=bitstr(sig2,n)

dist=1.d0*(hdist(str1,str2,N))


!print*,dist

distav(ik)=distav(ik)+dist
dist2(ik)=dist2(ik)+dist**2
enddo

do step=1,2
do ik=1,xnm

sig1=Ugen(ik,step)
str1=bitstr(sig1,n)
stepother=3-step

distmin=n*1.d0
do j=1,stcmax(stepother)
str2=bitstr(statelist(j,stepother),N)
dist=1.d0*(hdist(str1,str2,N))
if (dist<distmin) distmin=dist
enddo

kdist(ik,step)=kdist(ik,step)+distmin
kdist2(ik,step)=kdist2(ik,step)+distmin**2
enddo
enddo

enddo  !********************************* end of realizations

nmut=nmut/(1.d0*rea)
nmut2=nmut2/(1.d0*rea)
mutc=mutc/(1.d0*rea)
mutc2=mutc2/(1.d0*rea)
!smcfone2=smcfone2/(1.d0*rea)
!smcftot2=smcftot2/(1.d0*rea)

smcfone=smcfone/(1.d0*rea)
distav=distav/(1.d0*rea)
dist2=dist2/(1.d0*rea)
dist2=sqrt(dist2-(1.d0*distav)**2)/sqrt((1.d0*rea))
kdist=kdist/(1.d0*rea)
kdist2=kdist2/(1.d0*rea)

do step=1,2
do ik=0,n
 smcftot(ik,step)=smcftot(ik,step)/(1.d0*smcnrea(ik,step))
 smcftot2(ik,step)=smcftot2(ik,step)/(1.d0*smcnrea(ik,step))
 enddo
enddo

do step=1,2
do ik=1,xnm
!if (reacount(ik,step).ne.(0.d0)) smcftot(ik,step)=smcftot(ik,step)/(1.d0*reacount(ik,step))
mutc(ik,step)=mutc(ik,step)/(1.d0*reacount(ik,step))
kdist2(ik,step)=sqrt(kdist2(ik,step)-kdist(ik,step)**2)
nmut2(ik,step)=sqrt(nmut2(ik,step)-(nmut(ik,step))**2)
enddo

do ik=0,n
smcfone2(ik,step)=smcfone2(ik,step)/(1.d0*rea)
smcfone2(ik,step)=sqrt(smcfone2(ik,step)-smcfone(ik,step)**2)
!if (reacount(ik,step).ne.(0.d0)) smcftot2(ik,step)=smcftot2(ik,step)/(1.d0*reacount(ik,step))
smcftot2(ik,step)=sqrt(smcftot2(ik,step)-smcftot(ik,step)**2)

!mutc2(ik,step)=sqrt(mutc2(ik,step)-mutc(ik,step)**2)
enddo

enddo

smcfone2=smcfone2/sqrt(1.d0*rea)
smcftot2=smcftot2/sqrt(1.d0*rea)
nmut2=nmut2/sqrt(1.d0*rea)
kdist2=kdist2/sqrt(1.d0*rea)

!do ik=1,xnm 
!print*,xpts(ik),distav(ik),dist2(ik),reacount(ik)!nmut(ik),nmut2(ik),msig(ik),msig2(ik),rsig(ik),sxpts(ik),snmut(ik),snmut2(ik)
!enddo

do ik=0,n
write(10,*) ik,smcfone(ik,1),smcfone2(ik,1),smcfone(ik,2),smcfone2(ik,2),smcftot(ik,1),smcftot2(ik,1)&
,smcftot(ik,2),smcftot2(ik,2)
enddo

do ik=1,xnm
write(30,*) xpts(ik),nmut(ik,1), nmut2(ik,1), nmut(ik,2), nmut2(ik,2),kdist(ik,1),kdist2(ik,1),kdist(ik,2), kdist2(ik,2)
enddo

!do i=0,n 

!write(20,*) " "
do j=0,n
 ! cent=smcn(i,j,1)/(1.d0*smcnrea(i,1))
 !dent=smcn(i,j,2)/(1.d0*smcnrea(i,2))
!if (i==0) dent=0.d0
! if (i==n) cent=0.d0
!write(20,*) i,j,cent,dent
write(20,*) j,(1.d0*smcdtot(j,1))/(1.d0*smcdtotrea(1)),(1.d0*smcdtot(j,2))/(1.d0*smcdtotrea(2))
enddo

!enddo

do i=0,n 

write(20,*) " "
do j=0,n
  cent=smcn(i,j,1)/(1.d0*smcnrea(i,1))
 dent=smcn(i,j,2)/(1.d0*smcnrea(i,2))
if (i==0) dent=0.d0
 if (i==n) cent=0.d0
write(20,*) i,j,cent,dent
!write(20,*) j,(1.d0*smcdtot(j,1))/(1.d0*smcdtotrea(1)),(1.d0*smcdtot(j,2))/(1.d0*smcdtotrea(2))
enddo

enddo

end program main




  
  

  
