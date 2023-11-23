subroutine grid
use params
use mini
implicit none

double precision :: d_eps, put_eps0(neps), dummy
integer :: iter_eps
integer, parameter :: itmax_eps=5000
double precision, parameter :: tol_eps=1e-9

integer :: t

double precision, dimension(nz,nz) :: pi_z

double precision, dimension(neta,neta) :: pi_eta
double precision, dimension(neta) :: grid_eta

double precision, dimension(niota,niota) :: pi_iota
double precision, dimension(niota) :: grid_iota

integer :: id
integer :: jeta, idp

integer :: find

    do t=2,t1

        do ia=1,intbirth
        grid_a(t,ia)=amin+(zero-amin)*(ia-one)/(intbirth-one)
        end do
            
        do ia=intbirth,na
        grid_a(t,ia)=zero+(amax-zero)*(real(ia-intbirth)/real(na-intbirth))**1.5
        end do

    grid_a(t,intbirth)=zero
        
    end do
    
grid_a(1,:)=grid_a(2,:)

    do ik=1,nk
    grid_k(ik)=kmin+(kmax-kmin)*(((ik-one)/(nk-one))**1.5)
    end do

grid_k(1)=kmin
grid_k(nk)=kmax
    
find=0

    do ik=1,nk

        if ( grid_k(ik) .ge. k_e .and. find .eq. 0 ) then
        find=1
        grid_k(ik)=k_e
        intbirth_k=ik
        end if
        
    end do

    if ( rank .eq. 0 ) then
    print*,'intbirth_k',intbirth_k,grid_k(intbirth_k)
    end if

    if ( nz .eq. 1 ) then
    
    grid_z=1
    pi_z=1
    
    else
    
        call tauchen (nz,zero,zero,var_z,grid_z,pi_z)
    grid_z=exp(grid_z)

    end if

    if ( niota .eq. 1 ) then
    
    grid_iota=1
    pi_iota=1
    
    else
    
!         call tauchen (niota,rho_iota,mu_iota,var_iota,grid_iota,pi_iota)
!     grid_iota=exp(grid_iota)

    !covid shock: firm has 0 or 1 productivity
    grid_iota(1)=1
    grid_iota(2)=0
!     pi_iota(:,1)=1-shock         !fraction of essential firms
!     pi_iota(:,2)=shock         !fraction of non-essential firms

    pi_iota(1,1)=1
    pi_iota(1,2)=0
    
    pi_iota(2,1)=1-rho_shock         !fraction of non-essential firms
    pi_iota(2,2)=rho_shock         !fraction of non-essential firms
    
        if ( niota .ne. 2 ) then
        print*,'error in grid.f90'
            call exit
        end if
        
    end if
    
        call tauchen (neta,rho_eps,mu_eps,var_eps_t(1),grid_eta,pi_eta)
    grid_eta=exp(grid_eta)
       
put_eps0=zero

    do iiota=1,niota
    put_eps0((iiota-1)*neta*nz+1)=pi_iota(1,iiota)/sum(pi_iota(1,:))
    end do

! pi_iota=zero
! 
!     do iiota=1,niota
!     pi_iota(iiota,iiota)=one
!     end do   
   
id=0
    
    do iiota=1,niota
    do iz=1,nz
    do ieta=1,neta
    
    id=id+1
    grid_eps(id)=grid_eta(ieta)*grid_iota(iiota)
    id_z(id)=iz
    id_iota(id)=iiota
    
    idp=0

        do jiota=1,niota    
        do jz=1,nz    
        do jeta=1,neta
        idp=idp+1
        pi_eps(id,idp)=pi_eta(ieta,jeta)*pi_z(iz,jz)*pi_iota(iiota,jiota)
        end do
        end do
        end do
        
    end do
    end do    
    end do

d_eps=1
iter_eps=0
! put_eps0=one/neps

    do while ( iter_eps .le. itmax_eps .and. d_eps .ge. tol_eps )

    put_eps=matmul(put_eps0,pi_eps)
    d_eps=maxval(abs(put_eps0-put_eps))
    put_eps0=put_eps

    iter_eps=iter_eps+1

    end do

    if ( rank .eq. 0 .and. iter_eps .ge. itmax_eps ) then
    print*,'put_eps did not converge ',d_eps
        call exit
    end if

grid_eps_original=grid_eps    
grid_eps=grid_eps/(sum(grid_eps*put_eps))   

dummy=zero

    do ieps=1,neps
    iz=id_z(ieps)
    dummy=dummy+put_eps(ieps)*grid_z(iz)
    end do
        
grid_z=grid_z/dummy

    do t=1,t1
    grid_eps_t(t,:)=grid_eps
    end do

    if ( niota .eq. 2 ) then
    grid_eps_t(1,nz*neta+1:neps)=grid_eps_t(1,1:nz*neta)
    grid_eps_t(t1,nz*neta+1:neps)=grid_eps_t(t1,1:nz*neta)
    end if
    
!covid shock for firms    
!     do t=2,3        !number of no operation periods
!     do t=1,t1        !number of no operation periods
!     do ieps=1,neps
!     
!     iiota=id_iota(ieps)
!     
!         if ( iiota .eq. 2 ) then
!         grid_eps_t(t,ieps)=zero
!         end if
!         
!     end do
!     end do

    do t=1,t1
    put_eps_t(t,:)=put_eps
    end do

    do t=2,min(t1-1,30)
    
        if ( niota .ne. 2 ) then
        
        put_eps_t(t,:)=put_eps
        
        else if ( niota .eq. 2 ) then
        
        ieps=0
            
            do iz=1,nz
            do ieta=1,neta
            
            ieps=ieps+1
            
            put_eps_t(t,neta*nz+ieps)=put_eps(ieps)*shock*(rho_shock**(t-2.))
            put_eps_t(t,ieps)=put_eps(ieps)-put_eps_t(t,neta*nz+ieps)
            
            end do
            end do
            
        end if
        
    end do
    
    if ( rank .eq. 0 ) then
    
        open (unit=1,file='output_grid_eps.txt')
        write(1,*)'ieps grid_eps grid_z put_eps grid_eps_original cut_off_data'
            do ieps=1,neps
            iz=id_z(ieps)
            write(1,'(300F12.7)')real(ieps),grid_eps_t(1,ieps),grid_z(iz),put_eps(ieps),pi_eps(ieps,:),grid_eps_original(ieps),cut_off_data
            end do    
        close (1)
        
        open (unit=1,file='output_grid_k.txt')
        write(1,*)'ik k'
            do ik=1,nk
            write(1,'(300F12.7)')real(ik),grid_k(ik)
            end do    
        close (1)

        open (unit=1,file='output_grid_a.txt')
        write(1,*)'ia a_1 a_t1'
            do ia=1,na
            write(1,'(300F12.7)')real(ia),grid_a(1,ia),grid_a(t1,ia)
            end do    
        close (1)
        
        open (unit=1,file='output_entry_distribution.txt')
            do ieps=1,neps
            write(1,'(300F12.7)')grid_eps(ieps),put_eps_t(:,ieps)
            end do
        close(1)
        
    end if
    
end subroutine grid
