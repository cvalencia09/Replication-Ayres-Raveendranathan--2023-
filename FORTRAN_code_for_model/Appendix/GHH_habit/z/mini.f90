module mini
use params

contains

function func_adj (kp,k,d,time)
implicit none
double precision :: func_adj
double precision, intent(in) :: kp, k
integer(kind=1), intent(in) :: d
integer, intent(in) :: time
double precision :: x_i

x_i=kp-(d*(1-delta_k)+(1-d)*chi)*k
    
    if ( cap_adj .eq. 0 ) then  !no adjustment cost
    
    func_adj=zero+lambda*(x_i/k-delta_k)**2.

    else if ( cap_adj .eq. 1 ) then !costly reversibility of investment
    
        if ( x_i .ge. zero ) then
        func_adj=zero+lambda*(x_i/k-delta_k)**2.
        else if ( x_i .lt. zero ) then
        func_adj=-gamma*x_i+lambda*(x_i/k-delta_k)**2.
        end if
        
    else if ( cap_adj .eq. 2 ) then !non-convex

        if ( x_i .ge. gamma1*k .and. x_i .le. gamma2*k ) then
        func_adj=zero+lambda*(x_i/k-delta_k)**2.        
        else 
        func_adj=w_t(time)*f_k+lambda*(x_i/k-delta_k)**2.
        end if
        
    else
    
    print*,'error in func_adj, mini.f90'
        call exit
    
    end if

end function func_adj
                
subroutine vfn_i
implicit none
include 'mpif.h'

integer(kind=1), allocatable, dimension(:) :: g_i_o_mpi, g_binding_collateral_mpi
double precision, allocatable, dimension(:) :: g_i_a_mpi
integer, allocatable, dimension(:) :: g_i_a_integer_mpi
integer(kind=2), allocatable, dimension(:) :: g_i_k_mpi
double precision, allocatable, dimension(:) :: v_i_mpi, e_v_i_mpi

integer(kind=1), allocatable, dimension(:) :: g_i_o_recmpi, g_binding_collateral_recmpi
double precision, allocatable, dimension(:) :: g_i_a_recmpi
integer, allocatable, dimension(:) :: g_i_a_integer_recmpi
integer(kind=2), allocatable, dimension(:) :: g_i_k_recmpi
double precision, allocatable, dimension(:) :: v_i_recmpi, e_v_i_recmpi

double precision :: vtemp, ev, kp_un(neps)
integer :: find

integer :: iter_v
double precision :: d_v, dummy
double precision :: d_v_mpi, d_v_recmpi(nproc)

integer :: iter_re
double precision, dimension(neps,nk,na) :: y2
integer :: ja_temp

    call mpi_barrier(mpi_comm_world,ierr)

    allocate (g_i_o_recmpi(np*nproc),g_i_k_recmpi(np*nproc),g_i_a_recmpi(np*nproc),v_i_recmpi(np*nproc),e_v_i_recmpi(np*nproc),&
    g_binding_collateral_recmpi(np*nproc),g_i_a_integer_recmpi(np*nproc))

!initial guess for value function

    if ( read_v_i .eq. 0 ) then

    read_v_i=0
    
        allocate (v_i_mpi(np))
    
        do cnt=itop,iend
    
        ieps=zmap(cnt,1)
        ik=zmap(cnt,2)
        ia=zmap(cnt,3)
        iz=id_z(ieps)
        
        eps=grid_eps(ieps)
        k=grid_k(ik)
        a=grid_a(date,ia)
        kp=kmin
    
        l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
        
        !operate
        d_o=1
        profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))+&
        (1.+tau_x)*(1-delta_k)*k+a-func_adj(kp,k,d_o,date)-f_o_good
            
        !exit
        d_o=0

            if ( default_economy .eq. 0 ) then            
            vtemp=(1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date)                    
            else if ( default_economy .eq. 1 ) then
            vtemp=max((1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date),min(zero,(1.+tau_x)*k*chi-func_adj(kp,k,d_o,date)))                                    
            end if

            if ( profit .gt. zero ) then
            v_i_mpi(cnt-itop+1)=max(profit,vtemp)
            else if ( profit .le. zero ) then
            v_i_mpi(cnt-itop+1)=vtemp
            end if
                            
        end do

        call mpi_allgather (v_i_mpi,np,mpi_double_precision,v_i_recmpi,np,mpi_double_precision,mpi_comm_world,ierr)
    
        deallocate (v_i_mpi)
        
    g_i_o=1

    else if ( read_v_i .eq. 1 ) then
    
        do cnt=1,ni        
        ieps=zmap(cnt,1)
        ik=zmap(cnt,2)
        ia=zmap(cnt,3)        
        v_i_recmpi(cnt)=v_i(ieps,ik,ia)
        g_i_o_recmpi(cnt)=g_i_o(ieps,ik,ia)
        end do            

    end if
    
d_v=1
d_v_mpi=1
iter_v=0

    do while ( iter_v .le. itmax_v .and. d_v .ge. tol_v )

        !bond price schedule given operate/exit decisions    
        do ieps=1,neps
        do jk=1,nk
        do ja=1,na
            
            if ( ja .lt. intbirth ) then
            
                if ( default_economy .eq. 0 ) then
                
                g_q(ieps,jk,ja)=one/(one+r)-tau_a
                
                else if ( default_economy .eq. 1 ) then                
                
                g_q(ieps,jk,ja)=zero                
                
                    do jeps=1,neps
                    
                        if ( pi_eps(ieps,jeps) .gt. zero ) then
                        g_q(ieps,jk,ja)=g_q(ieps,jk,ja)+pi_eps(ieps,jeps)*&
                        (g_i_o(jeps,jk,ja)*(-grid_a(date,ja))+&
                        (1-g_i_o(jeps,jk,ja))*min(-grid_a(date,ja),max(zero,(1.+tau_x)*grid_k(jk)*chi-func_adj(zero,grid_k(jk),g_i_o(jeps,jk,ja),date))))
                        end if
                        
                    end do
                
                g_q(ieps,jk,ja)=g_q(ieps,jk,ja)/((one+r)*(-grid_a(date,ja)))-tau_a
                
                end if
                
            else if ( ja .ge. intbirth ) then
            
            g_q(ieps,jk,ja)=one/(one+r)
            
            end if
        
        end do
        end do
        end do
    
        !expectations
        allocate (e_v_i_mpi(np))
    
        do cnt=itop,iend
    
        ieps=zmap(cnt,1)
        jk=zmap(cnt,2)
        ja=zmap(cnt,3)

        e_v_i_mpi(cnt-itop+1)=zero
    
            do jeps=1,neps
                if ( pi_eps(ieps,jeps) .gt. zero ) then
                e_v_i_mpi(cnt-itop+1)=e_v_i_mpi(cnt-itop+1)+pi_eps(ieps,jeps)*v_i_recmpi(zmap_inv(jeps,jk,ja))
                end if
            end do

        end do

        call mpi_allgather (e_v_i_mpi,np,mpi_double_precision,e_v_i_recmpi,np,mpi_double_precision,mpi_comm_world,ierr)
    
        deallocate (e_v_i_mpi)

        !incumbent's problem
        allocate (g_i_o_mpi(np),g_i_k_mpi(np),g_i_a_mpi(np),v_i_mpi(np),g_binding_collateral_mpi(np),g_i_a_integer_mpi(np))
    
        do cnt=itop,iend
    
        ieps=zmap(cnt,1)
        ik=zmap(cnt,2)
        ia=zmap(cnt,3)
        iz=id_z(ieps)
        
        eps=grid_eps(ieps)
        k=grid_k(ik)
        a=grid_a(date,ia)

        g_binding_collateral_mpi(cnt-itop+1)=1

        find=0
    
            !unconstrained
            if ( grid_k(g_un_k(ieps,ik)) .lt. zero ) then
            print*,'g_un_k(ieps,ik) .lt. zero, vfn_i, mini.f90',g_un_k(ieps,ik)
                call exit
            end if
            
        jk=g_un_k(ieps,ik)
        d_o=g_un_o(ieps,ik)
        kp=grid_k(jk)
            
            if ( d_o .eq. 0 ) then  
                
            ja=intbirth    
            ap=grid_a(date,ja) 
            
                if ( kp .ne. kmin ) then
                print*,'kp .ne. kmin, vfn_i, mini.f90'
                    call exit
                end if

                if ( default_economy .eq. 0 ) then            
                vtemp=(1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date)                    
                else if ( default_economy .eq. 1 ) then
                vtemp=max((1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date),min(zero,(1.+tau_x)*k*chi-func_adj(kp,k,d_o,date)))                                                    
                end if
                                
                if ( find .eq. 0 ) then
        
                find=1
                g_i_o_mpi(cnt-itop+1)=d_o
                g_i_k_mpi(cnt-itop+1)=jk
                g_i_a_mpi(cnt-itop+1)=ap
                v_i_mpi(cnt-itop+1)=vtemp
                g_binding_collateral_mpi(cnt-itop+1)=0
                g_i_a_integer_mpi(cnt-itop+1)=-1
            
                else if ( find .eq. 1 ) then

                print*,'find .eq. 1, vfn_i, mini.f90'
                    call exit
                
                end if
                
            else if ( d_o .eq. 1 ) then

            l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))

                if ( isnan(l) ) then
                print*,'isnan l, vfn_i, mini.f90 ',l
                    call exit
                end if
            
            ap=a_un(ieps,ik)             
            profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))+&
            (1.+tau_x)*((1-delta_k)*k-kp)+a-ap/(1+r)-func_adj(kp,k,d_o,date)-f_o_good
            
                if ( profit .ge. zero ) then
!                 if ( profit .ge. zero .or. profit .lt. zero ) then

                exp_ap=e_v_i_recmpi(zmap_inv(ieps,jk,:))                       
                    call onedlin (grid_a(date,:),exp_ap,na,ap,ev)
                    
                vtemp=profit+(1-delta_f)*beta*ev
                
                    if ( find .eq. 0 ) then
                    
                    find=1
                    g_i_o_mpi(cnt-itop+1)=d_o
                    g_i_k_mpi(cnt-itop+1)=jk
                    g_i_a_mpi(cnt-itop+1)=ap
                    v_i_mpi(cnt-itop+1)=vtemp
                    g_binding_collateral_mpi(cnt-itop+1)=0
                    g_i_a_integer_mpi(cnt-itop+1)=-1
                    
                    else if ( find .eq. 1 ) then
                    
                    print*,'find .eq. 1, vfn_i, mini.f90'
                        call exit
                        
                    end if
                    
                else if ( profit .lt. zero ) then
                
                !exit
                d_o=0
                jk=1
                ja=intbirth            
                kp=grid_k(jk)
                ap=grid_a(date,ja)

                    if ( default_economy .eq. 0 ) then            
                    vtemp=(1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date)                    
                    else if ( default_economy .eq. 1 ) then
                    vtemp=max((1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date),min(zero,(1.+tau_x)*k*chi-func_adj(kp,k,d_o,date)))                                    
                    end if
                
                    if ( find .eq. 0 ) then
        
                    find=1
                    g_i_o_mpi(cnt-itop+1)=d_o
                    g_i_k_mpi(cnt-itop+1)=jk
                    g_i_a_mpi(cnt-itop+1)=ap
                    v_i_mpi(cnt-itop+1)=vtemp
                    g_i_a_integer_mpi(cnt-itop+1)=-1
        
                    else if ( find .eq. 1 ) then

                    print*,'find .eq. 1, vfn_i, mini.f90'
                        call exit
                
                    end if
                
                !operate
                d_o=1
                l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))

                    if ( isnan(l) ) then
                    print*,'isnan l, mini.f90 ',k
                        call exit
                    end if

                kpmin=kmin
                kpmax=kmax
                                    
                    !operate (constrained)            
                    do jk=1,nk
                
                    kp=grid_k(jk)

                    constant=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))+&
                    (1.+tau_x)*((1-delta_k)*k-kp)+a-func_adj(kp,k,d_o,date)-f_o_good 
                    
                        if ( constant .ge. zero ) then
                        
                        q=one/(one+r)
                        ap=constant/q
                        
                        else if ( constant .lt. zero ) then

                        q=one/(one+r)-tau_a
                        ap=constant/q
                        
                        end if
                        
                        if ( default_economy .eq. 0 .or. constant .ge. zero ) then                            
                         
                        ap=min(grid_a(date,na),ap)
                                                
                            if ( ap .ge. max(abar-theta*k,grid_a(date,1)) .and. kp .ge. kpmin .and. kp .le. kpmax ) then                                    

                            profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))+&
                            (1.+tau_x)*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,date)-f_o_good                            
                            
                            exp_ap=e_v_i_recmpi(zmap_inv(ieps,jk,:))                         
                                call onedlin (grid_a(date,:),exp_ap,na,ap,ev)
                    
                            vtemp=profit+(1-delta_f)*beta*ev
                            
                                if ( find .eq. 0 ) then

                                print*,'find .eq. 0, mini.f90'                                
                                    call exit
                                        
                                else if ( find .eq. 1 ) then
                    
                                    if ( vtemp .gt. v_i_mpi(cnt-itop+1) ) then
                                    g_i_o_mpi(cnt-itop+1)=d_o
                                    g_i_k_mpi(cnt-itop+1)=jk
                                    g_i_a_mpi(cnt-itop+1)=ap
                                    v_i_mpi(cnt-itop+1)=vtemp
                                    g_i_a_integer_mpi(cnt-itop+1)=-1
                                    end if
                        
                                end if
                                
                            g_binding_collateral_mpi(cnt-itop+1)=0
                        
                            end if
                        
                        else if ( default_economy .eq. 1 .and. constant .lt. zero ) then

                            do ja=1,intbirth
                            
                            q=g_q(ieps,jk,ja)
                            ap=grid_a(date,ja)
                            
                                if ( kp .ge. kpmin .and. kp .le. kpmax .and. constant-q*ap .ge. zero ) then                                    

                                profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))+&
                                (1.+tau_x)*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,date)-f_o_good                            
                        
                                exp_ap=e_v_i_recmpi(zmap_inv(ieps,jk,:))                         
                                    call onedlin (grid_a(date,:),exp_ap,na,ap,ev)
                    
                                vtemp=profit+(1-delta_f)*beta*ev
                            
                                    if ( find .eq. 0 ) then

                                    print*,'find .eq. 0, mini.f90'                                
                                        call exit
                                        
                                    else if ( find .eq. 1 ) then
                    
                                        if ( vtemp .gt. v_i_mpi(cnt-itop+1) ) then
                                        g_i_o_mpi(cnt-itop+1)=d_o
                                        g_i_k_mpi(cnt-itop+1)=jk
                                        g_i_a_mpi(cnt-itop+1)=ap
                                        v_i_mpi(cnt-itop+1)=vtemp
                                        g_i_a_integer_mpi(cnt-itop+1)=ja                                        
                                        end if
                        
                                    end if
                                
                                g_binding_collateral_mpi(cnt-itop+1)=0
                        
                                end if                            
                                
                            end do
                        
                        end if
    
                    end do
                                            
                d_o=g_i_o_mpi(cnt-itop+1)
                
                    if ( d_o .eq. 1 ) then
                    
                    jk=g_i_k_mpi(cnt-itop+1)                    
                    ap=g_i_a_mpi(cnt-itop+1)
                    kp=grid_k(jk)                        

                        if ( ap .ge. zero ) then
                                
                        q=one/(one+r)
                                
                        else if ( ap .lt. zero ) then
                                
                            if ( default_economy .eq. 0 ) then                                
                            q=one/(one+r)-tau_a
                            else if ( default_economy .eq. 1 ) then
                            ja=g_i_a_integer_mpi(cnt-itop+1)
                            
                                if ( ja .le. 0 ) then
                                print*,'ja .le. 0'
                                    call exit
                                end if
                                
                            q=g_q(ieps,jk,ja)
                            end if
                            
                        end if
                                
                    profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))&
                    +(1.+tau_x)*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,date)-f_o_good                         
                    
                        if ( profit .lt. -1e-10 ) then
                        print*,'profit .lt. -1e-10, vfn_i, mini.f90',profit
                        print*,'iter_v',iter_v
                        print*,'q',q
                        print*,'ap',ap                        
                            call exit
                        end if
                        
                    end if
                
                end if

            end if
            
            if ( g_binding_collateral_mpi(cnt-itop+1) .eq. 1 .and. g_i_o_mpi(cnt-itop+1) .eq. 1 ) then
            print*,'g_binding_collateral_mpi(cnt-itop+1) .eq. 1 .and. g_i_o_mpi(cnt-itop+1) .eq. 1, mini.f90'
                call exit
            end if
            
            if ( cnt .eq. itop ) then      
            
            d_v_mpi=max(abs(v_i_mpi(cnt-itop+1)-v_i_recmpi(zmap_inv(ieps,ik,ia)))&
            ,abs(g_i_a_mpi(cnt-itop+1)-g_i_a_recmpi(zmap_inv(ieps,ik,ia)))&
            ,real(abs(g_i_k_mpi(cnt-itop+1)-g_i_k_recmpi(zmap_inv(ieps,ik,ia))))&
            ,abs(real(g_i_o_mpi(cnt-itop+1)-g_i_o_recmpi(zmap_inv(ieps,ik,ia)))))
                        
            else if ( cnt .ne. itop ) then            

            dummy=max(abs(v_i_mpi(cnt-itop+1)-v_i_recmpi(zmap_inv(ieps,ik,ia)))&
            ,abs(g_i_a_mpi(cnt-itop+1)-g_i_a_recmpi(zmap_inv(ieps,ik,ia)))&
            ,real(abs(g_i_k_mpi(cnt-itop+1)-g_i_k_recmpi(zmap_inv(ieps,ik,ia))))&
            ,abs(real(g_i_o_mpi(cnt-itop+1)-g_i_o_recmpi(zmap_inv(ieps,ik,ia)))))
            
                if ( dummy .gt. d_v_mpi ) then
                d_v_mpi=dummy
                end if
                
            end if

            
        end do

        call mpi_allgather (g_i_o_mpi,np,mpi_integer1,g_i_o_recmpi,np,mpi_integer1,mpi_comm_world,ierr)    
        call mpi_allgather (g_i_k_mpi,np,mpi_integer2,g_i_k_recmpi,np,mpi_integer2,mpi_comm_world,ierr)    
        call mpi_allgather (g_i_a_mpi,np,mpi_double_precision,g_i_a_recmpi,np,mpi_double_precision,mpi_comm_world,ierr)        
        call mpi_allgather (v_i_mpi,np,mpi_double_precision,v_i_recmpi,np,mpi_double_precision,mpi_comm_world,ierr)
        call mpi_allgather (d_v_mpi,1,mpi_double_precision,d_v_recmpi,1,mpi_double_precision,mpi_comm_world,ierr)
        call mpi_allgather (g_binding_collateral_mpi,np,mpi_integer1,g_binding_collateral_recmpi,np,mpi_integer1,mpi_comm_world,ierr)    
        call mpi_allgather (g_i_a_integer_mpi,np,mpi_integer,g_i_a_integer_recmpi,np,mpi_integer,mpi_comm_world,ierr)        

        deallocate (g_i_o_mpi,g_i_k_mpi,g_i_a_mpi,v_i_mpi,g_binding_collateral_mpi,g_i_a_integer_mpi)

    d_v=maxval(d_v_recmpi,1)
    iter_v=iter_v+1
    
!         if ( rank .eq. 0 ) then
!         print*,'iter_v, d_v ',iter_v,d_v
!         end if
        
        !value function iteration without updating policy function
        do iter_re=1,15
        
            !expectations
            allocate (e_v_i_mpi(np))
    
            do cnt=itop,iend
    
            ieps=zmap(cnt,1)
            jk=zmap(cnt,2)
            ja=zmap(cnt,3)

            e_v_i_mpi(cnt-itop+1)=zero
    
                do jeps=1,neps
                
                    if ( pi_eps(ieps,jeps) .gt. zero ) then
                    e_v_i_mpi(cnt-itop+1)=e_v_i_mpi(cnt-itop+1)+pi_eps(ieps,jeps)*v_i_recmpi(zmap_inv(jeps,jk,ja))
                    end if
                    
                end do

            end do

            call mpi_allgather (e_v_i_mpi,np,mpi_double_precision,e_v_i_recmpi,np,mpi_double_precision,mpi_comm_world,ierr)
    
            deallocate (e_v_i_mpi)
        
            !update value function
            allocate (v_i_mpi(np))
        
            do cnt=itop,iend
    
            ieps=zmap(cnt,1)
            ik=zmap(cnt,2)
            ia=zmap(cnt,3)
            iz=id_z(ieps)
    
            eps=grid_eps(ieps)
            k=grid_k(ik)
            a=grid_a(date,ia)
        
            d_o=g_i_o_recmpi(zmap_inv(ieps,ik,ia))
            jk=g_i_k_recmpi(zmap_inv(ieps,ik,ia))
            ap=g_i_a_recmpi(zmap_inv(ieps,ik,ia))
            kp=grid_k(jk)
        
                if ( d_o .eq. 0 ) then

                    if ( default_economy .eq. 0 ) then            
                    v_i_mpi(cnt-itop+1)=(1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date)                    
                    else if ( default_economy .eq. 1 ) then
                    v_i_mpi(cnt-itop+1)=max((1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date),min(zero,(1.+tau_x)*k*chi-func_adj(kp,k,d_o,date)))                                    
                    end if
                
                else if ( d_o .eq. 1 ) then
            
                l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))

                    if ( ap .ge. zero ) then
                                
                    q=one/(one+r)
                                
                    else if ( ap .lt. zero ) then
                                
                        if ( default_economy .eq. 0 ) then                                
                        q=one/(one+r)-tau_a
                        else if ( default_economy .eq. 1 ) then
                        ja=g_i_a_integer_recmpi(zmap_inv(ieps,ik,ia))

                            if ( ja .le. 0 ) then
                            print*,'ja .le. 0'
                                call exit
                            end if
                        
                        q=g_q(ieps,jk,ja)
                        end if
                            
                    end if
                                
                profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))+&
                (1.+tau_x)*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,date)-f_o_good
                    
                exp_ap=e_v_i_recmpi(zmap_inv(ieps,jk,:))  
                    
                    call onedlin (grid_a(date,:),exp_ap,na,ap,ev)
                        
                v_i_mpi(cnt-itop+1)=profit+(1-delta_f)*beta*ev
            
                end if
                
            end do
        
            call mpi_allgather (v_i_mpi,np,mpi_double_precision,v_i_recmpi,np,mpi_double_precision,mpi_comm_world,ierr)
            
            deallocate (v_i_mpi)
        
        end do
        
    end do
    
    do cnt=1,ni
    
    ieps=zmap(cnt,1)
    ik=zmap(cnt,2)
    ia=zmap(cnt,3)
    
    v_i(ieps,ik,ia)=v_i_recmpi(zmap_inv(ieps,ik,ia))
    g_i_o(ieps,ik,ia)=g_i_o_recmpi(zmap_inv(ieps,ik,ia))
    g_i_k(ieps,ik,ia)=g_i_k_recmpi(zmap_inv(ieps,ik,ia))
    g_i_a(ieps,ik,ia)=g_i_a_recmpi(zmap_inv(ieps,ik,ia))
    g_binding_collateral(ieps,ik,ia)=g_binding_collateral_recmpi(zmap_inv(ieps,ik,ia))
    g_i_a_integer(ieps,ik,ia)=g_i_a_integer_recmpi(zmap_inv(ieps,ik,ia))

    end do

    if ( rank .eq. 0 ) then

        if ( d_v .ge. tol_v ) then
        print*,'d_v .ge. tol_v, mini.f90 ',d_v
!             call exit
        end if
        
        if ( date .eq. 1 ) then
        
            open (unit=1,file='output_vfn_pfn.txt')        
            write(1,*)'eps k a v kp ap d'
                do ieps=5,neps,5       
                do ik=1,nk
                do ia=1,na
                write(1,'(300F12.7)')grid_eps(ieps),grid_k(ik),grid_a(date,ia),&
                v_i(ieps,ik,ia),grid_k(g_i_k(ieps,ik,ia)),g_i_a(ieps,ik,ia),real(g_i_o(ieps,ik,ia))
                end do            
                end do
                end do
            close(1)
            
            open (unit=1,file='output_g_q.txt')
            write(1,*)'ap q iz ieta iiota'
                do ja=1,na
                ieps=0
                    do iiota=1,niota
                    do iz=1,nz
                    do ieta=1,neta
                    ieps=ieps+1
                    write(1,'(300F12.7)')grid_a(date,ja),g_q(ieps,intbirth_k,ja),real(iz),real(ieta),real(iiota)
                    end do
                    end do
                    end do
                end do
            close(1)
        
        end if
        
    end if

    deallocate (g_i_o_recmpi,g_i_k_recmpi,g_i_a_recmpi,v_i_recmpi,e_v_i_recmpi,g_binding_collateral_recmpi,g_i_a_integer_recmpi)

    call MPI_BCAST (alpha,1,mpi_double_precision,0,mpi_comm_world,ierr)

    if ( rank .eq. 0 ) then
    print*,'sum g_binding_collateral',sum(g_binding_collateral)
    end if
    
end subroutine vfn_i

subroutine vfn_un
implicit none
include 'mpif.h'

integer(kind=1), allocatable, dimension(:) :: g_un_o_mpi
integer(kind=2), allocatable, dimension(:) :: g_un_k_mpi
double precision, allocatable, dimension(:) :: v_un_mpi, e_v_un_mpi

integer(kind=1), allocatable, dimension(:) :: g_un_o_recmpi
integer(kind=2), allocatable, dimension(:) :: g_un_k_recmpi
double precision, allocatable, dimension(:) :: v_un_recmpi, e_v_un_recmpi

double precision :: vtemp, ev, vtemp0
integer :: find

integer :: iter_v
double precision :: d_v, dummy
double precision :: d_v_mpi, d_v_recmpi(nproc)

integer :: iter_re
double precision, dimension(neps,nk) :: y2
integer :: jk_temp

    call mpi_barrier(mpi_comm_world,ierr)

    allocate (g_un_o_recmpi(np_un*nproc),g_un_k_recmpi(np_un*nproc),v_un_recmpi(np_un*nproc),e_v_un_recmpi(np_un*nproc))
    
!initial guess for value function
    if ( read_v_un .eq. 0 ) then

    read_v_un=0
    
        allocate (v_un_mpi(np_un))
    
        do cnt=itop_un,iend_un
    
        ieps=zmap_un(cnt,1)
        ik=zmap_un(cnt,2)
        iz=id_z(ieps)

        eps=grid_eps(ieps)
        k=grid_k(ik)
    
        !exit
        d_o=0
        jk=1
        kp=grid_k(jk)
        v_un_mpi(cnt-itop_un+1)=(1.+tau_x)*k*chi-func_adj(kp,k,d_o,date)
            
        !operate
        d_o=1
        jk=1
        kp=grid_k(jk)
        l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
        profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))+&
        (1.+tau_x)*(1-delta_k)*k-func_adj(kp,k,d_o,date)-f_o_good
    
        v_un_mpi(cnt-itop_un+1)=max(profit,v_un_mpi(cnt-itop_un+1))
                
        end do

        call mpi_allgather (v_un_mpi,np_un,mpi_double_precision,v_un_recmpi,np_un,mpi_double_precision,mpi_comm_world,ierr)
    
        deallocate (v_un_mpi)

    else if ( read_v_un .eq. 1 ) then
    
        do cnt=1,ni_un        
        ieps=zmap_un(cnt,1)
        ik=zmap_un(cnt,2)
        v_un_recmpi(cnt)=v_un(ieps,ik)
        end do            

    end if

d_v=1
d_v_recmpi=0
v_un_recmpi=0
iter_v=0

    do while ( iter_v .le. 1500 .and. d_v .ge. tol_v )
    
        !expectations
        allocate (e_v_un_mpi(np_un))
    
        do cnt=itop_un,iend_un
    
        ieps=zmap_un(cnt,1)
        jk=zmap_un(cnt,2)

        e_v_un_mpi(cnt-itop_un+1)=zero
    
            do jeps=1,neps
                if ( pi_eps(ieps,jeps) .gt. zero ) then
                e_v_un_mpi(cnt-itop_un+1)=e_v_un_mpi(cnt-itop_un+1)+pi_eps(ieps,jeps)*v_un_recmpi(zmap_inv_un(jeps,jk))
                end if
            end do

        end do

        call mpi_allgather (e_v_un_mpi,np_un,mpi_double_precision,e_v_un_recmpi,np_un,mpi_double_precision,mpi_comm_world,ierr)
    
        deallocate (e_v_un_mpi)

        !incumbent's problem
        allocate (g_un_o_mpi(np_un),g_un_k_mpi(np_un),v_un_mpi(np_un))
        
        do cnt=itop_un,iend_un
    
        ieps=zmap_un(cnt,1)
        ik=zmap_un(cnt,2)
        iz=id_z(ieps)
        
        eps=grid_eps(ieps)
        k=grid_k(ik)

        find=0
                
        !exit
        d_o=0
        jk=1
        kp=grid_k(1)
        vtemp=(1.+tau_x)*k*chi-func_adj(kp,k,d_o,date)
        
            if ( find .eq. 0 ) then

            find=1
            g_un_o_mpi(cnt-itop_un+1)=d_o
            g_un_k_mpi(cnt-itop_un+1)=jk
            v_un_mpi(cnt-itop_un+1)=vtemp
            
            else if ( find .eq. 1 ) then

            print*,'find .eq. 1, vfn_un,mini.f90'
                call exit
                            
            end if
                
        !operate
        d_o=1    
        l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
        
            if ( isnan(l) ) then
            print*,'isnan l, vfn_un, mini.f90 ',k
                call exit
            end if

        kpmin=kmin
        kpmax=kmax                                          
        kp_temp=-one
        jk_temp=-1
        
            do jk=1,nk

            kp=grid_k(jk)
            
                if ( kp .ge. kpmin .and. kp .le. kpmax ) then
                            
                profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))+&
                (1.+tau_x)*((1-delta_k)*k-kp)-func_adj(kp,k,d_o,date)-f_o_good                        
                
                ev=e_v_un_recmpi(zmap_inv_un(ieps,jk))
                vtemp0=profit+(1-delta_f)*beta*ev
                
                    if ( kp_temp .lt. zero ) then
                
                    kp_temp=kp
                    jk_temp=jk
                    vtemp=vtemp0
                
                    else if ( kp_temp .ge. zero ) then
                
                        if ( vtemp0 .gt. vtemp ) then
                        kp_temp=kp
                        jk_temp=jk
                        vtemp=vtemp0
                        end if
                    
                    end if
            
                end if
                
            end do

            if ( kp_temp .ge. zero ) then
                        
            kp=kp_temp
            jk=jk_temp
            profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))+&
            (1.+tau_x)*((1-delta_k)*k-kp)-func_adj(kp,k,d_o,date)-f_o_good        

            ev=e_v_un_recmpi(zmap_inv_un(ieps,jk))
            vtemp=profit+(1-delta_f)*beta*ev
                
                if ( find .eq. 0 ) then

                print*,'find .eq. 0, vfn_un, mini.f90'
                    call exit
                    
                else if ( find .eq. 1 ) then
                    
                    if ( vtemp .gt. v_un_mpi(cnt-itop_un+1) ) then
                    g_un_o_mpi(cnt-itop_un+1)=d_o
                    g_un_k_mpi(cnt-itop_un+1)=jk
                    v_un_mpi(cnt-itop_un+1)=vtemp
                    end if
                        
                end if
                
            else if ( kp_temp .lt. zero ) then
            
            print*,'kp_temp .lt. zero, mini.f90'
                call exit

            end if
                                  
            if ( cnt .eq. itop_un ) then      
            
            d_v_mpi=max(abs(v_un_mpi(cnt-itop_un+1)-v_un_recmpi(zmap_inv_un(ieps,ik)))&
            ,real(abs(g_un_k_mpi(cnt-itop_un+1)-g_un_k_recmpi(zmap_inv_un(ieps,ik))))&
            ,real(abs(g_un_o_mpi(cnt-itop_un+1)-g_un_o_recmpi(zmap_inv_un(ieps,ik)))))
                        
            else if ( cnt .ne. itop_un ) then            

            dummy=max(abs(v_un_mpi(cnt-itop_un+1)-v_un_recmpi(zmap_inv_un(ieps,ik)))&
            ,real(abs(g_un_k_mpi(cnt-itop_un+1)-g_un_k_recmpi(zmap_inv_un(ieps,ik))))&
            ,real(abs(g_un_o_mpi(cnt-itop_un+1)-g_un_o_recmpi(zmap_inv_un(ieps,ik)))))
            
                if ( dummy .gt. d_v_mpi ) then
                d_v_mpi=dummy
                end if
                
            end if

        end do

        call mpi_allgather (g_un_o_mpi,np_un,mpi_integer1,g_un_o_recmpi,np_un,mpi_integer1,mpi_comm_world,ierr)    
        call mpi_allgather (g_un_k_mpi,np_un,mpi_integer2,g_un_k_recmpi,np_un,mpi_integer2,mpi_comm_world,ierr)    
        call mpi_allgather (v_un_mpi,np_un,mpi_double_precision,v_un_recmpi,np_un,mpi_double_precision,mpi_comm_world,ierr)
        call mpi_allgather (d_v_mpi,1,mpi_double_precision,d_v_recmpi,1,mpi_double_precision,mpi_comm_world,ierr)

        deallocate (g_un_o_mpi,g_un_k_mpi,v_un_mpi)

    d_v=maxval(d_v_recmpi,1)
    iter_v=iter_v+1    
        
    end do
    
    do cnt=1,ni_un
    
    ieps=zmap_un(cnt,1)
    ik=zmap_un(cnt,2)
    
    v_un(ieps,ik)=v_un_recmpi(zmap_inv_un(ieps,ik))
    g_un_o(ieps,ik)=g_un_o_recmpi(zmap_inv_un(ieps,ik))
    g_un_k(ieps,ik)=g_un_k_recmpi(zmap_inv_un(ieps,ik))

        if ( cnt .ne. zmap_inv_un(ieps,ik) ) then
        print*,'error in parallel, vfn_un, check uniform_parallel'
            call exit
        end if

    end do

    if ( rank .eq. 0 ) then

        if ( d_v .ge. tol_v ) then
        print*,'d_v .ge. tol_v in vfn iteration, vfn_un, mini.f90 ',d_v
!             call exit
        end if
        
    end if

    deallocate (g_un_o_recmpi,g_un_k_recmpi,v_un_recmpi,e_v_un_recmpi)
        
    if ( rank .eq. 0 ) then
    
        open(unit=1,file='output_pfn_unconstrained.txt')
            do ik=1,nk
            write(1,'(300F12.7)')grid_k(ik),(grid_k(g_un_k(ieps,ik)),ieps=5,neps,5)
            end do
        close(1)
            
        open(unit=1,file='output_profit_unconstrained.txt')
            do ik=1,nk
            write(1,'(300F12.7)')grid_k(ik),(((Z*grid_eps(ieps))**(1-alpha*nu))*(((grid_k(ik)**alpha)*&
                (((w/(((Z*grid_eps(ieps))**(1-alpha*nu))*(1-alpha)*nu*(grid_k(ik)**(alpha*nu))))**(one/((1-alpha)*nu-1)))**(1-alpha)))**nu)&
                -w*(((w/(((Z*grid_eps(ieps))**(1-alpha*nu))*(1-alpha)*nu*(grid_k(ik)**(alpha*nu))))**(one/((1-alpha)*nu-1)))+f_o*grid_z(id_z(ieps)))+&
                (1.+tau_x)*((1-delta_k)*grid_k(ik)-grid_k(g_un_k(ieps,ik)))-f_o_good,ieps=5,neps,5)
            end do
        close(1)

    end if
    
!solve for minimum saving rule
a_un=zero       !guess for minimum saving rule given ieps and ik
a_tilde1=zero
d_v=1
iter_v=0

    do while ( d_v .ge. tol_v .and. iter_v .le. 1500 )

        !minimum savings to be unconstrained in each state in current period
        do ieps=1,neps
        do ik=1,nk
    
        iz=id_z(ieps)
        eps=grid_eps(ieps)
        k=grid_k(ik)
        d_o=g_un_o(ieps,ik)
        jk=g_un_k(ieps,ik)
        kp=grid_k(jk)
                        
            if ( d_o .eq. 0 ) then

            a_tilde(ieps,ik)=zero
            
            else if ( d_o .eq. 1 ) then
            
            l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))                            
            
            a_tilde(ieps,ik)=a_un(ieps,ik)/(1+r)-((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)+w*(l+f_o*grid_z(iz))-&
            (1.+tau_x)*((1-delta_k)*k-kp)+func_adj(kp,k,d_o,date)+f_o_good
            
            end if

        end do
        end do
    
        !solve for savings to be unconstrained in next period
        do ieps=1,neps
        do ik=1,nk

        eps=grid_eps(ieps)
        k=grid_k(ik)
        d_o=g_un_o(ieps,ik)
        jk=g_un_k(ieps,ik)
        kp=grid_k(jk)
        
        find=0
    
            do jeps=1,neps
        
                if ( pi_eps(ieps,jeps) .gt. zero ) then
            
                    if ( find .eq. 0 ) then
                    
                    find=1
                    a_un1(ieps,ik)=a_tilde(jeps,jk)
                    
                    else if ( find .eq. 1 ) then
                
                    dummy=a_tilde(jeps,jk)
                    
                        if ( dummy .ge. a_un1(ieps,ik) ) then
                        a_un1(ieps,ik)=dummy
                        end if
                    
                    end if
        
                end if
            
            end do
        
        a_un1(ieps,ik)=max(a_un1(ieps,ik),zero)
        
        end do
        end do
        
    d_v=max(maxval(abs(a_un1-a_un)),maxval(abs(a_tilde1-a_tilde)))
    iter_v=iter_v+1
    
        if ( d_v .ge. tol_v ) then
        a_un=a_un1
        a_tilde1=a_tilde
        end if
        
    end do
    
    if ( rank .eq. 0 ) then

        if ( d_v .ge. tol_v ) then
        print*,'d_v .ge. tol_v in minimum savings rule, vfn_un, mini.f90 ',d_v
!            call exit
        end if
    
        open(unit=1,file='output_unconstrained_saving_rule.txt')
            do ik=1,nk
            write(1,'(300F13.7)')grid_k(ik),(a_un(ieps,ik),ieps=5,neps,5),(a_tilde(ieps,ik),ieps=5,neps,5)
            end do
        close(1)
                
    end if
    
    if ( maxval(a_un(1:neta*nz,:)) .gt. amax .and. rank .eq. 0 ) then
    print*,'increase amax, vfn_un, mini.f90',maxval(a_un(1:neta*nz,:)),amax
       call exit
    end if
    
    if ( rank .eq. 0 ) then
    print*,'maxval(a_un(1:neta*nz,:))',maxval(a_un(1:neta*nz,:))
    end if

    do ieps=1,neps
    do ik=1,nk
    a_un(ieps,ik)=min(amax,a_un(ieps,ik))
    end do
    end do
    
end subroutine vfn_un

SUBROUTINE rouwenhorst(n,rho,mu_z,var_e,Z,Pi)
IMPLICIT NONE
!z_{t+1}=(1-rho)*mu_z+rho*z_t+eps_t
!z_{t+1}-mu_z=rho*(z_t-mu_z)+eps_t
!y_{t+1}=rho*y_t+eps_t
!rho: persistence of AR(1) process
!mu_z: Mean of z
!var_e: Variance of epsilon
!n: number of finite states
DOUBLE PRECISION, INTENT(IN) :: rho,mu_z,var_e
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(OUT):: Z(n,1),Pi(n,n)
DOUBLE PRECISION :: p
DOUBLE PRECISION :: sig_z, psi
DOUBLE PRECISION :: Y(n,1)
INTEGER :: i,j,k
DOUBLE PRECISION :: T(n,n,n),T1(n,n,n),T2(n,n,n),T3(n,n,n),T4(n,n,n)
sig_z=SQRT(var_e/(1.-(rho**2.)))
psi=(SQRT(n-1.))*sig_z
p=(1.+rho)/2.

	DO i=1,n
	Y(i,1)=-psi+2.*psi*(i-1.)/(n-1.)
	Z(i,1)=Y(i,1)+mu_z
	ENDDO
	
	if ( n .eq. 1 ) then
	pi(1,1)=1
	z(1,1)=mu_z
	end if

	IF (n.EQ.2) THEN
	Pi(1,1)=p
	Pi(1,2)=1.-p
	Pi(2,1)=1.-p
	Pi(2,2)=p
	ENDIF
T=0.
T1=0.
T2=0.
T3=0.
T4=0.
	IF (n.GT.2) THEN
	T(2,1,1)=p
	T(2,1,2)=1.-p
	T(2,2,1)=1.-p
	T(2,2,2)=p
		DO i=3,n
			DO j=1,i-1
				DO k=1,i-1
				T1(i,j,k)=p*T(i-1,j,k)
				T2(i,j,k+1)=(1.-p)*T(i-1,j,k)
				T3(i,j+1,k)=(1.-p)*T(i-1,j,k)
				T4(i,j+1,k+1)=p*T(i-1,j,k)
				ENDDO		
			ENDDO
		T(i,:,:)=T1(i,:,:)+T2(i,:,:)+T3(i,:,:)+T4(i,:,:)
			DO j=2,i-1
			T(i,j,:)=T(i,j,:)/2.	
			ENDDO
		ENDDO
	
		Pi(:,:)=T(n,:,:)

	ENDIF


END SUBROUTINE rouwenhorst

subroutine uniform_parallel
implicit none
include 'mpif.h'

real, allocatable, dimension(:) :: rn
integer, dimension(3) :: zmaptemp0
integer :: cnt0

integer, dimension(2) :: zmaptemp0_un

! gridpoints per processor
np=(ni-1)/nproc+1
itop=rank*np+1
iend=MIN((rank+1)*np,ni)

    cnt0=0
    
    allocate (rn(ni-1))
	
    do ieps=1,neps            
    do ik=1,nk
    do ia=1,na
    
    cnt0=cnt0+1
    zmap(cnt0,1)=ieps
    zmap(cnt0,2)=ik
    zmap(cnt0,3)=ia
    
        if ( cnt0 .le. ni-1 ) then
            call random_number(rn(cnt0))
        end if
            
    end do
    end do
    end do
    
    do cnt=ni-1,1,-1
      
    zmaptemp0(:)=zmap(int(rn(cnt)*cnt)+1,:)
    zmap(int(rn(cnt)*cnt)+1,:)=zmap(cnt+1,:)
    zmap(cnt+1,:)=zmaptemp0(:)
    zmap_inv(zmap(cnt+1,1),zmap(cnt+1,2),zmap(cnt+1,3))=cnt+1
                 
    end do

    deallocate (rn)
     
    cnt=1
zmap_inv(zmap(cnt,1),zmap(cnt,2),zmap(cnt,3))=cnt

! gridpoints per processor
np_un=(ni_un-1)/nproc+1
itop_un=rank*np_un+1
iend_un=MIN((rank+1)*np_un,ni_un)

    cnt0=0
    
    allocate (rn(ni_un-1))
	
    do ieps=1,neps            
    do ik=1,nk
    
    cnt0=cnt0+1
    zmap_un(cnt0,1)=ieps
    zmap_un(cnt0,2)=ik
    
        if ( cnt0 .le. ni_un-1 ) then
            call random_number(rn(cnt0))
        end if
            
    end do
    end do
    
    do cnt=ni_un-1,1,-1
      
    zmaptemp0_un(:)=zmap_un(int(rn(cnt)*cnt)+1,:)
    zmap_un(int(rn(cnt)*cnt)+1,:)=zmap_un(cnt+1,:)
    zmap_un(cnt+1,:)=zmaptemp0_un(:)
    zmap_inv_un(zmap_un(cnt+1,1),zmap_un(cnt+1,2))=cnt+1
                 
    end do

    deallocate (rn)
     
    cnt=1
zmap_inv_un(zmap_un(cnt,1),zmap_un(cnt,2))=cnt

end subroutine uniform_parallel

FUNCTION brent(ax,bx,cx,tol,xmin,func)
IMPLICIT NONE

double precision, INTENT(IN) :: ax,bx,cx,tol
double precision, INTENT(OUT) :: xmin
double precision :: brent

INTEGER, PARAMETER :: ITMAX=1000
double precision , PARAMETER :: CGOLD=0.3819660 ,ZEPS=1.0e-3*epsilon(ax)
INTEGER :: iter
double precision  :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm

	interface 
		function func (x)
		double precision, intent(in) :: x
		double precision :: func
		end function func
	end interface 

	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.0
	fx=func(x)
	fv=fx
	fw=fx
	do iter=1,ITMAX
		xm=0.5*(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.0*tol1
		if (abs(x-xm) <= (tol2-0.5*(b-a))) then
			xmin=x
			brent=fx
			RETURN
		end if
		if (abs(e) > tol1) then
			r=(x-w)*(fx-fv)
			q=(x-v)*(fx-fw)
			p=(x-v)*q-(x-w)*r
			q=2.0*(q-r)
			if (q > 0.0) p=-p
			q=abs(q)
			etemp=e
			e=d
			if (abs(p) >= abs(0.5*q*etemp) .or. &
				p <= q*(a-x) .or. p >= q*(b-x)) then
				e=merge(a-x,b-x, x >= xm )
				d=CGOLD*e
			else
				d=p/q
				u=x+d
				if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
			end if
		else
			e=merge(a-x,b-x, x >= xm )
			d=CGOLD*e
		end if
		u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
		fu=func(u)
		if (fu <= fx) then
			if (u >= x) then
				a=x
			else
				b=x
			end if
			call shft(v,w,x,u)
			call shft(fv,fw,fx,fu)
		else
			if (u < x) then
				a=u
			else
				b=u
			end if
			if (fu <= fw .or. w == x) then
				v=w
				fv=fw
				w=u
				fw=fu
			else if (fu <= fv .or. v == x .or. v == w) then
				v=u
				fv=fu
			end if
		end if
	end do
!	call nrerror('brent: exceed maximum iterations')
END FUNCTION brent

SUBROUTINE shft(a,b,c,d)
DOUBLE PRECISION , INTENT(OUT) :: a
DOUBLE PRECISION , INTENT(INOUT) :: b,c
DOUBLE PRECISION , INTENT(IN) :: d
	a=b
	b=c
	c=d
END SUBROUTINE shft

subroutine onedlin (x,y,n,x1,y1)
implicit none
integer, intent(in) :: n
double precision, intent(in) :: x(n), y(n), x1
double precision, intent(out) :: y1
integer :: point
    
    if ( x1 .lt. x(1) ) then
    point=1
    y1=y(point)+(x1-x(point))*(y(point+1)-y(point))/(x(point+1)-x(point))
    else if ( x1 .gt. x(n) ) then
    point=n-1
    y1=y(point)+(x1-x(point))*((y(point+1)-y(point))/(x(point+1)-x(point)))
    
    else
    
    point=1
	do while ( point .lt. n )

		if ( x1 .ge. x(point) .and. x1 .le. x(point+1) ) then
		y1=y(point)+(x1-x(point))*(y(point+1)-y(point))/(x(point+1)-x(point))
		point=point+n+100		
		end if

	point=point+1
	end do

    end if

end subroutine onedlin

subroutine basefun (grid_x,npx,x,vals,inds) 
implicit none

! this subroutine returns the values and the indices of the two basis
! functions that are positive on a given x in the grid_x

double precision,intent(in) :: x
integer , intent(in):: npx
double precision, intent(in) :: grid_x (npx)
double precision, intent(out) ::vals(2)
integer ,intent(out) ::inds(2)
integer :: i
logical :: elem
	
	elem=.false.
	i=1
	do while ( (.not. elem) .and. (i<=npx-1) )
		elem= ((grid_x(i) <= x) .and. (grid_x(i+1) >= x))
		i=i+1 
	end do


	select case(elem)

	case(.true.)
		vals(2)=( x-grid_x(i-1) )/(grid_x(i)-grid_x(i-1))
		vals(1)=( grid_x(i)-x )/(grid_x(i)-grid_x(i-1))
		inds(2)=i
		inds(1)=i-1

	case(.false.)

	if (x>grid_x(npx)) then
	   vals(2)=( x-grid_x(npx-1) )/(grid_x(npx)-grid_x(npx-1))
	   vals(1)=( grid_x(npx)-x )/(grid_x(npx)-grid_x(npx-1))
	   inds(2)=npx
	   inds(1)=npx-1

	endif

    if (x<grid_x(1)) then
	   vals(2)=( x-grid_x(1) )/(grid_x(2)-grid_x(1))
	   vals(1)=( grid_x(2)-x )/(grid_x(2)-grid_x(1))
	   inds(2)=2
	   inds(1)=1
	endif

	end select

    if ( vals(1) .lt. zero .or. vals(2) .lt. zero ) then
    print*,'out of bound in basefun mini.f90 '
        call exit
    end if
    
end subroutine basefun

SUBROUTINE tauchen(n,rho,mu_z,var_e,Z,Pi)
!z_{t+1}=(1-rho)*mu_z+rho*z_t+eps_t
!z_{t+1}-mu_z=rho*(z_t-mu_z)+eps_t
!y_{t+1}=rho*y_t+eps_t
!rho: persistence of AR(1) process
!mu_z: Mean of z
!var_e: Variance of epsilon
!n: number of finite states
INTEGER, INTENT(IN) :: n
DOUBLE PRECISION, INTENT(IN) :: rho,mu_z,var_e
DOUBLE PRECISION, INTENT(OUT) :: Z(n,1),Pi(n,n)
DOUBLE PRECISION :: Y(n,1),sig_z,w
INTEGER, PARAMETER :: m=3
INTEGER :: i,j
DOUBLE PRECISION :: x,mu,sd,pr
double precision :: dummy(n/2+1)

mu=0.
sd=SQRT(var_e)

sig_z=SQRT(var_e/(1.-(rho**2.)))

! 	DO i=1,n	
! 	Y(i,1)=-m*sig_z+2.*m*sig_z*((i-1.)/(n-1.))**1.
! 	Z(i,1)=Y(i,1)+mu_z	
! 	ENDDO

	DO i=1,n/2+1	
! 	Y(i,1)=-m*sig_z+(log(cut_off_data)+m*sig_z)*((i-1.)/(n/2+1-1.))**(.5)	
! 	Y(i,1)=-m*sig_z+(log(cut_off_data)+m*sig_z)*((i-1.)/(n/2+1-1.))**(1./1.5)	
	Y(i,1)=log(cut_off_data)+(-m*sig_z-log(cut_off_data))*((i-1.)/(n/2+1-1.))**(2.5)	
	ENDDO

	dummy=Y(1:n/2+1,1)
	
    call sort(dummy,n/2+1)
	
Y(1:n/2+1,1)=dummy

	DO i=n/2+1,n	
	Y(i,1)=log(cut_off_data)+(m*sig_z-log(cut_off_data))*((i-n/2-1.)/(n-n/2-1.))**(2.5)	
	ENDDO	

	do i=1,n
	Z(i,1)=Y(i,1)+mu_z	
	end do
	
	DO i=1,n

	w=Y(2,1)-Y(1,1)
	x=Y(1,1)+w/2.-rho*Y(i,1)
	CALL Ncdf(x,mu,sd,pr)
	Pi(i,1)=pr

        w=Y(n,1)-Y(n-1,1)	
	x=Y(n,1)-w/2.-rho*Y(i,1)
	CALL Ncdf(x,mu,sd,pr)
	Pi(i,n)=1.-pr

		DO j=2,n-1

		w=Y(j+1,1)-Y(j,1)
		x=Y(j,1)+w/2.-rho*Y(i,1)
		CALL Ncdf(x,mu,sd,pr)
		Pi(i,j)=pr

		w=Y(j,1)-Y(j-1,1)		
		x=Y(j,1)-w/2.-rho*Y(i,1)
		CALL Ncdf(x,mu,sd,pr)
		Pi(i,j)=Pi(i,j)-pr
		
		ENDDO
		
	ENDDO
	
	do i=1,n
	Pi(i,:)=Pi(i,:)/sum(Pi(i,:))
	end do

END SUBROUTINE tauchen

SUBROUTINE Ncdf(x,mu,sd,pr)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x,mu,sd
DOUBLE PRECISION, INTENT(OUT) :: pr
DOUBLE PRECISION :: z
DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462643383279502884197

z=(x-mu)/sd

	IF (z .GT. 0.) THEN
	pr=.5+.5*SQRT(1.-EXP(-SQRT(PI/8.)*(z**2.)))
	ELSE
	pr=.5-.5*SQRT(1.-EXP(-SQRT(PI/8.)*(z**2.)))
	ENDIF
END SUBROUTINE Ncdf

SUBROUTINE spline(x,y,yp1,ypn,y2)
IMPLICIT NONE
DOUBLE PRECISION , DIMENSION(:), INTENT(IN) :: x,y
DOUBLE PRECISION , INTENT(IN) :: yp1,ypn
DOUBLE PRECISION , DIMENSION(:), INTENT(OUT) :: y2
INTEGER :: n
DOUBLE PRECISION , DIMENSION(size(x)) :: a,b,c,r

	n=size(x)
	c(1:n-1)=x(2:n)-x(1:n-1)
	r(1:n-1)=6.0*((y(2:n)-y(1:n-1))/c(1:n-1))
	r(2:n-1)=r(2:n-1)-r(1:n-2)
	a(2:n-1)=c(1:n-2)
	b(2:n-1)=2.0*(c(2:n-1)+a(2:n-1))
	b(1)=1.0
	b(n)=1.0
	if (yp1 > 0.99e30) then
		r(1)=0.0
		c(1)=0.0
	else
		r(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
		c(1)=0.5
	end if
	if (ypn > 0.99e30) then
		r(n)=0.0
		a(n)=0.0
	else
		r(n)=(-3.0/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
		a(n)=0.5
	end if
	call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
END SUBROUTINE spline

SUBROUTINE tridag(a,b,c,r,u)
IMPLICIT NONE
DOUBLE PRECISION , DIMENSION(:), INTENT(IN) :: a,b,c,r
DOUBLE PRECISION , DIMENSION(:), INTENT(OUT) :: u
DOUBLE PRECISION , DIMENSION(size(b)) :: gam
INTEGER :: n,j
DOUBLE PRECISION  :: bet
	n=size(a)+1
	bet=b(1)
	if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
	u(1)=r(1)/bet
	do j=2,n
		gam(j)=c(j-1)/bet
		bet=b(j)-a(j-1)*gam(j)
		if (bet == 0.0) &
			call nrerror('tridag_ser: Error at code stage 2')
		u(j)=(r(j)-a(j-1)*u(j-1))/bet
	end do
	do j=n-1,1,-1
		u(j)=u(j)-gam(j+1)*u(j+1)
	end do
END SUBROUTINE tridag

FUNCTION splint(xa,ya,y2a,x)
IMPLICIT NONE
DOUBLE PRECISION , DIMENSION(:), INTENT(IN) :: xa,ya,y2a
DOUBLE PRECISION , INTENT(IN) :: x
DOUBLE PRECISION  :: splint
INTEGER :: khi,klo,n
DOUBLE PRECISION :: a,b,h
	n=size(xa)
	klo=max(min(locate(xa,x),n-1),1)
	khi=klo+1
	h=xa(khi)-xa(klo)
	if (h == 0.0) call nrerror('bad xa input in splint')
	a=(xa(khi)-x)/h
	b=(x-xa(klo))/h
	splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0
END FUNCTION splint

SUBROUTINE nrerror(string)
CHARACTER(LEN=*), INTENT(IN) :: string
	write (*,*) 'nrerror: ',string
	STOP 'program terminated by nrerror'
END SUBROUTINE nrerror


FUNCTION locate(xx,x)
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: xx
DOUBLE PRECISION, INTENT(IN) :: x
INTEGER :: locate
INTEGER :: n,jl,jm,ju
LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	jl=0
	ju=n+1
	do
		if (ju-jl <= 1) exit
		jm=(ju+jl)/2
		if (ascnd .eqv. (x >= xx(jm))) then
			jl=jm
		else
			ju=jm
		end if
	end do
	if (x == xx(1)) then
		locate=1
	else if (x == xx(n)) then
		locate=n-1
	else
		locate=jl
	end if
END FUNCTION locate

subroutine sort(x,n)
implicit none
integer, intent(in) :: n
double precision, dimension(n,2), intent(inout) :: x

integer :: i, location
double precision :: y(n)

y(:)=x(:,1)

	do i=1,n-1
	location=findmin(y,i,n)
		call swap(x(i,:),x(location,:))
	y(:)=x(:,1)
	end do

end subroutine sort

function findmin (x,beg,en)
implicit none
double precision, dimension(:), intent(in) :: x
integer, intent(in) :: beg, en

double precision :: minimum
integer :: location, i
integer :: findmin

minimum=x(beg)
location=beg

	do i=beg+1,en
		if ( x(i) .le. minimum ) then
		minimum=x(i)
		location=i
		end if
	end do

findmin=location

end function findmin

subroutine swap(a,b)
implicit none
double precision, intent(inout) :: a(2), b(2)
double precision :: temp(2)
temp=a
a=b
b=temp
end subroutine swap

end module mini 


