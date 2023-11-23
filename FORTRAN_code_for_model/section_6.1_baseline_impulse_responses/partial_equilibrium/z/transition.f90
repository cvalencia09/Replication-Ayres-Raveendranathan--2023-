subroutine transition
use params
use mini
implicit none

include 'mpif.h'

integer(kind=1), allocatable, dimension(:) :: g_i_o_mpi, g_i_constrained_mpi
integer(kind=2), allocatable, dimension(:) :: g_i_k_mpi
double precision, allocatable, dimension(:) :: g_i_a_mpi
double precision, allocatable, dimension(:) :: v_i_mpi, e_v_i_mpi

integer(kind=1), allocatable, dimension(:) :: g_i_o_recmpi, g_i_constrained_recmpi
integer(kind=2), allocatable, dimension(:) :: g_i_k_recmpi
double precision, allocatable, dimension(:) :: g_i_a_recmpi
double precision, allocatable, dimension(:) :: v_i_recmpi, e_v_i_recmpi

double precision :: vtemp, ev
integer :: find

integer :: t

integer :: iter_bi
double precision :: d_bi
double precision :: wmin, wmax
double precision :: mumax, mumin

integer(kind=1), allocatable, dimension(:,:,:,:) :: g_i_o_t
double precision, allocatable, dimension(:,:,:,:) :: g_i_a_t, g_q_t
integer(kind=2), allocatable, dimension(:,:,:,:) :: g_i_k_t
integer(kind=1), allocatable, dimension(:,:) :: g_e_o_t
double precision, dimension(neps,nk,na) :: omega_ss0, omega_ss1

double precision, parameter :: weight=.999, weight_w=.975

double precision :: vals(2)
integer :: inds(2)

double precision, parameter :: tol_tr=.005
double precision :: itmax_tr=0
double precision :: d_tr

double precision, allocatable, dimension(:,:,:,:) :: omega_ss_age0, omega_ss_age1 
double precision :: dummy

double precision, dimension(t1) :: r_t1, w_t1
double precision, dimension(neps) :: v_entry

!unconstrained firm's problem
integer(kind=1), allocatable, dimension(:) :: g_un_o_mpi
integer(kind=2), allocatable, dimension(:) :: g_un_k_mpi
double precision, allocatable, dimension(:) :: v_un_mpi, e_v_un_mpi

integer(kind=1), allocatable, dimension(:) :: g_un_o_recmpi
integer(kind=2), allocatable, dimension(:) :: g_un_k_recmpi
double precision, allocatable, dimension(:) :: v_un_recmpi, e_v_un_recmpi
integer(kind=1), allocatable, dimension(:,:,:) :: g_un_o_t
integer(kind=2), allocatable, dimension(:,:,:) :: g_un_k_t

double precision, dimension(neps) :: v_entry_un, ap_entry_un
integer(kind=2), dimension(neps) :: kp_entry_un
integer(kind=1), dimension(neps) :: d_o_entry_un
double precision :: vtemp0

integer :: jk_temp

double precision :: tot_adj_cost

    if ( rank .eq. 0 ) then
    print*,' '
    print*,'transition '
    end if
    
    if ( rank .eq. 0 ) then
    
        allocate (g_i_o_t(t1,neps,nk,na),g_i_k_t(t1,neps,nk,na),g_i_a_t(t1,neps,nk,na),g_e_o_t(t1,neps),g_q_t(t1,neps,nk,na),g_un_k_t(t1,neps,nk),&
        g_un_o_t(t1,neps,nk))
        
    g_i_o_t(1,:,:,:)=g_i_o_initial
    g_i_k_t(1,:,:,:)=g_i_k_initial
    g_i_a_t(1,:,:,:)=g_i_a_initial
    g_e_o_t(1,:)=g_e_o_initial
    g_q_t(1,:,:,:)=g_q_initial
    g_un_k_t(1,:,:)=g_un_k_initial
    g_un_o_t(1,:,:)=g_un_o_initial
    
    g_i_o_t(t1,:,:,:)=g_i_o_final
    g_i_k_t(t1,:,:,:)=g_i_k_final
    g_i_a_t(t1,:,:,:)=g_i_a_final
    g_e_o_t(t1,:)=g_e_o_final
    g_q_t(t1,:,:,:)=g_q_final
    g_un_k_t(t1,:,:)=g_un_k_final
    g_un_o_t(t1,:,:)=g_un_o_final
    
        do ieps=1,neps
        do ik=1,nk
        do ia=1,na
        
            if ( g_i_o_t(1,ieps,ik,ia) .ne. g_i_o_t(t1,ieps,ik,ia) ) then
            print*,'g_i_o_t(1,ieps,ik,ia) .ne. g_i_o_t(t1,ieps,ik,ia)'
                call exit                
            end if
            
        end do
        end do
        end do
        
        do ieps=1,neps
        do ik=1,nk
        
            if ( abs(a_un_t(1,ieps,ik)-a_un_t(t1,ieps,ik)) .ge. 1e-4 ) then
            print*,'abs(a_un_t(1,ieps,ik)-a_un_t(t1,ieps,ik)) .ge. 1e-4'
            print*,a_un_t(1,ieps,ik),a_un_t(t1,ieps,ik)
                call exit
            end if

        end do
        end do

        do ieps=1,neps
        do ik=1,nk
        
            if ( g_un_k_t(1,ieps,ik) .ne. g_un_k_t(t1,ieps,ik) ) then
            print*,'g_un_k_t(1,ieps,ik) .ne. g_un_k_t(t1,ieps,ik)'
            print*,g_un_k_t(1,ieps,ik),g_un_k_t(t1,ieps,ik)
                call exit
            end if
            
        end do
        end do

        do ieps=1,neps
        do ik=1,nk
        
            if ( g_un_o_t(1,ieps,ik) .ne. g_un_o_t(t1,ieps,ik) ) then
            print*,'g_un_o_t(1,ieps,ik) .ne. g_un_o_t(t1,ieps,ik)'
                call exit
            end if
            
        end do
        end do
        
        do ieps=1,neps
        do ik=1,nk
        
            if ( abs(a_tilde_t(1,ieps,ik)-a_tilde_t(t1,ieps,ik)) .ge. 1e-4 ) then
            print*,'abs(a_tilde_t(1,ieps,ik)-a_tilde_t(t1,ieps,ik)) .ge. 1e-4'
            print*,a_tilde_t(1,ieps,ik),a_tilde_t(t1,ieps,ik)
                call exit
            end if
            
        end do
        end do

        do ieps=1,neps
        do ik=1,nk
        do ia=1,na
        
            if ( abs(g_i_a_t(1,ieps,ik,ia)-g_i_a_t(t1,ieps,ik,ia)) .ge. 1e-4 ) then
            print*,'abs(g_i_a_t(1,ieps,ik,ia)-g_i_a_t(t1,ieps,ik,ia)) .ge. 1e-4'
            end if
            
        end do
        end do
        end do
        
        do ieps=1,neps
        do ik=1,nk
        do ia=1,na
        
            if ( g_i_k_t(1,ieps,ik,ia) .ne. g_i_k_t(t1,ieps,ik,ia) ) then
            print*,'g_i_k_t(1,ieps,ik,ia) .ne. g_i_k_t(t1,ieps,ik,ia)'
                call exit                
            end if
            
        end do
        end do
        end do
 
        do ieps=1,neps
        
            if ( g_e_o_t(1,ieps) .ne. g_e_o_t(t1,ieps) ) then
            print*,'g_e_o_t(1,ieps) .ne. g_e_o_t(t1,ieps)'
                call exit                
            end if
            
        end do
        
        do ieps=1,neps
        do ik=1,nk
        do ia=1,na
        
            if ( abs(g_q_t(1,ieps,ik,ia)-g_q_t(t1,ieps,ik,ia)) .ge. 1e-4 ) then
            print*,'abs(g_q_t(1,ieps,ik,ia)-g_q_t(t1,ieps,ik,ia)) .ge. 1e-4'
            end if
            
        end do
        end do
        end do
 
    end if
    
!initial conditions
! b_t(2)=b_t(1)
r_t(2)=r_t(1)
firm_debt_t(2)=firm_debt_t(1)

r_t1(1)=r_t(1)
r_t1(t1)=r_t(t1)
r_t1(2)=r_t1(1)
w_t1(1)=w_t(1)
w_t1(t1)=w_t(t1)

!guess
r_t(3:t1-1)=r_t(1)

!          open (unit=1,file='input_r_t.txt')
!              do t=3,t1-1
!              read(1,*)r_t(t)
!              end do
!          close(1)

    if ( ex_entry .eq. 1 ) then

    w_t(2:t1-1)=w_t(1)
    
!          open (unit=1,file='input_w_t.txt')
!              do t=2,t1-1
!              read(1,*)w_t(t)
!             end do
!          close(1)

    end if
    
d_tr=1
iter_tr=0

    do while ( iter_tr .le. itmax_tr .and. d_tr .ge. tol_tr )
 
        !solve for sequence of c_t
        do t=t1-1,2,-1
        c_t(t)=c_t(t+1)/(beta*(1+r_t(t+1)))
        end do
        
        !solve for sequence of w_t
        allocate (g_i_o_recmpi(np*nproc),g_i_k_recmpi(np*nproc),g_i_a_recmpi(np*nproc),v_i_recmpi(np*nproc),e_v_i_recmpi(np*nproc))        
        allocate (g_un_o_recmpi(np_un*nproc),g_un_k_recmpi(np_un*nproc),v_un_recmpi(np_un*nproc),e_v_un_recmpi(np_un*nproc))

        allocate (g_i_constrained_recmpi(np*nproc))
        
        !value function estimates from fistst: actual firm        
        do cnt=1,ni     
        
        ieps=zmap(cnt,1)
        ik=zmap(cnt,2)
        ia=zmap(cnt,3)                
        
        v_i_recmpi(cnt)=v_i(ieps,ik,ia)                
        g_i_o_recmpi(cnt)=g_i_o(ieps,ik,ia)
        
            if ( rank .eq. 0 .and. cnt .ne. zmap_inv(ieps,ik,ia) ) then
            print*,'rank .eq. 0 .and. cnt .ne. zmap_inv(ieps,ik,ia), transition.f90'
                call exit
            end if
            
        end do

        do cnt=1,ni_un     
        
        ieps=zmap_un(cnt,1)
        ik=zmap_un(cnt,2)
        
        v_un_recmpi(cnt)=v_un(ieps,ik)                
        
            if ( rank .eq. 0 .and. cnt .ne. zmap_inv_un(ieps,ik) ) then
            print*,'rank .eq. 0 .and. cnt .ne. zmap_inv_un(ieps,ik), transition.f90'
                call exit
            end if
            
        end do
            
        !backward induction to solve for sequence of w
        if ( rank .eq. 0 ) then
        print*,' '
        print*,'backward induction'
        end if
        
        do t=t1-1,2,-1                                            
            
            !bond price schedule
            do ieps=1,neps
            do jk=1,nk
            do ja=1,na
    
                if ( ja .lt. intbirth ) then
                
                    if ( default_economy .eq. 0 ) then
                    g_q(ieps,jk,ja)=one/(one+r_t(min(t+1,t1-1)))-tau_a_t(t)
                    else if ( default_economy .eq. 1 ) then
                    
                    g_q(ieps,jk,ja)=zero                
                
                        do jeps=1,neps
                    
                            if ( pi_eps(ieps,jeps) .gt. zero ) then
                            d_o=g_i_o_recmpi(zmap_inv(jeps,jk,ja))                            
                            g_q(ieps,jk,ja)=g_q(ieps,jk,ja)+pi_eps(ieps,jeps)*&
                            (d_o*(-grid_a(t+1,ja))+&
                            (1-d_o)*min(-grid_a(t+1,ja),max(zero,(1.+tau_x_t(t+1))*grid_k(jk)*chi-func_adj(zero,grid_k(jk),d_o,t+1))))
                            end if
                        
                        end do
                
                    g_q(ieps,jk,ja)=g_q(ieps,jk,ja)/((one+r_t(min(t+1,t1-1)))*(-grid_a(t+1,ja)))-tau_a_t(t)
                
                    end if
                    
                else if ( ja .ge. intbirth ) then
                g_q(ieps,jk,ja)=one/(one+r_t(min(t+1,t1-1)))
                end if
            
            end do
            end do
            end do
        
            ! actual firm's problem: expectations
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

            !unconstrained firm's problem: expectations
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
            
            if ( ex_entry .eq. 0 ) then
            
            !entrant's problem                    
                open (unit=1,file='input_instst_w.txt')
                read(1,*)wmin,wmax
                close(1)
    
            wmin=wmin*.9
            wmax=wmax*1.1
        
            iter_bi=0
            d_bi=1            
            
                do while ( iter_bi .le. itmax_bi .and. abs(d_bi) .ge. tol_bi )
    
                w_t(t)=(wmin+wmax)/2
                            
                    !entrant: problem
                    do ieps=1,neps
        
                    ik=intbirth_k
                    ia=intbirth
                    iz=id_z(ieps)
                    
                    eps=grid_eps_t(t,ieps)
                    k=grid_k(ik)
                    a=grid_a(t,ia)

                    !entrant: unconstrained problem                              
                    find=0
                
                    !exit
                    d_o=0
                    jk=1
                    kp=grid_k(jk)
                    
                    vtemp=(1.+tau_x_t(t))*k*chi-func_adj(kp,k,d_o,t)
        
                        if ( find .eq. 0 ) then
                        find=1
                        d_o_entry_un(ieps)=d_o
                        v_entry_un(ieps)=vtemp
                        kp_entry_un(ieps)=jk                    
                        else if ( find .eq. 1 ) then
                        print*,'find .eq. 1, entrant: unconstrained problem, transition.f90'
                            call exit                            
                        end if
                                    
                    !operate
                    d_o=1    
                    l=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
        
                        if ( isnan(l) ) then
                        print*,'isnan l, entrant: unconstrained problem, transition.f90 ',k
                            call exit
                        end if

                    kpmin=kmin
                    kpmax=kmax
                      
                    kp_temp=-one
                    jk_temp=-1
        
                        do jk=1,nk

                        kp=grid_k(jk)
            
                            if ( kp .ge. kpmin .and. kp .le. kpmax ) then
                            
                            profit=((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w_t(t)*(l+f_o_t(t)*grid_z(iz))&
                            +(1.+tau_x_t(t))*((1-delta_k)*k-kp)-func_adj(kp,k,d_o,t)

                            ev=e_v_un_recmpi(zmap_inv_un(ieps,jk)) 
                            vtemp0=profit+(1-delta_f)*beta*ev*c_t(t)/c_t(t+1)
            
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
                        profit=((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w_t(t)*(l+f_o_t(t)*grid_z(iz))&
                        +(1.+tau_x_t(t))*((1-delta_k)*k-kp)-func_adj(kp,k,d_o,t)

                        ev=e_v_un_recmpi(zmap_inv_un(ieps,jk)) 
                        vtemp=profit+(1-delta_f)*beta*ev*c_t(t)/c_t(t+1)
                
                            if ( find .eq. 0 ) then

                            print*,'find .eq. 0, entrant: unconstrained problem, transition.f90'
                                call exit
                            
                            else if ( find .eq. 1 ) then
                    
                                if ( vtemp .gt. v_entry_un(ieps) ) then
                                d_o_entry_un(ieps)=d_o
                                v_entry_un(ieps)=vtemp
                                kp_entry_un(ieps)=jk
                                end if

                            end if

                        else if ( kp_temp .lt. zero ) then
                    
                        print*,'kp_temp .lt. zero, entrant: unconstrained problem, transition.f90'
                            call exit
                        
                        end if
                
                    !minimum savings to be unconstrained: entrant                
                    find=0
                    jk=kp_entry_un(ieps)                    
    
                        do jeps=1,neps
        
                            if ( pi_eps(ieps,jeps) .gt. zero ) then
            
                                if ( find .eq. 0 ) then
                    
                                find=1
                                ap_entry_un(ieps)=a_tilde_t(t+1,jeps,jk)
                    
                                else if ( find .eq. 1 ) then
                
                                dummy=a_tilde_t(t+1,jeps,jk)

                                    if ( dummy .ge. ap_entry_un(ieps) ) then
                                    ap_entry_un(ieps)=dummy
                                    end if
                    
                                end if
        
                            end if
            
                        end do
        
                    ap_entry_un(ieps)=max(ap_entry_un(ieps),zero)        
                
                    !entrant problem: actual
                    find=0
                
                    !behave like unconstrained firm
                    d_o=d_o_entry_un(ieps)
                    jk=kp_entry_un(ieps)
                    kp=grid_k(jk)
                    
                        if ( d_o .eq. 0 ) then
                    
                        ja=intbirth   
                        ap=grid_a(t+1,ja)
            
                            if ( kp .ne. kmin ) then
                            print*,'kp .ne. kmin, entrant problem: actual, transition.f90'
                                call exit
                            end if

                            if ( default_economy .eq. 0 ) then
                            v_entry(ieps)=(1.+tau_x_t(t))*k*chi+a-func_adj(kp,k,d_o,t)
                            else if ( default_economy .eq. 1 ) then
                            v_entry(ieps)=max((1.+tau_x_t(t))*k*chi+a-func_adj(kp,k,d_o,t),min(zero,(1.+tau_x_t(t))*k*chi-func_adj(kp,k,d_o,t)))
                            end if
                            
                        else if ( d_o .eq. 1 ) then
                
                        l=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
        
                            if ( isnan(l) ) then
                            print*,'isnan l, entrant problem: actual, transition.f90 ',l
                                call exit
                            end if
            
                        ap=ap_entry_un(ieps)

                            if ( ap .ge. zero ) then
                                
                            q=one/(one+r_t(min(t+1,t1-1)))
                                
                            else if ( ap .lt. zero ) then
                                
                                if ( default_economy .eq. 0 ) then                                
                                q=one/(one+r_t(min(t+1,t1-1)))-tau_a_t(t)
                                else if ( default_economy .eq. 1 ) then
                                    call onedlin (grid_a(t+1,:),g_q(ieps,jk,:),na,ap,q) 
                                end if

                            end if
                                        
                        profit=((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w_t(t)*(l+f_o_t(t)*grid_z(iz))&
                        +(1.+tau_x_t(t))*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,t)
                        
                            if ( profit .ge. zero ) then
!                             if ( profit .ge. zero .or. profit .lt. zero ) then
                        
                            exp_ap=e_v_i_recmpi(zmap_inv(ieps,jk,:))           
                        
                                call onedlin (grid_a(t+1,:),exp_ap,na,ap,ev)
                                                        
                            v_entry(ieps)=profit+(1-delta_f)*beta*ev*c_t(t)/c_t(t+1)
                        
                            else if ( profit .lt. zero ) then
                        
                            !behave like constrained firm
                            !exit
                            d_o=0
                            jk=1
                            ja=intbirth                
                            kp=grid_k(jk)
                            ap=grid_a(t+1,ja)
                            
                                if ( default_economy .eq. 0 ) then                            
                                vtemp=(1.+tau_x_t(t))*k*chi+a-func_adj(kp,k,d_o,t)
                                else if ( default_economy .eq. 1 ) then                            
                                vtemp=max((1.+tau_x_t(t))*k*chi+a-func_adj(kp,k,d_o,t),min(zero,(1.+tau_x_t(t))*k*chi-func_adj(kp,k,d_o,t)))
                                end if
                    
                                if ( find .eq. 0 ) then
            
                                find=1
                                v_entry(ieps)=vtemp
        
                                else if ( find .eq. 1 ) then

                                print*,'find .eq. 0, entrant problem: actual, transition.f90'
                                    call exit
                
                                end if
            
                            !operate
                            d_o=1    
                            l=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
                
                                if ( isnan(l) ) then
                                print*,'isnan l, entrant problem: actual, transition.f90'
                                    call exit
                                end if

                            kpmin=kmin
                            kpmax=kmax
                            
                                !operate (constrained)                
                                do jk=1,nk
           
                                kp=grid_k(jk)

                                constant=((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w_t(t)*(l+f_o_t(t)*grid_z(iz))+&
                                a+(1.+tau_x_t(t))*((1-delta_k)*k-kp)-func_adj(kp,k,d_o,t)

                                    if ( constant .ge. zero ) then

                                    q=one/(one+r_t(min(t+1,t1-1)))
                                    ap=constant/q
                        
                                    else if ( constant .lt. zero ) then
                            
                                        if ( default_economy .eq. 0 ) then                            
                                        q=one/(one+r_t(min(t+1,t1-1)))-tau_a_t(t)
                                        ap=constant/q
                                        else if ( default_economy .eq. 1 ) then
                                        print*,'fix default, transition.f90'
                                            call exit                        
                                        end if
                            
                                    end if

                                ap=min(grid_a(t+1,na),ap)                                
                                
                                    if ( ap .ge. max(abar_t(t)-theta_t(t)*k,grid_a(t+1,1)) .and. kp .ge. kpmin .and. kp .le. kpmax ) then
                            
                                        if ( ap .ge. zero ) then
                                        q=one/(one+r_t(min(t+1,t1-1)))
                                        else if ( ap .lt. zero ) then        
                                        
                                            if ( default_economy .eq. 0 ) then
                                            q=one/(one+r_t(min(t+1,t1-1)))-tau_a_t(t)
                                            else if ( default_economy .eq. 1 ) then                                            
                                                call onedlin (grid_a(t+1,:),g_q(ieps,jk,:),na,ap,q) 
                                            end if
                                            
                                        end if
                                
                                    profit=((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w_t(t)*(l+f_o_t(t)*grid_z(iz))&
                                    +(1.+tau_x_t(t))*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,t)
                                        
                                    exp_ap=e_v_i_recmpi(zmap_inv(ieps,jk,:)) 
                                    
                                            call onedlin (grid_a(t+1,:),exp_ap,na,ap,ev)
                                            
                                    vtemp=profit+(1-delta_f)*beta*ev*c_t(t)/c_t(t+1)
                        
                                        if ( find .eq. 0 ) then

                                        print*,'find .eq. 0, entrant problem: actual, transition.f90'                                
                                            call exit
                        
                                        else if ( find .eq. 1 ) then
                        
                                            if ( vtemp .gt. v_entry(ieps) ) then
                                            v_entry(ieps)=vtemp
                                            end if
                        
                                        end if
        
                                    end if    

                                end do                        
                                                
                            end if
                        
                        end if
                    
                    end do

                v_e=zero
                
                    do ieps=1,neps
                    v_e=v_e+put_eps(ieps)*v_entry(ieps)
                    end do
  
                    if ( v_e-w_t(t)*f_e_t(t)-f_c-k_e .ge. zero ) then
                    wmin=w_t(t)
                    else if ( v_e-w_t(t)*f_e_t(t)-f_c-k_e .lt. zero ) then
                    wmax=w_t(t)
                    end if
            
                d_bi=max(abs(v_e-w_t(t)*f_e_t(t)-f_c-k_e),abs(wmax-wmin))
                iter_bi=iter_bi+1
                
                end do
            
                if ( rank .eq. 0 .and. abs(d_bi) .ge. tol_bi ) then
                print*,'iterations exceeded in bisection for w_t, transition.f90 ',t,d_bi            
!                    call exit
                end if       
            
            end if
            
            !incumbent problem: unconstrained
            allocate (g_un_o_mpi(np_un),g_un_k_mpi(np_un),v_un_mpi(np_un))
        
            do cnt=itop_un,iend_un
    
            ieps=zmap_un(cnt,1)
            ik=zmap_un(cnt,2)
            iz=id_z(ieps)
            
            eps=grid_eps_t(t,ieps)
            k=grid_k(ik)

            find=0
        
            !exit
            d_o=0
            jk=1
            kp=grid_k(jk)
            vtemp=(1.+tau_x_t(t))*k*chi-func_adj(kp,k,d_o,t)
        
                if ( find .eq. 0 ) then
                find=1
                g_un_o_mpi(cnt-itop_un+1)=d_o
                g_un_k_mpi(cnt-itop_un+1)=jk
                v_un_mpi(cnt-itop_un+1)=vtemp                        
                else if ( find .eq. 1 ) then
                print*,'find .eq. 1, incumbent problem: unconstrained, transition.f90'
                    call exit                            
                end if        
                
            !operate
            d_o=1    
            l=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
        
                if ( isnan(l) ) then
                print*,'isnan l, incumbent problem: unconstrained, transition.f90 ',k
                    call exit
                end if

            kpmin=kmin
            kpmax=kmax
                                      
            kp_temp=-one
            jk_temp=-1
        
                do jk=1,nk

                kp=grid_k(jk)
            
                    if ( kp .ge. kpmin .and. kp .le. kpmax ) then
                            
                    profit=((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w_t(t)*(l+f_o_t(t)*grid_z(iz))&
                    +(1.+tau_x_t(t))*((1-delta_k)*k-kp)-func_adj(kp,k,d_o,t)
                                        
                    ev=e_v_un_recmpi(zmap_inv_un(ieps,jk)) 
                    vtemp0=profit+(1-delta_f)*beta*ev*c_t(t)/c_t(t+1)
            
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
                profit=((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w_t(t)*(l+f_o_t(t)*grid_z(iz))&
                +(1.+tau_x_t(t))*((1-delta_k)*k-kp)-func_adj(kp,k,d_o,t)
            
                ev=e_v_un_recmpi(zmap_inv_un(ieps,jk))  
                vtemp=profit+(1-delta_f)*beta*ev*c_t(t)/c_t(t+1)
                
                    if ( find .eq. 0 ) then

                    print*,'find .eq. 0, incumbent problem: unconstrained, transition.f90'
                        call exit
                        
                    else if ( find .eq. 1 ) then
                    
                        if ( vtemp .gt. v_un_mpi(cnt-itop_un+1) ) then
                        g_un_o_mpi(cnt-itop_un+1)=d_o
                        g_un_k_mpi(cnt-itop_un+1)=jk
                        v_un_mpi(cnt-itop_un+1)=vtemp
                        end if
                        
                    end if

                else if ( kp_temp .lt. zero ) then
                
                print*,'kp_temp .lt. zero, incumbent problem: unconstrained, transition.f90'
                    call exit
                    
                end if
                      
            end do            
 
            call mpi_allgather (g_un_o_mpi,np_un,mpi_integer1,g_un_o_recmpi,np_un,mpi_integer1,mpi_comm_world,ierr)    
            call mpi_allgather (g_un_k_mpi,np_un,mpi_integer2,g_un_k_recmpi,np_un,mpi_integer2,mpi_comm_world,ierr)    
            call mpi_allgather (v_un_mpi,np_un,mpi_double_precision,v_un_recmpi,np_un,mpi_double_precision,mpi_comm_world,ierr)

            deallocate (g_un_o_mpi,g_un_k_mpi,v_un_mpi)

            !minimum savings to be unconstrained next period
            do ieps=1,neps
            do ik=1,nk

            eps=grid_eps_t(t,ieps)
            k=grid_k(ik)
            jk=g_un_k_recmpi(zmap_inv_un(ieps,ik))
            d_o=g_un_o_recmpi(zmap_inv_un(ieps,ik))
        
            kp=grid_k(jk)
            find=0
    
                do jeps=1,neps
        
                    if ( pi_eps(ieps,jeps) .gt. zero ) then
            
                        if ( find .eq. 0 ) then
                    
                        find=1
                        a_un_t(t,ieps,ik)=a_tilde_t(t+1,jeps,jk)
                    
                        else if ( find .eq. 1 ) then
                
                        dummy=a_tilde_t(t+1,jeps,jk)

                            if ( dummy .ge. a_un_t(t,ieps,ik) ) then
                            a_un_t(t,ieps,ik)=dummy
                            end if
                    
                        end if
                
                    end if
            
                end do
        
            a_un_t(t,ieps,ik)=max(a_un_t(t,ieps,ik),zero)
            a_un_t(t,ieps,ik)=min(a_un_t(t,ieps,ik),amax)
        
            end do
            end do
                        
            !incumbent problem: actual
            allocate (g_i_o_mpi(np),g_i_k_mpi(np),g_i_a_mpi(np),v_i_mpi(np))
            allocate (g_i_constrained_mpi(np))

            do cnt=itop,iend
        
            ieps=zmap(cnt,1)
            ik=zmap(cnt,2)
            ia=zmap(cnt,3)
            iz=id_z(ieps)
            
            eps=grid_eps_t(t,ieps)
            k=grid_k(ik)
            a=grid_a(t,ia)

            find=0    

            !unconstrained firm's policy function
            jk=g_un_k_recmpi(zmap_inv_un(ieps,ik))
            d_o=g_un_o_recmpi(zmap_inv_un(ieps,ik))

            kp=grid_k(jk)
            
                if ( d_o .eq. 0 ) then
            
                ja=intbirth    
                ap=grid_a(t+1,ja) 
            
                    if ( kp .ne. kmin ) then
                    print*,'kp .ne. kmin, incumbent problem: actual, transition.f90'
                        call exit
                    end if

                    if ( default_economy .eq. 0 ) then
                    vtemp=(1.+tau_x_t(t))*k*chi+a-func_adj(kp,k,d_o,t)
                    else if ( default_economy .eq. 1 ) then
                    vtemp=max((1.+tau_x_t(t))*k*chi+a-func_adj(kp,k,d_o,t),min(zero,(1.+tau_x_t(t))*k*chi-func_adj(kp,k,d_o,t)))
                    end if
                    
                    if ( find .eq. 0 ) then
            
                    find=1
                    g_i_o_mpi(cnt-itop+1)=d_o
                    g_i_k_mpi(cnt-itop+1)=jk
                    g_i_a_mpi(cnt-itop+1)=ap
                    v_i_mpi(cnt-itop+1)=vtemp
                    g_i_constrained_mpi(cnt-itop+1)=0
                    
                    else if ( find .eq. 1 ) then

                    print*,'find .eq. 0, incumbent problem: actual, transition.f90'
                        call exit
            
                    end if            
            
                else if ( d_o .eq. 1 ) then

                l=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))

                    if ( isnan(l) ) then
                    print*,'isnan l, incumbent problem: actual, transition.f90 ',l
                        call exit
                    end if
                                                
                ap=a_un_t(t,ieps,ik)

                    if ( ap .ge. zero ) then
                    q=one/(one+r_t(min(t+1,t1-1)))
                    else if ( ap .lt. zero ) then                            
                    
                        if ( default_economy .eq. 0 ) then
                        q=one/(one+r_t(min(t+1,t1-1)))-tau_a_t(t)                        
                        else if ( default_economy .eq. 1 ) then
                            call onedlin (grid_a(t+1,:),g_q(ieps,jk,:),na,ap,q)
                        end if
                        
                    end if       

                profit=((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w_t(t)*(l+f_o_t(t)*grid_z(iz))&
                +(1.+tau_x_t(t))*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,t)
                    
                    if ( profit .ge. zero ) then
!                     if ( profit .ge. zero .or. profit .lt. zero ) then
                    
                    exp_ap=e_v_i_recmpi(zmap_inv(ieps,jk,:))      
                    
                            call onedlin (grid_a(t+1,:),exp_ap,na,ap,ev)
                    
                    vtemp=profit+(1-delta_f)*beta*ev*c_t(t)/c_t(t+1)
                    
                        if ( find .eq. 0 ) then

                        find=1
                        g_i_o_mpi(cnt-itop+1)=d_o
                        g_i_k_mpi(cnt-itop+1)=jk
                        g_i_a_mpi(cnt-itop+1)=ap
                        v_i_mpi(cnt-itop+1)=vtemp
                        g_i_constrained_mpi(cnt-itop+1)=0
                    
                        else if ( find .eq. 1 ) then
                        
                        print*,'find .eq. 1, incumbent problem: actual, transition.f90'                                
                            call exit
                        
                        end if
                                                    
                    else if ( profit .lt. zero ) then

                    g_i_constrained_mpi(cnt-itop+1)=1
                    
                    !exit
                    d_o=0
                    jk=1
                    ja=intbirth   
                    kp=grid_k(jk)
                    ap=grid_a(t+1,ja) 

                        if ( default_economy .eq. 0 ) then
                        vtemp=(1.+tau_x_t(t))*k*chi+a-func_adj(kp,k,d_o,t)
                        else if ( default_economy .eq. 1 ) then
                        vtemp=max((1.+tau_x_t(t))*k*chi+a-func_adj(kp,k,d_o,t),min(zero,(1.+tau_x_t(t))*k*chi-func_adj(kp,k,d_o,t)))
                        end if
                
                        if ( find .eq. 0 ) then
            
                        find=1
                        g_i_o_mpi(cnt-itop+1)=d_o
                        g_i_k_mpi(cnt-itop+1)=jk
                        g_i_a_mpi(cnt-itop+1)=ap
                        v_i_mpi(cnt-itop+1)=vtemp
        
                        else if ( find .eq. 1 ) then

                        print*,'find .eq. 0, incumbent problem: actual, transition.f90'
                            call exit
            
                        end if

                    !operate
                    d_o=1    
                    l=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
                
                        if ( isnan(l) ) then
                        print*,'isnan l, incumbent problem: actual, transition.f90'
                            call exit
                        end if
                    
                    kpmin=kmin
                    kpmax=kmax
                
                        !operate (constrained)                
                        do jk=1,nk
                    
                        kp=grid_k(jk)
                        
                        constant=((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w_t(t)*(l+f_o_t(t)*grid_z(iz))+&
                        a+(1.+tau_x_t(t))*((1-delta_k)*k-kp)-func_adj(kp,k,d_o,t)

                            if ( constant .ge. zero ) then
                        
                            q=one/(one+r_t(min(t+1,t1-1)))
                            ap=constant/q
                        
                            else if ( constant .lt. zero ) then
                            
                                if ( default_economy .eq. 0 ) then                            
                                q=one/(one+r_t(min(t+1,t1-1)))-tau_a_t(t)
                                ap=constant/q
                                else if ( default_economy .eq. 1 ) then
                                print*,'fix default, transition.f90'
                                    call exit                        
                                end if
                        
                            end if
                        
                        ap=min(grid_a(t,na),ap)
                        
                            if ( ap .ge. max(abar_t(t)-theta_t(t)*k,grid_a(t+1,1)) .and. kp .ge. kpmin .and. kp .le. kpmax ) then

                                if ( ap .ge. zero ) then
                                q=one/(one+r_t(min(t+1,t1-1)))
                                else if ( ap .lt. zero ) then                            
                                
                                    if ( default_economy .eq. 0 ) then
                                    q=one/(one+r_t(min(t+1,t1-1)))-tau_a_t(t)
                                    else if ( default_economy .eq. 1 ) then                                                                        
                                        call onedlin (grid_a(t+1,:),g_q(ieps,jk,:),na,ap,q)
                                    end if
                                    
                                end if       
                            
                            profit=((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w_t(t)*(l+f_o_t(t)*grid_z(iz))&
                            +(1.+tau_x_t(t))*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,t)
                                
                            exp_ap=e_v_i_recmpi(zmap_inv(ieps,jk,:))      
                            
                                    call onedlin (grid_a(t+1,:),exp_ap,na,ap,ev)
                    
                            vtemp=profit+(1-delta_f)*beta*ev*c_t(t)/c_t(t+1)
                    
                                if ( find .eq. 0 ) then

                                print*,'find .eq. 0, incumbent problem: actual, transition.f90'                                
                                    call exit
                                    
                                else if ( find .eq. 1 ) then
                        
                                    if ( vtemp .gt. v_i_mpi(cnt-itop+1) ) then
                                    g_i_o_mpi(cnt-itop+1)=d_o
                                    g_i_k_mpi(cnt-itop+1)=jk
                                    g_i_a_mpi(cnt-itop+1)=ap
                                    v_i_mpi(cnt-itop+1)=vtemp
                                    end if
                        
                                end if
    
                            end if
    
                        end do

                    end if
                    
                end if                
                
            end do

            call mpi_allgather (g_i_o_mpi,np,mpi_integer1,g_i_o_recmpi,np,mpi_integer1,mpi_comm_world,ierr)    
            call mpi_allgather (g_i_k_mpi,np,mpi_integer2,g_i_k_recmpi,np,mpi_integer2,mpi_comm_world,ierr)        
            call mpi_allgather (g_i_a_mpi,np,mpi_double_precision,g_i_a_recmpi,np,mpi_double_precision,mpi_comm_world,ierr)    
            call mpi_allgather (v_i_mpi,np,mpi_double_precision,v_i_recmpi,np,mpi_double_precision,mpi_comm_world,ierr)
            call mpi_allgather (g_i_constrained_mpi,np,mpi_integer1,g_i_constrained_recmpi,np,mpi_integer1,mpi_comm_world,ierr)    
        
            deallocate (g_i_o_mpi,g_i_k_mpi,g_i_a_mpi,v_i_mpi)
            deallocate (g_i_constrained_mpi)
                
            !minimum savings to be unconstrained in each state in current period
            do ieps=1,neps
            do ik=1,nk
    
            iz=id_z(ieps)
            eps=grid_eps_t(t,ieps)
            k=grid_k(ik)
            d_o=g_un_o_recmpi(zmap_inv_un(ieps,ik))
            jk=g_un_k_recmpi(zmap_inv_un(ieps,ik))
            kp=grid_k(jk)
            
                if ( d_o .eq. 0 ) then

                a_tilde_t(t,ieps,ik)=zero
            
                else if ( d_o .eq. 1 ) then
            
                l=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))                            
                
                a_tilde_t(t,ieps,ik)=a_un_t(t,ieps,ik)/(1+r_t(min(t+1,t1-1)))-&
                ((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)+w_t(t)*(l+f_o_t(t)*grid_z(iz))-&
                (1.+tau_x_t(t))*((1-delta_k)*k-kp)+func_adj(kp,k,d_o,t)
            
                end if

            end do
            end do
            
            if ( rank .eq. 0 ) then

                if ( ex_entry .eq. 0 ) then
                
                g_e_o_t(t,:)=1
            
                else if ( ex_entry .eq. 1 ) then
            
                    do ieps=1,neps

                        if ( v_i_recmpi(zmap_inv(ieps,intbirth_k,intbirth)) .ge. w_t(t)*f_e_t(t)+f_c+k_e ) then
                        g_e_o_t(t,ieps)=1
                        else if ( v_i_recmpi(zmap_inv(ieps,intbirth_k,intbirth)) .lt. w_t(t)*f_e_t(t)+f_c+k_e ) then
                        g_e_o_t(t,ieps)=0
                        end if
                    
                    end do

!             g_e_o_t(t,:)=g_e_o_t(1,:)
                    
                    if ( t .eq. t1-1 ) then
                    
                        open (unit=1,file='output_transition_entry_cutoff.txt')
                        write(1,*)'date d_e_o v_e cost mass' 
                        close(1)
                        
                    end if
                                        
                    open (unit=1,file='output_transition_entry_cutoff.txt',action='write',position='append')
                        do ieps=1,neps
                        write(1,'(300F12.7)')real(t),real(g_e_o_t(t,ieps)),v_i_recmpi(zmap_inv(ieps,intbirth_k,intbirth)),w_t(t)*f_e_t(t)+f_c+k_e,&
                        put_eps(ieps)
                        end do                    
                    close(1)            
                                
                end if
            
            print*,'t, w',t,w_t(t)

                do ieps=1,neps
                do ik=1,nk
                
                g_un_k_t(t,ieps,ik)=g_un_k_recmpi(zmap_inv_un(ieps,ik))
                
!                     if ( g_un_k_t(t,ieps,ik) .ne. g_un_k_t(t+1,ieps,ik) ) then
!                     print*,'g_un_k_t(t,ieps,ik) .ne. g_un_k_t(t+1,ieps,ik)',t
!                         call exit
!                     end if
                    
                end do 
                end do

                do ieps=1,neps
                do ik=1,nk
                
                g_un_o_t(t,ieps,ik)=g_un_o_recmpi(zmap_inv_un(ieps,ik))
                
!                     if ( g_un_o_t(t,ieps,ik) .ne. g_un_o_t(t+1,ieps,ik) ) then
!                     print*,'g_un_o_t(t,ieps,ik) .ne. g_un_o_t(t+1,ieps,ik)',t
!                         call exit
!                     end if
                    
                end do 
                end do
                
                do ieps=1,neps
                do ik=1,nk
                
!                     if ( abs(a_un_t(t,ieps,ik)-a_un_t(t+1,ieps,ik)) .ge. 1e-4 ) then
!                     print*,'abs(a_un_t(t,ieps,ik)-a_un_t(t+1,ieps,ik)) .ge. 1e-4',t
!                         call exit
!                     end if

                end do 
                end do

                do ieps=1,neps
                do ik=1,nk
                
!                     if ( abs(a_tilde_t(t,ieps,ik)-a_tilde_t(t+1,ieps,ik)) .ge. 1e-4 ) then
!                     print*,'abs(a_tilde_t(t,ieps,ik)-a_tilde_t(t+1,ieps,ik)) .ge. 1e-4',t
!                     print*,a_tilde_t(t,ieps,ik),a_tilde_t(t+1,ieps,ik)
!                         call exit
!                     end if

                end do
                end do
                
                do cnt=1,ni
    
                ieps=zmap(cnt,1)
                ik=zmap(cnt,2)
                ia=zmap(cnt,3)
    
                g_i_o_t(t,ieps,ik,ia)=g_i_o_recmpi(cnt)
                g_i_k_t(t,ieps,ik,ia)=g_i_k_recmpi(cnt)
                g_i_a_t(t,ieps,ik,ia)=g_i_a_recmpi(cnt)
                g_q_t(t,ieps,ik,ia)=g_q(ieps,ik,ia)
                
                    if ( rank .eq. 0 .and. cnt .ne. zmap_inv(ieps,ik,ia) ) then
                    print*,'error in uniform_parallel, transition.f90'
                        call exit
                    end if
                    
!                     if ( g_i_o_t(t,ieps,ik,ia) .ne. g_i_o_t(t+1,ieps,ik,ia) ) then
!                     print*,'g_i_o_t(t,ieps,ik,ia) .ne. g_i_o_t(t+1,ieps,ik,ia)',t
!                         call exit
!                     end if
!                     
!                     if ( abs(g_i_a_t(t,ieps,ik,ia)-g_i_a_t(t+1,ieps,ik,ia)) .ge. 1e-4 ) then
!                     print*,'abs(g_i_a_t(t,ieps,ik,ia)-g_i_a_t(t+1,ieps,ik,ia)) .ge. 1e-4',t
!                     print*,g_i_a_t(t,ieps,ik,ia),g_i_a_t(t+1,ieps,ik,ia)
!                         call exit
!                     end if
! 
!                     if ( abs(g_i_k_t(t,ieps,ik,ia)-g_i_k_t(t+1,ieps,ik,ia)) .gt. 0 ) then                    
!                     print*,'abs(g_i_k_t(t,ieps,ik,ia)-g_i_k_t(t+1,ieps,ik,ia)) .gt. 0',t
!                     print*,g_i_k_t(t,ieps,ik,ia),g_i_k_t(t+1,ieps,ik,ia)
!                         call exit
!                     end if
! 
!                     if ( abs(g_q_t(t,ieps,ik,ia)-g_q_t(t+1,ieps,ik,ia)) .ge. 1e-4 ) then
!                     print*,'abs(g_q_t(t,ieps,ik,ia)-g_q_t(t+1,ieps,ik,ia)) .ge. 1e-4',t
!                         call exit
!                     end if
                    
                end do

!                 do ieps=1,neps
!                 
!                     if ( g_e_o_t(t,ieps) .ne. g_e_o_t(t+1,ieps) ) then
!                     print*,'g_e_o_t(t,ieps) .ne. g_e_o_t(t+1,ieps)',t
!                         call exit
!                     end if
!                 
!                 end do

            end if

        end do
        
    deallocate (g_i_o_recmpi,g_i_k_recmpi,g_i_a_recmpi,v_i_recmpi,e_v_i_recmpi)
    deallocate (g_un_o_recmpi,g_un_k_recmpi,v_un_recmpi,e_v_un_recmpi)        
    deallocate (g_i_constrained_recmpi)
    
        !begin simulation 
        if ( rank .eq. 0 ) then

            if ( niota .ne. 2 ) then
            
            omega_ss0=omega_ss_initial
            
            else if ( niota .eq. 2 ) then
                        
            ieps=0
            
                do iz=1,nz
                do ieta=1,neta
                ieps=ieps+1
                omega_ss0(ieps,:,:)=(1-shock)*omega_ss_initial(ieps,:,:)
                omega_ss0(neta*nz+ieps,:,:)=shock*omega_ss_initial(ieps,:,:)            
                end do
                end do
                
            end if
    
            do t=2,t1-1
        
            g_q=g_q_t(t,:,:,:)
            
            !solve for sequence of h_t(t)
    !         h_t(t)=1-psi_t(t)*c_t(t)/w_t(t)
 
                if ( ex_entry .eq. 0 ) then

                h_t(t)=(w_t(t)*(1-tau_l_t(t))/(psi_t(t)*c_t(t)))**(one/phi)
                
                    open (unit=1,file='input_instst_mu.txt')
                    read(1,*)mumin,mumax
                    close(1)
      
                mumin=mumin*.1
                mumax=mumax*1.5

                else if ( ex_entry .eq. 1 ) then
 
                mu_t(t)=ex_mu
            
                end if
            
            d_bi=1
            iter_bi=0

                do while ( iter_bi .le. itmax_bi .and. abs(d_bi) .ge. tol_bi )

                    if ( ex_entry .eq. 0 ) then
                    mu_t(t)=(mumin+mumax)/2
                    end if
                    
                !dividends and labor demand
                div_t(t)=-mu_t(t)*(w_t(t)*f_e_t(t)+f_c+k_e)*sum(put_eps(:)*g_e_o_t(t,:))
                lab_dem=mu_t(t)*f_e_t(t)*sum(put_eps(:)*g_e_o_t(t,:))
    
                gdp_t(t)=mu_t(t)*(w_t(1)*f_e_t(t)+f_c)*sum(put_eps(:)*g_e_o_t(t,:))
                ngdp_t(t)=mu_t(t)*(w_t(t)*f_e_t(t)+f_c)*sum(put_eps(:)*g_e_o_t(t,:)) 
                inv_t(t)=mu_t(t)*(w_t(1)*f_e_t(t)+f_c+k_e)*sum(put_eps(:)*g_e_o_t(t,:))               
                ninv_t(t)=mu_t(t)*(w_t(t)*f_e_t(t)+f_c+k_e)*sum(put_eps(:)*g_e_o_t(t,:))   
                exits_t(t)=zero
                exits_binding_t(t)=zero        
                std_dev_growth_emp_t(t)=zero
                std_dev_growth_sales_t(t)=zero
                tot_adj_cost=zero
    
                    if ( t+1 .ne. t1 ) then
                    firm_debt_t(t+1)=zero                    
                    end if

                b_t(t)=zero
                
                    do ieps=1,neps
                    do ik=1,nk
                    do ia=1,na

                        if ( ( ik .eq. intbirth_k .and. ia .eq. intbirth ) .or. omega_ss0(ieps,ik,ia) .gt. zero )  then
                
                        iz=id_z(ieps)
                        eps=grid_eps_t(t,ieps)
                        k=grid_k(ik)
                        a=grid_a(t,ia)
        
                        d_o=g_i_o_t(t,ieps,ik,ia)
                        jk=g_i_k_t(t,ieps,ik,ia)
                        ap=g_i_a_t(t,ieps,ik,ia)
        
                        kp=grid_k(jk)
        
                            if ( d_o .eq. 1 ) then    
                    
                            l=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
                            
                                if ( ap .ge. zero ) then
                                q=one/(one+r_t(min(t+1,t1-1)))
                                else if ( ap .lt. zero ) then                   
                                
                                    if ( default_economy .eq. 0 ) then
                                    q=one/(one+r_t(min(t+1,t1-1)))-tau_a_t(t)
                                    else if ( default_economy .eq. 1 ) then                                    
                                        call onedlin (grid_a(t+1,:),g_q_t(t,ieps,jk,:),na,ap,q)
                                    end if
                                    
                                end if       

                            profit=((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w_t(t)*(l+f_o_t(t)*grid_z(iz))&
                            +(1.+tau_x_t(t))*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,t)
                                
                                if ( rank .eq. 0 .and. isnan(l) ) then
                                print*,'isnan l, transition.f90',k
                                    call exit
                                end if
            
                            else if ( d_o .eq. 0 ) then            

                            l=zero            
                            
                                if ( default_economy .eq. 0 ) then                                            
                                profit=(1.+tau_x_t(t))*k*chi+a-func_adj(kp,k,d_o,t)
                                else if ( default_economy .eq. 1 ) then                                 
                                profit=max((1.+tau_x_t(t))*k*chi+a-func_adj(kp,k,d_o,t),min(zero,(1.+tau_x_t(t))*k*chi-func_adj(kp,k,d_o,t)))                                
                                end if
                                        
                            end if
            
                            if ( ik .eq. intbirth_k .and. ia .eq. intbirth ) then         
                            
                                if ( cap_adj .ne. 2 ) then
                                lab_dem=lab_dem+(l+d_o*f_o_t(t)*grid_z(iz))*(omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))
                                else if ( cap_adj .eq. 2 ) then
                                lab_dem=lab_dem+(l+d_o*f_o_t(t)*grid_z(iz)+(func_adj(kp,k,d_o,t)-lambda*(kp/k-1)**2.)/w_t(t))*&
                                (omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))
                                end if
                                
                            div_t(t)=div_t(t)+profit*(omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))  
                            gdp_t(t)=gdp_t(t)+(((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))*&
                            (omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))
                            ngdp_t(t)=ngdp_t(t)+(((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))*&
                            (omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))           
                            
                                if ( cap_adj .ne. 2 ) then
                                inv_t(t)=inv_t(t)+&
                                (kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+func_adj(kp,k,d_o,t))*(omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))
                                else if ( cap_adj .eq. 2 ) then
                                inv_t(t)=inv_t(t)+&
                                (kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+(func_adj(kp,k,d_o,t)-lambda*(kp/k-1)**2.)*w_t(1)/w_t(t)&
                                +lambda*(kp/k-1)**2.)*(omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))
                                end if
                                
                            ninv_t(t)=ninv_t(t)+&
                            (kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+func_adj(kp,k,d_o,t))*(omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))
                            exits_t(t)=exits_t(t)+(1-d_o+d_o*delta_f)*(omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps)) 

                                if ( t+1 .ne. t1 ) then                                        
                                firm_debt_t(t+1)=firm_debt_t(t+1)+min(ap,zero)*(omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))
                                end if

                                if ( default_economy .eq. 0 ) then
                                b_t(t)=b_t(t)-a*(omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))
                                else if ( default_economy .eq. 1 ) then                                
                                b_t(t)=b_t(t)+(-a*d_o+(1-d_o)*min(-a,max(zero,(1.+tau_x_t(t))*k*chi-func_adj(kp,k,d_o,t))))*&
                                (omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))                                 
                                end if
                                
                            tot_adj_cost=tot_adj_cost+func_adj(kp,k,d_o,t)*(omega_ss0(ieps,ik,ia)+mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps))
                            
                            else                        
                            
                                if ( cap_adj .ne. 2 ) then                            
                                lab_dem=lab_dem+(l+d_o*f_o_t(t)*grid_z(iz))*omega_ss0(ieps,ik,ia)
                                else if ( cap_adj .eq. 2 ) then
                                lab_dem=lab_dem+(l+d_o*f_o_t(t)*grid_z(iz)+(func_adj(kp,k,d_o,t)-lambda*(kp/k-1)**2.)/w_t(t))*omega_ss0(ieps,ik,ia)
                                end if
                                
                            div_t(t)=div_t(t)+profit*omega_ss0(ieps,ik,ia)       
                            gdp_t(t)=gdp_t(t)+(((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))*&
                            omega_ss0(ieps,ik,ia)
                            ngdp_t(t)=ngdp_t(t)+(((Z_t(t)*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))*&
                            omega_ss0(ieps,ik,ia)        
                            
                                if ( cap_adj .ne. 2 ) then                            
                                inv_t(t)=inv_t(t)+&
                                (kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+func_adj(kp,k,d_o,t))*omega_ss0(ieps,ik,ia)                        
                                else if ( cap_adj .eq. 2 ) then
                                inv_t(t)=inv_t(t)+&
                                (kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+(func_adj(kp,k,d_o,t)-lambda*(kp/k-1)**2.)*w_t(1)/w_t(t)&
                                +lambda*(kp/k-1)**2.)*omega_ss0(ieps,ik,ia)                        
                                end if
                                
                            ninv_t(t)=ninv_t(t)+&
                            (kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+func_adj(kp,k,d_o,t))*omega_ss0(ieps,ik,ia)
                            exits_t(t)=exits_t(t)+(1-d_o+d_o*delta_f)*omega_ss0(ieps,ik,ia)    

                                if ( t+1 .ne. t1 ) then
                                firm_debt_t(t+1)=firm_debt_t(t+1)+min(ap,zero)*omega_ss0(ieps,ik,ia)    
                                end if
                                                            
                                if ( default_economy .eq. 0 ) then
                                b_t(t)=b_t(t)-a*omega_ss0(ieps,ik,ia)                                
                                else if ( default_economy .eq. 1 ) then                                
                                b_t(t)=b_t(t)+(-a*d_o+(1-d_o)*min(-a,max(zero,(1.+tau_x_t(t))*k*chi-func_adj(kp,k,d_o,t))))*omega_ss0(ieps,ik,ia)                                
                                end if
                                
                            tot_adj_cost=tot_adj_cost+func_adj(kp,k,d_o,t)*omega_ss0(ieps,ik,ia)
                                
                            end if
            
                        dummy=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
                                                    
                        end if
                    
                    end do
                    end do
                    end do        
                    
                    if ( ex_entry .eq. 0 ) then
                    
                        if ( lab_dem-h_t(t) .le. zero ) then
                        mumin=mu_t(t)
                        else if ( lab_dem-h_t(t) .gt. zero ) then
                        mumax=mu_t(t)
                        end if

                    d_bi=max(abs(lab_dem-h_t(t)),abs(mumax-mumin))

                    else if ( ex_entry .eq. 1 ) then
                    
                    d_bi=zero
                    
                    end if

                iter_bi=iter_bi+1
    
                end do

                if ( abs(d_bi) .ge. tol_bi ) then
                print*,'iterations exceeded in bisection for mu_t, transition.f90 ',t,d_bi
                end if

                if ( rank .eq. 0 ) then
                print*,'t, mu',t,mu_t(t)
                end if
                
            avg_firm_size_t(t)=h_t(t)/(sum(omega_ss0)+mu_t(t)*sum(put_eps(:)*g_e_o_t(t,:)))
            labor_share_t(t)=w_t(t)*h_t(t)/ngdp_t(t)
            
            firm_capital_t(t)=mu_t(t)*k_e*sum(put_eps(:)*g_e_o_t(t,:))

                do ik=1,nk
                firm_capital_t(t)=firm_capital_t(t)+grid_k(ik)*sum(omega_ss0(:,ik,:))
                end do
            
            !solve for next period distribution        
            omega_ss1=zero
    
                do ieps=1,neps
                do ik=1,nk
                do ia=1,na

                    if ( ( ik .eq. intbirth_k .and. ia .eq. intbirth ) .or. omega_ss0(ieps,ik,ia) .gt. zero )  then
            
                    d_o=g_i_o_t(t,ieps,ik,ia)
                    jk=g_i_k_t(t,ieps,ik,ia)
                
                        if ( d_o .eq. 1 ) then

                        ap=g_i_a_t(t,ieps,ik,ia)
                    
                            call basefun (grid_a(t+1,:),na,ap,vals,inds) 
                    
                            do jeps=1,neps
                            
                                if ( pi_eps(ieps,jeps) .gt. zero ) then
                    
                                    if ( ik .eq. intbirth_k .and. ia .eq. intbirth ) then
                                    omega_ss1(jeps,jk,inds(1))=omega_ss1(jeps,jk,inds(1))+&
                                    vals(1)*pi_eps(ieps,jeps)*(omega_ss0(ieps,ik,ia)+put_eps(ieps)*mu_t(t)*g_e_o_t(t,ieps))*(1-delta_f)
                                    omega_ss1(jeps,jk,inds(2))=omega_ss1(jeps,jk,inds(2))+&
                                    vals(2)*pi_eps(ieps,jeps)*(omega_ss0(ieps,ik,ia)+put_eps(ieps)*mu_t(t)*g_e_o_t(t,ieps)) *(1-delta_f)                           
                                    else
                                    omega_ss1(jeps,jk,inds(1))=omega_ss1(jeps,jk,inds(1))+&
                                    vals(1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)
                                    omega_ss1(jeps,jk,inds(2))=omega_ss1(jeps,jk,inds(2))+&
                                    vals(2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)                         
                                    end if
                            
                                end if
                                
                            end do
                
                        end if
            
                    end if
                
                end do
                end do
                end do
 
                if ( ex_entry .eq. 1 ) then
                h_t(t)=lab_dem
                end if
 
            mass_t(t)=sum(omega_ss0)            
            firms_constrained_t(t)=-1	!sum(omega_ss1(:,:,1))/(sum(omega_ss0)+mu_t(t))
            omega_ss0=omega_ss1                    

            end do
            
            do t=2,t1-1

                if ( ex_entry .eq. 0 ) then 
                
                w_t1(t)=w_t(t)

                else if ( ex_entry .eq. 1 ) then                             
                
                w_t1(t)=(h_t(t)**phi)*(psi_t(t)*c_t(t))/(1-tau_l_t(t))
                                
                    if ( isnan(w_t1(t1)) ) then
                    print*,'isnan(w_t1(t1)), transition.f90',t,h_t(t)
                        call exit
                    end if
                
                end if
            
            c_t(t)=w_t(t)*h_t(t)*(1-tau_l_t(t))+b_t(t)+div_t(t)-b_t(min(t+1,t1-1))/(one+r_t(min(t+1,t1-1)))&
            -tau_a_t(t)*firm_debt_t(t+1)+tau_l_t(t)*w_t(t)*h_t(t)+(ninv_t(t)-&
            sum(put_eps(:)*g_e_o_t(t,:))*mu_t(t)*(w_t(t)*f_e_t(t)+f_c+k_e)-tot_adj_cost)*tau_x_t(t)                                                           

            exdem_t(t)=abs(ngdp_t(t)-c_t(t)-ninv_t(t))
            
            end do
            
            do t=3,t1-1        
            r_t1(t)=c_t(t)/(beta*c_t(t-1))-one
            end do
                
            open (unit=1,file='check.txt')
            write(1,*)'t r r1 w w1'
                do t=1,t1
                write(1,'(300F12.7)')real(t),r_t(t),r_t1(t),w_t(t),w_t1(t)
                end do
            close(1)
                    
!         d_tr=max(maxval(abs(b_t+firm_debt_t)),maxval(abs(r_t-r_t1)))
        d_tr=max(maxval(abs(r_t-r_t1)),maxval(abs(w_t-w_t1)))
        iter_tr=iter_tr+1
    
        print*,'iter_tr, d_tr ',iter_tr,d_tr
!         print*,maxval(abs(b_t+firm_debt_t)),maxval(abs(r_t-r_t1))
!         print*,maxloc(abs(b_t+firm_debt_t)),maxloc(abs(r_t-r_t1))        

            if ( d_tr .ge. tol_tr ) then
            r_t=weight*r_t+(1-weight)*r_t1
!             w_t=.7*w_t+(1-.7)*w_t1
            w_t=weight_w*w_t+(1-weight_w)*w_t1
            end if    
        
            open (unit=1,file='output_transition_aggregates.txt')
            write(1,*)'t z psi theta tau_a f_e f_o tau_x tau_l C H B r w div entry mass rgdp rinv firms_constrained exit exit_b'
            write(1,*)'avg_firm_size labor_share std_dev_growth_emp std_dev_growth_sales debt capital ngdp ninv exdem'
                do t=1,t1
                write(1,'(300F12.7)')real(t),z_t(t),psi_t(t),theta_t(t),tau_a_t(t),f_e_t(t),f_o_t(t),tau_x_t(t),tau_l_t(t),&
                c_t(t),h_t(t),b_t(t),r_t(t),w_t(t),div_t(t),mu_t(t)*sum(put_eps(:)*g_e_o_t(t,:)*g_i_o_t(t,:,intbirth_k,intbirth)),mass_t(t),&                
                gdp_t(t),inv_t(t),firms_constrained_t(t),exits_t(t),exits_binding_t(t),avg_firm_size_t(t),labor_share_t(t),&
                std_dev_growth_emp_t(t),std_dev_growth_sales_t(t),firm_debt_t(t),firm_capital_t(t),ngdp_t(t),ninv_t(t),&
                exdem_t(t)
                end do
            close(1)
	
            open (unit=1,file='input_r_t.txt')
                 do t=3,t1-1
                 write(1,*)r_t(t)
                 end do
            close(1)

            if ( ex_entry .eq. 1 ) then
            
                open (unit=1,file='input_w_t.txt')
                    do t=2,t1-1
                    write(1,*)w_t(t)
                    end do
                close(1)

            end if
            
        end if

        call MPI_BCAST (d_tr,1,mpi_double_precision,0,mpi_comm_world,ierr)        
        call MPI_BCAST (iter_tr,1,mpi_integer,0,mpi_comm_world,ierr)        
        call MPI_BCAST (r_t,t1,mpi_double_precision,0,mpi_comm_world,ierr)
        call MPI_BCAST (w_t,t1,mpi_double_precision,0,mpi_comm_world,ierr)
        
    end do

    !end simulation
    
    if ( rank .eq. 0 ) then
    
    allocate (omega_ss_age0(120,neps,nk,na), omega_ss_age1(120,neps,nk,na))
    
    !compute employment by age    
    
        if ( niota .ne. 2 ) then
            
        omega_ss_age0=omega_ss_age_initial
            
        else if ( niota .eq. 2 ) then
                        
        ieps=0
            
            do iz=1,nz
            do ieta=1,neta
            ieps=ieps+1
            omega_ss_age0(:,ieps,:,:)=(1-shock)*omega_ss_age_initial(:,ieps,:,:)
            omega_ss_age0(:,neta*nz+ieps,:,:)=shock*omega_ss_age_initial(:,ieps,:,:)            
            end do
            end do
                
        end if    

        do t=2,t1-1

        omega_ss_age0(1,:,:,:)=zero
        omega_ss_age0(1,:,intbirth_k,intbirth)=mu_t(t)*put_eps(:)*g_e_o_t(t,:)
    
        !compute firms by age and employment by age
    
        emp_by_age_t(t,:)=zero
        exit_by_age_t(t,:)=zero
        exit_binding_by_age_t(t,:)=zero
    
            do tc=1,120
        
            firms_by_age_t(t,tc)=sum(omega_ss_age0(tc,:,:,:))
            
                do ieps=1,neps
                do ik=1,nk
                do ia=1,na
        
                iz=id_z(ieps)
                eps=grid_eps_t(t,ieps)
                k=grid_k(ik)
                a=grid_a(t,ia)
        
                d_o=g_i_o_t(t,ieps,ik,ia)
                jk=g_i_k_t(t,ieps,ik,ia)
                ap=g_i_a_t(t,ieps,ik,ia)
        
                kp=grid_k(jk)
        
                    if ( d_o .eq. 1 ) then                        
                    l=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))                
                    else if ( d_o .eq. 0 ) then                        
                    l=zero                        
                    end if
                    
                    if ( cap_adj .eq. 2 ) then
                    l=l+(func_adj(kp,k,d_o,t)-lambda*(kp/k-1)**2.)/w_t(t)
                    end if                        
                                    
                    if ( tc .eq. 1 ) then
                    l=l+f_e_t(t)
                    end if
                
                emp_by_age_t(t,tc)=emp_by_age_t(t,tc)+(l+d_o*f_o_t(t)*grid_z(iz))*omega_ss_age0(tc,ieps,ik,ia)
                exit_by_age_t(t,tc)=exit_by_age_t(t,tc)+(1-d_o+d_o*delta_f)*omega_ss_age0(tc,ieps,ik,ia)
            
                dummy=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
            
                end do
                end do
                end do
    
            end do
        
        !solve for next period distribution        
        omega_ss_age1=zero
    
            do tc=1,120
        
                if ( tc .eq. 120 ) then
                tcc=120
                else if ( tc .le. 119 ) then
                tcc=tc+1
                end if
            
                do ieps=1,neps
                do ik=1,nk
                do ia=1,na

                    if ( omega_ss_age0(tc,ieps,ik,ia) .gt. zero )  then
            
                    d_o=g_i_o_t(t,ieps,ik,ia)
                    jk=g_i_k_t(t,ieps,ik,ia)
                    
                        if ( d_o .eq. 1 ) then

                        ap=g_i_a_t(t,ieps,ik,ia)
                    
                            call basefun (grid_a(t+1,:),na,ap,vals,inds) 
                    
                            do jeps=1,neps
                            
                                if ( pi_eps(ieps,jeps) .gt. zero ) then     
                                omega_ss_age1(tcc,jeps,jk,inds(1))=omega_ss_age1(tcc,jeps,jk,inds(1))+&
                                vals(1)*pi_eps(ieps,jeps)*omega_ss_age0(tc,ieps,ik,ia)*(1-delta_f)
                                omega_ss_age1(tcc,jeps,jk,inds(2))=omega_ss_age1(tcc,jeps,jk,inds(2))+&
                                vals(2)*pi_eps(ieps,jeps)*omega_ss_age0(tc,ieps,ik,ia)*(1-delta_f)                            
                                end if
                                
                            end do
                
                        end if
            
                    end if
                
                end do
                end do
                end do
            
            end do
        
        omega_ss_age0=omega_ss_age1

        end do
                
    deallocate (omega_ss_age0, omega_ss_age1)
    
    !compute job creation and job destruction
        if ( niota .ne. 2 ) then
            
        omega_ss0=omega_ss_initial
            
        else if ( niota .eq. 2 ) then
                        
        ieps=0
            
            do iz=1,nz
            do ieta=1,neta
            ieps=ieps+1
            omega_ss0(ieps,:,:)=(1-shock)*omega_ss_initial(ieps,:,:)
            omega_ss0(neta*nz+ieps,:,:)=shock*omega_ss_initial(ieps,:,:)            
            end do
            end do
                
        end if
    
        do t=1,t1-2
    
        job_creation_t(t+1)=zero
        job_destruction_t(t+1)=zero

        omega_ss1=zero
    
            do ieps=1,neps
            do ik=1,nk
            do ia=1,na    
    
                if ( omega_ss0(ieps,ik,ia) .gt. zero )  then
        
                iz=id_z(ieps)
                eps=grid_eps_t(t,ieps)
                k=grid_k(ik)
                a=grid_a(t,ia)
        
                d_o=g_i_o_t(t,ieps,ik,ia)
                jk=g_i_k_t(t,ieps,ik,ia)
                ap=g_i_a_t(t,ieps,ik,ia)
        
                kp=grid_k(jk)
        
                    if ( d_o .eq. 1 ) then                    
                    l=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o_t(t)*grid_z(iz)                
                    else if ( d_o .eq. 0 ) then                        
                    l=zero                                        
                    end if

                    if ( cap_adj .eq. 2 ) then
                    l=l+(func_adj(kp,k,d_o,t)-lambda*(kp/k-1)**2.)/w_t(t)
                    end if                        
                    
                    if ( d_o .eq. 1 ) then
            
                        call basefun (grid_a(t+1,:),na,ap,vals,inds) 
                    
                        do jeps=1,neps
                        
                        jz=id_z(jeps)
                        
                            if ( pi_eps(ieps,jeps) .gt. zero ) then                    
                            
                            ja=inds(1)

                                if ( g_i_o_t(t+1,jeps,jk,ja) .eq. 1 ) then                    
                                lp=(w_t(t+1)/(((Z_t(t+1)*grid_eps_t(t+1,jeps))**(1-alpha*nu))&
                                *(1-alpha)*nu*(kp**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o_t(t+1)*grid_z(jz)         
                                else if ( g_i_o_t(t+1,jeps,jk,ja) .eq. 0 ) then                        
                                lp=zero                                        
                                end if

                                if ( cap_adj .eq. 2 ) then
                                lp=lp+(func_adj(grid_k(g_i_k_t(t+1,jeps,jk,ja)),kp,g_i_o_t(t+1,jeps,jk,ja),t+1)-&
                                lambda*(grid_k(g_i_k_t(t+1,jeps,jk,ja))/kp-1)**2.)/w_t(t+1)
                                end if                        
                                
                                if ( lp-l .ge. zero ) then
                                job_creation_t(t+1)=job_creation_t(t+1)+&
                                (lp-l)*vals(1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)
                                else if ( lp-l .lt. zero ) then
                                job_destruction_t(t+1)=job_destruction_t(t+1)+&
                                (l-lp)*vals(1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)
                                end if
                    
                            ja=inds(2)
    
                                if ( g_i_o_t(t+1,jeps,jk,ja) .eq. 1 ) then                    
                                lp=(w_t(t+1)/(((Z_t(t+1)*grid_eps_t(t+1,jeps))**(1-alpha*nu))&
                                *(1-alpha)*nu*(kp**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o_t(t+1)*grid_z(jz)                                   
                                else if ( g_i_o_t(t+1,jeps,jk,ja) .eq. 0 ) then                        
                                lp=zero                                        
                                end if

                                if ( cap_adj .eq. 2 ) then
                                lp=lp+(func_adj(grid_k(g_i_k_t(t+1,jeps,jk,ja)),kp,g_i_o_t(t+1,jeps,jk,ja),t+1)-&
                                lambda*(grid_k(g_i_k_t(t+1,jeps,jk,ja))/kp-1)**2.)/w_t(t+1)
                                end if                        
                                
                                if ( lp-l .ge. zero ) then
                                job_creation_t(t+1)=job_creation_t(t+1)+&
                                (lp-l)*vals(2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f) 
                                else if ( lp-l .lt. zero ) then
                                job_destruction_t(t+1)=job_destruction_t(t+1)+&
                                (l-lp)*vals(2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f) 
                                end if
                        
                            omega_ss1(jeps,jk,inds(1))=omega_ss1(jeps,jk,inds(1))+&
                            vals(1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)                            
                            omega_ss1(jeps,jk,inds(2))=omega_ss1(jeps,jk,inds(2))+&
                            vals(2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)                        
                            
                            end if
                            
                        end do

                    job_destruction_t(t+1)=job_destruction_t(t+1)+(l-zero)*omega_ss0(ieps,ik,ia)*delta_f
                    
                    else if ( d_o .eq. 0 ) then
            
                    lp=zero
                    job_destruction_t(t+1)=job_destruction_t(t+1)+(l-lp)*omega_ss0(ieps,ik,ia)
            
                    end if                
                
                end if
        
                if ( ik .eq. intbirth_k .and. ia .eq. intbirth ) then
        
                iz=id_z(ieps)
                eps=grid_eps_t(t,ieps)
                k=grid_k(ik)
                a=grid_a(t,ia)
        
                d_o=g_i_o_t(t,ieps,ik,ia)
                jk=g_i_k_t(t,ieps,ik,ia)
                ap=g_i_a_t(t,ieps,ik,ia)
        
                kp=grid_k(jk)
        
                    if ( d_o .eq. 1 ) then                    
                    l=(w_t(t)/(((Z_t(t)*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))&
                    +f_o_t(t)*grid_z(iz)+f_e_t(t)        
                    else if ( d_o .eq. 0 ) then                        
                    l=f_e_t(t)                                        
                    end if

                    if ( cap_adj .eq. 2 ) then
                    l=l+(func_adj(kp,k,d_o,t)-lambda*(kp/k-1)**2.)/w_t(t)
                    end if                        
                    
                    if ( d_o .eq. 1 ) then
            
                        call basefun (grid_a(t+1,:),na,ap,vals,inds) 
                    
                        do jeps=1,neps                    
                        
                        jz=id_z(jeps)
                        
                            if ( pi_eps(ieps,jeps) .gt. zero ) then

                            ja=inds(1)

                                if ( g_i_o_t(t+1,jeps,jk,ja) .eq. 1 ) then                    
                                lp=(w_t(t+1)/(((Z_t(t+1)*grid_eps_t(t+1,jeps))**(1-alpha*nu))&
                                *(1-alpha)*nu*(kp**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o_t(t+1)*grid_z(jz)         
                                else if ( g_i_o_t(t+1,jeps,jk,ja) .eq. 0 ) then                        
                                lp=zero                                        
                                end if

                                if ( cap_adj .eq. 2 ) then
                                lp=lp+(func_adj(grid_k(g_i_k_t(t+1,jeps,jk,ja)),kp,g_i_o_t(t+1,jeps,jk,ja),t+1)-&
                                lambda*(grid_k(g_i_k_t(t+1,jeps,jk,ja))/kp-1)**2.)/w_t(t+1)
                                end if                        
                                
                                if ( lp-l .ge. zero ) then
                                job_creation_t(t+1)=job_creation_t(t+1)+&
                                (lp-l)*vals(1)*pi_eps(ieps,jeps)*mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps)*(1-delta_f)
                                else if ( lp-l .lt. zero ) then
                                job_destruction_t(t+1)=job_destruction_t(t+1)+&
                                (l-lp)*vals(1)*pi_eps(ieps,jeps)*mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps)*(1-delta_f)
                                end if
                    
                            ja=inds(2)

                                if ( g_i_o_t(t+1,jeps,jk,ja) .eq. 1 ) then                    
                                lp=(w_t(t+1)/(((Z_t(t+1)*grid_eps_t(t+1,jeps))**(1-alpha*nu))&
                                *(1-alpha)*nu*(kp**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o_t(t+1)*grid_z(jz)         
                                else if ( g_i_o_t(t+1,jeps,jk,ja) .eq. 0 ) then                        
                                lp=zero                                        
                                end if
                                
                                if ( cap_adj .eq. 2 ) then
                                lp=lp+(func_adj(grid_k(g_i_k_t(t+1,jeps,jk,ja)),kp,g_i_o_t(t+1,jeps,jk,ja),t+1)-&
                                lambda*(grid_k(g_i_k_t(t+1,jeps,jk,ja))/kp-1)**2.)/w_t(t+1)
                                end if                        
    
                                if ( lp-l .ge. zero ) then
                                job_creation_t(t+1)=job_creation_t(t+1)+&
                                (lp-l)*vals(2)*pi_eps(ieps,jeps)*mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps)*(1-delta_f)
                                else if ( lp-l .lt. zero ) then
                                job_destruction_t(t+1)=job_destruction_t(t+1)+&
                                (l-lp)*vals(2)*pi_eps(ieps,jeps)*mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps)*(1-delta_f)
                                end if
                    
                            omega_ss1(jeps,jk,inds(1))=omega_ss1(jeps,jk,inds(1))+&
                            vals(1)*pi_eps(ieps,jeps)*mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps)*(1-delta_f)                            
                            omega_ss1(jeps,jk,inds(2))=omega_ss1(jeps,jk,inds(2))+&
                            vals(2)*pi_eps(ieps,jeps)*mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps)*(1-delta_f)                      
                            
                            end if
                            
                        end do

                    job_destruction_t(t+1)=job_destruction_t(t+1)+(l-zero)*mu_t(t)*put_eps(ieps)*g_e_o_t(t,ieps)*delta_f
                    
                    else if ( d_o .eq. 0 ) then
            
                    lp=zero
                    job_destruction_t(t+1)=job_destruction_t(t+1)+(l-lp)*mu_t(t)*g_e_o_t(t,ieps)*put_eps(ieps)
                    
                    end if 
                
                d_o=g_i_o_t(t+1,ieps,ik,ia)
                jk=g_i_k_t(t+1,ieps,ik,ia)
                ap=g_i_a_t(t+1,ieps,ik,ia)
                
                kp=grid_k(jk)
        
                    if ( d_o .eq. 1 ) then                    
                    l=(w_t(t+1)/(((Z_t(t+1)*grid_eps_t(t+1,ieps))**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))&
                    +f_o_t(t+1)*grid_z(iz)+f_e_t(t+1)        
                    else if ( d_o .eq. 0 ) then                        
                    l=f_e_t(t+1)                                        
                    end if              
                    
                    if ( cap_adj .eq. 2 ) then
!                     l=l+func_adj(kp,k,d_o,t+1)/w_t(t+1)
                    l=l+(func_adj(kp,k,d_o,t+1)-lambda*(kp/k-1)**2.)/w_t(t+1)
                    end if
                            
                job_creation_t(t+1)=job_creation_t(t+1)+l*mu_t(t+1)*put_eps(ieps)*g_e_o_t(t+1,ieps)
            
                end if
                    
            end do
            end do
            end do
    
        omega_ss0=omega_ss1
    
        end do
        
        open (unit=1,file='output_transition_firms_age.txt')
            do t=1,t1
            write(1,'(300F12.7)')real(t),firms_by_age_t(t,:),exit_by_age_t(t,:)
            end do
        close(1)    
    
        open (unit=1,file='output_transition_emp_age.txt')
            do t=1,t1
            write(1,'(300F12.7)')real(t),emp_by_age_t(t,:),job_creation_t(t),job_destruction_t(t)
            end do
        close(1)  
        
        if ( d_tr .ge. tol_tr ) then
        print*,'d_tr .ge. tol_tr, transition.f90 ',d_tr
        end if
        
    print*,'completed transition'
    
        deallocate (g_i_o_t, g_i_k_t, g_i_a_t, g_e_o_t, g_un_k_t, g_un_o_t)

    end if
    
end subroutine transition
