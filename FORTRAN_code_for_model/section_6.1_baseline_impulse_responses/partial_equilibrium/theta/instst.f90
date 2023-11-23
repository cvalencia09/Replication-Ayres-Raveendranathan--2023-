subroutine instst
use params
use mini
implicit none
include 'mpif.h'

integer :: iter_bi
double precision :: d_bi

double precision :: wmin, wmax
double precision :: mumin, mumax

double precision, dimension(neps,nk,na) :: omega_ss0, omega_ss1
integer :: iter_dist
double precision :: d_dist

double precision, dimension(neps,nk,na,2) :: vals
integer, dimension(neps,nk,na,2) :: inds

double precision :: sales0, sales1, surv
double precision :: avg_sales0, avg_sales1
double precision :: std_sales0, std_sales1
double precision :: dummy
double precision :: temp, temp_k
double precision :: tot_inv_wedge

integer :: cnt_i_k
double precision, allocatable, dimension(:,:) :: dist_i_k, dist_i_k1
double precision, parameter :: std_dev_i_k_data=.337

integer :: ik0, ieps0, ia0, iter_simulate
double precision, dimension(neps,nk,na) :: omega_ss2, omega_ss3

double precision, allocatable, dimension(:,:,:,:) :: omega_temp0, omega_temp1
double precision, allocatable, dimension(:) :: job_creation_mpi, job_destruction_mpi
double precision, allocatable, dimension(:) :: job_creation_recmpi, job_destruction_recmpi

date=1

z=z_t(date)
abar=abar_t(date)
tau_a=tau_a_t(date)
psi=psi_t(date)
f_e=f_e_t(date)
grid_eps=grid_eps_t(date,:)
f_o=f_o_t(date)
tau_x=tau_x_t(date)
tau_l=tau_l_t(date)
theta=theta_t(date)

r=(one/beta-1)

    if ( ex_entry .eq. 0 ) then
    
    !bisection for w   
        open (unit=1,file='input_instst_w.txt')
        read(1,*)wmin,wmax
        close(1)

    wmax=1.5*wmax
    wmin=.5*wmin

    iter_bi=0   
    d_bi=1

        if ( rank .eq. 0 ) then
        print*,'instst '
        print*,'iter_bi, w, d_bi '
        end if
    
        do while ( iter_bi .le. itmax_bi .and. abs(d_bi) .ge. tol_bi )
    
        w=(wmin+wmax)/2
        w_t(date)=w
    
            call vfn_un
            
            call vfn_i
        
        v_e=sum(put_eps(:)*v_i(:,intbirth_k,intbirth))
        
            if ( v_e-w*f_e-f_c-k_e .ge. zero ) then
            wmin=w
            else if ( v_e-w*f_e-f_c-k_e .lt. zero ) then
            wmax=w
            end if

        d_bi=max(abs(v_e-w*f_e-f_c-k_e),abs(wmax-wmin))
        iter_bi=iter_bi+1

            if ( rank .eq. 0 ) then
        
            write(6,'(300F12.7)')real(iter_bi),w,v_e,w*f_e+f_c+k_e
        
                open (unit=1,file='input_instst_w.txt')
                write(1,*)wmin,wmax
                close(1)
            
            end if
    
        end do
    
        if ( rank .eq. 0 .and. d_bi .ge. tol_bi ) then
        print*,'iterations exceeded at zero profit condition, instst.f90 ',d_bi
            call exit
        end if

    end if
        
    if ( ex_entry .eq. 0 ) then
    
    !bisection for mu    
     open (unit=1,file='input_instst_mu.txt')
     read(1,*)mumin,mumax
     close(1)

    mumin=mumin*.5
    mumax=mumax*2.0

    else if ( ex_entry .eq. 1 ) then

    mu=ex_mu
    
        open (unit=1,file='input_instst_w.txt')
        read(1,*)wmin,wmax
        close(1)    
        
!     wmax=1.2*wmax
!     wmin=.8*wmin
    
    end if
    
d_bi=1
iter_bi=0

    if ( rank .eq. 0 ) then
    
    print*,' '
    
        if ( ex_entry .eq. 0 ) then    
        print*,'iter_bi, mu, d_bi'
        else if ( ex_entry .eq. 1 ) then
        print*,'iter_bi, w, d_bi'
        end if
        
    end if

    allocate(dist_i_k(ni,2))

    do while ( iter_bi .le. itmax_bi .and. abs(d_bi) .ge. tol_bi )
        
        if ( ex_entry .eq. 0 ) then        
        
        mu=(mumin+mumax)/2
        g_e_o=1
        
        else if ( ex_entry .eq. 1 ) then        
        
        w=(wmin+wmax)/2
        w_t(date)=w
    
            call vfn_un
            
            call vfn_i        
        
            do ieps=1,neps

                if ( v_i(ieps,intbirth_k,intbirth) .ge. w*f_e+f_c+k_e ) then
                g_e_o(ieps)=1
                else if ( v_i(ieps,intbirth_k,intbirth) .lt. w*f_e+f_c+k_e ) then
                g_e_o(ieps)=0
                end if
                    
            end do

            if ( rank .eq. 0 ) then
            
                open(unit=1,file='output_instst_entry_cutoff.txt')
                write(1,*)'date d_e_o v_e cost mass' 
                    do ieps=1,neps
                    write(1,'(300F12.7)')real(date),real(g_e_o(ieps)),v_i(ieps,intbirth_k,intbirth),w*f_e+f_c+k_e,put_eps(ieps)
                    end do                    
                close(1)
                
            end if
            
        cut_off_model=grid_eps_original((neps-sum(g_e_o))+1)
                
        end if

        if ( ( ex_entry .eq. 0 .and. iter_bi .eq. 0 ) .or. ex_entry .eq. 1 ) then
        
            do ieps=1,neps
            do ik=1,nk
            do ia=1,na

            ap=g_i_a(ieps,ik,ia)
    
                call basefun (grid_a(date,:),na,ap,vals(ieps,ik,ia,:),inds(ieps,ik,ia,:))     
        
            end do
            end do
            end do
        
        end if
        
    !simulation
    omega_ss0=zero
    d_dist=1
    iter_dist=0

        do while ( iter_dist .le. itmax_dist .and. d_dist .ge. tol_dist )
    
        omega_ss1=zero

        avg_sales0=zero    
        avg_sales1=zero
        surv=zero
        
            do ieps=1,neps
            do ik=1,nk
            do ia=1,na

                if ( ( ik .eq. intbirth_k .and. ia .eq. intbirth ) .or. omega_ss0(ieps,ik,ia) .gt. zero )  then
            
                d_o=g_i_o(ieps,ik,ia)                                
                jk=g_i_k(ieps,ik,ia)

                    if ( d_o .eq. 1 ) then
            
                        do jeps=1,neps

                            if ( pi_eps(ieps,jeps) .gt. zero ) then
                        
                                if ( ik .eq. intbirth_k .and. ia .eq. intbirth ) then
                            
                                omega_ss1(jeps,jk,inds(ieps,ik,ia,1))=omega_ss1(jeps,jk,inds(ieps,ik,ia,1))+&
                                vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*(omega_ss0(ieps,ik,ia)+put_eps(ieps)*mu*g_e_o(ieps))*(1-delta_f)
                                omega_ss1(jeps,jk,inds(ieps,ik,ia,2))=omega_ss1(jeps,jk,inds(ieps,ik,ia,2))+&
                                vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*(omega_ss0(ieps,ik,ia)+put_eps(ieps)*mu*g_e_o(ieps))*(1-delta_f)                            

                                else
                                
                                omega_ss1(jeps,jk,inds(ieps,ik,ia,1))=omega_ss1(jeps,jk,inds(ieps,ik,ia,1))+&
                                vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)
                                omega_ss1(jeps,jk,inds(ieps,ik,ia,2))=omega_ss1(jeps,jk,inds(ieps,ik,ia,2))+&
                                vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)                            
                            
                                end if
                        
                            end if
                        
                        end do
                
                    end if
            
                end if
                
            end do
            end do
            end do
        
        iter_dist=iter_dist+1
        d_dist=maxval(abs(omega_ss0-omega_ss1))
        omega_ss0=omega_ss1

        end do
        
        if ( rank .eq. 0 .and. d_dist .ge. tol_dist ) then
        print*,'iterations exceeded in simulation, instst.f90 ',d_dist
!            call exit
        end if
    
    !dividends and labor demand
    div=-mu*(w*f_e+f_c+k_e)*sum(put_eps*g_e_o)
    lab_dem=mu*f_e*sum(put_eps*g_e_o)    
    gdp=mu*(w*f_e+f_c)*sum(put_eps*g_e_o)
    ngdp=mu*(w*f_e+f_c)*sum(put_eps*g_e_o)
    inv=mu*(w*f_e+f_c+k_e)*sum(put_eps*g_e_o)
    ninv=mu*(w*f_e+f_c+k_e)*sum(put_eps*g_e_o)
    firms_constrained=-1  !sum(omega_ss0(:,:,1))/(sum(omega_ss0)+mu)
    exits=zero
    exits_binding=zero
    std_dev_growth_emp=zero
    std_dev_growth_sales=zero
    avg_sales=zero
    var_sales=zero
    
    auto_cov_sales=zero
    std_sales0=zero
    std_sales1=zero

    b=zero
    neg_b=zero
    
    avg_i_k=zero
    std_dev_i_k=zero    
    i_k_inaction=zero
    i_k_pos=zero
    i_k_neg=zero
    i_k_pos_spike=zero
    i_k_neg_spike=zero    
    temp=zero
    temp_k=zero
    
    tot_inv_wedge=zero
    
    cnt_i_k=0    
    dist_i_k=0
    
!         if ( rank .eq. 0 ) then
!             open (unit=1,file='output_instst_i_k_dist.txt')
!             write(1,*)'k kp i_k mass'
!             close(1)
!         end if
    
        do ieps=1,neps
        do ik=1,nk
        do ia=1,na

            if ( ( ik .eq. intbirth_k .and. ia .eq. intbirth ) .or. omega_ss0(ieps,ik,ia) .gt. zero )  then
        
            eps=grid_eps(ieps)
            k=grid_k(ik)
            a=grid_a(date,ia)
            iz=id_z(ieps)
        
            d_o=g_i_o(ieps,ik,ia)
            jk=g_i_k(ieps,ik,ia)
            ap=g_i_a(ieps,ik,ia)        
            kp=grid_k(jk)
        
                if ( d_o .eq. 1 ) then                    
                
                l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
                
                    if ( ap .ge. zero ) then
                                
                    q=one/(one+r)
                                
                    else if ( ap .lt. zero ) then
                                
                        if ( default_economy .eq. 0 ) then                                
                        q=one/(one+r)-tau_a
                        else if ( default_economy .eq. 1 ) then
                            call onedlin (grid_a(date,:),g_q(ieps,jk,:),na,ap,q) 
                        end if
                            
                    end if
                                
                profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))&
                +(1.+tau_x)*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,date)            
                    
                    if ( rank .eq. 0 .and. isnan(l) ) then
                    print*,'isnan l, instst.f90',k
                        call exit
                    end if
                
                else if ( d_o .eq. 0 ) then                        
                
                l=zero         

                    if ( default_economy .eq. 0 ) then                            
                    profit=(1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date)                                        
                    else if ( default_economy .eq. 1 ) then
                    profit=max((1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date),min(zero,(1.+tau_x)*k*chi-func_adj(kp,k,d_o,date)))                                    
                    end if            
            
                end if
            
                if ( ik .eq. intbirth_k .and. ia .eq. intbirth ) then       

                    if ( cap_adj .ne. 2 ) then
                    lab_dem=lab_dem+(l+d_o*f_o*grid_z(iz))*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                    else if ( cap_adj .eq. 2 ) then
                    lab_dem=lab_dem+(l+d_o*f_o*grid_z(iz)+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)/w)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                    end if
                    
                div=div+profit*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))                       
                gdp=gdp+(((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))*&
                (omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                ngdp=ngdp+(((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))*&
                (omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))                   
                
                    if ( cap_adj .eq. 2 ) then
                    gdp=gdp+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                    ngdp=ngdp+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))                   
                    end if
                                        
                inv=inv+(kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+func_adj(kp,k,d_o,date))*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                ninv=ninv+(kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+func_adj(kp,k,d_o,date))*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))                
                exits=exits+(1-d_o+d_o*delta_f)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))   
                avg_sales=avg_sales+(((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))*&
                (omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps)) 
                var_sales=var_sales+((((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))**2)*&
                (omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps)) 
                
                    if ( default_economy .eq. 0 ) then                    
                    b=b-a*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps)) 
                    else if ( default_economy .eq. 1 ) then                    
                    b=b+(-a*d_o+(1-d_o)*min(-a,max(zero,(1.+tau_x)*k*chi-func_adj(kp,k,d_o,date))))*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps)) 
                    end if
                    
                neg_b=neg_b-min(ap,zero)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps)) 
                tot_inv_wedge=tot_inv_wedge+tau_x*(kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))

                    if ( d_o .eq. 1 .and. k .gt. zero ) then
                    
                    x_i=kp-(1-delta_k)*k*d_o-chi*k*(1-d_o)

                    temp=temp+(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))

                    cnt_i_k=cnt_i_k+1
                    dist_i_k(cnt_i_k,1)=x_i/k
                    dist_i_k(cnt_i_k,2)=omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps)                   
                    
!                         if ( rank .eq. 0 ) then                        
!                             open(unit=1,file='output_instst_i_k_dist.txt',action='write',position='append')        
!                             write(1,'(300F14.9)')k,kp,x_i/k,omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)
!                             close(1)                            
!                         end if
                        
                        if ( x_i/k .gt. -.01 .and. x_i/k .lt. .01 ) then                        
                        i_k_inaction=i_k_inaction+(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                        end if
                        
                        if ( x_i/k .ge. .01 ) then
                        i_k_pos=i_k_pos+(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                        end if
                        
                        if ( x_i/k .le. -.01 ) then
                        i_k_neg=i_k_neg+(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                        end if
                        
                        if ( x_i/k .gt. .2 ) then
                        i_k_pos_spike=i_k_pos_spike+(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                        end if
                        
                        if ( x_i/k .lt. -.2 ) then
                        i_k_neg_spike=i_k_neg_spike+(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                        end if
                        
                    end if
                    
                else      
                
                    if ( cap_adj .ne. 2 ) then     
                    lab_dem=lab_dem+(l+d_o*f_o*grid_z(iz))*omega_ss0(ieps,ik,ia)
                    else if ( cap_adj .eq. 2 ) then                
                    lab_dem=lab_dem+(l+d_o*f_o*grid_z(iz)+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)/w)*omega_ss0(ieps,ik,ia)
                    end if
                    
                div=div+profit*omega_ss0(ieps,ik,ia)       
                gdp=gdp+(((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))*&
                omega_ss0(ieps,ik,ia)                   
                ngdp=ngdp+(((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))*&
                omega_ss0(ieps,ik,ia)   

                    if ( cap_adj .eq. 2 ) then
                    gdp=gdp+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)*omega_ss0(ieps,ik,ia)
                    ngdp=ngdp+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)*omega_ss0(ieps,ik,ia)              
                    end if                
                
                inv=inv+(kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+func_adj(kp,k,d_o,date))*omega_ss0(ieps,ik,ia)                       
                ninv=ninv+(kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+func_adj(kp,k,d_o,date))*omega_ss0(ieps,ik,ia)                                                   
                exits=exits+(1-d_o+d_o*delta_f)*omega_ss0(ieps,ik,ia)
                avg_sales=avg_sales+(((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))*&
                omega_ss0(ieps,ik,ia) 
                var_sales=var_sales+((((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))**2)*&
                omega_ss0(ieps,ik,ia)

                    if ( default_economy .eq. 0 ) then                    
                    b=b-a*omega_ss0(ieps,ik,ia)
                    else if ( default_economy .eq. 1 ) then                    
                    b=b+(-a*d_o+(1-d_o)*min(-a,max(zero,(1.+tau_x)*k*chi-func_adj(kp,k,d_o,date))))*omega_ss0(ieps,ik,ia)
                    end if
                
                neg_b=neg_b-min(ap,zero)*omega_ss0(ieps,ik,ia)
                tot_inv_wedge=tot_inv_wedge+tau_x*(kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k)*omega_ss0(ieps,ik,ia)
                
                     if ( d_o .eq. 1 .and. k .gt. zero ) then

                    x_i=kp-(1-delta_k)*k*d_o-chi*k*(1-d_o)
                    temp=temp+omega_ss0(ieps,ik,ia)                    

                    cnt_i_k=cnt_i_k+1
                    dist_i_k(cnt_i_k,1)=x_i/k
                    dist_i_k(cnt_i_k,2)=omega_ss0(ieps,ik,ia)                

!                         if ( rank .eq. 0 ) then                        
!                             open(unit=1,file='output_instst_i_k_dist.txt',action='write',position='append')        
!                             write(1,'(300F14.9)')k,kp,x_i/k,omega_ss0(ieps,ik,ia)
!                             close(1)                            
!                         end if
                        
                        if ( x_i/k .gt. -.01 .and. x_i/k .lt. .01 ) then                        
                        i_k_inaction=i_k_inaction+omega_ss0(ieps,ik,ia)
                        end if
                        
                        if ( x_i/k .ge. .01 ) then
                        i_k_pos=i_k_pos+omega_ss0(ieps,ik,ia)
                        end if
                        
                        if ( x_i/k .le. -.01 ) then
                        i_k_neg=i_k_neg+omega_ss0(ieps,ik,ia)
                        end if
                        
                        if ( x_i/k .gt. .2 ) then
                        i_k_pos_spike=i_k_pos_spike+omega_ss0(ieps,ik,ia)
                        end if
                        
                        if ( x_i/k .lt. -.2 ) then
                        i_k_neg_spike=i_k_neg_spike+omega_ss0(ieps,ik,ia)
                        end if
                        
                    end if
                
                end if
                
            dummy=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
                
                if ( d_o .eq. 1 ) then
            
                sales0=(((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu))
            
                    do jeps=1,neps

                        if ( pi_eps(ieps,jeps) .gt. zero ) then
                    
                        lp=(w/(((Z*grid_eps(jeps))**(1-alpha*nu))*(1-alpha)*nu*(kp**(alpha*nu))))**(one/((1-alpha)*nu-1))         
                        sales1=(((Z*grid_eps(jeps))**(1-alpha*nu))*(((kp**alpha)*(lp**(1-alpha)))**nu))
                    
                            if ( ik .eq. intbirth_k .and. ia .eq. intbirth ) then                            
                            
                                if ( g_i_o(jeps,inds(ieps,ik,ia,1),ja) .eq. 1 ) then
                                std_sales0=std_sales0+((avg_sales0-sales0)**2)*&
                                vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*(omega_ss0(ieps,ik,ia)+put_eps(ieps)*mu*g_e_o(ieps))*(1-delta_f)
                            
                                std_sales1=std_sales1+((avg_sales1-sales1)**2)*&
                                vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*(omega_ss0(ieps,ik,ia)+put_eps(ieps)*mu*g_e_o(ieps))*(1-delta_f)
                            
                                auto_cov_sales=auto_cov_sales+(sales0-avg_sales0)*(sales1-avg_sales1)*&
                                vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*(omega_ss0(ieps,ik,ia)+put_eps(ieps)*mu*g_e_o(ieps))*(1-delta_f)
                                end if
                                
                                if ( g_i_o(jeps,inds(ieps,ik,ia,2),ja) .eq. 1 ) then
                                std_sales0=std_sales0+((avg_sales0-sales0)**2)*&
                                vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*(omega_ss0(ieps,ik,ia)+put_eps(ieps)*mu*g_e_o(ieps))*(1-delta_f)                                
                            
                                std_sales1=std_sales1+((avg_sales1-sales1)**2)*&
                                vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*(omega_ss0(ieps,ik,ia)+put_eps(ieps)*mu*g_e_o(ieps))*(1-delta_f)                                
                            
                                auto_cov_sales=auto_cov_sales+(sales0-avg_sales0)*(sales1-avg_sales1)*&
                                vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*(omega_ss0(ieps,ik,ia)+put_eps(ieps)*mu*g_e_o(ieps))*(1-delta_f)                                
                                end if
                            
                            else
                            
                                if ( g_i_o(jeps,inds(ieps,ik,ia,1),ja) .eq. 1 ) then
                                std_sales0=std_sales0+((avg_sales0-sales0)**2)*&
                                vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)
                            
                                std_sales1=std_sales1+((avg_sales1-sales1)**2)*&
                                vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)
                            
                                auto_cov_sales=auto_cov_sales+(sales0-avg_sales0)*(sales1-avg_sales1)*&
                                vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)
                                end if
                                
                                if ( g_i_o(jeps,inds(ieps,ik,ia,2),ja) .eq. 1 ) then
                                std_sales0=std_sales0+((avg_sales0-sales0)**2)*&
                                vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)                                
                            
                                std_sales1=std_sales1+((avg_sales1-sales1)**2)*&
                                vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)                                
                            
                                auto_cov_sales=auto_cov_sales+(sales0-avg_sales0)*(sales1-avg_sales1)*&
                                vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)                                
                                end if
                            
                            end if
                    
                        end if
                    
                    end do
                
                end if

            end if
            
        end do
        end do
        end do

!          if ( rank .eq. 0 ) then
!          print*,'avg_i_k',avg_i_k
!          print*,'std_dev_i_k',std_dev_i_k
!          print*,'i_k_inaction',i_k_inaction
!          print*,'i_k_pos',i_k_pos
!          print*,'i_k_neg',i_k_neg
!          print*,'i_k_pos_spike',i_k_pos_spike
!          print*,'i_k_neg_spike',i_k_neg_spike
!          end if
        
!     h=(1-psi*(r*b+div)/w)/(1+psi)
    
    !bisection for labor    
    d_bih=one
    iter_bih=0
    hmin=zero
    hmax=2.
    
        do while ( iter_bih .le. itmax_bih .and. abs(d_bih) .ge. tol_bih )
    
        h=(hmin+hmax)/2.

            if ( w*h+b*r/(1+r)+div+tau_a*neg_b+tot_inv_wedge-w*(1-tau_l)/(psi*(h**phi)) .le. zero ) then
            hmin=h
            else if ( w*h+b*r/(1+r)+div+tau_a*neg_b+tot_inv_wedge-w*(1-tau_l)/(psi*(h**phi)) .gt. zero ) then
            hmax=h
            end if
            
        d_bih=max(abs(hmax-hmin),abs(w*h+b*r/(1+r)+div+tau_a*neg_b+tot_inv_wedge-w*(1-tau_l)/(psi*(h**phi))))
        iter_bih=iter_bih+1
        
        end do
        
        if ( iter_bih .ge. itmax_bih ) then
        print*,'iterations exceeded labor bisection, instst'
            call exit
        end if

    c=w*h+b*r/(1+r)+div+tau_a*neg_b+tot_inv_wedge
    exdem=abs(ngdp-c-ninv)

        if ( ex_entry .eq. 0 ) then
        
            if ( lab_dem-h .le. zero ) then
            mumin=mu
            else if ( lab_dem-h .gt. zero ) then
            mumax=mu
            end if

        d_bi=max(abs(lab_dem-h),abs(mumax-mumin))

        else if ( ex_entry .eq. 1 ) then

            if ( lab_dem-h .gt. zero ) then
            wmin=w
            else if ( lab_dem-h .le. zero ) then
            wmax=w
            end if

        d_bi=abs(wmax-wmin)	!max(abs(lab_dem-h),abs(wmax-wmin))

        end if
        
    iter_bi=iter_bi+1
    
        if ( rank .eq. 0 ) then

            if ( ex_entry .eq. 0 ) then
    
            write(6,'(300F12.7)')real(iter_bi),mu,d_bi
            print*,'exdem, h',exdem,lab_dem,h
        
                open (unit=1,file='input_instst_mu.txt')
                write(1,*)mumin,mumax
                close(1)            
            
            else if ( ex_entry .eq. 1 ) then
            
            write(6,'(300F12.7)')real(iter_bi),w,d_bi
            print*,'exdem, h',exdem,lab_dem,h
            print*,'sum(g_e_o)',sum(g_e_o)
                
!                open (unit=1,file='input_instst_w.txt')
!                write(1,*)wmin,wmax
!                close(1)                    
                
            end if
            
!             open (unit=1,file='output_instst_dist.txt')
!                 do ik=1,nk
!                 do ia=1,na 
!             
!                     if ( ik .eq. intbirth_k .and. ia .eq. intbirth ) then
!                     write(1,'(300F12.7)')grid_k(ik),grid_a(date,ia),(sum(omega_ss0(:,ik,ia))+mu)/(sum(omega_ss0)+mu)
!                     else
!                     write(1,'(300F12.7)')grid_k(ik),grid_a(date,ia),sum(omega_ss0(:,ik,ia))/(sum(omega_ss0)+mu)
!                     end if
!                     
!                 end do
!                 end do
!             close(1)
            
        end if
    
    end do
    
    if ( rank .eq. 0 .and. abs(d_bi) .ge. tol_bi ) then
    print*,'iterations exceeded in labor market clearing, instst.f90 ',d_bi
        call exit
    end if  

    dist_i_k(:,2)=dist_i_k(:,2)/sum(dist_i_k(:,2))

    allocate(dist_i_k1(cnt_i_k,2))

    dist_i_k1=zero
    
        do cnt=1,cnt_i_k
        dist_i_k1(cnt,:)=dist_i_k(cnt,:)
        end do
        
!         call sort (dist_i_k1,cnt_i_k)

        do cnt=1,cnt_i_k

!             if ( sum(dist_i_k1(1:cnt,2)) .le. .95 ) then
            if ( abs(dist_i_k1(cnt,1)) .le. 3*std_dev_i_k_data ) then
            avg_i_k=avg_i_k+dist_i_k1(cnt,1)*dist_i_k1(cnt,2)
            std_dev_i_k=std_dev_i_k+(dist_i_k1(cnt,1)**2)*dist_i_k1(cnt,2)
            temp_k=temp_k+dist_i_k1(cnt,2)
            end if
            
        end do
        
    deallocate(dist_i_k,dist_i_k1)
        
    avg_i_k=avg_i_k/temp_k
    std_dev_i_k=std_dev_i_k/temp_k
    std_dev_i_k=sqrt(std_dev_i_k-avg_i_k**2)
    i_k_inaction=i_k_inaction/temp
    i_k_pos=i_k_pos/temp
    i_k_neg=i_k_neg/temp
    i_k_pos_spike=i_k_pos_spike/temp
    i_k_neg_spike=i_k_neg_spike/temp
    
avg_firm_size=h/(sum(omega_ss0)+mu*sum(put_eps*g_e_o))
labor_share=w*h/ngdp
avg_sales=avg_sales/(sum(omega_ss0)+mu*sum(put_eps*g_e_o))    
var_sales=var_sales/(sum(omega_ss0)+mu*sum(put_eps*g_e_o))    
var_sales=var_sales-avg_sales**2
std_dev_sales=sqrt(var_sales)

std_sales0=std_sales0/surv
std_sales1=std_sales1/surv
auto_cov_sales=auto_cov_sales/surv

std_sales0=sqrt(std_sales0)
std_sales1=sqrt(std_sales1)

auto_corr_sales=auto_cov_sales/(std_sales0*std_sales1)

firm_debt=zero

    do ia=1,intbirth
    firm_debt=firm_debt+grid_a(date,ia)*sum(omega_ss0(:,:,ia))
    end do

firm_capital=mu*k_e*sum(put_eps*g_e_o)
    
    do ik=1,nk
    firm_capital=firm_capital+grid_k(ik)*sum(omega_ss0(:,ik,:))
    end do
    
    if ( rank .eq. 0 ) then
    
        open (unit=1,file='output_instst_dist.txt')
        write(1,*)'k a mass'
            do ik=1,nk
            do ia=1,na 
            
                if ( ik .eq. intbirth_k .and. ia .eq. intbirth ) then
                write(1,'(300F12.7)')grid_k(ik),grid_a(date,ia),(sum(omega_ss0(:,ik,ia))+mu*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth)))&
                /(sum(omega_ss0)+mu*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth)))
                else
                write(1,'(300F12.7)')grid_k(ik),grid_a(date,ia),sum(omega_ss0(:,ik,ia))&
                /(sum(omega_ss0)+mu*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth)))
                end if
                
            end do
            end do
        close(1)
    
        if ( sum(omega_ss0(:,nk,:)) .gt. zero ) then
        print*,'increase kmax, instst.f90 ',sum(omega_ss0(:,nk,:))/(sum(omega_ss0)+mu*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth)))
        end if
    
        if ( sum(omega_ss0(:,:,1)) .gt. zero ) then
        print*,'decrease amin, instst.f90 ',sum(omega_ss0(:,:,1))/(sum(omega_ss0)+mu*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth)))
        end if
        
        open (unit=1,file='output_instst_emp_dist.txt')
        write(1,*)'size mass total_employment'
        
        ik=intbirth_k
        k=grid_k(ik)
        
            do ieps=1,neps            
            iz=id_z(ieps)
            eps=grid_eps(ieps)
            d_o=g_i_o(ieps,ik,intbirth)
            
                if ( g_e_o(ieps) .eq. 1 .and. d_o .eq. 1 ) then
                l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(iz)      
                write(1,'(300F12.7)')l,put_eps(ieps)*mu*g_e_o(ieps),l*put_eps(ieps)*mu*g_e_o(ieps)                                                                
                end if
                
            end do
        
            do ieps=1,neps            
            do ik=1,nk
            do ia=1,na
            
                if ( omega_ss0(ieps,ik,ia) .gt. zero ) then

                iz=id_z(ieps)
                k=grid_k(ik)
                eps=grid_eps(ieps)
                d_o=g_i_o(ieps,ik,ia)
            
                    if ( d_o .eq. 1 ) then
                    l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(iz)   
                    write(1,'(300F12.7)')l,omega_ss0(ieps,ik,ia),l*omega_ss0(ieps,ik,ia)                                                                                        
                    end if
                
                end if
                
            end do
            end do
            end do            
            
        close(1)
        
    print*,'firm debt/gdp',-firm_debt/gdp,neg_b/gdp
    print*,'outliers',1-temp_k
    
    end if

r_t(date)=r
w_t(date)=w
h_t(date)=h
    c_t(date)=c
exdem_t(date)=exdem
b_t(date)=b    
div_t(date)=div
mu_t(date)=mu
mass_t(date)=sum(omega_ss0)

gdp_t(date)=gdp
ngdp_t(date)=ngdp
inv_t(date)=inv
ninv_t(date)=ninv
firms_constrained_t(date)=firms_constrained
exits_t(date)=exits
exits_binding_t(date)=exits_binding
avg_firm_size_t(date)=avg_firm_size
labor_share_t(date)=labor_share
std_dev_growth_emp_t(date)=std_dev_growth_emp
std_dev_growth_sales_t(date)=std_dev_growth_sales
firm_debt_t(date)=firm_debt
firm_capital_t(date)=firm_capital

omega_ss_initial=omega_ss0

g_i_o_initial=g_i_o
g_i_k_initial=g_i_k
g_i_a_initial=g_i_a
g_e_o_initial=g_e_o
g_q_initial=g_q
g_un_k_initial=g_un_k
g_un_o_initial=g_un_o

a_un_t(date,:,:)=a_un
a_tilde_t(date,:,:)=a_tilde

!compute employment by age
emp_by_age=zero
firms_by_age=zero
exit_by_age=zero
exit_binding_by_age=zero

omega_ss0=zero
omega_ss0(:,intbirth_k,intbirth)=mu*put_eps*g_e_o

    
    do iter_dist=1,119       !age0, age1, age2, age,3, age4, age5

    firms_by_age(iter_dist)=sum(omega_ss0)
    
        if ( rank .eq. 0 ) then    
        omega_ss_age_initial(iter_dist,:,:,:)=omega_ss0
        end if
        
    omega_ss1=zero
        
        do ieps=1,neps
        do ik=1,nk
        do ia=1,na
    
            if ( omega_ss0(ieps,ik,ia) .gt. zero )  then
        
            iz=id_z(ieps)
            
            eps=grid_eps(ieps)
            k=grid_k(ik)
            a=grid_a(date,ia)
        
            d_o=g_i_o(ieps,ik,ia)
            jk=g_i_k(ieps,ik,ia)
            ap=g_i_a(ieps,ik,ia)
        
            kp=grid_k(jk)
        
                if ( d_o .eq. 1 ) then                    
                l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))                
                else if ( d_o .eq. 0 ) then                        
                l=zero     
                end if
                
                if ( iter_dist .eq. 1 ) then
                l=l+f_e
                end if                

                if ( cap_adj .eq. 2 ) then
                l=l+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)/w
                end if
                
            emp_by_age(iter_dist)=emp_by_age(iter_dist)+(l+d_o*f_o*grid_z(iz))*omega_ss0(ieps,ik,ia)
            exit_by_age(iter_dist)=exit_by_age(iter_dist)+(1-d_o+d_o*delta_f)*omega_ss0(ieps,ik,ia)
                
                if ( d_o .eq. 1 ) then
            
                    do jeps=1,neps       
                    
                        if ( pi_eps(ieps,jeps) .gt. zero ) then
                        omega_ss1(jeps,jk,inds(ieps,ik,ia,1))=omega_ss1(jeps,jk,inds(ieps,ik,ia,1))+&
                        vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)
                        omega_ss1(jeps,jk,inds(ieps,ik,ia,2))=omega_ss1(jeps,jk,inds(ieps,ik,ia,2))+&
                        vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)                                                    
                        endif
                    
                    end do
                
                end if
            
            end if
        
        end do
        end do
        end do

    omega_ss0=omega_ss1
        
    end do

firms_by_age(120)=mass_t(date)+mu_t(date)*sum(put_eps*g_e_o)-sum(firms_by_age(1:119))  
emp_by_age(120)=h_t(date)-sum(emp_by_age(1:119))
exit_by_age(120)=exits_t(date)-sum(exit_by_age(1:119))
exit_binding_by_age(120)=exits_binding_t(date)-sum(exit_binding_by_age(1:119))

firms_by_age_t(date,:)=firms_by_age(1:120) 
emp_by_age_t(date,:)=emp_by_age(1:120)
exit_by_age_t(date,:)=exit_by_age(1:120)
exit_binding_by_age_t(date,:)=exit_binding_by_age(1:120)

    if ( rank .eq. 0 ) then
    
        do ieps=1,neps
        do ik=1,nk
        do ia=1,na
    
            if ( ia .eq. intbirth .and. ik .eq. intbirth_k ) then
            omega_ss_age_initial(120,ieps,ik,ia)=omega_ss_initial(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps)-sum(omega_ss_age_initial(1:119,ieps,ik,ia))        
            else        
            omega_ss_age_initial(120,ieps,ik,ia)=omega_ss_initial(ieps,ik,ia)-sum(omega_ss_age_initial(1:119,ieps,ik,ia))
            end if
                    
        end do
        end do
        end do

    end if
    
!compute job creation and job destruction
omega_ss0=omega_ss_initial 
job_creation=zero
job_destruction=zero

    do ieps=1,neps
    do ik=1,nk
    do ia=1,na    
    
        if ( omega_ss0(ieps,ik,ia) .gt. zero )  then
        
        iz=id_z(ieps)
        eps=grid_eps(ieps)
        k=grid_k(ik)
        a=grid_a(date,ia)
        
        d_o=g_i_o(ieps,ik,ia)
        jk=g_i_k(ieps,ik,ia)
        ap=g_i_a(ieps,ik,ia)
        
        kp=grid_k(jk)
        
            if ( d_o .eq. 1 ) then                    
            l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(iz)                
            else if ( d_o .eq. 0 ) then                        
            l=zero                                        
            end if

            if ( cap_adj .eq. 2 ) then
            l=l+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)/w
            end if
            
            if ( d_o .eq. 1 ) then
            
                do jeps=1,neps
                jz=id_z(jeps)
                    if ( pi_eps(ieps,jeps) .gt. zero ) then
                
                    ja=inds(ieps,ik,ia,1)

                        if ( g_i_o(jeps,jk,ja) .eq. 1 ) then                    
                        lp=(w/(((Z*grid_eps(jeps))**(1-alpha*nu))*(1-alpha)*nu*(kp**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(jz)         
                        else if ( g_i_o(jeps,jk,ja) .eq. 0 ) then                        
                        lp=zero                                        
                        end if

                        if ( cap_adj .eq. 2 ) then
                        lp=lp+(func_adj(grid_k(g_i_k(jeps,jk,ja)),kp,g_i_o(jeps,jk,ja),date)-lambda*(grid_k(g_i_k(jeps,jk,ja))/kp-1)**2.)/w                                            
                        end if                        
                        
                        if ( lp-l .ge. zero ) then
                        job_creation=job_creation+&
                        (lp-l)*vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)
                        else if ( lp-l .lt. zero ) then
                        job_destruction=job_destruction+&
                        (l-lp)*vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)
                        end if
                    
                    ja=inds(ieps,ik,ia,2)

                        if ( g_i_o(jeps,jk,ja) .eq. 1 ) then                    
                        lp=(w/(((Z*grid_eps(jeps))**(1-alpha*nu))*(1-alpha)*nu*(kp**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(jz)                                   
                        else if ( g_i_o(jeps,jk,ja) .eq. 0 ) then                        
                        lp=zero                                        
                        end if
                        
                        if ( cap_adj .eq. 2 ) then
                        lp=lp+(func_adj(grid_k(g_i_k(jeps,jk,ja)),kp,g_i_o(jeps,jk,ja),date)-lambda*(grid_k(g_i_k(jeps,jk,ja))/kp-1)**2.)/w                                            
                        end if                                                
                
                        if ( lp-l .ge. zero ) then
                        job_creation=job_creation+&
                        (lp-l)*vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f) 
                        else if ( lp-l .lt. zero ) then
                        job_destruction=job_destruction+&
                        (l-lp)*vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f) 
                        end if
                
                    end if
                
                end do
                
            job_destruction=job_destruction+(l-zero)*omega_ss0(ieps,ik,ia)*delta_f
                
            else if ( d_o .eq. 0 ) then
            
            lp=zero
            job_destruction=job_destruction+(l-lp)*omega_ss0(ieps,ik,ia)
            
            end if                
                
        end if
        
        if ( ik .eq. intbirth_k .and. ia .eq. intbirth ) then
        
        iz=id_z(ieps)
        
        eps=grid_eps(ieps)
        k=grid_k(ik)
        a=grid_a(date,ia)
        
        d_o=g_i_o(ieps,ik,ia)
        jk=g_i_k(ieps,ik,ia)
        ap=g_i_a(ieps,ik,ia)
        
        kp=grid_k(jk)
        
            if ( d_o .eq. 1 ) then                    
            l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(iz)+f_e
            else if ( d_o .eq. 0 ) then                        
            l=f_e                                        
            end if
            
            if ( cap_adj .eq. 2 ) then
            l=l+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)/w
            end if                        
            
            if ( d_o .eq. 1 ) then
            
                do jeps=1,neps                    

                jz=id_z(jeps)
                
                    if ( pi_eps(ieps,jeps) .gt. zero ) then
                
                    ja=inds(ieps,ik,ia,1)

                        if ( g_i_o(jeps,jk,ja) .eq. 1 ) then                    
                        lp=(w/(((Z*grid_eps(jeps))**(1-alpha*nu))*(1-alpha)*nu*(kp**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(jz)                                   
                        else if ( g_i_o(jeps,jk,ja) .eq. 0 ) then                        
                        lp=zero                                        
                        end if
                        
                        if ( cap_adj .eq. 2 ) then
                        lp=lp+(func_adj(grid_k(g_i_k(jeps,jk,ja)),kp,g_i_o(jeps,jk,ja),date)-lambda*(grid_k(g_i_k(jeps,jk,ja))/kp-1)**2.)/w                                                                    
                        end if                        
                    
                        if ( lp-l .ge. zero ) then
                        job_creation=job_creation+&
                        (lp-l)*vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*mu*put_eps(ieps)*g_e_o(ieps)*(1-delta_f)
                        else if ( lp-l .lt. zero ) then
                        job_destruction=job_destruction+&
                        (l-lp)*vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*mu*put_eps(ieps)*g_e_o(ieps)*(1-delta_f)
                        end if
                    
                    ja=inds(ieps,ik,ia,2)

                        if ( g_i_o(jeps,jk,ja) .eq. 1 ) then                    
                        lp=(w/(((Z*grid_eps(jeps))**(1-alpha*nu))*(1-alpha)*nu*(kp**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(jz)                                  
                        else if ( g_i_o(jeps,jk,ja) .eq. 0 ) then                        
                        lp=zero                                        
                        end if

                        if ( cap_adj .eq. 2 ) then
                        lp=lp+(func_adj(grid_k(g_i_k(jeps,jk,ja)),kp,g_i_o(jeps,jk,ja),date)-lambda*(grid_k(g_i_k(jeps,jk,ja))/kp-1)**2.)/w                                            
                        end if                        
                        
                        if ( lp-l .ge. zero ) then
                        job_creation=job_creation+&
                        (lp-l)*vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*mu*put_eps(ieps)*g_e_o(ieps)*(1-delta_f)
                        else if ( lp-l .lt. zero ) then
                        job_destruction=job_destruction+&
                        (l-lp)*vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*mu*put_eps(ieps)*g_e_o(ieps)*(1-delta_f)
                        end if
                
                    end if
                
                end do
                
            job_destruction=job_destruction+(l-zero)*mu*put_eps(ieps)*g_e_o(ieps)*delta_f                    
            
            else if ( d_o .eq. 0 ) then
            
            lp=zero
            job_destruction=job_destruction+(l-lp)*mu*put_eps(ieps)*g_e_o(ieps)
                
            end if 
                        
        job_creation=job_creation+l*mu*put_eps(ieps)*g_e_o(ieps)
        
        end if
                    
    end do
    end do
    end do
    
job_creation_t(date)=job_creation
job_destruction_t(date)=job_destruction

!life cycle properties
emp_by_age=zero
firms_by_age=zero
exit_by_age=zero
exit_binding_by_age=zero
eps_by_age=zero
k_by_age=zero
a_by_age=zero
job_creation_by_age=zero
job_destruction_by_age=zero
div_by_age=zero
div_operating_by_age=zero

avg_log_emp_by_age=zero
std_dev_log_emp_by_age=zero
dummy_by_age=zero

omega_ss0=zero
omega_ss0(:,intbirth_k,intbirth)=mu*put_eps*g_e_o

    do iter_dist=1,max_age      !age0, age1, age2, age,3, age4, age5

        if ( rank .eq. 0 ) then
!        print*,'iter_dist',iter_dist
        end if

    firms_by_age(iter_dist)=sum(omega_ss0)

    omega_ss1=zero
        
        do ieps=1,neps
        do ik=1,nk
        do ia=1,na
    
            if ( omega_ss0(ieps,ik,ia) .gt. zero )  then
        
            iz=id_z(ieps)
            
            eps=grid_eps(ieps)
            k=grid_k(ik)
            a=grid_a(date,ia)        
        
            d_o=g_i_o(ieps,ik,ia)
            jk=g_i_k(ieps,ik,ia)
            ap=g_i_a(ieps,ik,ia)
        
            kp=grid_k(jk)
        
                if ( d_o .eq. 1 ) then                    
                
                l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
                
                    if ( ap .ge. zero ) then
                                
                    q=one/(one+r)
                                
                    else if ( ap .lt. zero ) then
                                
                        if ( default_economy .eq. 0 ) then                                
                        q=one/(one+r)-tau_a
                        else if ( default_economy .eq. 1 ) then
                            call onedlin (grid_a(date,:),g_q(ieps,jk,:),na,ap,q) 
                        end if
                            
                    end if
                                
                profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))&
                +(1.+tau_x)*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,date)            
                    
                    if ( rank .eq. 0 .and. isnan(l) ) then
                    print*,'isnan l, instst.f90',k
                        call exit
                    end if

                div_operating_by_age(iter_dist)=div_operating_by_age(iter_dist)+profit*omega_ss0(ieps,ik,ia)

                else if ( d_o .eq. 0 ) then                        

                l=zero                                          

                    if ( default_economy .eq. 0 ) then                                            
                    profit=(1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date)                    
                    else if ( default_economy .eq. 1 ) then
                    profit=max((1.+tau_x)*k*chi+a-func_adj(kp,k,d_o,date),min(zero,(1.+tau_x)*k*chi-func_adj(kp,k,d_o,date)))                                    
                    end if            
                    
                end if
                
            k_by_age(iter_dist)=k_by_age(iter_dist)+k*omega_ss0(ieps,ik,ia)
            a_by_age(iter_dist)=a_by_age(iter_dist)+a*omega_ss0(ieps,ik,ia)
            div_by_age(iter_dist)=div_by_age(iter_dist)+profit*omega_ss0(ieps,ik,ia)
        
                if ( d_o .eq. 1 ) then                    
                l=l+f_o*grid_z(iz)
                eps_by_age(iter_dist)=eps_by_age(iter_dist)+eps*omega_ss0(ieps,ik,ia)
                else if ( d_o .eq. 0 ) then                        
                l=zero     
                end if

                if ( cap_adj .eq. 2 ) then
                l=l+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)/w
                end if                        

                if ( iter_dist .eq. 1 .or. iter_dist .eq. 2 .or. iter_dist .eq. 3 .or. iter_dist .eq. 4 ) then
                job_creation_by_age(iter_dist)=job_creation_by_age(iter_dist)+l*omega_ss0(ieps,ik,ia)         
                end if
                
                if ( iter_dist .eq. 1 ) then                
!                 l=l+f_e                
                emp_by_age(iter_dist)=emp_by_age(iter_dist)+(l+f_e)*omega_ss0(ieps,ik,ia)
                else
                emp_by_age(iter_dist)=emp_by_age(iter_dist)+(l)*omega_ss0(ieps,ik,ia)                
                end if

                if ( d_o .eq. 1 ) then
                avg_log_emp_by_age(iter_dist)=avg_log_emp_by_age(iter_dist)+log(l)*omega_ss0(ieps,ik,ia)
                std_dev_log_emp_by_age(iter_dist)=std_dev_log_emp_by_age(iter_dist)+(log(l)**2)*omega_ss0(ieps,ik,ia)
                dummy_by_age(iter_dist)=dummy_by_age(iter_dist)+omega_ss0(ieps,ik,ia)
                end if
                
            exit_by_age(iter_dist)=exit_by_age(iter_dist)+(1-d_o+d_o*delta_f)*omega_ss0(ieps,ik,ia)
                       
                if ( d_o .eq. 1 ) then

                    !next period distribution                
                    do jeps=1,neps       
                    
                    jz=id_z(jeps)
                    
                        if ( pi_eps(ieps,jeps) .gt. zero ) then                        
                
                        ja=inds(ieps,ik,ia,1)
                        omega_ss1(jeps,jk,ja)=omega_ss1(jeps,jk,ja)+vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)      

                        ja=inds(ieps,ik,ia,2)                        
                        omega_ss1(jeps,jk,ja)=omega_ss1(jeps,jk,ja)+vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)                                   
                
                        end if
                    
                    end do
                
                end if
            
            end if
        
        end do
        end do
        end do

        if ( compute_job_flows .eq. 1 .and. iter_dist+4 .le. max_age ) then
            
		if ( rank .eq. 0 ) then
		print*,'iter_dist',iter_dist
		end if

        allocate (job_creation_mpi(np),job_destruction_mpi(np))        
        
        job_creation_mpi=zero
        job_destruction_mpi=zero
        
            do cnt=itop,iend
    
            ieps=zmap(cnt,1)
            ik=zmap(cnt,2)
            ia=zmap(cnt,3)
            iz=id_z(ieps)
        
                if ( omega_ss0(ieps,ik,ia) .gt. zero )  then
        
                iz=id_z(ieps)
            
                eps=grid_eps(ieps)
                k=grid_k(ik)
                a=grid_a(date,ia)        
        
                d_o=g_i_o(ieps,ik,ia)
                jk=g_i_k(ieps,ik,ia)
                ap=g_i_a(ieps,ik,ia)
        
                kp=grid_k(jk)
        
                    if ( d_o .eq. 1 ) then                    
                
                    l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(iz)
                
                    !simulate 4 times to compute job creation and job destruction
                    omega_ss2=zero
                    omega_ss2(ieps,ik,ia)=omega_ss0(ieps,ik,ia)
                
                        do iter_simulate=1,4

                        omega_ss3=zero
                    
                            do ieps0=1,neps
                            do ik0=1,nk
                            do ia0=1,na
                        
                                if ( omega_ss2(ieps0,ik0,ia0) .gt. zero ) then
                            
                                jk=g_i_k(ieps0,ik0,ia0)
                        
                                    if ( g_i_o(ieps0,ik0,ia0) .eq. 1 ) then
                            
                                        do jeps=1,neps       
                    
                                        jz=id_z(jeps)
                    
                                            if ( pi_eps(ieps,jeps) .gt. zero ) then                        
                        
                                            ja=inds(ieps0,ik0,ia0,1)                            
                                            omega_ss3(jeps,jk,ja)=omega_ss3(jeps,jk,ja)+vals(ieps0,ik0,ia0,1)*pi_eps(ieps,jeps)*omega_ss2(ieps0,ik0,ia0)*(1-delta_f)      

                                            ja=inds(ieps0,ik0,ia0,2)                        
                                            omega_ss3(jeps,jk,ja)=omega_ss3(jeps,jk,ja)+vals(ieps0,ik0,ia0,2)*pi_eps(ieps,jeps)*omega_ss2(ieps0,ik0,ia0)*(1-delta_f)                                                    
                        
                                            end if
                    
                                        end do            
                                                        
                                    else if ( g_i_o(ieps0,ik0,ia0) .eq. 0 ) then
                                
                                    lp=zero
                                
                                        if ( cap_adj .eq. 2 ) then
                                        print*,'this needs to be fixed, instst.f90'
                                            call exit                                
                                        end if                                  
                                
                                        if ( lp-l .ge. zero ) then
                                        job_creation_mpi(cnt-itop+1)=job_creation_mpi(cnt-itop+1)+(lp-l)*omega_ss2(ieps0,ik0,ia0)
                                        else if ( lp-l .lt. zero ) then
                                        job_destruction_mpi(cnt-itop+1)=job_destruction_mpi(cnt-itop+1)+(l-lp)*omega_ss2(ieps0,ik0,ia0)
                                        end if
                                
                                    end if
                            
                                end if
                            
                            end do
                            end do
                            end do
                    
                        omega_ss2=omega_ss3
                    
                        end do

                        do jeps=1,neps
                        do jk=1,nk
                        do ja=1,na
                    
                            if ( omega_ss2(jeps,jk,ja) .gt. zero ) then

                            jz=id_z(jeps)
			    kp=grid_k(jk)
                        
                                if ( g_i_o(jeps,jk,ja) .eq. 1 ) then                    
                                lp=(w/(((Z*grid_eps(jeps))**(1-alpha*nu))*(1-alpha)*nu*(kp**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(jz)         
                                else if ( g_i_o(jeps,jk,ja) .eq. 0 ) then                        
                                lp=zero                                        
                                end if

                                if ( cap_adj .eq. 2 ) then
                                lp=lp+(func_adj(grid_k(g_i_k(jeps,jk,ja)),kp,g_i_o(jeps,jk,ja),date)-lambda*(grid_k(g_i_k(jeps,jk,ja))/kp-1)**2.)/w                                  
                                end if   
                            
                                if ( lp-l .ge. zero ) then
                                job_creation_mpi(cnt-itop+1)=job_creation_mpi(cnt-itop+1)+(lp-l)*omega_ss2(jeps,jk,ja)
                                else if ( lp-l .lt. zero ) then
                                job_destruction_mpi(cnt-itop+1)=job_destruction_mpi(cnt-itop+1)+(l-lp)*omega_ss2(jeps,jk,ja)
                                end if
                                                        
                            end if
                        
                        end do
                        end do
                        end do                    

                    end if
        
                end if
            
            end do
        
            allocate (job_creation_recmpi(np*nproc),job_destruction_recmpi(np*nproc))
        
        job_creation_recmpi=zero
        job_destruction_recmpi=zero
        
            call mpi_allgather (job_creation_mpi,np,mpi_double_precision,job_creation_recmpi,np,mpi_double_precision,mpi_comm_world,ierr)
            call mpi_allgather (job_destruction_mpi,np,mpi_double_precision,job_destruction_recmpi,np,mpi_double_precision,mpi_comm_world,ierr)
              
        job_creation_by_age(iter_dist+4)=sum(job_creation_recmpi)
        job_destruction_by_age(iter_dist+4)=sum(job_destruction_recmpi)

            deallocate (job_creation_mpi,job_destruction_mpi)        
            deallocate (job_creation_recmpi,job_destruction_recmpi)
        
        end if
        
    avg_log_emp_by_age(iter_dist)=avg_log_emp_by_age(iter_dist)/dummy_by_age(iter_dist)
    std_dev_log_emp_by_age(iter_dist)=std_dev_log_emp_by_age(iter_dist)/dummy_by_age(iter_dist)
    std_dev_log_emp_by_age(iter_dist)=sqrt(std_dev_log_emp_by_age(iter_dist)-avg_log_emp_by_age(iter_dist)**2.)
        
    omega_ss0=omega_ss1
        
    end do
                        
    if ( rank .eq. 0 ) then
    
        open (unit=1,file='output_life_cycle.txt')
        write(1,*)'age firms emp exit eps k a exit_bind job_creat job_destr div div_ope' 
        write(1,*)'avg_log_emp std_dev_log_emp'
            do iter_dist=1,max_age
            write(1,'(300F12.7)')real(iter_dist),firms_by_age(iter_dist),emp_by_age(iter_dist),exit_by_age(iter_dist),&
            eps_by_age(iter_dist),k_by_age(iter_dist),a_by_age(iter_dist),exit_binding_by_age(iter_dist),&
            job_creation_by_age(iter_dist),job_destruction_by_age(iter_dist),div_by_age(iter_dist),div_operating_by_age(iter_dist),&
            avg_log_emp_by_age(iter_dist),std_dev_log_emp_by_age(iter_dist)
            end do
        close(1)
        
    end if

!investment moments after 30 years    
omega_ss0=zero
omega_ss0(:,intbirth_k,intbirth)=mu*put_eps*g_e_o
omega_ss0=omega_ss0+omega_ss_initial    

    do iter_dist=1,119
    
    omega_ss1=zero

        do ieps=1,neps
        do ik=1,nk
        do ia=1,na

            if ( omega_ss0(ieps,ik,ia) .gt. zero )  then
            
            d_o=g_i_o(ieps,ik,ia)                                
            jk=g_i_k(ieps,ik,ia)

                if ( d_o .eq. 1 ) then
            
                    do jeps=1,neps

                        if ( pi_eps(ieps,jeps) .gt. zero ) then                                
                        omega_ss1(jeps,jk,inds(ieps,ik,ia,1))=omega_ss1(jeps,jk,inds(ieps,ik,ia,1))+&
                        vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)
                        omega_ss1(jeps,jk,inds(ieps,ik,ia,2))=omega_ss1(jeps,jk,inds(ieps,ik,ia,2))+&
                        vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*omega_ss0(ieps,ik,ia)*(1-delta_f)                                                        
                        end if
                        
                    end do
                
                end if
            
            end if
                
        end do
        end do
        end do

    omega_ss0=omega_ss1

    end do

    if ( rank .eq. 0 ) then
    
    allocate (omega_temp0(nk,neps,nk,na),omega_temp1(nk,neps,nk,na))

    omega_temp0=zero

    do ieps=1,neps
    do ik=1,nk
    do ia=1,na        
    omega_temp0(ik,ieps,ik,ia)=omega_ss0(ieps,ik,ia)            
    end do
    end do
    end do
    
    do iter_dist=1,4
    
    omega_temp1=zero

        do ik0=1,nk
        do ieps=1,neps
        do ik=1,nk
        do ia=1,na

            if ( omega_temp0(ik0,ieps,ik,ia) .gt. zero )  then
            
            d_o=g_i_o(ieps,ik,ia)                                
            jk=g_i_k(ieps,ik,ia)

                if ( d_o .eq. 1 ) then
            
                    do jeps=1,neps

                        if ( pi_eps(ieps,jeps) .gt. zero ) then                                
                        omega_temp1(ik0,jeps,jk,inds(ieps,ik,ia,1))=omega_temp1(ik0,jeps,jk,inds(ieps,ik,ia,1))+&
                        vals(ieps,ik,ia,1)*pi_eps(ieps,jeps)*omega_temp0(ik0,ieps,ik,ia)*(1-delta_f)
                        omega_temp1(ik0,jeps,jk,inds(ieps,ik,ia,2))=omega_temp1(ik0,jeps,jk,inds(ieps,ik,ia,2))+&
                        vals(ieps,ik,ia,2)*pi_eps(ieps,jeps)*omega_temp0(ik0,ieps,ik,ia)*(1-delta_f)                                                        
                        end if
                        
                    end do
                
                end if
            
            end if
                
        end do
        end do
        end do
        end do
        
    omega_temp0=omega_temp1
    
    end do    
            
avg_i_k30=zero
std_dev_i_k30=zero     
i_k_inaction30=zero
i_k_pos30=zero
i_k_neg30=zero
i_k_pos_spike30=zero
i_k_neg_spike30=zero    
i_k_pos_high_spike30=zero

temp=zero
temp_k=zero
    
avg_i_k30check=zero
std_dev_i_k30check=zero     

    do ik0=1,nk
    do ieps=1,neps
    do ik=1,nk
    do ia=1,na

        if ( omega_temp0(ik0,ieps,ik,ia) .gt. zero )  then
        
!         eps=grid_eps(ieps)
!         k=grid_k(ik)
!         a=grid_a(date,ia)
!         
!         d_o=g_i_o(ieps,ik,ia)
!         jk=g_i_k(ieps,ik,ia)
!         ap=g_i_a(ieps,ik,ia)
!         
!         kp=grid_k(jk)

        d_o=1
        k=grid_k(ik0)
        kp=grid_k(ik)
        
            if ( d_o .eq. 1 .and. k .gt. zero ) then
                    
!            x_i=kp-(1-delta_k)*k*d_o-chi*k*(1-d_o)
            x_i=kp-((1-delta_k)**4)*k*d_o
            temp=temp+omega_temp0(ik0,ieps,ik,ia)
                    
                if ( abs(x_i) .le. 3*std_dev_i_k_data ) then
                temp_k=temp_k+omega_temp0(ik0,ieps,ik,ia)
                avg_i_k30=avg_i_k30+(x_i/k)*omega_temp0(ik0,ieps,ik,ia)
                std_dev_i_k30=std_dev_i_k30+((x_i/k)**2.)*omega_temp0(ik0,ieps,ik,ia)
                end if

            avg_i_k30check=avg_i_k30check+(x_i/k)*omega_temp0(ik0,ieps,ik,ia)
            std_dev_i_k30check=std_dev_i_k30check+((x_i/k)**2.)*omega_temp0(ik0,ieps,ik,ia)
                                            
                if ( x_i/k .gt. -.01 .and. x_i/k .lt. .01 ) then                        
                i_k_inaction30=i_k_inaction30+omega_temp0(ik0,ieps,ik,ia)
                end if
                        
                if ( x_i/k .ge. .01 ) then
                i_k_pos30=i_k_pos30+omega_temp0(ik0,ieps,ik,ia)
                end if
                    
                if ( x_i/k .le. -.01 ) then
                i_k_neg30=i_k_neg30+omega_temp0(ik0,ieps,ik,ia)
                end if
                        
                if ( x_i/k .gt. .2 ) then
                i_k_pos_spike30=i_k_pos_spike30+omega_temp0(ik0,ieps,ik,ia)
                end if
                        
                if ( x_i/k .lt. -.2 ) then
                i_k_neg_spike30=i_k_neg_spike30+omega_temp0(ik0,ieps,ik,ia)
                end if
            
                if ( x_i/k .gt. .8 ) then
                i_k_pos_high_spike30=i_k_pos_high_spike30+omega_temp0(ik0,ieps,ik,ia)
                end if

                if ( rank .eq. 0 ) then                        
!                	open(unit=1,file='output_instst_i_k_dist.txt',action='write',position='append')        
!                        write(1,'(300F14.9)')k,kp,x_i/k,omega_ss0(ieps,ik,ia)
!                        close(1)                            
                end if
            
            end if

        end if
            
    end do
    end do
    end do
    end do
    
    avg_i_k30=avg_i_k30/temp_k
    std_dev_i_k30=std_dev_i_k30/temp_k
    std_dev_i_k30=sqrt(std_dev_i_k30-avg_i_k30**2)
    i_k_inaction30=i_k_inaction30/temp
    i_k_pos30=i_k_pos30/temp
    i_k_neg30=i_k_neg30/temp
    i_k_pos_spike30=i_k_pos_spike30/temp
    i_k_neg_spike30=i_k_neg_spike30/temp

    i_k_pos_high_spike30=i_k_pos_high_spike30/temp

    avg_i_k30check=avg_i_k30check/temp
    std_dev_i_k30check=std_dev_i_k30check/temp
    std_dev_i_k30check=sqrt(std_dev_i_k30check-avg_i_k30check**2)

    deallocate (omega_temp0,omega_temp1)
    
    print*,'outliers_30',(temp-temp_k)/temp
    ! print*,'firms 30+',sum(omega_ss0)/(sum(omega_ss_initial)+mu_t(1)*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth))) 
    end if
            
    call MPI_BCAST (avg_i_k30,1,mpi_double_precision,0,mpi_comm_world,ierr)             
    call MPI_BCAST (std_dev_i_k30,1,mpi_double_precision,0,mpi_comm_world,ierr)             

    call MPI_BCAST (i_k_inaction30,1,mpi_double_precision,0,mpi_comm_world,ierr)             
    call MPI_BCAST (i_k_pos30,1,mpi_double_precision,0,mpi_comm_world,ierr)             
    call MPI_BCAST (i_k_neg30,1,mpi_double_precision,0,mpi_comm_world,ierr)             
    call MPI_BCAST (i_k_pos_spike30,1,mpi_double_precision,0,mpi_comm_world,ierr)             
    call MPI_BCAST (i_k_neg_spike30,1,mpi_double_precision,0,mpi_comm_world,ierr)             
    call MPI_BCAST (i_k_pos_high_spike30,1,mpi_double_precision,0,mpi_comm_world,ierr)             

    call MPI_BCAST (avg_i_k30check,1,mpi_double_precision,0,mpi_comm_world,ierr)             
    call MPI_BCAST (std_dev_i_k30check,1,mpi_double_precision,0,mpi_comm_world,ierr)             
    
end subroutine instst
