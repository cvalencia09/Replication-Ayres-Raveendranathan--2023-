subroutine fistst
use params
use mini
implicit none

integer :: iter_bi
double precision :: d_bi

double precision :: wmin, wmax
double precision :: mumin, mumax

double precision, dimension(neps,nk,na) :: omega_ss0, omega_ss1, omega_ss_final
integer :: iter_dist
double precision :: d_dist

double precision, dimension(neps,nk,na,2) :: vals
integer, dimension(neps,nk,na,2) :: inds

double precision :: dummy
double precision :: temp
double precision :: tot_inv_wedge

date=t1

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
put_eps=put_eps_t(date,:)

r=(one/beta-1)

    if ( ex_entry .eq. 0 ) then

    !bisection for w
        open (unit=1,file='input_instst_w.txt')
        read(1,*)wmin,wmax
        close(1)
    
    iter_bi=0
    d_bi=1

        if ( rank .eq. 0 ) then
        print*,' '
        print*,'fistst '
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
        
            write(6,'(300F12.7)')real(iter_bi),w,d_bi

                open (unit=1,file='input_fistst_w.txt')
                write(1,*)wmin,wmax
                close(1)

            end if
    
        end do
    
        if ( rank .eq. 0 .and. d_bi .ge. tol_bi ) then
        print*,'iterations exceeded at zero profit condition, fistst.f90 ',d_bi
        call exit
        end if

    end if
    
    if ( ex_entry .eq. 0 ) then
    
    !bisection for mu    
     open (unit=1,file='input_instst_mu.txt')
     read(1,*)mumin,mumax
     close(1)

!    mumin=mumin*.01
!    mumax=mumax*3.05

    else if ( ex_entry .eq. 1 ) then

    mu=ex_mu
    
        open (unit=1,file='input_instst_w.txt')
        read(1,*)wmin,wmax
        close(1)    
            
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

                open(unit=1,file='output_fistst_entry_cutoff.txt')
                write(1,*)'date d_e_o v_e cost mass' 
                    do ieps=1,neps
                    write(1,'(300F12.7)')real(date),real(g_e_o(ieps)),v_i(ieps,intbirth_k,intbirth),w*f_e+f_c+k_e,put_eps(ieps)
                    end do                    
                close(1)
                
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
        print*,'iterations exceeded in simulation, fistst.f90 ',d_dist
!             call exit
        end if
    
    !dividends and labor demand
    div=-mu*(w*f_e+f_c+k_e)*sum(put_eps*g_e_o)
    lab_dem=mu*f_e*sum(put_eps*g_e_o)    
    gdp=mu*(w_t(1)*f_e+f_c)*sum(put_eps*g_e_o)
    ngdp=mu*(w*f_e+f_c)*sum(put_eps*g_e_o)
    inv=mu*(w_t(1)*f_e+f_c+k_e)*sum(put_eps*g_e_o)    
    ninv=mu*(w*f_e+f_c+k_e)*sum(put_eps*g_e_o)
    firms_constrained=-1  !sum(omega_ss0(:,:,1))/(sum(omega_ss0)+mu)
    exits=zero
    exits_binding=zero    
    std_dev_growth_emp=zero
    std_dev_growth_sales=zero
    avg_sales=zero
    var_sales=zero

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
    
    tot_inv_wedge=zero
    collateral_binding=zero

    avg_prod_all_firms=zero
    avg_capital_all_firms=zero
    avg_net_assets_all_firms=zero
    avg_emp_all_firms=zero
    
    avg_prod_exits=zero
    avg_capital_exits=zero
    avg_net_assets_exits=zero
    avg_emp_exits=zero
    
    avg_spread=zero
    debtors=zero
    
        do ieps=1,neps
        do ik=1,nk
        do ia=1,na
    
            if ( ( ik .eq. intbirth_k .and. ia .eq. intbirth ) .or. omega_ss0(ieps,ik,ia) .gt. zero )  then
        
            eps=grid_eps(ieps)
            k=grid_k(ik)
            a=grid_a(date,ia)
            iz=id_z(ieps)
        
            d_o=g_i_o(ieps,ik,ia)
            ap=g_i_a(ieps,ik,ia)
            jk=g_i_k(ieps,ik,ia)
            
            kp=grid_k(jk)
        
                if ( d_o .eq. 1 ) then                    
                l=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))

                    if ( ap .ge. zero ) then
                                
                    q=one/(one+r)
                                
                    else if ( ap .lt. zero ) then
                                
                        if ( default_economy .eq. 0 ) then                                
                        q=one/(one+r)-tau_a
                        else if ( default_economy .eq. 1 ) then
                        ja=g_i_a_integer(ieps,ik,ia)
                        
                            if ( ja .le. 0 ) then
                            print*,'ja .le. 0'
                                call exit
                            end if
                        
                        q=g_q(ieps,jk,ja)
                        end if
                            
                    end if                
                                
                profit=((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-w*(l+f_o*grid_z(iz))&
                +(1.+tau_x)*((1-delta_k)*k-kp)+a-q*ap-func_adj(kp,k,d_o,date)-f_o_good            
                
                    if ( rank .eq. 0 .and. isnan(l) ) then
                    print*,'isnan l, fistst.f90',k
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
                gdp=gdp+(((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-f_o_good)*&
                (omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                ngdp=ngdp+(((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-f_o_good)*&
                (omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))                
                
                    if ( cap_adj .eq. 2 ) then
                    gdp=gdp+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)*(w_t(1)/w)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                    ngdp=ngdp+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))                   
                    end if                
                
                    if ( cap_adj .ne. 2 ) then                    
                    inv=inv+(kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+func_adj(kp,k,d_o,date))*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                    else if ( cap_adj .eq. 2 ) then                                
                    inv=inv+(kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+&
                    (func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)*w_t(1)/w+lambda*(kp/k-1)**2)*&
                    (omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))                   
                    end if
                
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
                collateral_binding=collateral_binding+g_binding_collateral(ieps,ik,ia)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                
                avg_prod_all_firms=avg_prod_all_firms+eps*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                avg_capital_all_firms=avg_capital_all_firms+k*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                avg_net_assets_all_firms=avg_net_assets_all_firms-a*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                avg_emp_all_firms=avg_emp_all_firms+((w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(iz))*&
                (omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                
                avg_prod_exits=avg_prod_exits+(1-d_o)*eps*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                avg_capital_exits=avg_capital_exits+(1-d_o)*k*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                avg_net_assets_exits=avg_net_assets_exits-(1-d_o)*a*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                avg_emp_exits=avg_emp_exits+(1-d_o)*((w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(iz))*&
                (omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                
                    if ( d_o .eq. 1 .and. k .gt. zero ) then
                    
                    x_i=kp-(1-delta_k)*k
                    
                    temp=temp+(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                    avg_i_k=avg_i_k+(x_i/k)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                    std_dev_i_k=std_dev_i_k+((x_i/k)**2.)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                        
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

                    if ( ap .lt. zero ) then
                    avg_spread=avg_spread+(one/q-1-r)*(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                    debtors=debtors+(omega_ss0(ieps,ik,ia)+mu*put_eps(ieps)*g_e_o(ieps))
                    end if
                    
                else          

                    if ( cap_adj .ne. 2 ) then     
                    lab_dem=lab_dem+(l+d_o*f_o*grid_z(iz))*omega_ss0(ieps,ik,ia)
                    else if ( cap_adj .eq. 2 ) then                
                    lab_dem=lab_dem+(l+d_o*f_o*grid_z(iz)+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)/w)*omega_ss0(ieps,ik,ia)
                    end if
                
                div=div+profit*omega_ss0(ieps,ik,ia)   
                gdp=gdp+(((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-f_o_good)*&
                omega_ss0(ieps,ik,ia)   
                ngdp=ngdp+(((Z*eps)**(1-alpha*nu))*(((k**alpha)*(l**(1-alpha)))**nu)-f_o_good)*&
                omega_ss0(ieps,ik,ia)                   

                    if ( cap_adj .eq. 2 ) then
                    gdp=gdp+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)*(w_t(1)/w)*omega_ss0(ieps,ik,ia)
                    ngdp=ngdp+(func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)*omega_ss0(ieps,ik,ia)                   
                    end if                
                
                    if ( cap_adj .ne. 2 ) then                    
                    inv=inv+(kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+func_adj(kp,k,d_o,date))*omega_ss0(ieps,ik,ia)
                    else if ( cap_adj .eq. 2 ) then                                
                    inv=inv+(kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+&
                    (func_adj(kp,k,d_o,date)-lambda*(kp/k-1)**2)*w_t(1)/w+lambda*(kp/k-1)**2)*&
                    omega_ss0(ieps,ik,ia)
                    end if

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
                collateral_binding=collateral_binding+g_binding_collateral(ieps,ik,ia)*omega_ss0(ieps,ik,ia)

                avg_prod_all_firms=avg_prod_all_firms+eps*omega_ss0(ieps,ik,ia)
                avg_capital_all_firms=avg_capital_all_firms+k*omega_ss0(ieps,ik,ia)
                avg_net_assets_all_firms=avg_net_assets_all_firms-a*omega_ss0(ieps,ik,ia)
                avg_emp_all_firms=avg_emp_all_firms+((w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(iz))*&
                omega_ss0(ieps,ik,ia)
                
                avg_prod_exits=avg_prod_exits+(1-d_o)*eps*omega_ss0(ieps,ik,ia)
                avg_capital_exits=avg_capital_exits+(1-d_o)*k*omega_ss0(ieps,ik,ia)
                avg_net_assets_exits=avg_net_assets_exits-(1-d_o)*a*omega_ss0(ieps,ik,ia)
                avg_emp_exits=avg_emp_exits+(1-d_o)*((w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))+f_o*grid_z(iz))*&
                omega_ss0(ieps,ik,ia)
                
                    if ( d_o .eq. 1 .and. k .gt. zero ) then
                    
                    x_i=kp-(1-delta_k)*k

                    temp=temp+omega_ss0(ieps,ik,ia)                    
                    avg_i_k=avg_i_k+(x_i/k)*omega_ss0(ieps,ik,ia)
                    std_dev_i_k=std_dev_i_k+((x_i/k)**2.)*omega_ss0(ieps,ik,ia)
                        
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

                    if ( ap .lt. zero ) then
                    avg_spread=avg_spread+(one/q-1-r)*omega_ss0(ieps,ik,ia)
                    debtors=debtors+omega_ss0(ieps,ik,ia)
                    end if
                    
                end if

            dummy=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
            
            end if
            
        end do
        end do
        end do
        
avg_prod_all_firms=avg_prod_all_firms/(sum(omega_ss0)+mu*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth)))
avg_capital_all_firms=avg_capital_all_firms/(sum(omega_ss0)+mu*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth)))
avg_net_assets_all_firms=avg_net_assets_all_firms/(sum(omega_ss0)+mu*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth)))
avg_emp_all_firms=avg_emp_all_firms/(sum(omega_ss0)+mu*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth)))
                
avg_prod_exits=avg_prod_exits/exits
avg_capital_exits=avg_capital_exits/exits
avg_net_assets_exits=avg_net_assets_exits/exits
avg_emp_exits=avg_emp_exits/exits

avg_spread=avg_spread/debtors

    avg_i_k=avg_i_k/temp
    std_dev_i_k=std_dev_i_k/temp
    std_dev_i_k=sqrt(std_dev_i_k-avg_i_k**2)
    i_k_inaction=i_k_inaction/temp
    i_k_pos=i_k_pos/temp
    i_k_neg=i_k_neg/temp
    i_k_pos_spike=i_k_pos_spike/temp
    i_k_neg_spike=i_k_neg_spike/temp
    
!     h=(1-psi*(r*b+div)/w)/(1+psi)      

    !bisection for labor    
    
        if ( GHH .eq. 0 ) then
        
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
            print*,'iterations exceeded labor bisection, fistst'
                call exit
            end if
        
        else if ( GHH .eq. 1 ) then
        
        h=(w*(1-tau_l)/psi)**(one/phi)
        
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

        d_bi=abs(wmax-wmin) !max(abs(lab_dem-h),abs(wmax-wmin))

        end if
        
    iter_bi=iter_bi+1
    
        if ( rank .eq. 0 ) then
        
            if ( ex_entry .eq. 0 ) then
    
            write(6,'(300F12.7)')real(iter_bi),mu,d_bi
            print*,'exdem, h',exdem,lab_dem,h
        
                open (unit=1,file='input_fistst_mu.txt')
                write(1,*)mumin,mumax
                close(1)            
            
            else if ( ex_entry .eq. 1 ) then
            
            write(6,'(300F12.7)')real(iter_bi),w,d_bi
            print*,'exdem, h',exdem,lab_dem,h
            print*,'sum(g_e_o)',sum(g_e_o)
        
                open (unit=1,file='input_fistst_w.txt')
                write(1,*)wmin,wmax
                close(1)

            end if
        
        end if
    
    end do
    
    if ( rank .eq. 0 .and. abs(d_bi) .ge. tol_bi ) then
    print*,'iterations exceeded in labor market clearing, fistst.f90 ',d_bi
        call exit
    end if    

avg_firm_size=h/(sum(omega_ss0)+mu*sum(put_eps*g_e_o))
labor_share=w*h/ngdp
avg_sales=avg_sales/(sum(omega_ss0)+mu*sum(put_eps*g_e_o))    
var_sales=var_sales/(sum(omega_ss0)+mu*sum(put_eps*g_e_o))    
var_sales=var_sales-avg_sales**2
std_dev_sales=sqrt(var_sales)
    
firm_debt=zero

    do ia=1,intbirth
    firm_debt=firm_debt+grid_a(date,ia)*sum(omega_ss0(:,:,ia))
    end do

firm_capital=mu*k_e*sum(put_eps*g_e_o)
    
    do ik=1,nk
    firm_capital=firm_capital+grid_k(ik)*sum(omega_ss0(:,ik,:))
    end do
    
    if ( rank .eq. 0 ) then
    
    print*,'accidental_defaults',collateral_binding/(sum(omega_ss0)+mu*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth)))
    
        open (unit=1,file='output_fistst_dist.txt')
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
        print*,'increase kmax, fistst.f90 ',sum(omega_ss0(:,nk,:))/(sum(omega_ss0)+mu*sum(put_eps*g_e_o))
        end if
        
        if ( sum(omega_ss0(:,:,1)) .gt. zero ) then
        print*,'decrease amin, fistst.f90 ',sum(omega_ss0(:,:,1))/(sum(omega_ss0)+mu*sum(put_eps*g_e_o))
        end if
    
    print*,'firm debt/gdp',-firm_debt/gdp,neg_b/gdp
    
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
collateral_binding_t(date)=collateral_binding
lab_dem_t(date)=lab_dem

avg_prod_all_firms_t(date)=avg_prod_all_firms
avg_capital_all_firms_t(date)=avg_capital_all_firms
avg_net_assets_all_firms_t(date)=avg_net_assets_all_firms
avg_emp_all_firms_t(date)=avg_emp_all_firms
                
avg_prod_exits_t(date)=avg_prod_exits
avg_capital_exits_t(date)=avg_capital_exits
avg_net_assets_exits_t(date)=avg_net_assets_exits
avg_emp_exits_t(date)=avg_emp_exits
avg_spread_t(date)=avg_spread

omega_ss_final=omega_ss0

g_i_o_final=g_i_o
g_i_k_final=g_i_k
g_i_a_final=g_i_a
g_e_o_final=g_e_o
g_q_final=g_q
g_un_k_final=g_un_k
g_un_o_final=g_un_o
g_binding_collateral_final=g_binding_collateral
g_i_a_integer_final=g_i_a_integer

a_un_t(date,:,:)=a_un
a_tilde_t(date,:,:)=a_tilde

!compute employment by age
emp_by_age=zero
firms_by_age=zero
exit_by_age=zero
exit_binding_by_age=zero
investment_by_age=zero
firms_collateral_binding_by_age=zero

omega_ss0=zero
omega_ss0(:,intbirth_k,intbirth)=mu*put_eps*g_e_o

    do iter_dist=1,119       !age0, age1, age2, age,3, age4, age5

    firms_by_age(iter_dist)=sum(omega_ss0)
    
    omega_ss1=zero
        
        do ieps=1,neps
        do ik=1,nk
        do ia=1,na
    
            if ( omega_ss0(ieps,ik,ia) .gt. zero )  then
        
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
            investment_by_age(iter_dist)=investment_by_age(iter_dist)+(kp-d_o*(1-delta_k)*k-(1-d_o)*chi*k+func_adj(kp,k,d_o,date))*omega_ss0(ieps,ik,ia)
            firms_collateral_binding_by_age(iter_dist)=firms_collateral_binding_by_age(iter_dist)+g_binding_collateral(ieps,ik,ia)*omega_ss0(ieps,ik,ia)
            
!             dummy=(w/(((Z*eps)**(1-alpha*nu))*(1-alpha)*nu*(k**(alpha*nu))))**(one/((1-alpha)*nu-1))
!  
!                 if ( a .ge. zero ) then
!                 
!                     if ( ((Z*eps)**(1-alpha*nu))*(((k**alpha)*(dummy**(1-alpha)))**nu)-w*(dummy+f_o)+(1.+tau_x)*(1-delta_k)*k+&
!                     (1+r)*a-(1-delta_f)*abar+theta*k .lt. zero ) then
!                             
!                         if ( d_o .eq. 0 ) then
!                         exit_binding_by_age(iter_dist)=exit_binding_by_age(iter_dist)+omega_ss0(ieps,ik,ia)
!                         else if ( d_o .ne. 0 ) then
!                         print*,'d_o .ne. 0, instst.f90 '
!                             call exit                    
!                         end if                                    
!                 
!                     end if
!                     
!                 else if ( a .lt. zero ) then
! 
!                     if ( ((Z*eps)**(1-alpha*nu))*(((k**alpha)*(dummy**(1-alpha)))**nu)-w*(dummy+f_o)+(1.+tau_x)*(1-delta_k)*k+&
!                     (1+r+tau_a)*a-(1-delta_f)*abar+theta*k .lt. zero ) then
!                             
!                         if ( d_o .eq. 0 ) then
!                         exit_binding_by_age(iter_dist)=exit_binding_by_age(iter_dist)+omega_ss0(ieps,ik,ia)
!                         else if ( d_o .ne. 0 ) then
!                         print*,'d_o .ne. 0, instst.f90 '
!                             call exit                    
!                         end if                                    
!                 
!                     end if
!                     
!                 end if
                
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

firms_by_age(120)=mass_t(date)+mu_t(date)*sum(put_eps*g_e_o)-sum(firms_by_age(1:119))  
emp_by_age(120)=h_t(date)-sum(emp_by_age(1:119))
exit_by_age(120)=exits_t(date)-sum(exit_by_age(1:119))
exit_binding_by_age(120)=exits_binding_t(date)-sum(exit_binding_by_age(1:119))
investment_by_age(120)=ninv-mu*(w*f_e+f_c+k_e)*sum(put_eps*g_e_o)-sum(investment_by_age(1:119))
firms_collateral_binding_by_age(120)=collateral_binding-sum(firms_collateral_binding_by_age(1:119))

firms_by_age_t(date,:)=firms_by_age(1:120) 
emp_by_age_t(date,:)=emp_by_age(1:120)
exit_by_age_t(date,:)=exit_by_age(1:120)
exit_binding_by_age_t(date,:)=exit_binding_by_age(1:120)
investment_by_age_t(date,:)=investment_by_age(1:120)
firms_collateral_binding_by_age_t(date,:)=firms_collateral_binding_by_age(1:120)

!compute job creation and job destruction
omega_ss0=omega_ss_final
job_creation=zero
job_destruction=zero

    do ieps=1,neps
    do ik=1,nk
    do ia=1,na    
    
        if ( omega_ss0(ieps,ik,ia) .gt. zero )  then
        
        eps=grid_eps(ieps)
        k=grid_k(ik)
        a=grid_a(date,ia)
        iz=id_z(ieps)
        
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
                
                    if ( pi_eps(ieps,jeps) .gt. zero ) then

                    jz=id_z(jeps)
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
        
        eps=grid_eps(ieps)
        k=grid_k(ik)
        a=grid_a(date,ia)
        iz=id_z(ieps)
        
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

end subroutine fistst
