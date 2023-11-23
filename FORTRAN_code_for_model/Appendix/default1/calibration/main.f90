program main
use params
use mini
implicit none
include 'mpif.h'

!mpif90 params.f90 mini.f90 grid.f90 instst.f90 fistst.f90 transition.f90 main.f90 -o ./main
!mpiexec -np 8 main.exe

!calibration variables
integer :: iter_cal
double precision, parameter :: pow=.3
double precision :: d_cal
double precision :: entrants, entrant_relative_size

double precision, parameter :: firm_debt_gdp_data=1.839557811, emp_data=one, lab_share_data=0.645285572
double precision, parameter :: std_dev_i_k_data=.337, five_year_survival_data=0.462291371, entrant_size_data=0.8088830
double precision, parameter :: entry_data=0.108798926, i_k_inaction_data=.081, i_k_spike_data=.204
double precision, parameter :: i_k_pos_spike_data=.186, i_k_pos_high_spike_data=.018
double precision, parameter :: cv_i_k_data=.337/.122, entrant_relative_size_data=0.5919066
double precision, parameter :: firms_21_plus_data=0.231510632, firms_5_firms0_5_data=0.118036885
double precision, parameter :: std_dev_log_emp_0=.9, std_dev_log_emp_19=1.3
double precision :: theta1, psi1, alpha1, var_eps1, f_o1, k_e1        

double precision, parameter :: avg_spread_data=1.0234
integer :: iter_shock
double precision, parameter :: gdp_drop_data1=0.959449549, gdp_drop_data2=(gdp_drop_data1+one)/2.

!transition variables
integer :: t

call MPI_INIT (ierr)
call MPI_COMM_SIZE (mpi_comm_world, nproc, ierr)
call MPI_COMM_RANK (mpi_comm_world, rank, ierr)

    if ( habit .eq. 1 .and. GHH .eq. 0 ) then
    print*,'habit .eq. 1 .and. GHH .eq. 0'
        call exit
    end if
        
	if ( rank .eq. 0 ) then
	allocate (omega_ss_age_initial(120,neps,nk,na))
	end if
    
read_v_i=0
read_v_un=0

call uniform_parallel

!initial guesses
! abar=-3.25987617384804

! theta=.93   !collateral constraint
! psi=1.75
! alpha=0.327
! var_eps=0.18
! f_o=.14
! k_e=1        !entry capital

!adjustment cost parameters
gamma=.05
gamma1=-.07
gamma2=.07
f_k=.1

lambda=one

var_z=.05
ex_mu=.1
var_iota=.3
cut_off_data=1.4437959
chi=0.492325089

    open (unit=1,file='input_parameters.txt')
    read(1,*)theta,psi,alpha,var_eps,f_o,k_e,gamma,lambda,ex_mu,var_iota,cut_off_model,f_o_good,chi
    close(1)
    
! lambda=zero
    
f_o_good=zero
cut_off_data=cut_off_model 

    !var_eps=var_eps/4
!f_o=f_o/4
!ex_mu=ex_mu/4
!theta=theta/4

!theta=1
!psi=1.8
!alpha=.32
!var_eps=.1
!f_o=.15
!k_e=1.5
!gamma=.03

!f_e_t=one   !normalization

abar_t=zero

theta_t=theta
psi_t=psi
var_eps_t=var_eps
f_o_t=f_o
f_e_t=one
psi_c_t=one

!shocks
z_t=one
tau_a_t=1.013**.25-1
tau_x_t=zero
tau_l_t=zero

!calibration
    
    do iter_cal=1,1

        call grid
      
        call instst

     entrants=sum(firms_by_age(2:5))          !mu_t(1)*sum(put_eps*g_e_o*g_i_o(:,intbirth_k,intbirth))
     entrant_relative_size=((emp_by_age(1)-f_e*mu*sum(put_eps*g_e_o))/firms_by_age(2))/(emp_by_age(6)/firms_by_age(7))
     firms_21_plus=sum(firms_by_age(86:max_age))/sum(firms_by_age(2:max_age))   !sum(firms_by_age(23:max_age))/sum(firms_by_age(2:max_age))
     
     d_cal=max(abs(firm_debt_gdp_data+firm_debt/ngdp),&
     abs(emp_data-h),&
     abs(lab_share_data-labor_share_t(1)),&
     abs(std_dev_i_k30check-std_dev_i_k_data),&
!     abs(std_dev_log_emp_19-std_dev_log_emp_by_age(20)),&
!     abs(std_dev_log_emp_0-std_dev_log_emp_by_age(1)),&     
     abs((sum(firms_by_age(22:25))/sum(firms_by_age(2:25)))-firms_5_firms0_5_data),&
     abs(entry_data-entrants/(firms_by_age(2)+sum(omega_ss_initial))),&
     abs(i_k_inaction_data-i_k_inaction30),&
     abs(i_k_pos_high_spike_data-i_k_pos_high_spike30),&
     abs(firms_21_plus-firms_21_plus_data))
!     abs(entrant_relative_size_data-entrant_relative_size))    
 
     k_e1=k_e*(firm_debt_gdp_data/(-firm_debt/ngdp))**(pow)
     psi1=psi*(h/emp_data)**pow
     alpha1=alpha*(labor_share_t(1)/lab_share_data)**pow
     var_eps1=var_eps*(std_dev_i_k_data/std_dev_i_k30check)**pow
     f_o1=f_o!*(entry_data/(entrants/(firms_by_age(2)+sum(omega_ss_initial))))**(pow/2)
     gamma=gamma!*(i_k_inaction_data/i_k_inaction30)**(pow/2)
     lambda=lambda!*(i_k_pos_high_spike30/i_k_pos_high_spike_data)**(pow/4)
     ex_mu=ex_mu!*(firms_21_plus/firms_21_plus_data)**(pow/4)
     chi=min(chi*(((1.+avg_spread)**4)/avg_spread_data)**.7,1-delta_k)
     
         if ( rank .eq. 0 ) then
         
         print*,''
         print*,'iter_cal ',iter_cal,d_cal
         print*,'k_e firm debt to GDP'
         write(6,'(300F14.7)')k_e,k_e1,firm_debt_gdp_data,-firm_debt/ngdp
         print*,'psi employment'   
         write(6,'(300F14.7)')psi,psi1,emp_data,h
         print*,'alpha labor share'
         write(6,'(300F14.7)')alpha,alpha1,lab_share_data,labor_share_t(1)
         print*,'var_eps std_dev_i_k '
         write(6,'(300F14.7)')var_eps,var_eps1,std_dev_i_k_data,std_dev_i_k30check              
         print*,'f_o entry_rate'
         write(6,'(300F14.7)')f_o,f_o1,entry_data,entrants/(firms_by_age(2)+sum(omega_ss_initial))
         print*,'gamma i_k_inaction'
         write(6,'(300F14.7)')gamma,gamma,i_k_inaction_data,i_k_inaction30
         print*,'lambda i_k_pos_high_spike'
         write(6,'(300F14.7)')lambda,lambda,i_k_pos_high_spike_data,i_k_pos_high_spike30
         print*,'mu firms_21_plus'
         write(6,'(300F14.7)')ex_mu,ex_mu,firms_21_plus_data,firms_21_plus    
         print*,'chi spread'
         write(6,'(300F14.7)')chi,chi,avg_spread_data,(1.+avg_spread)**4             
         print*,'cut_off_model',cut_off_model

!          print*,'avg_i_k',avg_i_k,avg_i_k30,avg_i_k30check
!          print*,'std_dev_i_k',std_dev_i_k,std_dev_i_k30,std_dev_i_k30check
!          print*,'i_k_inaction',i_k_inaction,i_k_inaction30
!          print*,'i_k_pos',i_k_pos,i_k_pos30
!          print*,'i_k_neg',i_k_neg,i_k_neg30
!          print*,'i_k_pos_spike',i_k_pos_spike,i_k_pos_spike30
!          print*,'i_k_neg_spike',i_k_neg_spike,i_k_neg_spike30
         
         print*,''                
             
         end if
             
     theta=theta1
     psi=psi1
     alpha=alpha1
     var_eps=var_eps1
     k_e=k_e1
     f_o=f_o1
     
     theta_t=theta
     psi_t=psi
     var_eps_t=var_eps
     f_o_t=f_o
     f_e_t=f_e    
 
         if ( rank .eq. 0 ) then	        
         
            if ( isnan(theta) .or. isnan(psi) .or. isnan(alpha) .or. isnan(var_eps) .or. isnan(f_o) .or. isnan(k_e) .or. isnan(gamma) &
            .or. isnan(lambda) .or. isnan(ex_mu) .or. isnan(var_iota) ) then
            print*,'isnan main.f90'
                call exit
            else          
            
!                open (unit=1,file='input_parameters.txt')
!                write(1,*)theta,psi,alpha,var_eps,f_o,k_e,gamma,lambda,ex_mu,var_iota,cut_off_model,f_o_good,chi
!                close(1)
                
            end if
            
         end if

	    call MPI_BCAST (k_e,1,mpi_double_precision,0,mpi_comm_world,ierr)        

    end do

! experiment
! z_t(2)=.95*z_t(1)
!       
!     do t=3,90
!     z_t(t)=.88*z_t(t-1)+.12*z_t(1)
!     end do
     
! theta_t(2)=.6*theta_t(1)
!      
!     do t=3,90
!     theta_t(t)=.95*theta_t(t-1)+.05*theta_t(1)
!     end do
        
!abar_t(2)=.75*abar_t(1)
     
!     do t=3,30
!     abar_t(t)=.65*abar_t(t-1)+.35*abar_t(1)
!     end do
    
!tau_a_t(2)=4*tau_a_t(1)
    
!    do t=3,30
!    tau_a_t(t)=.65*tau_a_t(t-1)+.35*tau_a_t(1)
!    end do


!     if ( rank .eq. 0 ) then
!         open(unit=4,file='output_comp_stats.txt',status='replace')
!         write(4,*)'var_eps gamma f_k gamma1 gamma2 avg_i_k std_dev_i_k i_k_inaction i_k_pos i_k_neg i_k_pos_spike i_k_neg_spike'
!         close(4)
!     end if
    
! read_v_i=0
! read_v_un=0
    
!     do iter_cal=1,10
! 
! !     var_eps=0.147161353264860/2.+(0.147161353264860*1.5-0.147161353264860/2.)*(iter_cal-1.)/(10.-1.)
! !     var_eps_t=var_eps
!     
! !     gamma=0.11+(0.11-0.01)*(iter_cal-1.)/(10.-1.)
! !     f_k=0.15+(0.15-0.05)*(iter_cal-1.)/(10.-1.)
! !     gamma1=-0.1+(-0.01+0.1)*(iter_cal-1.)/(10.-1.)
!     gamma2=0.01+(0.1-0.01)*(iter_cal-1.)/(10.-1.)
!     
!         call grid
!             
!         call instst
!             
!         if ( rank .eq. 0 ) then
!             open(unit=4,file='output_comp_stats.txt',action='write',position='append')  
!             write(4,'(300F14.9)')var_eps,gamma,f_k,gamma1,gamma2,avg_i_k,std_dev_i_k,i_k_inaction,i_k_pos,i_k_neg,i_k_pos_spike,i_k_neg_spike
!             close(4)
!         end if
!             
!     end do

!      call grid
           
! read_v_i=0
! read_v_un=0
!              
!      call instst
! 
! read_v_i=0
! read_v_un=0
!       
!      call fistst
!    
!      call transition

! experiment
! shock=.95*z_t(1)
! rho_shock=.88
! 
!       open (unit=1,file='input_shock.txt')
!       read(1,*)shock,rho_shock
!       close(1)
! 
!     do iter_shock=1,1
! 
!     z_t(2)=shock
!       
!         do t=3,90
!         z_t(t)=rho_shock*z_t(t-1)+(1-rho_shock)*z_t(1)
!         end do
! 
!         call grid
!     
! 	if ( iter_shock .eq. 1 ) then
! 
! 	read_v_i=0
!         read_v_un=0
!                 
! 	        call instst
!       
!         	call fistst
! 
! 	end if
!    

!    do t=1,t1
        
!        if ( habit .eq. 1 .and. psi_c_t(t) .ne. one ) then
!        print*,'habit .eq. 1 .and. psi_c_t(t) .ne. one'
!            call exit
!        end if
        
!    end do


!         call transition
! 
!     shock=shock*(gdp_drop_data1/(gdp_t(2)/gdp_t(1)))**.6
!     rho_shock=rho_shock*((gdp_t(15)/gdp_t(1))/gdp_drop_data2)**.6
!     
!         if ( rank .eq. 0 ) then
!         
!         print*,'iter_shock',iter_shock
!         print*,'shock gdp_drop_data1 gdp_t(2)/gdp_t(1)'
!         write(6,'(300F12.4)')shock,gdp_drop_data1,gdp_t(2)/gdp_t(1)
! 
!         print*,'rho_shock gdp_drop_data2 gdp_t(5)/gdp_t(1)'
!         write(6,'(300F12.4)')rho_shock,gdp_drop_data2,gdp_t(15)/gdp_t(1)
!         
! !           open (unit=1,file='input_shock.txt')
! !           write(1,*)shock,rho_shock
! !           close(1)
! 
!         end if
! 
!         call MPI_BCAST (shock,1,mpi_double_precision,0,mpi_comm_world,ierr)        
!         call MPI_BCAST (rho_shock,1,mpi_double_precision,0,mpi_comm_world,ierr)        
!     
!     end do     
     
	if ( rank .eq. 0 ) then
	deallocate (omega_ss_age_initial)
	end if

call MPI_FINALIZE(ierr)

end program main

