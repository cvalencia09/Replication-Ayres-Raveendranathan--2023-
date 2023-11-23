module params
implicit none

!mpi variables
integer :: ierr, rank, nproc

!parameters
integer, parameter :: usespline=0, cap_adj=1, GHH=0, Habit=0        !0=no adjustment cost, 1=costly reversibility, 2=non-convex
integer, parameter :: ex_entry=1, default_economy=0, compute_job_flows=0
double precision, parameter :: one=1, zero=0
double precision, parameter :: beta=.99
double precision, parameter :: nu=0.83625, delta_k=1-(1-.06)**.25, delta_f=0, chi=(1-delta_k)
double precision, parameter :: phi=.5
double precision, parameter :: lambda_c=0

double precision :: ex_mu

double precision :: lambda

double precision, parameter :: rho_eps=.859**.25, mu_eps=zero
double precision, parameter :: rho_iota=zero, mu_iota=zero
double precision :: var_iota
double precision, parameter :: f_c=zero

double precision :: z, abar, tau_a, psi, f_e, f_o, var_eps, tau_x, tau_l, k_e, theta, var_z, f_o_good
double precision :: gamma, gamma1, gamma2, f_k      !adjustment cost parameters
double precision :: alpha

integer, parameter :: t1=120

!state space
integer, parameter :: nz=1, neta=31, niota=1
integer, parameter :: neps=neta*nz*niota, nk=350, na=125
integer :: id_z(neps), id_iota(neps)
double precision :: grid_eps(neps), grid_k(nk), grid_a(t1,na), grid_z(nz), grid_eps_original(neps)
double precision :: eps, k, a, q
integer :: ieps, ik, ia, iz, ieta, iiota
integer :: jeps, jk, ja, jz, jiota
double precision, dimension(neps,neps) :: pi_eps
double precision, dimension(neps) :: put_eps
double precision, parameter :: kmin=1e-13, kmax=400, amax=450, amin=-99
integer, parameter :: intbirth=na/2+1
integer :: intbirth_k

!incumbent's value function and policy functions
double precision, dimension(neps,nk,na) :: v_i
integer(kind=1), dimension(neps,nk,na) :: g_i_o, g_binding_collateral   !g_binding_collateral captures the fact that firm cannot repay at all
double precision, dimension(neps,nk,na) :: g_i_a
integer, dimension(neps,nk,na) :: g_i_a_integer
integer(kind=2), dimension(neps,nk,na) :: g_i_k
double precision, dimension(neps,nk,na) :: g_q

!entrant's policy function
integer(kind=1), dimension(neps) :: g_e_o

!initial distribution
double precision, dimension(neps,nk,na) :: omega_ss_initial
integer(kind=1), dimension(neps,nk,na) :: g_i_o_initial, g_binding_collateral_initial
double precision, dimension(neps,nk,na) :: g_i_a_initial
integer, dimension(neps,nk,na) :: g_i_a_integer_initial
integer(kind=2), dimension(neps,nk,na) :: g_i_k_initial
integer(kind=1), dimension(neps) :: g_e_o_initial
double precision, dimension(neps,nk,na) :: g_q_initial

integer(kind=1), dimension(neps,nk,na) :: g_i_o_final, g_binding_collateral_final
double precision, dimension(neps,nk,na) :: g_i_a_final
integer, dimension(neps,nk,na) :: g_i_a_integer_final
integer(kind=2), dimension(neps,nk,na) :: g_i_k_final
integer(kind=1), dimension(neps) :: g_e_o_final
double precision, dimension(neps,nk,na) :: g_q_final

!value of entry
double precision :: v_e

!consumer's policy functions
double precision :: c, h, b, exdem, neg_b

!storing variables
double precision :: l, kp, ap, profit, kpmin, kpmax, x_i, constant
integer(kind=1) :: d_o
double precision, dimension(na) :: exp_ap

!prices and aggregates
double precision :: w, div, mu, lab_dem, r

!statistics
double precision :: gdp, inv, firms_constrained, exits, avg_firm_size, labor_share
double precision :: std_dev_growth_emp, std_dev_growth_sales, exits_binding
double precision :: firm_debt, firm_capital
double precision :: ngdp, ninv
double precision :: avg_i_k, std_dev_i_k, i_k_inaction, i_k_pos, i_k_neg, i_k_pos_spike, i_k_neg_spike
double precision :: avg_i_k30, std_dev_i_k30, i_k_inaction30, i_k_pos30, i_k_neg30, i_k_pos_spike30, i_k_neg_spike30, &
i_k_pos_high_spike30
double precision :: avg_i_k30check, std_dev_i_k30check


!for mpi
integer, parameter :: ni=neps*nk*na
integer :: cnt, np, itop, iend
integer :: zmap(ni,3)
integer, dimension(neps,nk,na) :: zmap_inv

double precision, dimension(t1) :: c_t, h_t, b_t, exdem_t, lab_dem_t
double precision, dimension(t1) :: w_t, div_t, mu_t, r_t

double precision, dimension(t1) :: z_t, abar_t, tau_a_t, psi_t, f_e_t, var_eps_t, tau_x_t, tau_l_t, theta_t, psi_c_t
double precision, dimension(t1) :: mass_t, f_o_t
double precision, dimension(t1,neps) :: grid_eps_t

!statistics
double precision, dimension(t1) :: gdp_t, inv_t, firms_constrained_t, exits_t, avg_firm_size_t, labor_share_t
double precision, dimension(t1) :: std_dev_growth_emp_t, std_dev_growth_sales_t, exits_binding_t
double precision, dimension(t1) :: firm_debt_t, firm_capital_t
double precision, dimension(t1) :: ngdp_t, ninv_t

double precision, parameter :: tol_bi=1e-6
double precision, parameter :: tol_v=1e-6
double precision, parameter :: tol_dist=1e-9
double precision, parameter :: tol_brent=1e-9

integer, parameter :: itmax_bi=100
integer, parameter :: itmax_v=100
integer, parameter :: itmax_dist=2500

!variable to specify if value function guess should be read (1) or not (0)
integer :: read_v_i

integer :: date, iter_tr

! age 0, 1, 2, 3, 4, 5, rest
double precision, dimension(t1,120) :: emp_by_age_t, firms_by_age_t, exit_by_age_t, exit_binding_by_age_t, investment_by_age_t, &
firms_collateral_binding_by_age_t
double precision, allocatable, dimension(:,:,:,:) :: omega_ss_age_initial
integer :: tc, tcc

double precision :: job_creation, job_destruction, lp
double precision, dimension(t1) :: job_creation_t, job_destruction_t

double precision :: avg_sales, var_sales, std_dev_sales
double precision :: auto_cov_sales, auto_corr_sales

integer, parameter :: max_age=600
double precision, dimension(max_age) :: emp_by_age, firms_by_age, exit_by_age, exit_binding_by_age, div_by_age, div_operating_by_age, &
firms_collateral_binding_by_age
double precision, dimension(max_age) :: eps_by_age, k_by_age, a_by_age, job_creation_by_age, job_destruction_by_age
double precision, dimension(max_age) :: avg_log_emp_by_age, std_dev_log_emp_by_age, dummy_by_age, investment_by_age

double precision :: yp1, ypn
double precision, dimension(nk) :: y2temp

!variables for labor bisection
double precision :: d_bih, hmin, hmax
integer :: iter_bih
double precision, parameter :: tol_bih=1e-9
integer, parameter :: itmax_bih=500

!variables for unconstrained firm problem
integer, parameter :: ni_un=neps*nk
integer :: np_un, itop_un, iend_un
integer :: zmap_un(ni_un,2)
integer, dimension(neps,nk) :: zmap_inv_un

integer :: read_v_un
double precision, dimension(neps,nk) :: v_un
integer(kind=1), dimension(neps,nk) :: g_un_o, g_un_o_initial, g_un_o_final
integer(kind=2), dimension(neps,nk) :: g_un_k, g_un_k_initial, g_un_k_final
double precision :: kp_temp
double precision, dimension(neps,nk) :: a_un, a_un1, a_tilde, a_tilde1    
double precision, dimension(t1,neps,nk) :: a_un_t, a_tilde_t

double precision :: firms_21_plus

double precision :: cut_off_model, cut_off_data
double precision :: shock, rho_shock
double precision, dimension(t1,neps) :: put_eps_t


double precision :: collateral_binding, collateral_binding_t(t1)

double precision :: avg_prod_all_firms, avg_capital_all_firms, avg_net_assets_all_firms, avg_emp_all_firms
double precision :: avg_prod_exits, avg_capital_exits, avg_net_assets_exits, avg_emp_exits
double precision, dimension(t1) :: avg_prod_all_firms_t, avg_capital_all_firms_t, avg_net_assets_all_firms_t, avg_emp_all_firms_t
double precision, dimension(t1) :: avg_prod_exits_t, avg_capital_exits_t, avg_net_assets_exits_t, avg_emp_exits_t

double precision, dimension(t1) :: output_loss_entry_t, output_gain_entry_t, output_loss_exit_t, output_gain_exit_t, &
output_loss_incumbent_t, output_gain_incumbent_t
double precision, dimension(t1) :: employment_loss_entry_t, employment_gain_entry_t, employment_loss_exit_t, employment_gain_exit_t, &
employment_loss_incumbent_t, employment_gain_incumbent_t
double precision, dimension(t1) :: investment_loss_entry_t, investment_gain_entry_t, investment_loss_exit_t, investment_gain_exit_t, &
investment_loss_incumbent_t, investment_gain_incumbent_t
double precision, dimension(t1) :: debt_loss_entry_t, debt_gain_entry_t, debt_loss_exit_t, debt_gain_exit_t, &
debt_loss_incumbent_t, debt_gain_incumbent_t

end module params
