Notes
 - Software used: FORTRAN with Intel parallel processing, Excel, and Matlab
 - Operating system used: Windows 10
 - General order in which the programs need to run (more details below)
	- Run FORTRAN code to generate model results
		- The 'test.sh' file in each sub-folder with the code contains the order in which programs need to be run
	- Copy results from FORTRAN output to excel and perform calculations
	- Import data from excel to matlab to generate figures in paper 
 - Expected computation time
	- With 32 processors, a steady state converges in less than 1 hour
	- With 32 processors, a transition converges in 24-36 hours


Description of files and folders

'Excel_files_with_data_and_model_output'
	- 'Appendix' contains excel files with calculations to generate results for figures in Online Appendix
		- '1980_recession' contains excel files with calculations for 1980s Double-Dip recession exercise
		- 'default1' contains excel files with calculations for model with default (calibration 1)
		- 'default2' contains excel files with calculations for model with default (calibration 2)
		- 'GHH' contains excel files with calculations for model GHH preferences
		- 'GHH_habit' contains excel files with calculations for model GHH + Habit formation preferences
		- 'larger_spread' contains excel files with calculations for model with a large spread
		- 'operating_cost_final_good' contains excel files with calculations for model with fixed operating costs in units of final good
		- 'perfectly_elastic_entry' contains excel files with calculations for model perfectly elastic mass of potential entrants
		- 'random_fixed_cost1' contains excel files with calculations for model random fixed operating costs (calibration 1)
		- 'random_fixed_cost2' contains excel files with calculations for model random fixed operating costs (calibration 1)

	- 'Section_5' contains excel files with calculations to generate results for figures in Section 5

	- 'Section_6.1' contains excel files with calculations to generate results for figures in Section 6.1 (Aggregate shocks, firm entry, and exit)
		- 'baseline' contains excel files with calculations to generate results for figures 5-7
		- 'partial_equilibrium' contains excel files with calculations to generate results for figure 7

	- 'Section_6.2' contains excel files with calculations to generate results for figures in Section 6.2 (Great Recession)

	- 'Section_6.3' contains excel files with calculations to generate results for figures in Section 6.3 (COVID-19 lockdown)

	- 'Entry_and_exit_BDS' contains data on firm entry and exit used from the BLS

	- 'GDP_consumption_investment_hours_firm_debt' contains data on GDP, consumption, investment, total hours, & firm debt

	- 'Productivity_deviations_BLS' contains data on productivity deviations used in the model for the Great Recession exercise (Section 6.2) and 1980s Double-Dip recession exercise (Online Appendix)


'FORTRAN_code_for_model' contains FORTRAN code used to generate model results

	- 'Appendix'
		- '1980_recession' contains code to generate results for Online Appendix A.3 (1980-Double-Dip recession exercise)
		- 'default1' contains code to generate results for Online Appendix A.5.1 (productivity shock with default - calibration 1)
		- 'default2' contains code to generate results for Online Appendix A.5.1 (productivity shock with default - calibration 2)
		- 'GHH'	contains code to generate results for Online Appendix A.5.2 with GHH preferences
			- 'calibration' contains code to calibrate model with GHH preferences
			- 'theta' contains code to compute transition with credit shock with GHH preferences
			- 'z' contains code to compute transition with productivity shock with GHH preferences
		- 'GHH_habit' contains code to generate results for Online Appendix A.5.2 with GHH + Habit formation preferences
			- 'theta'  contains code to compute transition with credit shock with GHH + Habit formation preferences
			- 'z'  contains code to compute transition with productivity shock with GHH + Habit formation preferences
		- 'large_spread' contains code to generate results for Online Appendix A.5.2 with a larger spread
			- 'calibration' contains code to calibrate model with a larger spread
			- 'theta' contains code to compute transition with credit shock with a larger spread
			- 'z' contains code to compute transition with productivity shock with a larger spread
		- 'operating_cost_final_good' contains code to generate results for Online Appendix A.5.2 with a fixed operating costs in unit of final good
			- 'calibration' contains code to calibrate model with a fixed operating costs in unit of final good
			- 'theta' contains code to compute transition with credit shock with fixed operating costs in unit of final good
			- 'z' contains code to compute transition with productivity shock with fixed operating costs in unit of final good
		- 'perfectly_elastic_entry' contains code to generate results for Online Appendix A.5.2 with perfectly elasitc mass of potential entrants
			- 'calibration' contains code to calibrate model with perfectly elasitc mass of potential entrants
			- 'theta' contains code to compute transition with credit shock with perfectly elastic mass of potential entrants
			- 'z' contains code to compute transition with credit shock with perfectly elastic mass of potential entrants
		- 'random_fixed_cost1' contains code to generate results for Online Appendix A.5.2 with random fixed operating cost (calibration 1)
			- 'calibration' contains code to calibrate model with random fixed operating cost (calibration 1)
			- 'theta' contains code to compute transition with credit shock with random fixed operating cost (calibration 1)
			- 'z' contains code to compute transition with productivity shock with random fixed operating cost (calibration 1)
		- 'random_fixed_cost2' contains code to generate results for Online Appendix A.5.2 with random fixed operating cost (calibration 2)
			- 'calibration' contains code to calibrate model with random fixed operating cost (calibration 2)
			- 'theta' contains code to compute transition with credit shock with random fixed operating cost (calibration 2)
			- 'z' contains code to compute transition with productivity shock with random fixed operating cost (calibration 2)

	- 'section_4_baseline_calibration' contains code to calibrate baseline model in Section 4 and results for Section 5 (except figures 4d-4e)

	- 'section_5_model_properties' contains code to generate results for Section 5 (figures 4d-4e) and Online Appendix
		- 'job_flows' contains code to compute job flows (job creation and job destruction) for Figures 4d-4e
			- 'benchmark' contains code to compute job flows (job creation and job destruction) for baseline
			- 'no_financial_friction' contains code to compute job flows for baseline without financial friction
			- 'no_irreversible_investment' contains code to compute job flows for baseline without irreversible investment
			- 'no_quadratic_adjustment' contains code to compute job flows for baseline without quadratic adjustment cost 
		- 'removing_frictions' contains code to generate results for Online Appendix A.1
			- 'no_financial_friction' contains code to compute baseline without financial frictions
			- 'no_irreversible_investment' contains code to compute baseline without irreversible investment
			- 'no_quadratic_adjustment' contains code to compute baseline without quadratic adjustment cost 

	- 'section_6.1_baseline_impulse_responses' contains code to generate results for Section 6.1 (Aggregate Shocks, firm entry, and exit)
		- 'baseline' contains code to generate results for impulse responses in baseline model in general equilibrium
			- 'theta' contains code to compute transition with credit shock in baseline model in general equilibrium
			- 'z' contains code to compute transition with productivity shock in baseline model in general equilibrium
		- 'partial_equilibrium' contains code to generate results for impulse responses in baseline model in partial equilibrium
			- 'theta' contains code to compute transition with credit shock in baseline model in partial equilibrium
			- 'z' contains code to compute transition with productivity shock in baseline model in partial equilibrium

	- 'section_6.2_Great_Recession' contains code to generate results for Section 6.2 (Great Recession)
		- 'great_recession_decomposition' contains code to decompose changes in output & hours due to differences in entry & exit (figure 12)
		- 'great_recession_theta' contains code to compute transition with only shock to credit in baseline model
		- 'great_recession_z' contains code to compute transition with only shock to productivity in baseline model
		- 'great_recession_z_and_theta' contains code to compute transition with both shocks to productivity and credit in baseline model

	- 'section_6.3_COVID' contains code to generate results for Section 6.3 (COVID-19 lockdown)
		- 'decomposition' contains code to decompose changes in output & hours due to differences in entry & exit 
		- 'no_operation' contains code to compute transition with only shock to firm operation in baseline model
		- 'no_operation_psi' contains code to compute transition with both shocks to firm operation and labor disutility in baseline model
		- 'psi' contains code to compute transition with only shock to labor disutility in baseline model


'Matlab_code_for_graphs' contains Matlab code used to generate figures

	- 'appendix' contains code to generate figures in Online Appendix
		- '1980_recession' contains code for figures 17 and 18
		- 'default' contains code for figure 20
		- 'life_cycle_properties_without_frictions' contains code for figure 15
		- 'other_model_specifications' contains code for figure 21
		- 'perfectly_elastic_entry' contains code for figure 16

	- 'section_2_data' contains code for figures 1 and 2 in Section 2

	- 'section_5_model_properties' contains code for figures 3 and 4 in Section 5

	- 'section_6.1_baseline_impulse_responses' contains code for figures 5-7 in Section 6.1 and figure 19 in Online Appendix
		- 'baseline' contains code for figures 5-6 in Section 6.1 and figure 19 in Online Appendix
		- 'partial_equilibrium' contains code for figure 7 in Section 6.1

	- 'section_6.2_Great_Recession' contains code for figures 8-12 in Section 6.2

	- 'section_6.3_COVID' contains code for figures 13-14 in Section 6.3

