#########################################
#										#
# Stefano Petrucci - s.petrucci@cern.ch	#
#										#
#########################################

#############################################
#											#
# USAGE: python this-script.py yyyy_mm_dd	#
#											#
#############################################

# This script reads every file with .dat extension inside the folder

import sys
import glob
import os
import re
import ntpath

#file_list = glob.glob('/home/petruccs-local/scripts/'+sys.argv[1]+'/*.dat')
file_list = glob.glob('./'+sys.argv[1]+'/*.dat')
initial_time = -1
final_time = 0
temp_final_time = 0
values_dict = {}
allow_dict = {}
# root dictionary to create the header for the root function "ReadFile". This scrip uses the filename (without the file extension) as dictionary key 
root_branch_dict = {
					'bi_31832_bctdc'			: 'bctdc3/F:bctdc3_yerrors',
					'bi_41435_bctdc'			: 'bctdc4/F:bctdc4_yerrors',
					'bi_51895_bctdc'			: 'bctdc51895/F:bctdc51895_yerrors',
					'bi_51897_bctdc'			: 'bctdc51897/F:bctdc51897_yerrors',
					'bi_31450_bctfr'			: 'bctfr3/F:bctfr3_yerrors:bctfr3_fast_intensity',
					'bi_51895_bctfr'			: 'bctfr5/F:bctfr5_yerrors:bctfr5_fast_intensity',
					'bshv_51793'				: 'bshv_51793_controllers/F:bshv_51793_lvdt:bshv_51793_res',
					'sps_bpmcol_box1_hor_up'	: 'bpm_box1_hor_up/F',
					'sps_bpmcol_box1_hor_down'	: 'bpm_box1_hor_down/F',
					'tacw_51998'				: 'tacw_51998_pos/F',
					'tal_52196_lin1'			: 'tal_52196_lin1_controllers/F:tal_52196_lin1_lvdt:tal_52196_lin1_res',
					'tal_52196_lin2'			: 'tal_52196_lin2_controllers/F:tal_52196_lin2_lvdt:tal_52196_lin2_res',
					'tal_52196_lin3'			: 'tal_52196_lin3_controllers/F:tal_52196_lin3_lvdt:tal_52196_lin3_res',
					'tcpc_51595_bp'				: 'tcpc_51595_bp_controllers/F:tcpc_51595_bp_lvdt:tcpc_51595_bp_res',
					'tcpc_51595_rot'			: 'tcpc_51595_rot_controllers/F:tcpc_51595_rot_lvdt:tcpc_51595_rot_res',
					'tcpc_51795_bp'				: 'tcpc_51795_bp_controllers/F:tcpc_51795_bp_lvdt:tcpc_51795_bp_res',
					'tcpc_51795_lin'			: 'tcpc_51795_lin_controllers/F:tcpc_51795_lin_lvdt:tcpc_51795_lin_res',
					'tcpc_51795_rot'			: 'tcpc_51795_rot_controllers/F:tcpc_51795_rot_lvdt:tcpc_51795_rot_res',
					'tcsm_51934'				: 'tcsm_51934_lvdt_ld/F:tcsm_51934_lvdt_lu:tcsm_51934_lvdt_rd:tcsm_51934_lvdt_ru',
					'tcxhw_51651_h1'			: 'tcxhw_51651_h1_controllers/F:tcxhw_51651_h1_lvdt:tcxhw_51651_h1_res',
					'tcxhw_51651_h2'			: 'tcxhw_51651_h2_controllers/F:tcxhw_51651_h2_lvdt:tcxhw_51651_h2_res',
					'tecs_51797'				: 'tecs_51797_a_position_ctrl/F:tecs_51797_a_position_lvdt:tecs_51797_b_position_ctrl:tecs_51797_b_position_lvdt:tecs_51797_a_angle_ctrl:tecs_51797_a_angle_lvdt:tecs_51797_b_angle_ctrl:tecs_51797_b_angle_lvdt',
					'tecs_51797_lindw'			: 'tecs_51797_lindw_controllers/F:tecs_51797_lindw_lvdt:tecs_51797_lindw_res',
					'tecs_51797_linup'			: 'tecs_51797_linup_controllers/F:tecs_51797_linup_lvdt:tecs_51797_linup_res',
					'tecs_51799'				: 'tecs_51799_a_position_ctrl/F:tecs_51799_a_position_lvdt:tecs_51799_b_position_ctrl:tecs_51799_b_position_lvdt:tecs_51799_a_angle_ctrl:tecs_51799_a_angle_lvdt:tecs_51799_b_angle_ctrl:tecs_51799_b_angle_lvdt',
					'tecs_51799_lindw'			: 'tecs_51799_lindw_controllers/F:tecs_51799_lindw_lvdt:tecs_51799_lindw_res',
					'tecs_51799_linup'			: 'tecs_51799_linup_controllers/F:tecs_51799_linup_lvdt:tecs_51799_linup_res',
					'tqcd_201271'				: 'tqcd_201271_controllers/F:tqcd_201271_lvdt:tqcd_201271_res',
					'tqcd_51794'				: 'tqcd_51794_controllers/F:tqcd_51794_lvdt:tqcd_51794_res',
					'tqcd_51991'				: 'tqcd_51991_controllers/F:tqcd_51991_lvdt:tqcd_51991_res',
					'xrph_51937_h1'				: 'xrph_51937_h1_controllers/F:xrph_51937_h1_lvdt:xrph_51937_h1_res',
					'xrph_51937_h2'				: 'xrph_51937_h2_controllers/F:xrph_51937_h2_lvdt:xrph_51937_h2_res',
					'xrph_52202_h1'				: 'xrph_52202_h1_controllers/F:xrph_52202_h1_lvdt:xrph_52202_h1_res',
					'xrph_52202_h2'				: 'xrph_52202_h2_controllers/F:xrph_52202_h2_lvdt:xrph_52202_h2_res',
					'A_counters'				: 'A_counters/F:A_counters_yerrors',
					'A_frev_counters'			: 'A_frev_counters/F:A_frev_counters_yerrors',
					'AA_counters'				: 'AA_counters/F:AA_counters_yerrors',
					'AA_frev_counters'			: 'AA_frev_counters/F:AA_frev_counters_yerrors',
					'AD_frev_counters'			: 'AD_frev_counters/F:AD_frev_counters_yerrors',
					'AD_I_counters'				: 'AD_I_counters/F:AD_I_counters_yerrors',
					'AD_I_frev_counters'		: 'AD_I_frev_counters/F:AD_I_frev_counters_yerrors',
					'AH_AI_counters'			: 'AH_AI_counters/F:AH_AI_counters_yerrors',
					'AH_AI_frev_counters'		: 'AH_AI_frev_counters/F:AH_AI_frev_counters_yerrors',
					'AH_frev_counters'			: 'AH_frev_counters/F:AH_frev_counters_yerrors',
					'AI_frev_counters'			: 'AI_frev_counters/F:AI_frev_counters_yerrors',
					'C_D_counters'				: 'C_D_counters/F:C_D_counters_yerrors',
					'C_D_frev_counters'			: 'C_D_frev_counters/F:C_D_frev_counters_yerrors',
					'C_frev_counters'			: 'C_frev_counters/F:C_frev_counters_yerrors',
					'D_frev_counters'			: 'D_frev_counters/F:D_frev_counters_yerrors',
					'E_frev_counters'			: 'E_frev_counters/F:E_frev_counters_yerrors',
					'E_N_counters'				: 'E_N_counters/F:E_N_counters_yerrors',
					'E_N_frev_counters'			: 'E_N_frev_counters/F:E_N_frev_counters_yerrors',
					'F_counters'				: 'F_counters/F',
					'frev_counters'				: 'frev_counters/F:frev_counters_yerrors',
					'frev1_counters'			: 'frev1_counters/F:frev1_counters_yerrors',
					'G_counters'				: 'G_counters/F:G_counters_yerrors',
					'G_frev_counters'			: 'G_frev_counters/F:G_frev_counters_yerrors',
					'G_M_counters'				: 'G_M_counters/F:G_M_counters_yerrors',
					'G_M_frev_counters'			: 'G_M_frev_counters/F:G_M_frev_counters_yerrors',
					'H_frev_counters'			: 'H_frev_counters/F:H_frev_counters_yerrors',
					'I_frev_counters'			: 'I_frev_counters/F:I_frev_counters_yerrors',
					'J_frev_counters'			: 'J_frev_counters/F:J_frev_counters_yerrors',
					'J_K_counters'				: 'J_K_counters/F:J_K_counters_yerrors',
					'J_K_frev_counters'			: 'J_K_frev_counters/F:J_K_frev_counters_yerrors',
					'K_frev_counters'			: 'K_frev_counters/F:K_frev_counters_yerrors',
					'M_counters'				: 'M_counters/F:M_counters_yerrors',
					'M_frev_counters'			: 'M_frev_counters/F:M_frev_counters_yerrors',
					'N_frev_counters'			: 'N_frev_counters/F:N_frev_counters',
					'Z_counters'				: 'Z_counters/F:Z_counters_yerrors',
					'Z_frev_counters'			: 'Z_frev_counters/F:Z_frev_counters_yerrors',
					'Z_H_frev_counters'			: 'Z_H_frev_counters/F:Z_H_frev_counters_yerrors',
					'CpFM1_counters'			: 'CpFM1_counters/F',
					'CpFM2_counters'			: 'CpFM2_counters/F',
					'CpFM1_CpFM2_counters'		: 'CpFM1_CpFM2_counters/F',
					'Trigger_counters'			: 'CpFM_trigger_counters/F',
					'medipix'					: 'medipix0_counters/F:medipix1_counters',
					'cpfm_tt20'					: 'cpfm_tt20_yerrors/F:cpfm_tt20_counters',
					'lss2_bump'					: 'lss2_bump/F',
					'tecs_51652'				: 'tecs_51652_a_position_ctrl/F:tecs_51652_a_position_lvdt:tecs_51652_b_position_ctrl:tecs_51652_b_position_lvdt:tecs_51652_a_angle_ctrl:tecs_51652_a_angle_lvdt:tecs_51652_b_angle_ctrl:tecs_51652_b_angle_lvdt',
					'tecs_51652_lindw'			: 'tecs_51652_lindw_controllers/F:tecs_51652_lindw_lvdt:tecs_51652_lindw_res',
					'tecs_51652_linup'			: 'tecs_51652_linup_controllers/F:tecs_51652_linup_lvdt:tecs_51652_linup_res'
}

# finding the initial time and the final time reading from the entire set of file
j = 0
for i in range(0, len(file_list)):
	with open(file_list[i]) as data_file:
		print file_list[i]
		j = 0
		for row in data_file:
			if float(row.split(' ')[0]) > 0 and initial_time < float(row.split(' ')[0]) and j == 0:
				initial_time = int(float(row.split(' ')[0]))
				j = 1
			if float(row.split(' ')[0]) < 0:
				j = 0
			temp_final_time = int(float(row.split(' ')[0]))
		if temp_final_time < final_time and i > 0:
			final_time = temp_final_time
		if i == 0:
			final_time = temp_final_time
print "--> initial_time = ", initial_time, " | final_time = ", final_time
my_list = []
# filling the addresses with empty strings
for k in range(initial_time, final_time + 1):
	values_dict[k] = []
	my_list = [str(k)]
	values_dict[k].extend(my_list)
	allow_dict[k] = 1

old_time = initial_time - 1
old_data = ''
data_len = 0
header = ['t/I']
missing_values = 0
for i in range(0, len(file_list)):
	with open(file_list[i]) as data_file:
		print re.sub(sys.argv[1]+'_', '', os.path.splitext(ntpath.basename(file_list[i]))[0])		
		header.append(root_branch_dict[re.sub(sys.argv[1]+'_', '', os.path.splitext(ntpath.basename(file_list[i]))[0])])		
		for row in data_file:
			time = int(round(float(row.split(' ')[0])))
			data = row.split(' ')
			data[-1] = data[-1].strip()
			data.pop(0)
			if time < initial_time:
				old_data = data
			if (time >= initial_time and time <= final_time and time != old_time and allow_dict[time]):
				if time == old_time + 1:
					values_dict[time].extend(data)					
				if time != old_time + 1:
					for jj in range(old_time + 1, time + 1):
						if allow_dict[jj]:
							values_dict[jj].extend(old_data)
							missing_values += 1
							allow_dict[jj] = 0
				old_time = time
				old_data = data
				allow_dict[time] = 0
				missing_values = 0
		for k in range(initial_time, final_time):
			allow_dict[k] = 1

with open(os.path.dirname(file_list[i])+'/aligned_data.txt', 'w') as output_file:
	output_file.write(':'.join(header))
	output_file.write('\n')

	for j in range(initial_time, final_time):
		output_file.write(' '.join(values_dict[j]))
		output_file.write('\n')
