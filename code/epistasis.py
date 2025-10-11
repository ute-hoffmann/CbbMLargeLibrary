# adapted from https://onlinelibrary.wiley.com/doi/full/10.1002/pro.2876

# read in single exchanges that exist

import pandas as pd

fitness_data = pd.read_csv("../WeightingStrategies_fitnessCalc/result_lfcSE_Pearson_Weighted_sameWeight.tsv", sep="\t")
satu_combi = pd.read_csv("../input/combi_satur_mutatios.tsv", names = ["sgRNA_target", "category"], sep="\t")
fitness_data = fitness_data.merge(satu_combi, on="sgRNA_target", how="left")

# read in fitness values
conditions = set(fitness_data["condition"])
dict_cond = {}
for cond in conditions:
	subset_fitness_data = fitness_data[fitness_data["condition"]==cond]
	
	single_amino_exchanges = subset_fitness_data[subset_fitness_data["category"]=="combiANDsatur"][["sgRNA_target", "norm"]].drop_duplicates()
	single_amino_exchanges["norm"] = single_amino_exchanges["norm"]
	
	combinatorial_variants = subset_fitness_data[subset_fitness_data["category"]=="combinatorial"][["sgRNA_target", "norm"]].drop_duplicates()
	
	dict_cond[cond] = {"satu":single_amino_exchanges, "combi": combinatorial_variants}
	
# somehow create network, systematic interrogation

with open("../lfcSE_weighted_outputImages/epistasis/table_values.tsv", "w") as f:
	f.write("Condition\tExchange\tExchange_value_and_ratio_in_WT\tVariant\tVariant_value\tVariant_minusOne\tVariant_minusOne_value\tVariant_ratio\n")
def create_variant_and_extract_value(list_of_values):
	variant_string = ""
	for pos in list_of_values:
		if pos:
			if variant_string:
				variant_string += "," + pos
			else:
				variant_string = pos
		try:
			value_float = float(combinatorial_variants[combinatorial_variants["sgRNA_target"]==variant_string]["norm"])
		except:
			value_float = None
			continue
	return(variant_string, value_float)
	
def write_values_in_file(condition_name, subst, subst_WT_value, variant, variant_value, variant_minus_one, variant_minus_one_value):
	# if one of the variant values is negative or close to 0 (<0.05), set it to 0.1 to avoid negative variant ratios; note: K214H has normed fitness values of approx. 0.2 to 0.3 in different conditions
	if variant_value < 0.1:
		variant_value = 0.1
	if variant_minus_one_value < 0.1:
		variant_minus_one_value = 0.1
		
	# calculate ratio 
	variant_improve_ratio = variant_value/variant_minus_one_value
	
	with open("../lfcSE_weighted_outputImages/epistasis/table_values.tsv", "a") as f:
		f.write(condition_name + "\t" + subst + "\t" + str(subst_WT_value) + "\t" + variant + "\t" + str(variant_value) + "\t" + variant_minus_one + "\t" + str(variant_minus_one_value) + "\t" + str(variant_improve_ratio) + "\n")
	return
	
for cond in dict_cond.keys():
	print(cond)
	single_amino_exchanges = dict_cond[cond]["satu"]
	combinatorial_variants = dict_cond[cond]["combi"]
	number_variants = 0
	for H141 in ["", "V", "L"]:
		if H141:
			H141var = "H141" + H141
			try:
				H141val = float(single_amino_exchanges[single_amino_exchanges["sgRNA_target"]==H141var]["norm"])
			except:
				H141val = None
				print("no value obtained for " + H141var)
		else: 
			H141var = ""
			H141val = None
		for M154 in ["", "E", "D", "K", "S", "A"]:
			if M154:
				M154var = "M154" + M154
				try:
					M154val = float(single_amino_exchanges[single_amino_exchanges["sgRNA_target"]==M154var]["norm"])
				except:
					M154val = None
					print("no value obtained for " + M154var)
			else:
				M154var = ""
				M154val = None
			for K261 in ["", "D", "A", "F", "E"]:
				if K261:
					K261var = "K261" + K261
					try:
						K261val = float(single_amino_exchanges[single_amino_exchanges["sgRNA_target"]==K261var]["norm"])
					except:
						K261val = None
						print("no value obtained for " + K261var)
				else:
					K261var = ""
					K261val = None
				for H265 in ["", "K", "A", "R", "E"]:
					if H265:	
						H265var = "H265" + H265
						try:
							H265val = float(single_amino_exchanges[single_amino_exchanges["sgRNA_target"]==H265var]["norm"])
						except:
							H265val = None
							print("no value obtained for " + H265var)
					else:
						H265var = ""
						H265val = None
					for V337 in ["", "A", "S"]:
						if V337:
							V337var = "V337" + V337
							try:
								V337val = float(single_amino_exchanges[single_amino_exchanges["sgRNA_target"]==V337var]["norm"])
							except:
								V337val = None
								print("no value obtained for " + V337var)
						else:
							V337var = ""
							V337val = None
						for Q406 in ["", "E", "D"]:
							if Q406:
								Q406var = "Q406" + Q406
								try:
									Q406val = float(single_amino_exchanges[single_amino_exchanges["sgRNA_target"]==Q406var]["norm"])
								except:
									Q406val = None
									print("no value obtained for " + Q406var)
							else:
								Q406var = ""
								Q406val = None
							for S446 in ["", "A", "I", "T"]:
								if S446:
									S446var = "S446" + S446
									try:
										S446val = float(single_amino_exchanges[single_amino_exchanges["sgRNA_target"]==S446var]["norm"])
									except:
										S446val = None
										print("no value obtained for " + S446var)
								else:
									S446var = ""
									S446val = None
								string_value = create_variant_and_extract_value([H141var, M154var, K261var, H265var, V337var, Q406var, S446var])
								try:
									a = string_value[1]
									number_variants += 1
								except:
									continue
								if number_variants % 1000 == 0:
									print(number_variants)
								if H141val:
									try:
										string_value_minusOne = create_variant_and_extract_value([M154var, K261var, H265var, V337var, Q406var, S446var])
										write_values_in_file(cond, H141var, H141val, string_value[0], string_value[1], string_value_minusOne[0], string_value_minusOne[1])
									except:
										continue
								if M154val:
									try:
										string_value_minusOne = create_variant_and_extract_value([H141var, K261var, H265var, V337var, Q406var, S446var])
										write_values_in_file(cond, M154var, M154val, string_value[0], string_value[1], string_value_minusOne[0], string_value_minusOne[1])
									except:
										continue
								if K261val:
									try:
										string_value_minusOne = create_variant_and_extract_value([H141var, M154var, H265var, V337var, Q406var, S446var])
										write_values_in_file(cond, K261var, K261val, string_value[0], string_value[1], string_value_minusOne[0], string_value_minusOne[1])
									except:
										continue
								if H265val:
									try:
										string_value_minusOne = create_variant_and_extract_value([H141var, M154var, K261var, V337var, Q406var, S446var])
										write_values_in_file(cond, H265var, H265val, string_value[0], string_value[1], string_value_minusOne[0], string_value_minusOne[1])
									except:
										continue
								if V337val:
									try:
										string_value_minusOne = create_variant_and_extract_value([H141var, M154var, K261var, H265var, Q406var, S446var])
										write_values_in_file(cond, V337var, V337val, string_value[0], string_value[1], string_value_minusOne[0], string_value_minusOne[1])
									except:	
										continue
								if Q406val:
									try:
										string_value_minusOne = create_variant_and_extract_value([H141var, M154var, K261var, H265var, V337var, S446var])
										write_values_in_file(cond, Q406var, Q406val, string_value[0], string_value[1], string_value_minusOne[0], string_value_minusOne[1])
									except:
										continue
								if S446val:
									try:
										string_value_minusOne = create_variant_and_extract_value([H141var, M154var, K261var, H265var, V337var, Q406var])
										write_values_in_file(cond, S446var, S446val, string_value[0], string_value[1], string_value_minusOne[0], string_value_minusOne[1])
									except:
										continue
