# load original sequence
cbbM_seq = ""
with open("Gallionella_Rubisco_I.txt") as f:
	for line in f:
		cbbM_seq = line

seq_to_name_dict = {}
with open("mutant_code.tsv") as f:
	header = True
	for line in f:
		if header:
			header = False
			continue
		seq_to_name_dict[line.split("\t")[0]] = line.strip("\n").split("\t")[1]

with open("MoCHI_predictions.tsv", "w") as pred_f:
	pred_f.write("Condition\tname\tprediction_mean\tprediction_std\n")
for condition in ["CLN2", "CLO2", "LD"]:
	list_written = []
	with open(condition + "/mochi_project/task_1/predictions/predicted_phenotypes_supp.txt") as f:
		header = True
		for line in f:
			if header: 
				header = False
				continue
			split_line = line.split("\t")
			seq = split_line[0]
			name = seq_to_name_dict[seq]
			if name in list_written:
				continue
			prediction_mean = split_line[17]
			prediction_std = split_line[18]
			if condition == "CLN2":
				cond = "CL_N2"
			elif condition == "CLO2":
				cond = "CL_O2"
			elif condition == "LD":
				cond = "LD"
			with open("MoCHI_predictions.tsv", "a") as pred_f:
				pred_f.write(cond + "\t" + name + "\t" + prediction_mean + "\t" + prediction_std + "\n")
				list_written.append(name)
				

