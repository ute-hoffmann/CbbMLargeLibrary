aa_to_maxASA = {}
with open("Tien2013_theoretical.txt") as f:
	for line in f:
		aa_to_maxASA[line[0]] = line.strip("\n").split("\t")[1]

resPos_to_aa = {}
with open("Gallionella_CbbM-I_noStrep.txt") as f:
	for line in f:
		for i in range(len(line.strip("\n"))):
			resPos_to_aa[str(i+1)] = line.strip("\n")[i]

res_dict = {}
aa_relSAS = []
aa_Strep = []
abs_SAS = []
calc_relSAS = []
with open("measure_sasa_solventAccesibleSurface_wholeProtein.defattr") as f:
	for line in f:
		split_line = line.strip("\n").split("\t")
		if split_line[0] != "":
			continue
		else:
			number = split_line[1].split(":")[1]
			if number == "471":
				break
			res_dict[number] = split_line[2]
			aa_pos = resPos_to_aa[number] + str(number)
			Strep_pos = resPos_to_aa[number] + str(int(number)+14)
			aa_relSAS.append(aa_pos)
			aa_Strep.append(Strep_pos)
			abs_SAS.append(split_line[2])
			calc_relSAS.append(float(split_line[2])/float(aa_to_maxASA[resPos_to_aa[number]]))
			
with open("relSAS_Gallionella_CbbM-I.tsv", "w") as f:
	f.write("aa_pos\tbaseAA\tabsSAS\trelSAS\n")
	for i in range(len(aa_relSAS)):
		f.write(aa_relSAS[i] + "\t" + aa_Strep[i] + "\t" + str(abs_SAS[i]) + "\t" + str(calc_relSAS[i]) + "\n")
