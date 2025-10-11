resPos_to_aa = {}
with open("Gallionella_CbbM-I_noStrep.txt") as f:
	for line in f:
		for i in range(len(line.strip("\n"))):
			resPos_to_aa[i+1] = line.strip("\n")[i]

list_aa = []
Strep_aa = []
with open("dimer_interface.tsv") as f:
	for line in f:
		split_line = line.strip("\n").split(" ")	
		if split_line[0] != "residue":
			continue
		else:
			number = int(split_line[2].split(":")[1])
			aa_pos = resPos_to_aa[number] + str(number)
			Strep_pos = resPos_to_aa[number] + str(number+14)
			list_aa.append(aa_pos)
			Strep_aa.append(Strep_pos)
			
with open("dimer_aa.tsv", "w") as f:
	f.write("aa_pos\tbaseAA\tdimer\n")
	for i in range(len(list_aa)):
		f.write(list_aa[i] + "\t" + Strep_aa[i] + "\tyes\n")
