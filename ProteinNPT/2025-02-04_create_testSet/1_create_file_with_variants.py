dict_variants = {}
with open("list_variants.txt") as comb_variants:
	mutation_counter = 0
	for line in comb_variants:
		dict_variants[mutation_counter] = line.strip("\n").split("\t")
		mutation_counter += 1

def add_mutations(list_of_created_variants, new_variants):
	list_with_new_variants = []
	new_variants.append("")
	for item in list_of_created_variants:
		for variant in new_variants:
			if item != "" and variant != "":
				list_with_new_variants.append(item + ":" + variant)
			elif item == "" and variant != "": 
				list_with_new_variants.append(variant)
			else:
				list_with_new_variants.append(item)
	list_of_created_variants.append(list_with_new_variants)
	return list_with_new_variants

list_variants = [""]
for mut in range(mutation_counter):
	list_variants = add_mutations(list_variants[:], dict_variants[mut][:])[:]
	
print("in total " + str(len(set(list_variants))) + " variants")

with open("list_combinations.txt", "w") as table_variants:		
	for var in set(list_variants):
		table_variants.write(var + "\n")
		
