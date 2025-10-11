with open("barcode_CDS_1occ.tsv") as f:
	for line in f:
		barcode = line.split("\t")[0]
		if barcode == "TGATGAGGACTCGTTACCTG" or barcode == "ATTACCGTTGATCCGACATG" or barcode == "ACCTGCAATGGAGAACCACC":
			print(line.split("\t")[1])
