# Create directory
mkdir -p results/2023-04-20_inSilicoEvolution_I
mkdir -p results/2023-04-20_inSilicoEvolution_I/replication_3mutations_selectedPos

LOGFILE=results/2023-04-20_inSilicoEvolution_I/2023-04-24_selectedPos_logfile.txt

mkdir -p results/2023-04-20_inSilicoEvolution_I/replication_3mutations_selectedPos

date > $LOGFILE 2>&1 # Log start time

#conda init bash
#conda activate deep_sequence
for i in {0..99}; do
  for T in 0.03; do
      # Evolve based on difference of candidate and best sequence scores
      Outfile="results/2023-04-20_inSilicoEvolution_I/replication_3mutations_selectedPos/${i}_${T}_Difference.tab"
      python source/2023-04-24_DeepSequence_inSilicoEvolution_I_selectedPos.py \
        input/Gallionella_Rubisco_I.txt \
        $Outfile \
        -s 200 -t 2 -T $T \
        -a2m input/TARGET_b0.3_I_.a2m \
        -DS GallionellaRubisco_I \
        -o results/2023-04-20_inSilicoEvolution_I >> $LOGFILE 2>&1

  done
done
date >> $LOGFILE 2>&1 # Log stop time
