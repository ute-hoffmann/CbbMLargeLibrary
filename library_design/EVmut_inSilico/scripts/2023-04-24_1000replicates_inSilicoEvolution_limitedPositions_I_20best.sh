# Create directory
mkdir -p results/I_replication_4mutations_limitedPos

# Iterate over 201 replicates
for i in {0..999}; do
  # Iterate over one temperature
  for T in 0.03; do
      # Evolve based on difference of candidate and best sequence scores
      Outfile="results/I_replication_4mutations_limitedPos/${i}_${T}_Difference.tab"
      python scripts/2023-01-13_in_silico_evolution_zeroShot_EVcouplings_limitedPositions.py \
        -s 500 -t 3 -T $T \
        -EV input/I_ab368e89ebe9421898840b5c014db5b3_TARGET_b0.3.model \
        -i input/Gallionella_Rubisco_I.txt \
        -o $Outfile -posInt input/2023-04-24_I_positions_of_interest.txt
  done
done

