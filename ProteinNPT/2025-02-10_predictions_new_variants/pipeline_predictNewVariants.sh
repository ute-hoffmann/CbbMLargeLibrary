source config.sh  # To get $proteinnpt_data_path
#conda activate proteinnpt_env

source $proteinnpt_data_path/proteinnpt_env/bin/activate # Uncomment if using python venv instead of conda env

export assay_data_location="/home/admin/ProteinNPT_cbbM/2025-02-04_create_testSet/2025-02-09_predict_new_variants.csv"
export MSA_location="/home/admin/ProteinNPT_cbbM/input/TARGET_b0.3_I_.a2m"
export target_seq="MWSHPQFEKGSGSGSDQSNRYANLNLKEEDLIKNGKHLLVAYKLIPAKGHGFLEVAAHVAAESSTGTNVEVSTTDDFTRGVDALVYEIDETAFGDDPVKGGGLFKVAYPVELFDPNLTDGTYNISHMWSLILGNNQGMGDHQGLRMLDFLVPEMMVRKFDGPSANISNLWKVLGRPETDGGYIAGTIIKPKLGLRPEPFAKACYDFWLGGDFIKNDEPQANQPFCPMEVVMPKVAEAMDRAQQATGQAKLFSANITADYYKEMIHRGDFVLETFAKYNSASHVAFLVDGFVTGPAGVTTARREFPDTFLHFHRAGHGAVTSYKSPMGMDPLCYMKLVRLMGASGMHTGTMGYGKMEGHGKETVLAYMLERDECQGPYFYQKWYGMKATTPIISGGMNALRLPGFFQNLGHGNVINTCGGGAFGHIDSPAAGGISLRQAYDCWKSGSDPIEYAKTHKEFARAFESFPKDGDKLFAGWREKLGVHK"

python pipeline.py \
    --proteinnpt_data_location ${proteinnpt_data_path} \
    --assay_data_location ${assay_data_location} \
    --MSA_location ${MSA_location} \
    --target_seq ${target_seq} \
    --target_config_location ${target_config_location_multi_objectives_3} \
    --fold_variable_name train_test_split \
    --test_fold_index 1
