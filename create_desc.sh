python pdbbind.py -B 12 -F 'csv' -O CASF13_refine_sybyl_B12.csv -P /home/dat/WORK/DB/PDBbind/v2013-refined/ -D sybyl -I /home/dat/WORK/DB/PDBbind/v2013-refined/INDEX_refined_data.2013
python pdbbind.py -B 12 -F 'csv' -O CASF13_refine_elements_B12.csv -P /home/dat/WORK/DB/PDBbind/v2013-refined/ -D elements -I /home/dat/WORK/DB/PDBbind/v2013-refined/INDEX_refined_data.2013
python pdbbind.py -B 4.5 -F 'csv' -O CASF13_refine_credo.csv -P /home/dat/WORK/DB/PDBbind/v2013-refined/ -D credo -I /home/dat/WORK/DB/PDBbind/v2013-refined/INDEX_refined_data.2013

python pdbbind.py -B 12 -F 'csv' -O CASF13_test_sybyl_B12.csv -P /home/dat/WORK/DB/PDBbind/v2013-core/ -D sybyl -I /home/dat/WORK/DB/PDBbind/v2013-core/INDEX_core_data.2013
python pdbbind.py -B 12 -F 'csv' -O CASF13_test_elements_B12.csv -P /home/dat/WORK/DB/PDBbind/v2013-core/ -D elements -I /home/dat/WORK/DB/PDBbind/v2013-core/INDEX_core_data.2013
python pdbbind.py -B 4.5 -F 'csv' -O CASF13_test_credo.csv -P /home/dat/WORK/DB/PDBbind/v2013-core/ -D credo -I /home/dat/WORK/DB/PDBbind/v2013-core/INDEX_core_data.2013

python pdbbind.py -B 12 -F 'csv' -O CASF12_training_sybyl.csv -P /home/dat/WORK/DB/PDBbind/v2012-refined/ -D sybyl -I /home/dat/WORK/DB/PDBbind/v2012/INDEX_refined_data.2012
python pdbbind.py -B 12 -F 'csv' -O CASF12_training_elements.csv -P /home/dat/WORK/DB/PDBbind/v2012-refined/ -D elements -I /home/dat/WORK/DB/PDBbind/v2012/INDEX_refined_data.2012
python pdbbind.py -F 'csv' -O CASF12_training_elements_B0.csv -P /home/dat/WORK/DB/PDBbind/v2012-refined/ -D elements -I /home/dat/WORK/DB/PDBbind/v2012/INDEX_refined_data.2012
python pdbbind.py -B 4.5 -F 'csv' -O CASF12_training_credo.csv -P /home/dat/WORK/DB/PDBbind/v2012-refined/ -D credo -I /home/dat/WORK/DB/PDBbind/v2012/INDEX_refined_data.2012

python pdbbind.py -B 12 -F 'csv' -O CASF12_test_sybyl.csv -P /home/dat/WORK/DB/PDBbind/v2012-core/ -D sybyl -I /home/dat/WORK/DB/PDBbind/v2012/INDEX_core_data.2012
python pdbbind.py -B 12 -F 'csv' -O CASF12_test_elements.csv -P /home/dat/WORK/DB/PDBbind/v2012-core/ -D elements -I /home/dat/WORK/DB/PDBbind/v2012/INDEX_core_data.2012
python pdbbind.py -B 4.5 -F 'csv' -O CASF12_test_credo.csv -P /home/dat/WORK/DB/PDBbind/v2012-core/ -D credo -I /home/dat/WORK/DB/PDBbind/v2012/INDEX_core_data.2012

python pdbbind.py -B 2 -O CASF07_test_elements.csv -D elements -P /home/dat/WORK/DB/PDBbind/v2007/ -I /home/dat/WORK/DB/PDBbind/v2007/INDEX.2007.core.data
python pdbbind.py --debug -D elements -O CASF13_test_elements.csv -P /home/dat/WORK/DB/PDBbind/v2013-core/  -I /home/dat/WORK/DB/PDBbind/v2013-core/INDEX_core_data.2013
python pdbbind.py --debug -D elements -O CASF12_test_elements.csv -P /home/dat/WORK/DB/PDBbind/v2012-core/ -I /home/dat/WORK/DB/PDBbind/v2012/INDEX_core_data.2012
