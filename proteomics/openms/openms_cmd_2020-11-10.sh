base_directory=~/proteomic_data
mzml_files=proteomic_data/wiff_files/mzml_peakpicker
pwd

echo $mzml_files

FileMerger -in *run1*.mzML* -out ../fractions_merged_run1_2020_11_02.mzML -threads 4
FileMerger -in *run2*.mzML* -out ../fractions_merged_run2_2020_11_02.mzML -threads 4

FileFilter -in wiff_files/mzml_peakpicker/fractions_merged_run1_2020_11_02.mzML -sort -out openms_cmd_outputs/filefilter_rtsort_run1.mzML -threads 4
FileFilter -in wiff_files/mzml_peakpicker/fractions_merged_run2_2020_11_02.mzML -sort -out openms_cmd_outputs/filefilter_rtsort_run2.mzML -threads 4

PercolatorAdapter -in openms_cmd_outputs/psmfeature_run1.idXML -out openms_cmd_outputs/percolator_run1.idXML -percolator_executable ~/Development/OpenMS/THIRDPARTY/Linux/64bit/Percolator/percolator -score_type q-value -enzyme trypsinp -threads 8 -verbose 5
PercolatorAdapter -in openms_cmd_outputs/psmfeature_run2.idXML -out openms_cmd_outputs/percolator_run2.idXML -percolator_executable ~/Development/OpenMS/THIRDPARTY/Linux/64bit/Percolator/percolator -score_type q-value -enzyme trypsinp -threads 8 -verbose 5

IDFilter -in openms_cmd_outputs/percolator_run1.idXML -out openms_cmd_outputs/idfilter1_run1.idXML -score:pep 0.05 -threads 4
IDFilter -in openms_cmd_outputs/percolator_run2.idXML -out openms_cmd_outputs/idfilter1_run2.idXML -score:pep 0.05 -threads 4

IsobaricAnalyzer -type itraq4plex -in openms_cmd_outputs/filefilter_rtsort_run1.mzML -out openms_cmd_outputs/isoanalse_run1.consensusXML -threads 6
IsobaricAnalyzer -type itraq4plex -in openms_cmd_outputs/filefilter_rtsort_run2.mzML -out openms_cmd_outputs/isoanalse_run2.consensusXML -threads 6

IDMapper -id openms_cmd_outputs/idfilter1_run1.idXML -in openms_cmd_outputs/isoanalse_run1.consensusXML -out openms_cmd_outputs/idmap_run1.consensusXML -rt_tolerance 0.1 -mz_reference precursor -feature:use_centroid_mz false -threads 6
IDMapper -id openms_cmd_outputs/idfilter1_run2.idXML -in openms_cmd_outputs/isoanalse_run2.consensusXML -out openms_cmd_outputs/idmap_run2.consensusXML -rt_tolerance 0.1 -mz_reference precursor -feature:use_centroid_mz false -threads 6

FileMerger -in openms_cmd_outputs/idmap_run1.consensusXML openms_cmd_outputs/idmap_run2.consensusXML -out openms_cmd_outputs/idmap_merged.consensusXML -annotate_file_origin -threads 4

Epifany -in openms_cmd_outputs/idmap_merged.consensusXML -out openms_cmd_outputs/epifany_2020-11-03.consensusXML -greedy_group_resolution remove_proteins_wo_evidence -algorithm:keep_best_PSM_only false

IDFilter -in openms_cmd_outputs/epifany_2020-11-03.consensusXML -out openms_cmd_outputs/idfilter2.consensusXML -score:protgroup 0.05 -remove_decoys -threads 4

IDConflictResolver -in openms_cmd_outputs/idfilter2.consensusXML -out openms_cmd_outputs/idconflictresolve_merged.consensusXML -threads 4

MSstatsConverter -in openms_cmd_outputs/idconflictresolve_merged.consensusXML -in_design experimental_design/fractions_merged_both_runs_2020-11-03.tsv -method ISO -out openms_cmd_outputs/workflow_results/msstatstmt_2020-11-03.csv

PercolatorAdapter -in psmfeature_example.mzid -out percadap_test.idXML -percolator_executable /home/mateus/Development/OpenMS/THIRDPARTY/Linux/64bit/Percolator/percolator -score_type pep -threads 8

MSstatsConverter -in proteomic_data/filefilter_run2.consensusXML -in_design proteomic_data/knime_ExperimentalDesign_fractions_merged.tsv -out proteomic_data/knime_outputs/run2_mstats_cmd.csv -method ISO

# Note that the docker damon must be running before docker can be used
# If it isn't running already it can be started with:
sudo systemctl start docker

sudo docker run -it --rm -v /home/mateus/proteomic_data:/proteomic_data proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /proteomic_data/Mateus_iTRAQ_run1_Fr01_4ul.wiff --filter "peakPicking true 1-" -o /proteomic_data/mzml_peakpicker

DecoyDatabase -enzyme Trypsin -in ~/proteomic_data/uniprot-proteome_UP000005640_human_reference_2020-10-01.fasta ~/proteomic_data/P00761_pig_trypsin_2020-10-01.fasta -out ~/proteomic_data/test_decoy_database.fasta

MSGFPlusAdapter -in /home/mateus/proteomic_data/openms_cmd_outputs/filefilter_rtsort_run1.mzML -out openms_cmd_outputs/msgfplusadapter_search_run1_2020_11_02.idXML -executable ~/Development/OpenMS/THIRDPARTY/All/MSGFPlus/MSGFPlus.jar -database ~/proteomic_data/test_decoy_database.fasta -fixed_modifications "Methylthio (C)" "iTRAQ4plex (N-term)" -variable_modifications "Oxidation (M)" -precursor_mass_tolerance 20 -enzyme Trypsin/P -protocol iTRAQ -instrument high_res -java_memory 15000 -threads 8

MSGFPlusAdapter -in /home/mateus/proteomic_data/openms_cmd_outputs/filefilter_rtsort_run2.mzML -out openms_cmd_outputs/msgfplusadapter_search_run2_2020_11_02.idXML -executable ~/Development/OpenMS/THIRDPARTY/All/MSGFPlus/MSGFPlus.jar -database ~/proteomic_data/test_decoy_database.fasta -fixed_modifications "Methylthio (C)" "iTRAQ4plex (N-term)" -variable_modifications "Oxidation (M)" -precursor_mass_tolerance 20 -enzyme Trypsin/P -protocol iTRAQ -instrument high_res -java_memory 15000 -threads 8

PeptideIndexer -in openms_cmd_outputs/msgfplusadapter_search_run1_2020_11_02.idXML -out openms_cmd_outputs/pepindex_msgfplus_run1.idXML
PeptideIndexer -in openms_cmd_outputs/msgfplusadapter_search_run2_2020_11_02.idXML -out openms_cmd_outputs/pepindex_msgfplus_run2.idXML

PSMFeatureExtractor -in openms_cmd_outputs/pepindex_msgfplus_run1.idXML -out openms_cmd_outputs/psmfeature_run1.idXML
PSMFeatureExtractor -in openms_cmd_outputs/pepindex_msgfplus_run2.idXML -out openms_cmd_outputs/psmfeature_run2.idXML

CometAdapter -in ~/Documents/openms_example_data/Example_Data/B1.mzML -out openms_example_data.idXML -comet_executable ~/Development/OpenMS/THIRDPARTY/Linux/64bit/Comet/comet.exe -database ~/Documents/openms_example_data/Example_Data/iPRG2016_shuff.fasta -precursor_mass_tolerance 20 -enzyme Trypsin -missed_cleavages 1 -fragment_mass_tolerance 0.01 -fixed_modifications "Carbamidomethyl (C)" -variable_modifications "Oxidation (M)" -threads 8

IDMerger -in openms_example_data.idXML -out idmerger_example_test.idXML -annotate_file_origin -merge_proteins_add_PSMs

PeptideIndexer -in idmerger_example_test.idXML -out pepindex_example.idXML

PSMFeatureExtractor -in pepindex_example.idXML -out psmfeature_example.mzid
