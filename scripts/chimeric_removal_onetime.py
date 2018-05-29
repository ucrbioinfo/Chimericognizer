import os.path
from shutil import copyfile

from chimeric import chimeric
from chimeric_cutting import cutting
from chimeric_cutting_ref import cutting_ref

def chimeric_removal_onetime(REFALIGNER, DIGEST, chimeric_sub_dir, optmap_list, optmap_type_list, fasta_file_input, fasta_file_output, num_threads, threshold_1, threshold_2, threshold_3, threshold_4, threshold_5, threshold_6):
    refaligner_dir = chimeric_sub_dir + "/refaligner"
    refaligner_prefix = "input_contigs"
    fasta_file_in_refaligner = refaligner_dir + "/" + refaligner_prefix + ".fasta"

    if not os.path.isdir(refaligner_dir):
        os.makedirs(refaligner_dir)
    copyfile(fasta_file_input, fasta_file_in_refaligner)
    
    out_list = []
    out2_list = []
    for i in range(len(optmap_list)):
        optmap = optmap_list[i]
        optmap_type = optmap_type_list[i]
        silicomap = refaligner_dir + "/" + refaligner_prefix + "_" + optmap_type + ".cmap"
        out = refaligner_dir + "/" + refaligner_prefix + "_" + optmap_type + "_BNG_VS_seq"
        out2 = refaligner_dir + "/" + refaligner_prefix + "_" + optmap_type
        command = "perl " + DIGEST + " -v -i " + fasta_file_in_refaligner + " -e " + optmap_type + " -m 0 -M 0"
        os.system(command)
        command = REFALIGNER + " -f -ref " + optmap + " -i " + silicomap + " -o " + out + " -endoutlier 1e-3 -outlier 1e-4 -extend 1 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 2.9 -resSD 0.75 -mres 1.2 -A 5 -biaswt 0 -M 3 -Mfast 0 -maxmem 256 -maxthreads " + str(num_threads) + " -deltaX 9 -deltaY 9 -xmapchim 14 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-15 -stdout -stderr -xmaplen -indel" 
        os.system(command)
        out_list.append(out)
        out2_list.append(out2)

    chimeric(out_list, out2_list, chimeric_sub_dir, threshold_1, threshold_2, threshold_3, threshold_4, threshold_5, threshold_6)
    cutting(out2_list, fasta_file_in_refaligner, fasta_file_output, chimeric_sub_dir)
    for i in range(len(optmap_list)):
        optmap = optmap_list[i]
        cols = optmap.split('/')
        new_optmap = chimeric_sub_dir + "/nochi_optmap_" + str(i) + "_" + cols[len(cols)-1]
        cutting_ref(i, optmap, new_optmap, chimeric_sub_dir)
    

