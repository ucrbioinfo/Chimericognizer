import os.path
from shutil import copyfile

from merge_inputs import merge_inputs
from chname_fasta import chname_fasta
from fasta_long_seqs import fasta_long_seqs
from chimeric_removal_onetime import chimeric_removal_onetime

def abs_path(path):
    cols = path.split('/')
    if cols[0] == '~':
        user = os.path.expanduser('~') 
        return user + path[1:]
    if cols[0] == '.':
        current = os.path.abspath(os.path.join(os.getcwd(), "."))
        return current + path[1:]
    if cols[0] == '..':
        parent = os.path.abspath(os.path.join(os.getcwd(), ".."))
        return parent + path[2:]
    return path

def output_paras(fasta_list_file, output_dir, optmap_list_file, num_threads, REFALIGNER, threshold_1, threshold_2, threshold_3, threshold_4, threshold_5, threshold_6):
    fo = file(output_dir+"/parameters.log", 'w')
    fo.write("fasta list files:\t"+fasta_list_file+"\n")
    fo.write("output directory:\t"+output_dir+"\n")
    fo.write("optical map list file:\t"+optmap_list_file+"\n")
    fo.write("number of threads:\t"+str(num_threads)+"\n")
    fo.write("REFALIGNER tool:\t"+REFALIGNER+"\n")
    fo.write("threshold_1:\t"+str(threshold_1)+"\n")
    fo.write("threshold_2:\t"+str(threshold_2)+"\n")
    fo.write("threshold_3\t"+str(threshold_3)+"\n")
    fo.write("threshold_4\t"+str(threshold_4)+"\n")
    fo.write("threshold_5\t"+str(threshold_5)+"\n")
    fo.write("threshold_6\t"+str(threshold_6)+"\n")
    fo.close()

def get_optmap_info(optmap_list_file, output_dir):
    fo = file(output_dir+"/optmap_info.log", 'w')
    optmap_list = []
    optmap_type_list = []
    for line in open(optmap_list_file):
        line = line.strip()
        if len(line) == 0:
            continue
        cols = line.split()
        optmap_type = cols[0]
        optmap = abs_path(cols[1])
        fo.write(optmap_type+"\t"+optmap+"\n")
        optmap_list.append(optmap)
        optmap_type_list.append(optmap_type)
    return optmap_list, optmap_type_list

def chimeric_removal(fasta_list_file, output_dir, DIGEST, REFALIGNER, optmap_list_file, num_threads, threshold_1, threshold_2, threshold_3, threshold_4, threshold_5, threshold_6):

    #tools

#    DIGEST = "./tools/fa2cmap_multi.pl"

    #intermediate files unrelated to each iteration
    input_merged_file = output_dir + "/input_merged.fasta" 
    input_merged_file_chname = output_dir + "/input_merged_chname.fasta"
    input_merged_file_chname_long = output_dir + "/input_merged_chname_long.fasta"
    final_nochi_file = output_dir + "/nochi_all.fasta"


    #check existence of output_dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    #pre-operations
    output_paras(fasta_list_file, output_dir, optmap_list_file, num_threads, REFALIGNER, threshold_1, threshold_2, threshold_3, threshold_4, threshold_5, threshold_6) 
    optmap_list, optmap_type_list = get_optmap_info(optmap_list_file, output_dir)
    merge_inputs(output_dir, fasta_list_file, input_merged_file)
    chname_fasta(input_merged_file, input_merged_file_chname)
    fasta_long_seqs(input_merged_file_chname, input_merged_file_chname_long) 

    #
    chimeric_removal_onetime(REFALIGNER, DIGEST, output_dir, optmap_list, optmap_type_list, input_merged_file_chname_long, final_nochi_file, num_threads, threshold_1, threshold_2, threshold_3, threshold_4, threshold_5, threshold_6)


    





