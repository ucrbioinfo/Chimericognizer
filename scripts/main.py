#!/usr/bin/python

import sys, getopt
import os.path

from chimeric_removal import chimeric_removal

def paras_in_file(paras_file):
    if not os.path.isfile(paras_file):    
        print "ERROR! Parameter file " + paras_file + " doesn't exist!"
        exit()
    argv = []
    for line in open(paras_file):
        line = line.strip()
        cols = line.split()
        argv += cols
    return argv   


def main():

    #default input and output
    paras_file = "./paras.txt"
    fasta_list_file = "./input_files_list.txt"
    output_dir = "./output_dir"
    optmap_list_file = "./optmap_files_list.txt"
 
    #default tools 
    REFALIGNER = "RefAligner" 

    #default parameters
    num_threads = 32
    threshold_1 = 1.5
    threshold_2 = 1
    threshold_3 = 25
    threshold_4 = 50000
    threshold_5 = 50000
    threshold_6 = 80000

    #obtaining parameters
    paras_string = "f:i:o:m:p:x:y:a:b:d:e:h:r:"
    opts, args = getopt.getopt(sys.argv[1:], paras_string)
    for op, value in opts:
        if op == "-f":
            paras_file = value
            argv_in_file = paras_in_file(paras_file)
            opts_in_file, args_in_file = getopt.getopt(argv_in_file, paras_string) 
            for fop, fvalue in opts_in_file:
                if fop == "-i":
                    fasta_list_file = fvalue
                elif fop == "-o":
                    output_dir = fvalue
                elif fop == "-m":
                    optmap_list_file = fvalue
                elif fop == "-p":
                    num_threads = int(fvalue)
                elif fop == "-x":
                    REFALIGNER = fvalue
                elif fop == "-y":
                    DIGEST = fvalue
                elif fop == "-a":
                    threshold_1 = float(fvalue)
                elif fop == "-b":
                    threshold_2 = float(fvalue)
                elif fop == "-d":
                    threshold_3 = float(fvalue)
                elif fop == "-e":
                    threshold_4 = float(fvalue)
                elif fop == "-h":
                    threshold_5 = float(fvalue)
                elif fop == "-r":
                    threshold_6 = float(fvalue)
        elif op == "-i":
            fasta_list_file = value
        elif op == "-o":
            output_dir = value
        elif op == "-m":
            optmap_list_file = value
        elif op == "-p":
            num_threads = int(value)
        elif op == "-x":
            REFALIGNER = value
        elif op == "-y":
            DIGEST = value
        elif op == "-a":
            threshold_1 = float(value)
        elif op == "-b":
            threshold_2 = float(value)
        elif op == "-d":
            threshold_3 = float(value)
        elif op == "-e":
            threshold_4 = float(value)
        elif op == "-h":
            threshold_5 = float(value)
        elif op == "-r":
            threshold_6 = float(value)

    chimeric_removal(paras_file, fasta_list_file, output_dir, DIGEST, REFALIGNER,  optmap_list_file, num_threads, threshold_1, threshold_2, threshold_3, threshold_4, threshold_5, threshold_6)        
if __name__ == "__main__":
    main()
