import sys
import os.path

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

def output_input(output_dir, fname_list):
    
    fo = file(output_dir+"/input.log", 'w')
    fos = file(output_dir+"/input_seqs.log", 'w')
    num_seqs = 0
    for fname in fname_list:
        begin = num_seqs + 1
        for line in open(fname):
            line = line.strip()
            if line == "":
                continue
            if line[0] == '>':
                seqname = line[1:]
                num_seqs += 1
                fos.write(str(num_seqs)+"\t"+seqname+"\t"+fname+"\n")
        end = num_seqs
        fo.write(fname+"\t"+str(begin)+"\t"+str(end)+"\n")
        
    fo.close()
    fos.close()

    

def merge_inputs(output_dir, fasta_list_file, merged_fasta_file):
    if not os.path.isfile(fasta_list_file):    
        print "ERROR! Inputs list file " + fasta_list_file + " doesn't exist!"
        exit()
    fname_list = set([])
    for line in open(fasta_list_file):
        line = line.strip()
        if line == "":
            continue
        fname_list.add(abs_path(line))
    output_input(output_dir, fname_list)
    fout = file(merged_fasta_file, 'w')
    for fn in fname_list:
        print "Reading file: ", fn
        f = open(fn, 'r')
        text = f.read()
        text = text.strip()
        fout.write(text+"\n")
        f.close()
    fout.close()

if __name__ == "__main__":
    merge_inputs(output_dir, sys.argv[1], sys.argv[2])


