import sys

def read_keyfile(myfile2, sname_qry):
    key_file = myfile2 + "_key.txt"
    with open(key_file) as f_key:
        for i in range(0, 4): # 4 header lines
            f_key.readline()  
        for line in f_key:
            line = line.strip()
            cols = line.split('\t')
            qry_id = int(cols[0])
            seq_name = cols[1]
            sname_qry[seq_name] = qry_id
    f_key.close()


def cutting(myfile2_list, fasta_file, new_fasta_file, output_dir):
    print '---------------cutting seqs-------------------'
    # read key_file
    sname_qry = {}
    for i in range(len(myfile2_list)):
        myfile2 = myfile2_list[i]
        read_keyfile(myfile2, sname_qry)

    # read qry_cuts_file
    qry_cuts_file = output_dir + "/qry_cuts.txt"
    qry_cuts = {}
    with open(qry_cuts_file) as f_qcuts:
        f_qcuts.readline()
        for line in f_qcuts:
            line = line.strip()
            cols = line.split('\t')
            qry_id = int(cols[0])
            positions = cols[1]
            cols = positions.split()
            cuts = []
            for c in cols:
                cuts.append(int(c))
            qry_cuts[qry_id] = cuts
    f_qcuts.close()

    # read fasta file
    seqs_list = {}
    for line in open(fasta_file):
        line = line.strip()
        if line[0] == '>':
            seq_name = line[1:]
            seqs_list[seq_name] = []
            current_name = seq_name
        else:
            seqs_list[current_name].append(line)
    

    seqs = {}
    for name in seqs_list:
        seq_list = seqs_list[name]
        whole_seq = ''.join(seq_list)
        seqs[name] = whole_seq

    # cutting
    f = file(output_dir+"/cuts_info.log", 'w')
    new_seqs = {}
    for name in seqs:
        if name not in sname_qry:
            new_seqs[name] = seqs[name]
            continue
        qry = sname_qry[name]
        if qry not in qry_cuts:
            new_seqs[name] = seqs[name]
            continue
        f.write(name+"\n")
        print name
        cuts = qry_cuts[qry]
        seq = seqs[name]
        subseq_list = []
        lastpos = 0
        for i in range(0, len(cuts)):
            pos = cuts[i]
            subseq = seq[lastpos:pos]
            newname = name+" |chimeric_"+str(lastpos+1)+"_to_"+str(pos)
            f.write(newname+"\n")
            print newname
            new_seqs[newname] = subseq
            lastpos = pos
        subseq = seq[lastpos:]
        newname = name+" |chimeric_"+str(lastpos+1)+"_to_"+str(len(seq))
        new_seqs[newname] = subseq
        f.write(newname+"\n")
        f.write("\n")
        print 
    f.close()

    # write new fasta file
    fout = file(new_fasta_file, 'w')
    for name in new_seqs:
        new_seq = new_seqs[name]
        fout.write(">"+name+"\n")
        fout.write(new_seq+"\n")
    fout.close()
    print '---------------Ends-------------------'
    

if __name__ == "__main__":
    cutting(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])



