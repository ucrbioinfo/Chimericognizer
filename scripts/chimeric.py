#!/usr/bin/python
import csv
import sys

from sets import Set

from collections import defaultdict

import collections
import math

from chimeric_alignment import Alignment 
from vertex_weights import vertex_weights
import vertex_cover_exhaust
import vertex_cover_aprx

def copy_alms(old_alms, removed):
    current_alms = {}
    for ref in old_alms:
        current_alms[ref] = {}
        for qry in old_alms[ref]:
            x = old_alms[ref][qry]
            if removed[ref, qry] == False:
                current_alms[ref][qry] = x
    return current_alms

def search_in_sorted_list(L, v):
    smaller_num = 0
    larger_num = 0
    for i in range(len(L)):
        if L[i] < v:
            smaller_num += 1
        elif L[i] > v:
            larger_num += 1
    return (smaller_num, larger_num)

def markers_in_qry_left_overhang(qry_markers, x):
    qry = x.qry
    pos = x.qrystartpos
    ori = x.orientation
    pos_list = qry_markers[qry]
    (smaller_num, larger_num) = search_in_sorted_list(pos_list, pos)
    if ori == '+':
        return smaller_num
    else:
        return larger_num

def markers_in_qry_right_overhang(qry_markers, x):
    qry = x.qry
    pos = x.qryendpos
    ori = x.orientation
    pos_list = qry_markers[qry]
    (smaller_num, larger_num) = search_in_sorted_list(pos_list, pos)
    if ori == '+':
        return larger_num
    else:
        return smaller_num


def DFS_ugraph(ug, visited, start, subgraph):
    visited[start] = True
    subgraph[start] = ug[start]
    preds = ug[start]
    for v in preds:
        if visited[v] == True:
            continue
        DFS_ugraph(ug, visited, v, subgraph)

def pre_process(optmap_i, myfile, myfile2, output_dir, min_confidence, minrefoverhang, minqryoverhang):
    header_lines = 10
    header = []

    all_alms = {} # stores all the Alignments for all groups, all_groups[ref] should contain molecule ref
    qualify_alms = {} # only keep one alignment(the one with highest confidence) for each contig in one molecule
    removed = {} # removed[ref,qry] == True means alignment for (ref, qry) is already removed

    # collecting alignments and store in all_groups
    print '---------------read .xmap file-------------------'
    with open(myfile+'.xmap', 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for i in range(header_lines): # 10 lines of header
            header.append(csvreader.next()) # save them
        # read the first non-header line
        while True:
            try:
                row = csvreader.next()
                x = Alignment(int(row[1]),int(row[2]),float(row[3]),float(row[4]),float(row[5]),
                                  float(row[6]),row[7],float(row[8]),row[9],float(row[10]),
                                  float(row[11]),int(row[12]),row[13])
                if x.ref not in all_alms:
                    all_alms[x.ref] = [x]
                else:
                    all_alms[x.ref].append(x)
            except StopIteration:
                break
    num_all_alms = 0
    for ref in all_alms:
        num_all_alms += len(all_alms[ref])
    print "In total, the number of alignments collected is ", num_all_alms
    
    
    # only keep one alignment(the one with highest confidence) for each contig in one molecule    
    for ref in all_alms:
        group = all_alms[ref]
        qry_bestx = {}
        for x in group:
            if x.qry not in qry_bestx:
                qry_bestx[x.qry] = x
            else:
                if x.confidence > qry_bestx[x.qry].confidence:
                    qry_bestx[x.qry] = x

        qualify_alms[ref] = {}
        for qry in qry_bestx:
            qualify_alms[ref][qry] = qry_bestx[qry]

    num_qualify_alms = 0
    for ref in qualify_alms:
        num_qualify_alms += len(qualify_alms[ref])
    print "In total, the number of alignments in qualify_alms is ", num_qualify_alms

    # initialize removed array
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            removed[ref,qry] = False
    print '---------------END-------------------'
    
    
    # remove low confidence alignments
    print '---------------Remove low quality alignments---------------'
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            x = qualify_alms[ref][qry]
            if x.confidence < min_confidence:
                removed[ref, qry] = True
                print 'alignment (', ref, ',', qry, ') is low quality and removed'
    num_alms = 0
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            if removed[ref, qry] == False:
                num_alms += 1
    print "After removing low confidence alignments, the number of alignments is ", num_alms
    print '---------------End---------------'
   
    print '---------------scaling-------------------'
    # calculating scaling
    qry_len = {}
    with open(myfile2+'_key.txt') as f_key:
        for i in range(0, 4): # 4 header lines
            f_key.readline()  
        for line in f_key:
            line = line.strip()
            cols = line.split('\t')
            qry_id = int(cols[0])
            seq_len = int(cols[2])
            qry_len[qry_id] = seq_len
    scaling = 0
    num = 0
    with open(myfile+'_q.cmap') as f_q:
        for i in range(0, 11): # 11 header lines
            f_q.readline()
        for line in f_q:
            line = line.strip()
            cols = line.split('\t')
            qry_id = int(cols[0])
            appr_len = float(cols[1])
            seq_len = qry_len[qry_id]
            scaling += appr_len/seq_len
            num += 1
    scaling /= num # scaling=1.02258059775

    # use scaling to adjsut coordinates
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            x = qualify_alms[ref][qry]
            x.qrystartpos /= scaling
            x.qryendpos /= scaling
            x.qrylen /= scaling
            x.refstartpos /= scaling
            x.refendpos /= scaling
            x.reflen /= scaling          
    print '---------------END-------------------'

    # find the reference-based coordinates for each contig
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            x = qualify_alms[ref][qry]
            if (x.orientation == '+'):
                x.qry_left_overlen = x.qrystartpos
                x.qry_right_overlen = x.qrylen - x.qryendpos
            else:
                x.qry_left_overlen = x.qrylen - x.qrystartpos
                x.qry_right_overlen = x.qryendpos
            x.start = x.refstartpos - x.qry_left_overlen
            x.end = x.refendpos + x.qry_right_overlen 
            x.ref_left_overlen = x.refstartpos
            x.ref_right_overlen = x.reflen - x.refendpos
            if (x.orientation == '+'):
                x.refstart = x.qrystartpos - x.ref_left_overlen
                x.refend = x.qryendpos + x.ref_right_overlen
            else:
                x.refstart = x.qryendpos - x.ref_right_overlen
                x.refend = x.qrystartpos + x.ref_left_overlen


    # read old optmap
    qry_markers = {}
    with open(myfile+'_q.cmap') as f_q:
        for i in range(11): # 10 lines of header
            header_line = f_q.readline()
        for line in f_q:
            line = line.strip()
            cols = line.split('\t')
            CMapId = int(cols[0])
            ContigLength = float(cols[1])
            NumSites = int(cols[2])       
            SiteID = int(cols[3])
            LabelChannel = cols[4]
            Position = float(cols[5])
            if LabelChannel == "0":
                continue
            if CMapId not in qry_markers:
                qry_markers[CMapId] = []
            Position /= scaling
            qry_markers[CMapId].append(Position)
    for CMapId in qry_markers:
        qry_markers[CMapId].sort()
    f_q.close()


    print '---------------candidate cutting sites-------------------'
    fpair = file(output_dir+"/chimeric_pairs_"+str(optmap_i)+".log", 'w')
    fpair.write("ref_id\tref_pos\tqry_id\tqry_pos\n")
    chimeric_pairs = []
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            x = qualify_alms[ref][qry]
            if (x.confidence > min_confidence):
                ref_left_overlen = x.refstartpos
                ref_right_overlen = x.reflen - x.refendpos
                if (x.qry_left_overlen > minqryoverhang and ref_left_overlen > minrefoverhang and markers_in_qry_left_overhang(qry_markers, x) > 0): 
                    chimeric_pairs.append((x.ref, x.refstartpos, x.qry, x.qrystartpos))
                    print (x.ref, x.refstartpos, x.qry, x.qrystartpos), "is a pair of candidate cutting sites"
                    fpair.write(str(x.ref)+"\t"+str(x.refstartpos)+"\t"+str(x.qry)+"\t"+str(x.qrystartpos)+"\n")
                if (x.qry_right_overlen > minqryoverhang and ref_right_overlen > minrefoverhang and markers_in_qry_right_overhang(qry_markers, x) > 0):
                    chimeric_pairs.append((x.ref, x.refendpos, x.qry, x.qryendpos))
                    print (x.ref, x.refendpos, x.qry, x.qryendpos), "is a pair of candidate cutting sites"
                    fpair.write(str(x.ref)+"\t"+str(x.refendpos)+"\t"+str(x.qry)+"\t"+str(x.qryendpos)+"\n")
    fpair.close()
    print '---------------END-------------------'

    current_alms = copy_alms(qualify_alms, removed)

    return chimeric_pairs, scaling, current_alms


def build_graph(vertices, vid_vertex, vertex_vid, edges, chimeric_pairs, optmap_i, num):
    for i in range(0, len(chimeric_pairs)):
        (ref, ref_pos, qry, qry_pos) = chimeric_pairs[i]
        ref_id = 'R_'+str(optmap_i)+"_"+str(ref)
        qry_id = 'Q'+str(qry)
        if ref_id not in vertices:
            vertices[ref_id] = set([ref_pos])
            vid_vertex[num] = (ref_id, ref_pos)
            vertex_vid[ref_id, ref_pos] = num
            vid_1 = num
            num += 1
        else: 
            found = False
            for p in vertices[ref_id]:
                if abs(p-ref_pos) <= 30000: # if the distance of 2 candidate cutting sites is smaller than 30,000bp, merge them
                    found = True
                    ref_pos_found = p
                    break
            if found == False:
                vertices[ref_id].add(ref_pos)
                vid_vertex[num] = (ref_id, ref_pos)
                vertex_vid[ref_id, ref_pos] = num
                vid_1 = num
                num += 1
            else:
                vid_1 = vertex_vid[ref_id, ref_pos_found]
        if qry_id not in vertices:
            vertices[qry_id] = set([qry_pos])
            vid_vertex[num] = (qry_id, qry_pos)
            vertex_vid[qry_id, qry_pos] = num
            vid_2 = num
            num += 1
        else:
            found = False
            for p in vertices[qry_id]:
                if abs(p-qry_pos) <= 30000: # if the distance of 2 candidate cutting sites is smaller than 30,000bp, merge them
                    found = True
                    qry_pos_found = p
                    break
            if found == False:
                vertices[qry_id].add(qry_pos)
                vid_vertex[num] = (qry_id, qry_pos)
                vertex_vid[qry_id, qry_pos] = num
                vid_2 = num
                num += 1
            else:
                vid_2 = vertex_vid[qry_id, qry_pos_found]
        vid_min = min(vid_1, vid_2)
        vid_max = max(vid_1, vid_2)
        edges.add((vid_min, vid_max))
    return num


def chimeric(myfile_list, myfile2_list, output_dir, ref_quality, qry_quality, max_product, min_confidence, minrefoverhang, minqryoverhang):
    
    # discard alignments below min_confidence
    # alignment overhangs above this number of bps are considered chimeric
#    minrefoverhang = 50000
#    minqryoverhang = 50000

    all_alms_list = []
    scaling_list = []
    chimeric_pairs_list = []
    for i in range(len(myfile_list)):
        myfile = myfile_list[i]
        myfile2 = myfile2_list[i]
        chimeric_pairs, scaling, alms = pre_process(i, myfile, myfile2, output_dir, min_confidence, minrefoverhang, minqryoverhang)
        scaling_list.append(scaling)
        chimeric_pairs_list.append(chimeric_pairs)
        all_alms_list.append(alms) 


    print '---------------build graph-------------------'
    vertices = {}
    vid_vertex = {}
    vertex_vid = {}
    edges = set([])
    num = 0

    for optmap_i in range(len(chimeric_pairs_list)):
        chimeric_pairs = chimeric_pairs_list[optmap_i]
        newnum = build_graph(vertices, vid_vertex, vertex_vid, edges, chimeric_pairs, optmap_i, num)
        num = newnum


    # store the undirected graph in ugraph
    ugraph = {}
    for vid in vid_vertex:
        ugraph[vid] = set([])
    for e in edges:   
        (vid_min, vid_max) = e
        ugraph[vid_min].add(vid_max)
        ugraph[vid_max].add(vid_min)

    # calculate weight for each vertex
    vid_weight = vertex_weights(all_alms_list, vid_vertex, ugraph, minrefoverhang, minqryoverhang, ref_quality, qry_quality, max_product)
    for vid in vid_weight:
        weight = vid_weight[vid]
        print vid_vertex[vid], weight

    fv = file(output_dir+"/vertices.log", 'w')
    fv.write("vid\t(seq_id, position)\tweights\n")
    for vid in vid_vertex:
        print vid, vid_vertex[vid] 
        fv.write(str(vid)+"\t")
        fv.write("( "+vid_vertex[vid][0]+", "+str(vid_vertex[vid][1])+" )\t") 
        fv.write(str(vid_weight[vid]))
        fv.write("\n")
    fv.close()

    fe = file(output_dir+"/edges.log", 'w')
    fe.write("vid_1\tvid_2\n")
    for v1, v2 in edges:
        print [v1, v2] 
        fe.write(str(v1)+"\t"+str(v2)+"\n")
    fe.close()
    print "In total, the number of vertices in graph is", len(vid_vertex)
    print "In total, the number of edges in graph is", len(edges)

    print '---------------END-------------------'
    
    print '---------------vertex cover-------------------'
    # divide ugraph into connected component
    visited = {}
    vid_list = []
    for vid in vid_vertex:
        visited[vid] = False
        vid_list.append(vid)
    vid_list.sort()

    ugraphs = []
    while True:
        found = False
        for vid in vid_list:
            if visited[vid] == False:
                found = True
                start = vid
                break
        if found == True: # there's still unvisited vertix in ugraph
            oneugraph = {}
            DFS_ugraph(ugraph, visited, start, oneugraph)
            ugraphs.append(oneugraph)
        else:
            break

    # call vertex cover algorithm
    whole_cover = set([])
    for ug in ugraphs:
        sub_vid_weight = {}
        for vid in ug:
            sub_vid_weight[vid] = vid_weight[vid]
        cover = vertex_cover_exhaust.VC(ug, sub_vid_weight)
        for v in cover:
            whole_cover.add(v)

    fc = file(output_dir+"/cover.log", 'w')
    fc.write("vid\t(seq_id, position)\tweights\n")
    for vid in whole_cover:
        print vid, vid_vertex[vid]
        fc.write(str(vid)+"\t")
        fc.write("( "+vid_vertex[vid][0]+", "+str(vid_vertex[vid][1])+" )\t") 
        fc.write(str(vid_weight[vid]))
        fc.write("\n")
    fc.close()
    print '---------------END-------------------'

    print '---------------optmap scaling back-------------------'

    # use scaling to adjsut coordinates
    ref_cuts_list = [] # list of cuts on ref
    qry_cuts_list = [] # list of cuts on qry
    for vid in whole_cover:
        (vname, vpos) = vid_vertex[vid]
        if vname[0] == 'R':
            cols = vname.split('_')
            optmap_i = int(cols[1])
            vpos = vpos * scaling_list[optmap_i]
            ref_cuts_list.append((vname, vpos))
        elif vname[0] == 'Q':
            vpos = int(math.ceil(float(vpos)))        
            qry_cuts_list.append((vname, vpos))
    print '---------------END-------------------'    
    
    print '---------------store cut information-------------------'
    qry_cuts = {}
    for vname, vpos in qry_cuts_list:
        qry = int(vname[1:])
        if qry not in qry_cuts:
            qry_cuts[qry] = []
        qry_cuts[qry].append(vpos)
    for qry in qry_cuts:
        qry_cuts[qry].sort()

    qry_cuts_file = output_dir + "/qry_cuts.txt"
    f_qcuts = file(qry_cuts_file, 'w')
    f_qcuts.write("qry_id\tcutting positions\n")
    for qry in qry_cuts:
        f_qcuts.write(str(qry)+"\t")
        for vpos in qry_cuts[qry]:
            f_qcuts.write(str(vpos)+" ")
        f_qcuts.write("\n")
    f_qcuts.close()

    ref_cuts = {}
    for vname, vpos in ref_cuts_list:
        cols = vname.split('_')
        optmap_i = int(cols[1])
        ref = int(cols[2])
        if optmap_i not in ref_cuts:
            ref_cuts[optmap_i] = {}
        if ref not in ref_cuts[optmap_i]:
            ref_cuts[optmap_i][ref] = []
        ref_cuts[optmap_i][ref].append(vpos)
    for optmap_i in ref_cuts:
        for ref in ref_cuts[optmap_i]:
            ref_cuts[optmap_i][ref].sort()

    for optmap_i in ref_cuts:
        ref_cuts_file = output_dir + "/ref_cuts_"+str(optmap_i)+".txt"
        f_rcuts = file(ref_cuts_file, 'w')
        f_rcuts.write("ref_id\tcutting positions\n")
        for ref in ref_cuts[optmap_i]:
            f_rcuts.write(str(ref)+"\t")
            for vpos in ref_cuts[optmap_i][ref]:
                f_rcuts.write(str(vpos)+" ")
            f_rcuts.write("\n")
        f_rcuts.close()
    
    print '---------------END-------------------'

if __name__ == "__main__":
    chimeric(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5]),float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8]), float(sys.argv[9]))



