
def ref_weight(alms, optmap_i, ref, vpos, edges, qry_quality, max_product, minrefoverhang, minqryoverhang):
    coverage = 0
    for qry in alms[optmap_i][ref]:
        x = alms[optmap_i][ref][qry]
        if x.start < vpos-minqryoverhang and x.end > vpos+minqryoverhang:
            if x.refstartpos < vpos - 30000 and x.refendpos < vpos - 30000: # alignment is on the left of this cutting site, so that it's not related
                continue
            if x.refstartpos > vpos + 30000 and x.refendpos > vpos + 30000: # alignment is on the right of this cutting site, so that it's not related
                continue
            coverage += 1
    weight = coverage / (1 + len(edges)*float(qry_quality)/max_product)
    print weight
    return weight

def qry_weight(alms_qry, qry, vpos, edges, ref_quality, max_product, minrefoverhang, minqryoverhang):
    coverage = 0
    for optmap_i in alms_qry[qry]:
        for ref in alms_qry[qry][optmap_i]:
            x = alms_qry[qry][optmap_i][ref]
            if x.refstart < vpos-minrefoverhang and x.refend > vpos+minrefoverhang:
                if x.qrystartpos < vpos - 30000 and x.qryendpos < vpos - 30000: # alignment is on the left of this cutting site, so that it's not related
                    continue
                if x.qrystartpos > vpos + 30000 and x.qryendpos > vpos + 30000: # alignment is on the right of this cutting site, so that it's not related
                    continue
                coverage += 1
    weight = coverage / (1 + len(edges)*float(ref_quality)/max_product)

    return weight
        
def vertex_weights(alms, vid_vertex, ugraph, minrefoverhang, minqryoverhang, ref_quality, qry_quality):
#    ref_quality = 1.5
#    qry_quality = 1
#    max_product = 10
    alms_qry = {} # a data struture for searching alignments by qry
    for optmap_i in range(len(alms)):
        for ref in alms[optmap_i]:
            for qry in alms[optmap_i][ref]:
                if qry not in alms_qry:
                    alms_qry[qry] = {}
                if optmap_i not in alms_qry[qry]:
                    alms_qry[qry][optmap_i] = {}
                alms_qry[qry][optmap_i][ref] = alms[optmap_i][ref][qry]

    max_product = len(alms) * ref_quality + len(alms_qry) * qry_quality
    vid_weight = {}
    for vid in vid_vertex:
        vname, vpos = vid_vertex[vid]
        if vname[0]=='R':
            cols = vname.split('_')
            optmap_i = int(cols[1])
            ref = int(cols[2])
            weight = ref_weight(alms, optmap_i, ref, vpos, ugraph[vid], qry_quality, max_product, minrefoverhang, minqryoverhang)
        else:
            weight = qry_weight(alms_qry, int(vname[1:]), vpos, ugraph[vid], ref_quality, max_product, minrefoverhang, minqryoverhang)
        vid_weight[vid] = weight     
    
    return vid_weight


