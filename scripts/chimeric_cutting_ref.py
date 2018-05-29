import sys

def search_in_sorted_list(array, value):
    for i in range(0, len(array)):
        if value < array[i]:
            return i
    return len(array)

def cutting_ref(optmap_i, old_optmap, new_optmap, output_dir):
    print '---------------cutting seqs-------------------'

    # read qry_cuts_file
    ref_cuts_file = output_dir + "/ref_cuts_" + str(optmap_i) + ".txt"
    ref_cuts = {}
    with open(ref_cuts_file) as f_rcuts:
        f_rcuts.readline()
        for line in f_rcuts:
            line = line.strip()
            cols = line.split('\t')
            ref_id = int(cols[0])
            positions = cols[1]
            cols = positions.split()
            cuts = []
            for c in cols:
                cuts.append(float(c))
            ref_cuts[ref_id] = cuts
    f_rcuts.close()

    # read old optmap
    last_markers = {}
    markers = {}
    maps_len = {}
    max_CMapId = 0
    headers = []
    with open(old_optmap) as f_map:
        for line in f_map:
            line = line.strip()
            if line == "":
                continue
            if line[0] == '#':
                headers.append(line)
                continue
            cols = line.split('\t')
            CMapId = int(cols[0])
            ContigLength = float(cols[1])
            NumSites = int(cols[2])       
            SiteID = int(cols[3])
            LabelChannel = cols[4]
            Position = float(cols[5])
            rest = ""
            if max_CMapId < CMapId:
                max_CMapId = CMapId
            for i in range(6, len(cols)):
                rest = rest + cols[i] + "\t"
            new_cols = (CMapId, ContigLength, NumSites, SiteID, LabelChannel, Position, rest)
            if LabelChannel == "0":
                maps_len[CMapId] = ContigLength
                last_markers[CMapId] = new_cols
            else: # =="1"
                if CMapId not in markers:
                    markers[CMapId] = []
                markers[CMapId].append(new_cols)
    # divide markers of old map into new maps
    new_markers = {}
    new_last_markers = {}
    new_maps_len = {}
    new_CMapId = max_CMapId + 1
    stable_markers_set = set([])
    newref_oldref = {}
    for ref in markers:
        if ref not in ref_cuts: 
            new_markers[ref] = markers[ref]
            stable_markers_set.add(ref)
            continue
        cuts = ref_cuts[ref]
        new_maps_len_list = [-1 for x in range(len(cuts)+1)]
        for i in range(len(cuts)):
            if i == 0:
                new_maps_len_list[i] = cuts[i]
            else:
                new_maps_len_list[i] = cuts[i] - cuts[i-1]
            new_maps_len_list[len(cuts)] = maps_len[ref] - cuts[len(cuts)-1]
        new_maps = [[] for x in range(len(cuts)+1)]
        for i in range(0, len(markers[ref])):
            (CMapId, ContigLength, NumSites, SiteID, LabelChannel, Position, rest) = markers[ref][i]
            index = search_in_sorted_list(cuts, Position)                
            if index == 0:
                Position = Position
            else:
                Position -= cuts[index-1]
            new_maps[index].append((CMapId, ContigLength, NumSites, SiteID, LabelChannel, Position, rest))
        for i in range(len(new_maps)):
            new_CMapId += 1
            new_markers[new_CMapId] = new_maps[i]
            new_maps_len[new_CMapId] = new_maps_len_list[i]
            new_last_markers[new_CMapId] = last_markers[ref]
            newref_oldref[new_CMapId] = ref

    fo = file(output_dir+"/optmap_cut_" + str(optmap_i) + ".log", 'w')
    fo.write("new_CMapId\told_CMapId\n")
    for newref in newref_oldref:
        oldref = newref_oldref[newref]
        fo.write(str(newref)+"\t"+str(oldref)+"\n")
    fo.close()

    # change markers 
    for ref in new_markers:
        if ref in stable_markers_set:
            (CMapId, ContigLength, NumSites, SiteID, LabelChannel, Position, rest) = last_markers[ref]
            new_markers[ref].append((CMapId, ContigLength, NumSites, SiteID, LabelChannel, Position, rest))
            continue
        marker_num = len(new_markers[ref])
        for i in range(len(new_markers[ref])):
            (CMapId, ContigLength, NumSites, SiteID, LabelChannel, Position, rest) = new_markers[ref][i]    
            new_CMapId = ref
            new_ContigLength = new_maps_len[ref]
            new_NumSites = marker_num
            new_SiteID = i+1
            new_markers[ref][i] = (new_CMapId, new_ContigLength, new_NumSites, new_SiteID, LabelChannel, Position, rest)
        # add last marker
        (CMapId, ContigLength, NumSites, SiteID, LabelChannel, Position, rest) = new_last_markers[ref]
        new_CMapId = ref
        new_ContigLength = new_maps_len[ref]
        new_NumSites = marker_num
        new_SiteID = marker_num + 1
        Position = new_maps_len[ref]
        new_markers[ref].append((new_CMapId, new_ContigLength, new_NumSites, new_SiteID, LabelChannel, Position, rest))
        
    # write into file
    fo = file(new_optmap, 'w')
    for i in range(len(headers)):
        line = headers[i]
#        if i == 7:
#            line = line + " changed to " + str(len(new_markers))
        fo.write(line+"\n") 
    for ref in new_markers:
        for i in range(len(new_markers[ref])):
            (CMapId, ContigLength, NumSites, SiteID, LabelChannel, Position, rest) = new_markers[ref][i]
            fo.write(str(CMapId)+"\t"+str(ContigLength)+"\t"+str(NumSites)+"\t"+str(SiteID)+"\t"+LabelChannel+"\t"+str(Position)+"\t"+rest+"\n")
    fo.close()    

if __name__ == "__main__":
    cutting_ref(sys.argv[1], sys.argv[2], sys.argv[3])



