import sys
from zeo.netstorage import AtomNetwork
from zeo.cycle import compute_centroid_4cycles, compute_face_centers

atmnet = AtomNetwork.read_from_CSSR("MgO.cssr", rad_file="MgO.rad")
compute_face_centers(atmnet)
sys.exit()

vornet = atmnet.perform_voronoi_decomposition()
oh_interstitials = compute_centroid_4cycles(vornet)

ids_set = set()
i = 0
while i < len(oh_interstitials):
    node_ids = frozenset(oh_interstitials[i]['ids'])
    if node_ids not in ids_set:
        ids_set.add(node_ids)
        i = i+1
    else:
        oh_interstitials.pop(i)
        
print len(oh_interstitials)
for inter in oh_interstitials:
    node_ids = inter['ids']
    print node_ids
    coord = inter['coords']
    print coord.x, coord.y, coord.z


