from zeo.netstorage import AtomNetwork

atmnet = AtomNetwork.read_from_CSSR("MgO.cssr", rad_file="MgO.rad")
a = atmnet.perform_voronoi_decomposition()
#vornet,edge_centers, face_centers = atmnet.perform_voronoi_decomposition()
#print 'edge_centers', edge_centers
#print 'face_centers', face_centers
print len(a)

print a
