from zeo.high_accuracy import high_accuracy_atmnet
from zeo.netstorage import AtomNetwork

atmnet = AtomNetwork.read_from_CSSR("MgO.cssr", rad_file="MgO.rad")
atmnet.write_to_XYZ("orig_mgo.xyz", False, True)
vornet,fcs = atmnet.perform_voronoi_decomposition()
vornet.write_to_XYZ("orig_mgo_voro.xyz", 0)
vornet.analyze_writeto_XYZ('orig_mgo', 0.4, atmnet)
high_accuracy_atmnet(atmnet, "DEF")
vornet,fcs = atmnet.perform_voronoi_decomposition()
vornet.write_to_XYZ("test_high.xyz", 0)
vornet.analyze_writeto_XYZ('mgo_high', 0.4, atmnet)
atmnet.write_to_CIF("highacc_MgO.cif")
atmnet.calculate_free_sphere_parameters('MgO.res')
#
#atmnet = AtomNetwork.read_from_CIF("mgo.cif", rad_file="MgO.rad")
#vornet = atmnet.perform_voronoi_decomposition()
#vornet.write_to_XYZ("mgo.xyz", 0)
#vornet.analyze_writeto_XYZ('mgo1', 0.4, atmnet)
#high_accuracy_atmnet(atmnet, "DEF")
#vornet = atmnet.perform_voronoi_decomposition()
#vornet.write_to_XYZ("mgo_high.xyz", 0)
#vornet.analyze_writeto_XYZ('mgo_high', 0.4, atmnet)
#atmnet.write_to_CIF("highacc_mgo.cif")
#atmnet.calculate_free_sphere_parameters('mgo.res')
#
#
#atmnet = AtomNetwork.read_from_CIF("mgo.cif", rad_file="MgO.rad")
#high_accuracy_atmnet(atmnet, "LOW")
#atmnet.write_to_XYZ("mgo_ha_lowset.xyz", False, True)
#vornet = atmnet.perform_voronoi_decomposition()
#vornet.write_to_XYZ("mgo_ha_low.xyz", 0)
#vornet.analyze_writeto_XYZ('mgo_ha_low', 0.4, atmnet)
