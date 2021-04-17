from zeo.netstorage import AtomNetwork
from zeo.area_volume import volume, surface_area

atmnet = AtomNetwork.read_from_CSSR("MgO_vac1.cssr", rad_file="MgO.rad")
vol_str = volume(atmnet, 0.1, 0.05, 20000)
lines = vol_str.split('\n')
for line in lines:
    if "Number_of_pockets" in line:
        print '---------'
        print line
        print '---------'
vol_str, ha_atmnet = volume(atmnet, 0.1, 0.05, 7000, True)
lines = vol_str.split('\n')
for line in lines:
    if "Number_of_pockets" in line:
        print '---------'
        print line
        print '---------'
#print "--------"
#print vol_str
#print "--------"
sa_str = surface_area(atmnet, 0.1, 0.05, 7000, False)
lines = sa_str.split('\n')
for line in lines:
    if "Number_of_pockets" in line:
        print '---------'
        print line.split()
        print '---------'
#print "--------"
#print sa_str
#print "--------"

