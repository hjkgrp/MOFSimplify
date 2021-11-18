from __future__ import print_function
import numpy as np
import re
import os
from ciftemplate2graph import isvert
from itertools import combinations
import datetime
import networkx as nx
from bbcif_properties import iscoord, isbond

def nn(string):
	return re.sub('[^a-zA-Z]','', string)

def nl(string):
	return re.sub('[^0-9]','', string)

def PBC3DF(c1, c2):
    
    diffa = c1[0] - c2[0]
    diffb = c1[1] - c2[1]
    diffc = c1[2] - c2[2]

    if diffa > 0.5:
        c2[0] = c2[0] + 1.0
    elif diffa < -0.5:
        c2[0] = c2[0] - 1.0
    
    if diffb > 0.5:
        c2[1] = c2[1] + 1.0
    elif diffb < -0.5:
        c2[1] = c2[1] - 1.0
 
    if diffc > 0.5:
        c2[2] = c2[2] + 1.0
    elif diffc < -0.5:
        c2[2] = c2[2] - 1.0
    
    return c2

def PBC3DF_sym(vec1, vec2):

	dX,dY,dZ = vec1 - vec2
			
	if dX > 0.5:
		s1 = 1 + 5 
		ndX = dX - 1.0
	elif dX < -0.5:
		s1 = -1 + 5
		ndX = dX + 1.0
	else:
		s1 = 0 + 5
		ndX = dX
				
	if dY > 0.5:
		s2 = 1 + 5
		ndY = dY - 1.0
	elif dY < -0.5:
		s2 = -1 + 5
		ndY = dY + 1.0
	else:
		s2 = 0 + 5
		ndY = dY
	
	if dZ > 0.5:
		s3 = 1 + 5
		ndZ = dZ - 1.0
	elif dZ < -0.5:
		s3 = -1 + 5
		ndZ = dZ + 1.0
	else:
		s3 = 0 + 5
		ndZ = dZ

	if str(s1) + str(s2) + str(s3) == '555':
		sym = '.'
	else:
		sym = '1_' + str(s1) + str(s2) + str(s3)

	return np.array([ndX,ndY,ndZ]), sym

def write_check_cif(template, placed_nodes, placed_edges, g, sp, sc_unit_cell):

	sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma = sp
	q = 0

	tpath = os.join('templates', template)

	with open(tpath, 'r') as tcif:

		tcif = tcif.read()
		tcif = filter(None, tcif.split('\n'))

	cpath = os.path.join('check_cifs', str(g) + '_check_scaled_placed_' + template)

	with open(cpath, 'w') as check:
		for line in tcif:
			s = line.split()
			if not isvert(s):
				if '_cell_length_a' in line:
					check.write('_cell_length_a   ' + str(sc_a))
				elif '_cell_length_b' in line:
					check.write('_cell_length_b   ' + str(sc_b))
				elif '_cell_length_c' in line:
					check.write('_cell_length_c   ' + str(sc_c))
				elif '_cell_angle_alpha' in line:
					check.write('_cell_angle_alpha   ' + str(sc_alpha))
				elif '_cell_angle_beta' in line:
					check.write('_cell_angle_beta   ' + str(sc_beta))
				elif '_cell_angle_gamma' in line:
					check.write('_cell_angle_gamma   ' + str(sc_gamma))
				else:
					check.write(line)
				check.write('\n')
			else:
				for n in placed_edges:

					q += 1
					name = re.sub('[0-9]','',n[0])
					if name == 'X':
						name = 'C'
					index = name + str(q)
					vec = np.array(list(map(float,[n[1],n[2],n[3]])))
					v = np.dot(np.linalg.inv(sc_unit_cell), vec)
					check.write('{:>5}{:>5}{:>20}{:>20}{:>20}{:>12}{:>8}{:>8}'.format(index,name,v[0],v[1],v[2],'0.00000','Uiso','1.00'))
					check.write('\n')

				for n in placed_nodes:

					q += 1
					name = re.sub('[0-9]','',n[0])
					if name == 'X':
						name = 'C'
					index = name + str(q)
					vec = np.array(list(map(float,[n[1],n[2],n[3]])))
					v = np.dot(np.linalg.inv(sc_unit_cell), vec)
					check.write('{:>5}{:>5}{:>20}{:>20}{:>20}{:>12}{:>8}{:>8}'.format(index,name,v[0],v[1],v[2],'0.00000','Uiso','1.00'))
					check.write('\n')

				break

def distance_search_bond(placed_all, bonds_all, sc_unit_cell, tol):

	fixed_bonds = []
	used_bonds = []
	fixed_bonds_append = fixed_bonds.append
	used_bonds_append = used_bonds.append
	
	for l in bonds_all:
		fixed_bonds_append([l[0],l[1],l[2],'.',l[4]])
		used_bonds_append((l[0],l[1]))

	connection_points = [line for line in placed_all if re.sub('[0-9]','',line[5]) == 'X']
	nbcount = 0
			
	for i in range(len(connection_points)):

		ielem = connection_points[i][0]
		ivec = np.dot(np.linalg.inv(sc_unit_cell), np.array([float(q) for q in connection_points[i][1:4]]))
		ibbid = int(connection_points[i][6])

		for j in range(i + 1, len(connection_points)):

			jelem = connection_points[j][0]
			jbbid = int(connection_points[j][6])

			if (ielem, jelem) not in used_bonds and (jelem, ielem) not in used_bonds:

				jvec = np.dot(np.linalg.inv(sc_unit_cell), list(map(float, connection_points[j][1:4])))
				DV, sym = PBC3DF_sym(ivec,jvec)
				dist = np.linalg.norm(np.dot(sc_unit_cell, DV))
	
				if dist < tol and ibbid != jbbid:
					nbcount += 1
					fixed_bonds_append([ielem, jelem, dist, sym, 'S'])
					break
				else:
					continue

	return fixed_bonds, nbcount

def bond_connected_components(placed_all, bonds_all, sc_unit_cell, max_length, bond_tol, nconnections_list, num_possible_XX_bonds):

	G = nx.Graph()

	for n in placed_all:
		isX = (re.sub('[0-9]','',n[-3]) == 'X')
		G.add_node(n[0], coords=np.array(list(map(float,(n[1:4])))), occ=n[4], isX=isX, bbcode=int(n[6]), nconnect=1, bbtype=n[-1])
	for l in bonds_all:
		G.add_edge(l[0], l[1], length=l[2], sym=l[3], ty=l[4], is_new_XX=False, order=(l[0],l[1]))

	for line in nconnections_list:
		node,nc = line
		G.nodes[node]['nconnect'] -= 1
		G.nodes[node]['nconnect'] += nc

	previous_degrees = dict((node,G.degree(node)) for node in G.nodes())

	ccs = list(nx.connected_components(G))
	count = 0
	bb_tol = max_length + 0.50*max_length
	print('distance search tolerance is', np.round(bb_tol,3), 'Angstroms')

	for connect_comp0, connect_comp1 in combinations(ccs,2):

		xname0 = [n for n in connect_comp0 if G.nodes[n]['isX']]
		xvecs0 = [np.dot(np.linalg.inv(sc_unit_cell),G.nodes[n]['coords']) for n in connect_comp0 if G.nodes[n]['isX']]
		com0 = np.average(xvecs0, axis=0)
		type0 = list(set([G.nodes[n]['bbtype'] for n in connect_comp0]))

		xname1 = [n for n in connect_comp1 if G.nodes[n]['isX']]
		xvecs1 = [np.dot(np.linalg.inv(sc_unit_cell),G.nodes[n]['coords']) for n in connect_comp1 if G.nodes[n]['isX']]
		com1 = np.average(xvecs1, axis=0)
		type1 = list(set([G.nodes[n]['bbtype'] for n in connect_comp1]))

		if len(type0) > 1 or len(type1) > 1:
			print(type0, type1)
			raise ValueError('building block indicated as both node and edge type')

		if len(xname0) == 0 or len(xname1) == 0:
			raise ValueError('There are connected components with no connection site atoms')

		type0 = type0[0]
		type1 = type1[0]
		# don't want to try any edge-edge bonds, this will always result in incorrect bonding that must be corrected later
		if type0 == 'edge' and type1 == 'edge':
			continue

		com_dist = np.linalg.norm(np.dot(sc_unit_cell, com0 - PBC3DF(com0,com1)))

		if com_dist <= bb_tol:

			min_dist = (100.0, 'no.3', 'the', 'larch')

			for xv0,xn0 in zip(xvecs0, xname0):
				for xv1,xn1 in zip(xvecs1, xname1):

					DV, sym = PBC3DF_sym(xv0,xv1)
					dist = np.linalg.norm(np.dot(sc_unit_cell, DV))

					if dist < min_dist[0]:
						min_dist = (dist, xn0, xn1, sym)

			if min_dist[0] < bond_tol:
				count += 1
				G.add_edge(min_dist[1], min_dist[2], length=min_dist[0], sym=min_dist[3], ty='S', is_new_XX=True, order=(min_dist[1],min_dist[2]))

	# attempt to bond any connection sites left without X neighbors
	# by forming bonds (distance search) within the X atoms withoug any X neighbors
	connection_nodes = [(n,data) for n,data in G.nodes(data=True) if data['isX']]
	no_bond_connection_sites = []
	
	for node,data in connection_nodes:

		bbcode0 = data['bbcode']
		nconnections = data['nconnect']
		nbors = list(G.neighbors(node))
		X_nbors = [nbor for nbor in nbors if (G.nodes[nbor]['isX'] and G.nodes[nbor]['bbcode'] != bbcode0)]

		if len(X_nbors) < nconnections:
			no_bond_connection_sites.append(node)

	for n0,n1 in combinations(no_bond_connection_sites, 2):
		
		vec0 = np.dot(np.linalg.inv(sc_unit_cell), G.nodes[n0]['coords'])
		vec1 = np.dot(np.linalg.inv(sc_unit_cell), G.nodes[n1]['coords'])
	
		bbcode0 = G.nodes[n0]['bbcode']
		bbcode1 = G.nodes[n1]['bbcode']
	
		bbtype0 = G.nodes[n0]['bbtype']
		bbtype1 = G.nodes[n1]['bbtype']
	
		if bbcode0 == bbcode1:
			continue
		if bbtype0 == 'edge' and bbtype1 == 'edge':
			continue
	
		DV, sym = PBC3DF_sym(vec0,vec1)
		dist = np.linalg.norm(np.dot(sc_unit_cell, DV))
	
		if dist < bond_tol:
			count += 1
			G.add_edge(n0, n1, length=dist, sym=sym, ty='S', is_new_XX=True, order=(n0,n1))

	# remove any extra bonds by keeping only bonds to the closest N X-neighbors
	# where N is the number of possible connections for each X
	for node,data in connection_nodes:
		
		vec0 = np.dot(np.linalg.inv(sc_unit_cell), G.nodes[node]['coords'])
		bbcode0 = data['bbcode']
		bbtype0 = data['bbtype']
	
		nbors = list(G.neighbors(node))
		X_nbors = [nbor for nbor in nbors if (G.nodes[nbor]['isX'] and G.nodes[nbor]['bbcode'] != bbcode0)]
		nconnections = data['nconnect']
	
		if len(X_nbors) > nconnections:
	
			nbor_dists = []
			for nbor in X_nbors:
				vec1 = np.dot(np.linalg.inv(sc_unit_cell), G.nodes[nbor]['coords'])
				dist = np.linalg.norm(np.dot(sc_unit_cell, vec0 - PBC3DF(vec0,vec1)))
				nbor_dists.append((dist,nbor))
			nbor_dists.sort(key=lambda x:x[0])
	
			for nbor in nbor_dists[nconnections:]:
				G.remove_edge(node, nbor[1])
				count -= 1

	wrong_connection_nodes = []
	wrong_connection_nodes_append = wrong_connection_nodes.append

	for node,data in connection_nodes:

		nconnect = data['nconnect']
		nbors = list(G.neighbors(node))
		bbcode0 = data['bbcode']

		X_nbors = [nbor for nbor in nbors if (G.nodes[nbor]['isX'] and G.nodes[nbor]['bbcode'] != bbcode0)]
		if len(X_nbors) != nconnect:
			wrong_connection_nodes_append(node)

	XX_bond_count = 0
	for e0,e1,data in G.edges(data=True):
		if data['is_new_XX']:
			XX_bond_count += 1
	
	if XX_bond_count != count:
		raise ValueError('bond counts do not match, there is a problem with the connection site bonding algorithm')

	bond_check_passed = True

	for node,data in G.nodes(data=True):

		isX = data['isX']
		degree = G.degree(node)
		previous_degree = previous_degrees[node]
		nconnections = data['nconnect']

		if isX and degree != previous_degree + nconnections:
			print('There are connection sites with too many/not enough bonds formed.')
			bond_check_passed = False
			break
		elif not isX and degree != previous_degree:
			print('The degree of a non-X atom has changed after XX bond formation.')
			bond_check_passed = False
			break

	if len(wrong_connection_nodes) > 0:
		print('There are connection sites with too few or too many bonds')
		bond_check_passed = False
	if XX_bond_count < num_possible_XX_bonds:
		print('The number of XX bonds formed', XX_bond_count, 'is less than the correct number', num_possible_XX_bonds)
		print('Try increasing the bond formation distance tolerance.')
		bond_check_passed = False
	if XX_bond_count > num_possible_XX_bonds:
		print('The number of XX bonds formed', XX_bond_count, 'is more than the correct number', num_possible_XX_bonds)
		print('Check your inputs for incorrect X atoms.')
		bond_check_passed = False

	fixed_bonds = []
	fixed_bonds_append = fixed_bonds.append

	for edge in G.edges(data=True):
		
		edict = edge[2]
		ty = edict['ty']
		leng = edict['length']
		sy = edict['sym']
		order = edict['order']
		fixed_bonds_append([order[0], order[1], leng, sy, ty])

	return fixed_bonds, count, bond_check_passed

def fix_bond_sym(bonds_all,placed_all,sc_unit_cell):
	
	coords_dict = dict((l[0],np.dot(np.linalg.inv(sc_unit_cell), list(map(float, l[1:4])))) for l in placed_all)

	fixed_bonds = []
	fixed_bonds_append = fixed_bonds.append
	for l in bonds_all:
		vec1 = coords_dict[l[0]]
		vec2 = coords_dict[l[1]]

		dist,sym = PBC3DF_sym(vec1,vec2)

		fixed_bonds_append([l[0],l[1],l[2],sym,l[4]])

	return fixed_bonds

def write_cif(placed_all, fixed_bonds, scaled_params, sc_unit_cell, cifname, charges, wrap_coords=True):

	sc_a,sc_b,sc_c,sc_alpha,sc_beta,sc_gamma = scaled_params

	opath = os.path.join('output_cifs', cifname)
	
	with open(opath, 'w') as out:
		out.write('data_' + cifname[0:-4] + '\n')
		out.write('_audit_creation_date              ' + datetime.datetime.today().strftime('%Y-%m-%d') + '\n')
		out.write("_audit_creation_method            'tobacco_3.0'" + '\n')
		out.write("_symmetry_space_group_name_H-M    'P1'" + '\n')
		out.write('_symmetry_Int_Tables_number       1' + '\n')
		out.write('_symmetry_cell_setting            triclinic' + '\n')
		out.write('loop_' + '\n')
		out.write('_symmetry_equiv_pos_as_xyz' + '\n')
		out.write('  x,y,z' + '\n')
		out.write('_cell_length_a                    ' + str(sc_a) + '\n')
		out.write('_cell_length_b                    ' + str(sc_b) + '\n')
		out.write('_cell_length_c                    ' + str(sc_c) + '\n')
		out.write('_cell_angle_alpha                 ' + str(sc_alpha) + '\n')
		out.write('_cell_angle_beta                  ' + str(sc_beta) + '\n')
		out.write('_cell_angle_gamma                 ' + str(sc_gamma) + '\n')
		out.write('loop_' + '\n')
		out.write('_atom_site_label' + '\n')
		out.write('_atom_site_type_symbol' + '\n')
		out.write('_atom_site_fract_x' + '\n')
		out.write('_atom_site_fract_y' + '\n')
		out.write('_atom_site_fract_z' + '\n')
		if charges:
			out.write('_atom_site_charge' + '\n')

		for l in placed_all:

			vec = list(map(float, l[1:4]))
			cvec = np.dot(np.linalg.inv(sc_unit_cell), vec)

			if wrap_coords:
				cvec = np.mod(cvec, 1) # makes sure that all fractional coordinates are in [0,1]

			if charges:
				out.write('{:7} {:>4} {:>15} {:>15} {:>15} {:>15}'.format(l[0], re.sub('[0-9]','',l[0]), "%.10f" % np.round(cvec[0],10), "%.10f" % np.round(cvec[1],10), "%.10f" % np.round(cvec[2],10), l[4]))
				out.write('\n')
			else:
				out.write('{:7} {:>4} {:>15} {:>15} {:>15}'.format(l[0], re.sub('[0-9]','',l[0]), "%.10f" % np.round(cvec[0],10), "%.10f" % np.round(cvec[1],10), "%.10f" % np.round(cvec[2],10)))
				out.write('\n')

		out.write('loop_' + '\n')
		out.write('_geom_bond_atom_site_label_1' + '\n')
		out.write('_geom_bond_atom_site_label_2' + '\n')
		out.write('_geom_bond_distance' + '\n')
		out.write('_geom_bond_site_symmetry_2' + '\n')
		out.write('_ccdc_geom_bond_type' + '\n')

		for e in fixed_bonds:
			out.write('{:7} {:>7} {:>5} {:>7} {:>3}'.format(e[0], e[1], "%.3f" % float(e[2]), e[3], e[4]))
			out.write('\n')

def cif_read(filename, charges=False):

	with open(filename,'r') as f:
		f = f.read()
		f = filter(None, f.split('\n'))

	names = []
	elems = []
	fcoords = []
	charge_list = []
	bonds = []

	for line in f:
		s = line.split()
		if '_cell_length_a' in line:
			a = s[1]
		if '_cell_length_b' in line:
			b = s[1]
		if '_cell_length_c' in line:
			c = s[1]
		if '_cell_angle_alpha' in line:
			alpha = s[1]
		if '_cell_angle_beta' in line:
			beta = s[1]
		if '_cell_angle_gamma' in line:
			gamma = s[1]
		if iscoord(s):
			
			names.append(s[0])
			elems.append(s[1])
			
			fvec = np.array([np.round(float(v),8) for v in s[2:5]])
			for dim in range(len(fvec)):
				if fvec[dim] < 0.0:
					fvec[dim] += 1.0
				elif fvec[dim] > 1.0:
					fvec[dim] -= 1.0

			fcoords.append(fvec)

			if charges:
				charge_list.append(float(s[-1]))
			else:
				charge_list.append(0.0)

		if isbond(s):
			bonds.append((s[0],s[1],s[2],s[3],s[4]))

	pi = np.pi
	a,b,c,alpha,beta,gamma = map(float,(a,b,c,alpha,beta,gamma))
	ax = a
	ay = 0.0
	az = 0.0
	bx = b * np.cos(gamma * pi / 180.0)
	by = b * np.sin(gamma * pi / 180.0)
	bz = 0.0
	cx = c * np.cos(beta * pi / 180.0)
	cy = (c * b * np.cos(alpha * pi /180.0) - bx * cx) / by
	cz = (c ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
	unit_cell = np.asarray([[ax,ay,az],[bx,by,bz],[cx,cy,cz]]).T
	
	ccoords = []

	for l in fcoords:
		vec = l
		vec = np.dot(unit_cell,vec)
		ccoords.append(vec)

	fcoords = np.asarray(fcoords)
	ccoords = np.asarray(ccoords)
	charges = np.asarray(charges)

	return elems, names, ccoords, fcoords, charge_list, bonds, (a,b,c,alpha,beta,gamma), unit_cell

def merge_catenated_cifs(comb, charges):
	
	all_read = [cif_read(cif) for cif in comb]
	natoms = [len(c[0]) for c in all_read]
	a,b,c,alpha,beta,gamma = all_read[0][6]
	
	if len(set(natoms)) != 1:
		raise ValueError('The number of atoms is not the same for the cifs:', comb)

	cifname = comb[0][0:-9] + '.cif'
	with open(cifname, 'w') as out:

		out.write('data_' + cifname[0:-4] + '\n')
		out.write('_audit_creation_date              ' + datetime.datetime.today().strftime('%Y-%m-%d') + '\n')
		out.write("_audit_creation_method            'tobacco_3.0'" + '\n')
		out.write("_symmetry_space_group_name_H-M    'P1'" + '\n')
		out.write('_symmetry_Int_Tables_number       1' + '\n')
		out.write('_symmetry_cell_setting            triclinic' + '\n')
		out.write('loop_' + '\n')
		out.write('_symmetry_equiv_pos_as_xyz' + '\n')
		out.write('  x,y,z' + '\n')
		out.write('_cell_length_a                    ' + str(a) + '\n')
		out.write('_cell_length_b                    ' + str(b) + '\n')
		out.write('_cell_length_c                    ' + str(c) + '\n')
		out.write('_cell_angle_alpha                 ' + str(alpha) + '\n')
		out.write('_cell_angle_beta                  ' + str(beta) + '\n')
		out.write('_cell_angle_gamma                 ' + str(gamma) + '\n')
		out.write('loop_' + '\n')
		out.write('_atom_site_label' + '\n')
		out.write('_atom_site_type_symbol' + '\n')
		out.write('_atom_site_fract_x' + '\n')
		out.write('_atom_site_fract_y' + '\n')
		out.write('_atom_site_fract_z' + '\n')
		if charges:
			out.write('_atom_site_charge' + '\n')

		count = 0
		for cif in all_read:

			elems, names, ccoords, fcoords, charge_list, bonds, (a,b,c,alpha,beta,gamma), unit_cell = cif
			shift = count*natoms[0]
			count += 1
			names = [nn(name) + str(int(nl(name)) + shift) for name in names]

			if charges:
				for name, coord, chg in zip(names, fcoords, charge_list):
					out.write('{:7} {:>4} {:>15.8f} {:>15.8f} {:>15.8f} {:>15.8f}'.format(name, nn(name), coord[0], coord[1], coord[2], chg))
					out.write('\n')
			else:
				for name, coord in zip(names, fcoords):
					out.write('{:7} {:>4} {:>15.8f} {:>15.8f} {:>15.8f}'.format(name, nn(name), coord[0], coord[1], coord[2]))
					out.write('\n')

		out.write('loop_' + '\n')
		out.write('_geom_bond_atom_site_label_1' + '\n')
		out.write('_geom_bond_atom_site_label_2' + '\n')
		out.write('_geom_bond_distance' + '\n')
		out.write('_geom_bond_site_symmetry_2' + '\n')
		out.write('_ccdc_geom_bond_type' + '\n')

		count = 0
		for cif in all_read:

			elems, names, ccoords, fcoords, charge_list, bonds, (a,b,c,alpha,beta,gamma), unit_cell = cif
			shift = count*natoms[0]
			count += 1
			bonds = [(nn(b[0]) + str(int(nl(b[0])) + shift), nn(b[1]) + str(int(nl(b[1])) + shift), b[2], b[3], b[4]) for b in bonds]
			
			for bond in bonds:
				out.write('{:7} {:>7} {:>5} {:>7} {:>3}'.format(*bond))
				out.write('\n')


