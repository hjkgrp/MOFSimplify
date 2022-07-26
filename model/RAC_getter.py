import sys
from molSimplify.Informatics.MOF.MOF_descriptors import get_primitive, get_MOF_descriptors;

def main():

	# user command line inputs
	cif_folder = sys.argv[1]
	name = sys.argv[2] # name of the MOF
	RACs_folder = sys.argv[3]

	# result log
	f = open(RACs_folder + 'RAC_getter_log.txt', 'w')

	try:
		# makes the linkers and sbus folders
		full_names, full_descriptors = get_MOF_descriptors(cif_folder + name + '_primitive.cif', 3, path = RACs_folder, xyzpath = RACs_folder + name + '.xyz');
	except ValueError:
		f.write('FAILED')
		f.close()
		return 'FAILED'
	except NotImplementedError:
		f.write('FAILED')
		f.close()
		return 'FAILED'
	except AssertionError:
		f.write('FAILED')
		f.close()
		return 'FAILED'

	if (len(full_names) <= 1) and (len(full_descriptors) <= 1): # this is a featurization check from MOF_descriptors.py
		f.write('FAILED')
		f.close()
		return 'FAILED'


if __name__ == "__main__":
	main()
