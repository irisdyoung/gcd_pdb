from __future__ import division
import math

# use gcd_pdb to represent models hierarchically and select elements in common
from gcd_pdb import hierarchical_pdb, gcd_pdb1_pdb2
from pdb_parsing_tools import get_record_from_record_dict, write_record_with_B

def mean(list_of_values):
  return sum(list_of_values)/len(list_of_values)

def get_shift(resi1, hier1, resi2, hier2, between='all'):
  if between == 'CA':
    atoms1 = [r for r in resi1.values() if r['atom'] == "CA"]
    atoms2 = [r for r in resi2.values() if r['atom'] == "CA"]
  else: # shift between average of all atom positions
    atoms1 = resi1.values()
    atoms2 = resi2.values()
  shift_vec = []
  for index in (0,1,2): # x,y,z
    records1 = [hier1.parsed_pdb[i[j]] for i in atoms1 for j in range(len(i))]
    records2 = [hier2.parsed_pdb[i[j]] for i in atoms2 for j in range(len(i))]
    position1 = mean([r['xyz'][index] for r in records1])
    position2 = mean([r['xyz'][index] for r in records2])
    shift_vec.append(abs(position2 - position1)) # magnitude, unsigned
  shift = math.sqrt(shift_vec[0]**2 + shift_vec[1]**2 + shift_vec[2]**2)
  return shift

def visualize_shifts(hier_pdb1, hier_pdb2, out_filename):
  # calculate shifts between matching residues in the two models
  # don't bother constructing an object with this info if we're just writing to file
  #
  # shifts_by_residue = OrderedDict((chain, OrderedDict((residue,
  #                       get_shift(hier_pdb1.records_by_atom[chain][residue],
  #                                 hier_pdb2.records_by_atom[chain][residue],
  #                                 between='all'))
  #                       for residue in hier_pdb1.records_by_atom[chain]))
  #                       for chain in hier_pdb1.records_by_atom)
  #
  # encode shifts of hier_pdb2 relative to hier_pdb1 in the B factor column
  # of hier_pdb1 and write to a specified new filename
  with open(out_filename, 'w') as out_pdb:
    hier_pdb1.write_header(out_pdb)
    record_number = 0
    for chain in hier_pdb1.records_by_atom:
      for resi in hier_pdb1.records_by_atom[chain]:
        shift = get_shift(hier_pdb1.records_by_atom[chain][resi], hier_pdb1,
                          hier_pdb2.records_by_atom[chain][resi], hier_pdb2,
                          between='all')
        for atom in hier_pdb1.records_by_atom[chain][resi]:
          for record_idx in hier_pdb1.records_by_atom[chain][resi][atom]:
            record_number += 1
            record = get_record_from_record_dict(hier_pdb1.parsed_pdb[record_idx], record_number)
            out_pdb.write(write_record_with_B(record, shift))
  print "shifts written to B factor column of", out_filename

if __name__ == "__main__":
  import sys
  pdb1_filename, pdb2_filename = sys.argv[1:3]
  pdb1, pdb2 = gcd_pdb1_pdb2(pdb1_filename, pdb2_filename)
  # superposed_pdb2 = superpose_hierarchical_pdbs(pdb1, pdb2)
  # presume pdbs have already been superposed as best they can
  visualize_shifts(pdb1, pdb2, pdb1_filename + "_relativeto_" + pdb2_filename)
  # vis.write_pymol_script()
  # print path to script?