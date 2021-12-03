from pdb_util import get_interatomic_distance
from gcd_pdb import read_pdb
from pdb_parsing_tools import get_resname, get_atom, isatom

# rename atoms of a particular residue according to a pair of templates
def rename_atoms_of_selected_residue(
  pdbfile, resname, template_pdb_start, template_pdb_target, newfilename):
  # first step is to construct the mapping from start to target template
  # for this we need to get the closest atom in template_pdb_target to each
  # atom in template_pdb_start. Assume templates are aligned.
  _, resis, ligands, solvent, ions, _ = read_pdb(template_pdb_start)
  records_start = [r for r in (resis + ligands + solvent + ions) if r['resname'] == resname]
  _, resis, ligands, solvent, ions, _ = read_pdb(template_pdb_target)
  records_target = [r for r in (resis + ligands + solvent + ions) if r['resname'] == resname]
  distance_matrix = []
  for rtarget in records_target:
    matrix_row = []
    for rstart in records_start:
      matrix_row.append(get_interatomic_distance(rtarget['xyz'], rstart['xyz']))
    distance_matrix.append(matrix_row)
  match_indices = [row.index(min(row)) for row in distance_matrix]
  records_match = [records_start[i] for i in match_indices]
  lookup = {}
  for i in xrange(len(records_match)):
    rtarget = records_target[i]
    rmatch = records_match[i]
    lookup[rmatch['atom']] = rtarget['atom']
    print 'replacing all instances of %s with %s' % (rmatch['atom'], rtarget['atom'])
  def update_record(record):
    new_atom = lookup[get_atom(record)]
    new_record = record[:12] + ("% 4s" % new_atom) + record[16:]
    return new_record
  with open(pdbfile, 'r') as oldfile:
    with open(newfilename, 'w') as newfile:
      count = 0
      for record in oldfile.readlines():
        if isatom(record) and get_resname(record) == resname.strip():
          newfile.write(update_record(record))
          count += 1
        else:
          newfile.write(record)
      print 'updated %i atom names' % count
      print 'updated file written to %s' % newfilename

if __name__ == "__main__":
  import sys
  rename_atoms_of_selected_residue(*sys.argv[1:6])
