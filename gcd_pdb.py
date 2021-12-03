# gcd_pdb: tool for comparing the parts in common (gcd = greatest common
# denominator) among two or more pdb files.
# Author: Iris Young
# Modified: 2021.12.03

from pdb_parsing_tools import *
from collections import OrderedDict
import sys
from residues_util import residues, solvent_molecules, metal_ions

debug = False

def redirect_print(print_process, args, target=None):
  if target is not None:
    orig_out = sys.stdout
    out = open(target, "wb")
    sys.stdout = out
    print_process(*args)
    sys.stdout = orig_out
  else:
    print_process(*args)

def read_pdb(pdb):
  if debug:
    print("ATTEMPTING TO READ A PDB")
  header = []
  parsed = []
  ligand = []
  ions = []
  solvent = []
  b_total = 0
  n_atoms = 0
  with open(pdb, "rb") as reader:
    while True:
      try:
        record = reader.next()
        if (record[0:6] == "HETATM" or record[0:4] == "ATOM"):
          atom_dict = {
          'atom':get_atom(record),
          'resname':get_resname(record),
          'chain':get_chain(record),
          'resid':get_resid(record),
          'xyz':get_xyz(record),
          'occ':get_occ(record),
          'b':get_b(record),
          'elem':record[76:78].strip(" ")}
          # if record[17:20].strip(" ") in (amino_acids + ordered_solvent): # XXX temporary
          if record[17:20].strip(" ") in residues:
            parsed.append(atom_dict)
          # elif record[17:20].strip(" ") in disordered_solvent:
          elif record[17:20].strip(" ") in solvent_molecules:
            solvent.append(atom_dict)
          elif record[17:20].strip(" ") in metal_ions:
            ions.append(atom_dict)
          else:
            ligand.append(atom_dict)
          b_total += get_b(record)
          n_atoms += 1
        elif record[0:4] == "LINK" or record[0:6] == "SSBOND" or record[0:6] == "CRYST1":
          header.append(record)
      except StopIteration:
        break
    b_scale = b_total / n_atoms
  if debug:
    print("READ A PDB")
  return header, parsed, ligand, solvent, ions, b_scale

class hierarchical_pdb(object):
  def __init__(self, parsed_pdb_and_header):
    if debug:
      print("ATTEMPTING TO ASSEMBLE A PDB OBJECT")
    self.header = parsed_pdb_and_header[0]
    self.parsed_pdb = parsed_pdb_and_header[1]
    self.ligand = parsed_pdb_and_header[2]
    self.solvent = parsed_pdb_and_header[3]
    self.ions = parsed_pdb_and_header[4]
    self.b_scale = parsed_pdb_and_header[5]
    # self.chain_set = set()
    # self.records_by_chain = dict()
    # self.residues_sets = dict()
    # self.records_by_residue = dict()
    # self.atoms_sets = dict()
    # self.records_by_atom = dict()
    if debug:
      print("...SORTING CHAINS")
    self.get_chains()
    if debug:
      print("...SORTING RESIDUES")
    self.get_residues()
    if debug:
      print("...SORTING ATOMS")
    self.get_atoms()
    if debug:
      print("DONE")
  def get_chains(self):
    self.records_by_chain = OrderedDict()
    self.chain_set = set()
    for idx in range(len(self.parsed_pdb)):
      record = self.parsed_pdb[idx]
      if record['chain'] in self.records_by_chain.keys():
        self.records_by_chain[record['chain']].append(idx)
      else:
        self.records_by_chain[record['chain']] = [idx]
        self.chain_set.add(record['chain'])
  def get_residues(self):
    self.residues_sets = OrderedDict((chain, set()) for chain in self.records_by_chain.keys())
    self.records_by_residue = OrderedDict((chain, OrderedDict()) for chain in self.records_by_chain.keys())
    for chain in self.records_by_chain.keys():
      for record_idx in self.records_by_chain[chain]:
        record = self.parsed_pdb[record_idx]
        # residue = "%03d%3s" % (record['resid'], record['resname'])
        residue = "%4d%3s" % (record['resid'], record['resname'])
        if residue in self.records_by_residue[chain].keys():
          self.records_by_residue[chain][residue].append(record_idx)
        else:
          self.records_by_residue[chain][residue] = [record_idx]
          self.residues_sets[chain].add(residue)
  def get_atoms(self):
    self.atoms_sets = OrderedDict((chain, OrderedDict((residue, set())
                      for residue in self.records_by_residue[chain].keys()))
                      for chain in self.records_by_chain.keys())
    self.records_by_atom = OrderedDict((chain, OrderedDict((residue, OrderedDict())
                           for residue in self.records_by_residue[chain].keys()))
                           for chain in self.records_by_chain.keys())
    for chain in self.records_by_chain.keys():
      for residue in self.records_by_residue[chain].keys():
        for record_idx in self.records_by_residue[chain][residue]:
          record = self.parsed_pdb[record_idx]
          if record['atom'] in self.records_by_atom[chain][residue].keys():
            self.records_by_atom[chain][residue][record['atom']].append(record_idx)
          else:
            self.records_by_atom[chain][residue][record['atom']] = [record_idx]
            self.atoms_sets[chain][residue].add(record['atom'])
  def write_chain(self, chain, outfile):
    for record in self.header:
      print(record[:-2])
      record_number = 0
    with open(outfile, "wb") as out:
      for residue in self.records_by_residue[chain].keys():
        for atom in self.records_by_atom[chain][residue].keys():
          record_idx = self.records_by_atom[chain][residue][atom][0]
          record = self.parsed_pdb[record_idx]
          record_number += 1
          record_num_str = get_record_num_str(None, record_number)
          ## XXX FIXME shouldn't use % #s. Remove all spaces and specify exact lengths.
          out.write("""ATOM  %5s %4s % 3s %-2s% 4s   % 8.3f% 8.3f% 8.3f % 4.2f %5.2f          % 2s\n""" % \
            (record_num_str, record['atom'], record['resname'], record['chain'], record['resid'],
             record['xyz'][0], record['xyz'][1], record['xyz'][2], record['occ'], record['b'], record['elem']))
  def write_header(self, open_outfile):
    if debug:
      print("ATTEMPTING TO WRITE HEADER ONLY TO OPEN FILE")
    for record in self.header:
      open_outfile.write(record[:-2])
  def print_pdb(self):
    if debug:
      print("ATTEMPTING TO PRINT HIERARCHICAL PDB")
    for record in self.header:
      print(record[:-2])
    record_number = 0
    for chain in self.records_by_chain.keys():
      if debug:
        print(10*"_")
        print("CHAIN " + chain)
      for residue in self.records_by_residue[chain].keys():
        if debug:
          print(10*"_")
          print("RESIDUE " + residue)
        # for record_idx in self.records_by_residue[chain][residue]:
        for atom in self.records_by_atom[chain][residue].keys():
          record_idx = self.records_by_atom[chain][residue][atom][0]
          record = self.parsed_pdb[record_idx]
          record_number += 1
          record_num_str = get_record_num_str(None, record_number)
          print("""ATOM  %5s %4s % 3s %-2s% 4s   % 8.3f% 8.3f% 8.3f % 4.2f %5.2f          % 2s""" % \
            (record_num_str, record['atom'], record['resname'], record['chain'], record['resid'],
              record['xyz'][0], record['xyz'][1], record['xyz'][2], record['occ'], record['b'], record['elem']))
    for record in self.ligand:
      record_number += 1
      record_num_str = get_record_num_str(None, record_number)
      print("""ATOM  %5s %4s % 3s %-2s% 4s   % 8.3f% 8.3f% 8.3f % 4.2f %5.2f          % 2s""" % \
        (record_num_str, record['atom'], record['resname'], record['chain'], record['resid'],
          record['xyz'][0], record['xyz'][1], record['xyz'][2], record['occ'], record['b'], record['elem']))
    if debug:
      print("END OF HIERARCHICAL PDB")
  def get_pdb_summary(self):
    if debug:
      print("ATTEMPTING TO COMPILE PDB SUMMARY")
    summary = {}
    for chain in self.records_by_chain.keys():
      current = []
      for residue in self.records_by_residue[chain].keys():
        current.append(residue)
      summary[chain] = current
    if debug:
      for chain in summary.keys():
        print(" ".join(summary[chain]))# e.g. A 12 ASN 13 LEU 14 TRP 15 GLU 16 ARG...

class selection_pdb(hierarchical_pdb):
  def __init__(self, hierarchical_pdb, selection_condition):
    """Use selection_condition to keep or discard each atom record in hierarchical_pdb,
    and use the selected atoms to create a new hierarchical_pdb."""
    self.header = hierarchical_pdb.header
    self.parsed_pdb = []
    self.ligand = []
    self.b_scale = hierarchical_pdb.b_scale
    for record in hierarchical_pdb.parsed_pdb:
      if selection_condition(record):
        self.parsed_pdb.append(record)
    for record in hierarchical_pdb.ligand:
      if selection_condition(record):
        self.ligand.append(record)
    self.get_chains()
    self.get_residues()
    self.get_atoms()

class gcd_pdb(hierarchical_pdb):
  def __init__(self, hierarchical_pdb_list):
    if len(hierarchical_pdb_list) == 0:
      raise Exception("Please supply at least one readable pdb")
    elif len(hierarchical_pdb_list) == 1:
      self.pdb_list = []
      if debug:
        print("TAKING GCD OF ONLY ONE PDB")
    else:
      self.pdb_list = hierarchical_pdb_list[1:]
      if debug:
        print("ATTEMPTING TO TRIM PDB")
    pdb0 = hierarchical_pdb_list[0]
    self.header = pdb0.header
    self.parsed_pdb = pdb0.parsed_pdb
    self.ligand = pdb0.ligand
    self.b_scale = pdb0.b_scale
    self.chain_set = pdb0.chain_set
    self.records_by_chain = pdb0.records_by_chain
    self.residues_sets = pdb0.residues_sets
    self.records_by_residue = pdb0.records_by_residue
    self.atoms_sets = pdb0.atoms_sets
    self.records_by_atom = pdb0.records_by_atom
    if debug:
      print("...SORTING CHAINS")
    self.get_chains()
    if debug:
      print("...SORTING RESIDUES")
    self.get_residues()
    if debug:
      print("...SORTING ATOMS")
    self.get_atoms()
    if debug:
      print("DONE")
  def get_chains(self):
    for pdb in self.pdb_list:
      self.chain_set = self.chain_set.intersection(pdb.chain_set)
      self.records_by_chain = OrderedDict((chain, self.records_by_chain[chain])
                              for chain in sorted(self.chain_set))
  def get_residues(self):
    for pdb in self.pdb_list:
      self.residues_sets = OrderedDict((chain, self.residues_sets[chain].intersection(pdb.residues_sets[chain]))
                           for chain in self.records_by_chain.keys())
      self.records_by_residue = OrderedDict((chain, OrderedDict((residue, self.records_by_residue[chain][residue])
                                for residue in sorted(self.residues_sets[chain])))
                                for chain in self.records_by_chain.keys())
  def get_atoms(self):
    for pdb in self.pdb_list:
      self.atoms_sets = OrderedDict((chain, OrderedDict((residue,
                        self.atoms_sets[chain][residue].intersection(pdb.atoms_sets[chain][residue]))
                        for residue in self.records_by_residue[chain].keys()))
                        for chain in self.records_by_chain.keys())
      self.records_by_atom = OrderedDict((chain, OrderedDict((residue, OrderedDict((atom,
                             self.records_by_atom[chain][residue][atom])
                                for atom in sorted(self.atoms_sets[chain][residue])))
                                for residue in self.records_by_residue[chain].keys()))
                                for chain in self.records_by_chain.keys())
  def print_with_cla_hem_pheo_oex(self):
    if debug:
      print("ATTEMPTING TO PRINT HIERARCHICAL PDB")
    for record in self.header:
      print(record[:-2])
    record_number = 0
    for chain in self.records_by_chain.keys():
      if debug:
        print(10*"_")
        print("CHAIN " + chain)
      for residue in self.records_by_residue[chain].keys():
        if debug:
          print(10*"_")
          print("RESIDUE " + residue)
        for record_idx in self.records_by_residue[chain][residue]:
          record = self.parsed_pdb[record_idx]
          record_number += 1
          record_num_str = get_record_num_str(None, record_number)
          print("""ATOM  %5s %4s % 3s %-2s% 4s   % 8.3f% 8.3f% 8.3f % 4.2f %5.2f          % 2s""" % \
            (record_num_str, record['atom'], record['resname'], record['chain'], record['resid'],
              record['xyz'][0], record['xyz'][1], record['xyz'][2], record['occ'], record['b'], record['elem']))
    if debug:
      print(10*"_")
      print("CHLOROPHYLLS, HEMES, PHEOPHYTINS AND OECS")
    for record in self.ligand:
      record_number += 1
      record_num_str = get_record_num_str(None, record_number)
      print("""ATOM  %5s %4s % 3s %-2s% 4s   % 8.3f% 8.3f% 8.3f % 4.2f %5.2f          % 2s""" % \
        (record_num_str, record['atom'], record['resname'], record['chain'], record['resid'],
          record['xyz'][0], record['xyz'][1], record['xyz'][2], record['occ'], record['b'], record['elem']))
    if debug:
      print("END OF HIERARCHICAL PDB")

class gcd_pdb_diffb(gcd_pdb):
  def __init__(self, gcd_pdb_list):
    self.gcd_list = gcd_pdb_list
    self.header = gcd_pdb_list[0].header
    self.parsed_pdb = gcd_pdb_list[0].parsed_pdb
    self.ligand = gcd_pdb_list[0].ligand
    self.chain_set = gcd_pdb_list[0].chain_set
    self.records_by_chain = gcd_pdb_list[0].records_by_chain
    self.residues_sets = gcd_pdb_list[0].residues_sets
    self.records_by_residue = gcd_pdb_list[0].records_by_residue
    self.atoms_sets = gcd_pdb_list[0].atoms_sets
    self.records_by_atom = gcd_pdb_list[0].records_by_atom
  def get_diff_b_factor_pdb(self):
    self.diff_b_pdb = []
    b_scales = [gcd.b_scale for gcd in self.gcd_list]
    for chain in self.records_by_chain.keys():
      for residue in self.records_by_residue[chain].keys():
        for atom in self.records_by_atom[chain][residue].keys():
          record_indices = [(gcd, gcd.records_by_atom[chain][residue][atom][0]) for gcd in self.gcd_list]
          records = [gcd.parsed_pdb[index] for (gcd, index) in record_indices]
          scaled_b_factors = [records[idx]['b']/b_scales[idx] for idx in range(len(records))]
          min_b = min(scaled_b_factors)
          max_b = max(scaled_b_factors)
          diff_b_record = self.parsed_pdb[self.records_by_atom[chain][residue][atom][0]]
          diff_b_record['b'] = max_b - min_b
          self.diff_b_pdb.append(diff_b_record)
  def print_diff_b_factor_pdb(self):
    record_number = 0
    for record in self.header:
      print(record[:-2])
    for record in self.diff_b_pdb:
      record_number += 1
      record_num_str = get_record_num_str(None, record_number)
      print("""ATOM  %5s %4s % 3s %-2s% 4s   % 8.3f% 8.3f% 8.3f % 4.2f %5.2f          % 2s""" % \
        (record_num_str, record['atom'], record['resname'], record['chain'], record['resid'],
          record['xyz'][0], record['xyz'][1], record['xyz'][2], record['occ'], record['b'], record['elem']))

def gcd_pdb1_pdb2(pdb1, pdb2, ligands=False, outname=None, outname2=None):
  hier_pdb1 = hierarchical_pdb(read_pdb(pdb1))
  hier_pdb2 = hierarchical_pdb(read_pdb(pdb2))
  gcd_pdb1 = gcd_pdb([hier_pdb1, hier_pdb2])
  gcd_pdb2 = gcd_pdb([hier_pdb2, hier_pdb1])
  if outname is not None:
    gcd_name = outname
  else:
    gcd_name = "gcd_%s_%s" % (pdb1.split(".pdb")[0], pdb2)
  with open(gcd_name, "wb") as gcd1:
    if ligands:
      redirect_print(gcd_pdb1.print_with_cla_hem_pheo_oex, [], target=gcd_name)
    else:
      redirect_print(gcd_pdb1.print_pdb, [], target=gcd_name)
  if outname2 is not None:
    gcd_name = outname2
  else:
    gcd_name = "gcd_%s_%s" % (pdb2.split(".pdb")[0], pdb1)
  with open(gcd_name, "wb") as gcd2:
    if ligands:
      redirect_print(gcd_pdb2.print_with_cla_hem_pheo_oex, [], target=gcd_name)
    else:
      redirect_print(gcd_pdb2.print_pdb, [], target=gcd_name)
  return (gcd_pdb1, gcd_pdb2)

# needs to be updated from scitbx to numpy/scipy
# def superpose_hierarchical_pdbs(fixed_hier, moving_hier):
#   from scitbx.array_family import flex
#   from scitbx.math.superpose import least_squares_fit
#   fixed, moving = flex.vec3_double(), flex.vec3_double()
#   moving_records = []
#   for tup in ((fixed_hier, fixed), (moving_hier, moving)):
#     hier, coord = tup
#     for chain in sorted(fixed_hier.chain_set):
#       for resi in sorted(fixed_hier.residues_sets[chain]):
#         for atom in sorted(fixed_hier.atoms_sets[chain][resi]):
#           record_idx = hier.records_by_atom[chain][resi][atom][0]
#           rec = hier.parsed_pdb[record_idx]
#           xyz = [rec['xyz']]
#           coord.append(xyz)
#           if hier is moving_hier:
#             moving_records.append(rec)
#   lsq_soln = least_squares_fit(fixed, moving)
#   superposed_coords = lsq_soln.other_sites_best_fit()
#   superposed_records = []
#   for idx in range(len(moving_records)):
#     old = moving_records[idx]
#     new = old
#     new['xyz'] = superposed_coords[idx]
#     superposed_records.append(new)
#   super_hier = hierarchical_pdb(
#     (moving_hier.header,
#      superposed_records,
#      moving_hier.ligand,
#      moving_hier.b_scale))
#   return super_hier
