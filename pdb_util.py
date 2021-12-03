#!/usr/bin/env python
#
# gcd_pdb_cmd.py
#
from pdb_parsing_tools import *
from gcd_pdb import read_pdb, hierarchical_pdb
from residues_util import modified_all

def phenix_fetch(pdbid):
  import asyncio # not tested with this yet
  assert pdbid.isalnum() and len(pdbid) == 4
  process = asyncio.create_subprocess_shell("phenix.fetch_pdb " + pdbid)
  process()
  return pdbid + ".pdb"

def print_ligand_list(pdb):
  ligand_records = read_pdb(pdb)[2]
  ligands = set()
  for r in ligand_records:
    rn = r['resname']
    ch = r['chain']
    ri = r['resid']
    ligands.add(rn + "_" + ch + "_" + str(ri))
  # ligand_set = set([l.split("_")[0] for l in ligands])
  lig_count = dict()
  for l in ligands:
    if l.split("_")[0] in lig_count.keys():
      lig_count[l.split("_")[0]] += 1
    else:
      lig_count[l.split("_")[0]] = 1
  total_ligand_count = 0
  for l in lig_count:
    print(l, lig_count[l])
    total_ligand_count += lig_count[l]
  print("total ligands:", total_ligand_count)

def filter_conf(old_filename, new_filename, conformer, remove=True):
  if remove:
    condition = lambda record: get_conf(record) == conformer
  else:
    condition = lambda record: get_conf(record) != conformer
  orig = open(old_filename, "rb")
  filtered = open(new_filename, "wb")
  while True:
    try:
      record = orig.next()
      if isatom(record) and condition(record):
        continue
      else:
        filtered.write(record)
    except StopIteration:
      break
  orig.close()
  filtered.close()
  return new_filename

def filter_pdb(old_filename, new_filename, only_chains=[], keep_resnames=[], keep_resnameids=[], keep_reschainids=[]):
  """Extract the atom records from old_filename and write them to new_filename according to some filters.
  Filter the entire PDB first by only_* arguments, if nonempty (produce intersection). Then, if ANY keep_*
  parameter is nonempty, include only records matching ANY keep_* parameter in the final result.

  old_filename: path to a pdb file
  new_filename: path to a filtered pdb file to be written
  only_chains: keep ONLY records in the specified chains (may be a list or a string of single-character chains)
  keep_resnames: keep all matching residues by name (CLA, OEX, TYR)
  keep_resnameids: keep all matching residues by name/id (OEX 601)
  keep_reschainids: keep all matching residues by full record (TYR A 161)
  """
  # filter by chain, if applicable
  if len(only_chains) > 0:
    def check_chain(record):
      ch = record[21:23].strip(" ")
      if get_chain(record) in only_chains:
        return record
  else:
    def check_chain(record):
      return record

  # filter by keep_* params, if applicable
  keeps = keep_resnames + keep_resnameids + keep_reschainids
  if len(keeps) > 0:
    def check_keeps(record):
      rn = get_resname(record)
      ch = get_chain(record)
      ri = get_resid(record)
      rci = get_res_chain_id(record)
      if rn in keep_resnames:
        return record
      if rci[:6].strip(" ") in keep_resnameids:
        return record
      if rci in keep_reschainids:
        return record
  else:
    def check_keeps(record):
      return record

  # add a step to make sure everything is inside the unit cell
  orig = open(old_filename, "rb")
  filtered = open(new_filename, "wb")

  while True:
    try:
      record = orig.next()
      if record[:6] == 'CRYST1':
        cell_a, cell_b, cell_c = get_ucell(record)
        filtered.write(record)
        break
    except StopIteration:
      raise Exception("No crystal record found in %s" % old_filename)

  # def modulo_unit_cell(record):
  #   x, y, z = get_xyz(record)
  #   if x < 0:
  #     x += cell_a
  #   if y < 0:
  #     y += cell_b
  #   if z < 0:
  #     z += cell_c
  #   return write_xyz(record, (x, y, z))

  # apply during a single pass through the pdb
  while True:
    try:
      record = orig.next()
      if isatom(record) and check_chain(record) and check_keeps(record):
        # filtered.write(modulo_unit_cell(record))
        filtered.write(record)
    except StopIteration:
      break

  orig.close()
  filtered.close()
  return new_filename

def average(list_of_xyz_tuples):
  x, y, z = zip(*list_of_xyz_tuples)
  avg_x = sum(x)/len(x)
  avg_y = sum(y)/len(y)
  avg_z = sum(z)/len(z)
  return (avg_x, avg_y, avg_z)

def put_centers(pdb):
  orig = open(pdb, "rb")
  with_centers_name = pdb.split(".pdb")[0] + "_CTR.pdb"
  with_centers = open(with_centers_name, "wb")
  resnames_template = {
  "CLA": ["NA", "NB", "NC", "ND"],
  "LBA": ["NA", "NB", "NC", "ND"],
  "PHO": ["NA", "NB", "NC", "ND"],
  "HEM": ["NA", "NB", "NC", "ND"],
  "OEX": ["MN1", "MN2", "MN3", "MN4", "CA1"],
  "FE":  ["FE"],
  "FE2": ["FE"],
  "PL9": ["C%d" % x for x in range(1, 7)],
  }
  reschainids_template = {
  "TYR D 160": ["OH"],
  "TYR A 161": ["OH"],
  "TYR d 160": ["OH"],
  "TYR a 161": ["OH"],
  }
  data = {}
  pheo_centers = []
  max_record_num = 0
  while True:
    try:
      record = orig.next()
      if isatom(record):
        with_centers.write(write_record(record, record_number=max_record_num))
        max_record_num += 1
      elif isaniso(record):
        with_centers.write(write_record(record, record_number=max_record_num - 1))
      elif isconnect(record):
        pass
      else:
        with_centers.write(record)
      if isatom(record):
        resname = get_resname(record)
        reschainid = get_res_chain_id(record)
        atom = get_atom(record)
        if (resname in resnames_template) or (reschainid in reschainids_template):
          if reschainid not in data:
            if reschainid in reschainids_template:
              data[reschainid] = {atom:None for atom in reschainids_template[reschainid]}
            else:
              data[reschainid] = {atom:None for atom in resnames_template[resname]}
          if atom in data[reschainid]:
            data[reschainid][atom] = get_xyz(record)
          if not (None in data[reschainid].values()):
            center = average(data[reschainid].values())
            with_centers.write(write_center(record, center, record_number=max_record_num))
            max_record_num += 1
            if resname == "PHO":
              pheo_centers.append(center)
              if len(pheo_centers) == 2:
                chain = get_chain(record)
                pcr_reschainid = "PHO %s 000" % ("A" if chain.isupper() else "a")
                with_centers.write(write_center(record, center, record_number=max_record_num,
                  atom="PCR", reschainid=pcr_reschainid))
                max_record_num += 1
            del data[reschainid]
    except StopIteration:
      break
  orig.close()
  with_centers.close()
  return with_centers_name

def get_interatomic_distance(xyz1, xyz2):
  import math
  diffs_sq = [(xyz2[i] - xyz1[i])**2 for i in range(3)]
  distance = math.sqrt(sum(diffs_sq))
  return distance

def get_closest_waters(pdb, chain=None, resname=None, resid=None, atom=None,
  water_chain=None, distance=5, max_waters=1):
  target_positions = []
  waters_distances = {}
  waters = []
  def match_target(record):
    if isatom(record) and \
    (chain is None or get_chain(record) == chain) and \
    (resname is None or get_resname(record) == resname) and \
    (resid is None or get_resid(record) == int(resid)) and \
    (atom is None or get_atom(record) == atom):
      target_positions.append(get_xyz(record))
  def match_water(record):
    if isatom(record) and \
    get_resname(record) in "HOH WAT OOO".split() and \
    (water_chain is None or get_chain(record) == water_chain):
      position = get_xyz(record)
      for target_position in target_positions:
        this_distance = get_interatomic_distance(position, target_position)
        if this_distance <= distance:
          return this_distance
  print("Finding waters within %4.2f Angstroms of target %s" % (distance,
    " ".join([i for i in (chain, resname, resid, atom) if i is not None])))
  model = open(pdb, "rb")
  for record in model:
    if match_target(record):
      target_positions.append(get_xyz(record))
  model.close()
  model = open(pdb, "rb")
  for record in model:
    dist = match_water(record)
    if dist:
      waters_distances[dist] = (get_chain(record),
                                get_resname(record),
                                str(get_resid(record)),
                                get_atom(record))
  dists_available = sorted(waters_distances.keys())
  for i in range(min(max_waters, len(dists_available))):
    d = dists_available[i]
    w = waters_distances[d]
    del waters_distances[d]
    waters.append(w)
    print(" ".join(w))
  return waters

def get_pdb_stats(pdb):
  hier = hierarchical_pdb(read_pdb(pdb))
  resis = set([d['resname'] for d in hier.parsed_pdb])
  print("modified residues:")
  mods = set([" ".join([rec['chain'],str(rec['resid']),rec['resname']])
    for rec in hier.parsed_pdb if rec['resname'] in modified_all])
  print("\n".join(sorted(mods)))
  print("\n")
  # chains with protein residues in them
  num_ch = len(hier.chain_set)
  # protein residues and atoms
  num_prot = 0
  for chain in hier.chain_set:
    for resn in hier.records_by_residue[chain]:
      num_prot += 1
  num_prot_atoms = len(hier.parsed_pdb)
  num_prot_non_H = len([rec for rec in hier.parsed_pdb if rec['elem'] != 'H'])
  total_B_prot = sum([rec['b'] for rec in hier.parsed_pdb])
  # ions
  ions = set([" ".join([rec['chain'],str(rec['resid']),rec['resname']]) for rec in hier.ions])
  num_ions = len(ions)
  # ligand molecules and atoms
  ligands = set([" ".join([rec['chain'],str(rec['resid']),rec['resname']]) for rec in hier.ligand])
  num_lig = len(ligands)
  num_lig_atoms = len(hier.ligand)
  num_lig_non_H = len([rec for rec in hier.ligand if rec['elem'] != 'H'])
  total_B_lig = sum([rec['b'] for rec in hier.ligand])
  # solvent molecules and atoms
  solvents = set([" ".join([rec['chain'],str(rec['resid']),rec['resname']]) for rec in hier.solvent])
  num_solv = len(solvents)
  num_solv_atoms = len(hier.solvent)
  num_solv_non_H = len([rec for rec in hier.solvent if rec['elem'] != 'H'])
  total_B_solv = sum([rec['b'] for rec in hier.solvent])
  # water molecules and atoms
  waters = set([" ".join([rec['chain'],str(rec['resid']),rec['resname']]) for rec in hier.solvent \
    if rec['resname'] in ['HOH', 'WAT', 'H2O', 'OOO']])
  num_wat = len(waters)
  num_wat_atoms = len([rec for rec in hier.solvent \
    if rec['resname'] in ['HOH', 'WAT', 'H2O', 'OOO']])
  num_wat_non_H = len([rec for rec in hier.solvent \
    if rec['resname'] in ['HOH', 'WAT', 'H2O', 'OOO'] and rec['elem'] != 'H'])
  # synthesize
  # avg_B = (total_B_prot + total_B_lig + total_B_solv) / (num_prot_atoms + num_lig_atoms + num_solv_atoms)
  # avg_B = (total_B_prot + total_B_lig) / (num_prot_atoms + num_lig_atoms)
  try:
    avg_B = total_B_prot / num_prot_atoms
  except ZeroDivisionError:
    avg_B = 0
  print("Number of chains:", num_ch)
  print()
  print("Number of protein residues:", num_prot)
  print("Number of protein atoms:", num_prot_atoms)
  print("Number of non-H protein atoms:", num_prot_non_H)
  print()
  print("Number of ions:", num_ions)
  print()
  print("Number of ligands:", num_lig)
  print("Number of ligand atoms:", num_lig_atoms)
  print("Number of non-H ligand atoms:", num_lig_non_H)
  print()
  print("Number of solvent molecules:", num_solv)
  print("Number of solvent atoms:", num_solv_atoms)
  print("Number of non-H solvent atoms:", num_solv_non_H)
  print()
  print("Number of water molecules:", num_wat)
  print("Numebr of water atoms:", num_wat_atoms)
  print("Number of non-H water atoms:", num_wat_non_H)
  print()
  print("Number of total atoms:", num_prot_atoms + num_lig_atoms + num_solv_atoms)
  print("Number of non-H atoms:", num_prot_non_H + num_lig_non_H + num_solv_non_H)
  print("Average B-factor: %5.2f" % avg_B)
  print()

def get_chains(pdb, monomer="both", exclude_chains=""):
  if monomer == "both":
    filter_chain = lambda chain: chain.isalnum()
  elif monomer == 1:
    filter_chain = lambda chain: chain.isupper()
  elif monomer == 2:
    filter_chain = lambda chain: chain.islower()
  else:
    raise Exception("could not interpret monomer keyword")
  chains = set()
  with open(pdb, "rb") as reader:
    while True:
      try:
        line = reader.next()
        if line[:6] in ["HETATM", "ATOM  "]:
          ch = line[21:23].strip(" ")
          if filter_chain(ch) and ch not in exclude_chains:
            chains.add(ch)
      except IndexError:
        continue
      except StopIteration:
        break
  return sorted(chains)

if __name__ == "__main__":
  import sys
  if sys.argv[1] == "print_ligand_list":
    if len(sys.argv) == 3:
      print_ligand_list(sys.argv[2])
    else:
      print("Usage: print_ligand_list model.pdb")
  elif sys.argv[1] == "get_pdb_stats":
    if len(sys.argv) == 3:
      get_pdb_stats(sys.argv[2])
    else:
      print("Usage: get_pdb_stats model.pdb")
  else:
    print("Usage: command *args")
