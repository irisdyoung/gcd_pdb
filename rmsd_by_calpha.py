from __future__ import division
from six.moves import range

from libtbx import easy_run
from gcd_pdb import gcd_pdb1_pdb2
import os

def get_chains(pdb):
  chain_set = set()
  with open(pdb, "rb") as reader:
    while True:
      try:
        record = reader.next()
        if (record[0:6] == "HETATM" or record[0:4] == "ATOM"):
          chain_set.add(record[21:23].strip(" "))
      except StopIteration:
        break
  return list(sorted(chain_set))

def rmsd_by_calpha(pdb1, pdb2, sele_fixed=None, sele_moving=None, name1=None, name2=None):
  if name1 is None:
    name1 = os.basename(pdb1).split(".pdb")[0]
  if name2 is None:
    name2 = os.basename(pdb2).split(".pdb")[0]
  if sele_fixed is None:
    sele_fixed = "name CA"
  else:
    sele_fixed = "name CA and " + sele_fixed
  if sele_moving is None:
    sele_moving = "name CA"
  else:
    sele_moving = "name CA and " + sele_moving
  command = "phenix.superpose_pdbs %s %s " % (pdb1, pdb2) \
  + "selection_fixed=\"%s\" selection_moving=\"%s\" " % (sele_fixed, sele_moving) \
  + "> superpose_%s_%s.log" % (name1, name2)
  easy_run.fully_buffered(command=command).raise_if_errors()
  return_string = "Superposition of %s %s with %s %s:\n" % (name1, sele_fixed, name2, sele_moving)
  with open("superpose_%s_%s.log" % (name1, name2), "rb") as log:
    for line in log:
      if "RMSD between fixed and moving atoms (start)" in line:
        return_string += "START RMSD: %s\t" % line.strip("\n").split(" ")[-1]
      elif "RMSD between fixed and moving atoms (final)" in line:
        return_string += "FINAL RMSD: %s\n\n" % line.strip("\n").split(" ")[-1]
  return return_string

def rmsd_by_calpha_vs_Suga(pdbs, names):
  with open("rmsd_by_calpha.out", "wb") as out:
    easy_run.fully_buffered(command = "phenix.fetch_pdb 4UB6").raise_if_errors()
    easy_run.fully_buffered(command = "phenix.fetch_pdb 4UB8").raise_if_errors()
    for idx in range(len(pdbs)):
      pdb = pdbs[idx]
      name = names[idx]
      for suga in "4UB6.pdb", "4UB8.pdb":
        gcd_pdb1_pdb2(pdb, suga, ligands=False, outname="gcd_refined.pdb")
        gcd_pdb1_pdb2(suga, pdb, ligands=False, outname="gcd_suga.pdb")
        print "Calculating rmsd of %s with %s (full dimer)" % (name, suga.split(".pdb")[0])
        return_string = rmsd_by_calpha("gcd_refined.pdb", "gcd_suga.pdb", name1=name, name2=suga.split(".pdb")[0])
        print return_string
        out.write(return_string)
        for chain in get_chains("gcd_refined.pdb"):
          print "Calculating rmsd of %s with %s (chain %s)" % (name, suga.split(".pdb")[0], chain)
          return_string = rmsd_by_calpha("gcd_refined.pdb", "gcd_suga.pdb",
                          sele_fixed="chain %s" % chain, sele_moving="chain %s" % chain,
                          name1=name, name2=suga.split(".pdb")[0])
          print return_string
          out.write(return_string)
  os.remove("gcd_suga.pdb_fitted.pdb")
  os.remove("gcd_refined.pdb")
  os.remove("gcd_suga.pdb")
  print "rmsd results written to rmsd_by_calpha.out"

if __name__ == "__main__":
  import sys
  pdb_list = sys.argv[1:]
  if "--withnames" in pdb_list:
    pdb_list.remove("--withnames")
    names = pdb_list[0::2]
    pdbs = pdb_list[1::2]
  else:
    names = map(lambda name: name.split(".pdb")[0], pdb_list)
    pdbs = pdb_list
  print "Generating rmsds for %d models against Suga models" % len(pdbs)
  rmsd_by_calpha_vs_Suga(pdbs, names)
