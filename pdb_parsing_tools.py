from hybrid_36 import hy36encode, hy36decode

# reading pdbs:
def get_chain(record):
  return record[20:22].strip(" ")

def get_resname(record):
  return record[17:20].strip(" ")

def get_resid(record):
  return int(record[22:26].strip(" "))

def get_res_chain_id(record):
  # the spacing is mission-critical!!
  return "%3s %-2s%-3s" % (\
    get_resname(record),
    get_chain(record),
    get_resid(record))

def get_atom(record):
  return record[12:16].strip(" ")

def get_xyz(record):
  x = float(record[30:38].strip(" "))
  y = float(record[38:46].strip(" "))
  z = float(record[46:54].strip(" "))
  return (x, y, z)

def get_occ(record):
  return float(record[56:60].strip(" "))

def get_b(record):
  return float(record[61:66].strip(" "))

def get_conf(record):
  return record[16:17].strip(" ")

def get_ucell(cryst_record):
  assert cryst_record[:6] == 'CRYST1'
  cell_a = float(cryst_record[7:15].strip(" "))
  cell_b = float(cryst_record[16:24].strip(" "))
  cell_c = float(cryst_record[25:33].strip(" "))
  return (cell_a, cell_b, cell_c)

def get_res_chain_ids_from_link(record):
  # the spacing is mission-critical!!
  rci1 = "%3s %-2s%-3s" % (\
    record[17:20].strip(" "),
    record[20:22].strip(" "),
    record[22:26].strip(" ")
    )
  rci2 = "%3s %-2s%-3s" % (\
    record[47:51].strip(" "),
    record[51:52].strip(" "),
    record[52:55].strip(" ")
    )
  return (rci1, rci2)

def get_res_chain_ids_from_disulfide(record):
  # the spacing is mission-critical!!
  rci1 = "%3s %-2s%-3s" % (\
    record[11:15].strip(" "),
    record[15:16].strip(" "),
    record[16:22].strip(" ")
    )
  rci2 = "%3s %-2s%-3s" % (\
    record[25:29].strip(" "),
    record[29:30].strip(" "),
    record[30:36].strip(" ")
    )
  return (rci1, rci2)

def isatom(record):
  return record[:6] in ["HETATM", "ATOM  "]

def isaniso(record):
  return record[:6] == "ANISOU"

def isdisulfide(record):
  return record[:6] == "SSBOND"

def isconnect(record):
  return record[:6] == "CONECT"

def islink(record):
  return record[:6] == "LINK  "

def isremark(record):
  return record[:6] == "REMARK"

def isucell(record):
  return record[:6] in ["CRYST1", "SCALE1", "SCALE2", "SCALE3"]

def isheader(record):
  return record[:6] in ["REMARK", "CRYST1", "SCALE1", "SCALE2", "SCALE3"]

# writing pdbs:
def get_record_num_str(record, override_record_number=None):
  if override_record_number is None:
    return record[6:11]
  else:
    return hy36encode(5, override_record_number) # width, value

def get_record_num_int(record):
  return hy36decode(5, record[6:11]) # width, value

def get_record_from_record_dict(record_dict, record_number):
  return """ATOM  %5s %4s % 3s %-2s% 4s   % 8.3f% 8.3f% 8.3f % 4.2f %5.2f          % 2s\n""" % \
  (get_record_num_str(None, override_record_number=record_number),
   record_dict['atom'], record_dict['resname'], record_dict['chain'],
   record_dict['resid'], record_dict['xyz'][0], record_dict['xyz'][1],
   record_dict['xyz'][2], record_dict['occ'], record_dict['b'], record_dict['elem'])

def write_record(record, record_number=None):
  record_num_str = get_record_num_str(record, override_record_number=record_number)
  new_record = record[:6] + record_num_str + record[11:]
  assert len(record) == len(new_record)
  return new_record

def write_record_with_B(record, B, record_number=None):
  record_num_str = get_record_num_str(record, override_record_number=record_number)
  new_record = record[:6] + record_num_str + record[11:60] + ("%6.2f" % B) + record[66:]
  if not (len(record) == len(new_record)):
    import pdb; pdb.set_trace()
  assert len(record) == len(new_record)
  return new_record

def write_record_as_conformer(record, occ, conf, record_number=None):
  record_num_str = get_record_num_str(record, override_record_number=record_number)
  new_record = record[:6] + record_num_str + record[11:16] + conf + record[17:56] + ("%4.2f" % occ) + record[60:]
  if not (len(record) == len(new_record)):
    import pdb; pdb.set_trace()
  assert len(record) == len(new_record)
  return new_record

def write_xyz(record, xyz_tuple, record_number=None):
  record_num_str = get_record_num_str(record, override_record_number=record_number)
  new_record = record[:6] + record_num_str + record[11:31] + ("%7.3f %7.3f %7.3f" % xyz_tuple) + record[54:]
  assert len(record) == len(new_record)
  return new_record

def write_center(template_record, xyz_tuple, record_number=None, atom="CTR", reschainid=None):
  record_num_str = get_record_num_str(template_record, override_record_number=record_number)
  if reschainid is None:
    reschainid = get_res_chain_id(template_record)
  ctr_record = template_record[:6] + ("%5d" % record_number) + template_record[11] + ("%4s" % atom) + " " + \
    reschainid + template_record[26:30] + ("%8.3f%8.3f%8.3f" % xyz_tuple) + template_record[54:]
  assert len(template_record) == len(ctr_record)
  return ctr_record
