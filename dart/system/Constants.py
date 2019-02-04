# General Constants
blank = ''

# Nomenclature
AAres3 = ['--- ', 'ALA ', 'ASN ', 'ASP ', 'ARG ', 'CYS ', 'GLN ', 'GLU ', 'GLY ', 'HIS ', 'ILE ', 'LEU ', 'LYS ',
          'MET ', 'PRO ', 'PHE ', 'SER ', 'THR ', 'TRP ', 'TYR ', 'VAL ', 'UNK ']
AAres1 = ['  - ', '  A ', '  N ', '  D ', '  R ', '  C ', '  Q ', '  E ', '  G ', '  H ', '  I ', '  L ', '  K ',
          '  M ', '  P ', '  F ', '  S ', '  T ', '  W ', '  Y ', '  V ', '  X ']
NAres3 = ['--- ', 'CYT ', 'THY ', 'GUA ', 'ADE ', 'URI ', 'CYT ', 'THY ', 'GUA ', 'ADE ', 'URI ', 'CYT ', 'GUA ',
          'ADE ']
NAres1 = ['  - ', '  C ', '  T ', '  G ', '  A ', '  U ', ' DC ', ' DT ', ' DG ', ' DA ', ' RU ', ' RC ', ' RG ',
          ' RA ']
IUPAC = [' O5*', ' C5*', ' C4*', ' O4*', ' C3*', ' O3*', ' C2*', ' C1*', ' O2*', ' C5M', ' OP1', ' OP2']
CNS = [" O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " C1'", " O2'", " C7 ", ' O1P', ' O2P']
BASEPAIRS = ['shear', 'stretch', 'stagger', 'buckle', 'proptw', 'opening']
BASEPAIR_STEPS = ['shift', 'slide', 'rise', 'tilt', 'roll', 'twist']
baselist3 = ['ADE', 'THY', 'GUA', 'CYT', 'URI']  # FIX
pairlist1_t = ['A', 'T', 'G', 'C', 'U']  # FIX
pairlist1_c = ['T', 'A', 'C', 'G', 'A']  # FIX
BASELIST3_T = ['ADE', 'THY', 'GUA', 'CYT', 'URI']
BASELIST1_T = ['A', 'T', 'G', 'C', 'U']
BASELIST3_C = ['THY', 'ADE', 'CYT', 'GUA', 'ADE']
BASELIST1_C = ['T', 'A', 'C', 'G', 'A']
PROTTHREE = ['ASN', 'ASP', 'ARG', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PRO', 'PHE', 'SER', 'TRP', 'TYR',
             'VAL', 'UNK']
PROTONE = ['N', 'D', 'R', 'Q', 'E', 'H', 'I', 'L', 'K', 'M', 'P', 'F', 'S', 'W', 'Y', 'V', 'X']
PAIRSTEPS = ['TT/AA', 'AA/TT', 'AT/AT', 'TA/TA', 'GG/CC', 'CC/GG', 'GC/GC', 'CG/CG', 'AC/GT', 'AG/CT', 'TC/GA', 'TG/CA',
             'GA/TC', 'GT/AC', 'CA/TG', 'CT/AG']

# Automodel scaling factors
SHEARSCALE = 1.0
STRETCHSCALE = 1.0
STAGGERSCALE = 1.0
BUCKLESCALE = 1.0
PROPTWSCALE = 1.0
OPENINGSCALE = 1.0
SHIFTSCALE = 0.8
SLIDESCALE = 0.8
RISESCALE = 0.0
TILTSCALE = 0.8
ROLLSCALE = 0.8
TWISTSCALE = 0.8

# Server related constants
MAXMB = 10000.0
MAXMODELS = 100
