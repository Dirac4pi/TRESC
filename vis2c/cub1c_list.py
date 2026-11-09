"""Read GTO Molden orbitals and report explicit SCPA screening candidates.

SCPA composition weights are coefficient squares grouped by atom and shell.
Optional PySCF or an explicit AO-overlap matrix provides true spatial matches
and population-defined Pipek-Mezey scores. No SCF, orbital localization, graphics
or cube generation is performed in inspection mode.
"""

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import sys

SYMBOLS = ('X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe '
            'Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In '
            'Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf '
            'Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am '
            'Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()
S_BLOCK = {1, 2, 3, 4, 11, 12, 19, 20, 37, 38, 55, 56, 87, 88}
D_BLOCK = set(range(21, 31)) | set(range(39, 49)) | set(range(72, 81)) | set(range(104, 113))
F_REGION = set(range(57, 72)) | set(range(89, 104))
L_ORDER = 'spdfg'


@dataclass
class Atom:
  number: int
  source_number: int
  symbol: str
  z: int
  coordinates: tuple


@dataclass
class Orbital:
  file_number: int
  number: int
  spin: str
  energy: float
  occupation: float
  symmetry: str
  weights: object
  coefficients: object = None
  metrics: object = None

  @property
  def identifier(self):
    return f'{self.number}{"A" if self.spin == "Alpha" else "B"}'


@dataclass
class Model:
  path: Path
  sha256: str
  atoms: list
  groups: list
  orbitals: list
  basis_count: int
  dimensions: dict
  coordinate_unit: str
  warnings: list
  ao_groups: object = None
  metric_info: object = None


def _clean(line):
  return line.split('#', 1)[0].split('!', 1)[0].strip()


def _float(text):
  value = float(text.replace('D', 'E').replace('d', 'e'))
  if not math.isfinite(value):
    raise ValueError('Nonfinite numerical value in Molden input.')
  return value


def resolve_molden(filename):
  path = Path(filename).expanduser()
  candidates = [path] if path.name.endswith(('.molden.input', '.molden')) else [
    Path(str(path) + '.molden.input'), Path(str(path) + '.molden'),
  ]
  for candidate in candidates:
    if candidate.is_file():
      return candidate.resolve()
  raise FileNotFoundError('No Molden file found: ' + ', '.join(map(str, candidates)))


def _sections(lines):
  sections = []
  for i, line in enumerate(lines):
    match = re.match(r'^\s*\[([^\]]+)\](.*)$', line)
    if match:
      sections.append((match[1].strip().lower(), match[2].strip(), i))
  result = {}
  for i, (name, suffix, start) in enumerate(sections):
    end = sections[i+1][2] if i+1 < len(sections) else len(lines)
    if name in ('atoms', 'gto', 'mo'):
      if name in result:
        raise ValueError(f'Multiple [{name}] sections are ambiguous; provide one structure/orbital set.')
      result[name] = (suffix, lines[start+1:end])
  if 'molden format' not in [s[0] for s in sections]:
    raise ValueError('Missing [Molden Format] header.')
  for required in ('atoms', 'gto', 'mo'):
    if required not in result:
      raise ValueError(f'Analysis requires [{required}]; only real GTO Molden input is supported.')
  declared = {'d': set(), 'f': set(), 'g': set()}
  implicit_f = False
  for name, _, _ in sections:
    tag = name.upper().replace(' ', '')
    if re.fullmatch(r'(?:(?:5D|6D|7F|10F|9G|15G))+', tag):
      for count, letter in re.findall(r'(\d+)([DFG])', tag):
        declared[letter.lower()].add(int(count))
      implicit_f |= tag == '5D'
  if any(len(values) > 1 for values in declared.values()):
    raise ValueError('Conflicting spherical/Cartesian shell declarations.')
  dims = {'s': 1, 'p': 3, 'd': 6, 'f': 7 if implicit_f else 10, 'g': 15}
  for letter, values in declared.items():
    if values:
      dims[letter] = next(iter(values))
  return result, dims


def _atoms(section):
  unit, lines = section
  unit = unit.strip().lower()
  if unit not in ('angs', 'au'):
    raise ValueError('[Atoms] must declare Angs or AU; coordinate units are not guessed.')
  atoms, ids = [], set()
  for line in lines:
    fields = _clean(line).split()
    if not fields:
      continue
    if len(fields) != 6:
      raise ValueError('Each [Atoms] record must have label, number, Z, x, y, z.')
    serial, z = int(fields[1]), int(fields[2])
    if serial < 1 or serial in ids or not 1 <= z <= 118:
      raise ValueError('Atom numbers must be positive/unique and Z must be in 1..118.')
    ids.add(serial)
    atoms.append(Atom(len(atoms)+1, serial, SYMBOLS[z], z, tuple(_float(v) for v in fields[3:])))
  if not atoms:
    raise ValueError('No atoms were found.')
  return atoms, 'angstrom' if unit == 'angs' else 'bohr'


def _basis(section, atoms, dims):
  lines = [_clean(line) for line in section[1] if _clean(line)]
  centers = {a.source_number: a.number for a in atoms}
  groups, indices, ao_group = [], {}, []
  i, center, shells = 0, None, 0
  while i < len(lines):
    fields = lines[i].split()
    if re.fullmatch(r'\d+', fields[0]):
      if len(fields) > 2 or (len(fields) == 2 and fields[1] != '0'):
        raise ValueError('Invalid GTO atom-center header.')
      serial = int(fields[0])
      if serial not in centers:
        raise ValueError(f'GTO center {serial} is not an [Atoms] source number.')
      center = centers[serial]
      i += 1
      continue
    shell = fields[0].lower()
    if center is None or shell not in ('s', 'p', 'd', 'f', 'g', 'sp') or len(fields) not in (2, 3):
      raise ValueError(f'Unsupported/invalid GTO shell header: {lines[i]}')
    count = int(fields[1])
    if count < 1:
      raise ValueError('A contracted shell must contain at least one primitive.')
    if len(fields) == 3:
      _float(fields[2])
    for j in range(count):
      if i+1+j >= len(lines):
        raise ValueError('Truncated GTO primitive list.')
      primitive = lines[i+1+j].split()
      expected = 3 if shell == 'sp' else 2
      if len(primitive) != expected or _float(primitive[0]) <= 0:
        raise ValueError('Unsupported primitive/contraction layout in [GTO].')
      for value in primitive[1:]:
        _float(value)
    for letter in ('s', 'p') if shell == 'sp' else (shell,):
      key = (center, letter)
      if key not in indices:
        indices[key] = len(groups)
        groups.append(key)
      ao_group.extend([indices[key]] * dims[letter])
    shells += 1
    i += count+1
  if not ao_group:
    raise ValueError('No supported GTO basis functions were found.')
  return groups, ao_group, shells


def _orbitals(section, ao_group, group_count, allow_sparse):
  import numpy as np
  orbitals, counters, metadata, coefficients = [], {'Alpha': 0, 'Beta': 0}, {}, {}
  basis_count = len(ao_group)
  missing_count = 0

  def finish():
    nonlocal metadata, coefficients, missing_count
    if not metadata and not coefficients:
      return
    if not {'ene', 'spin', 'occup'}.issubset(metadata) or not coefficients:
      raise ValueError('Each MO needs Ene, Spin, Occup and real AO coefficients.')
    if len(coefficients) != basis_count:
      if not allow_sparse:
        raise ValueError(
          f'MO {len(orbitals)+1} supplies {len(coefficients)}/{basis_count} AO entries. '
          'Check shell declarations; use --allow-sparse only for verified omitted-zero coefficients.'
        )
      missing_count += 1
    values = np.zeros(basis_count)
    for number, value in coefficients.items():
      values[number-1] = value
    scale = float(np.max(np.abs(values)))
    if scale == 0:
      raise ValueError('An all-zero MO cannot have normalized SCPA weights.')
    weights = (values / scale) ** 2
    weights /= weights.sum()
    grouped = np.bincount(np.asarray(ao_group), weights=weights, minlength=group_count)
    spin = metadata['spin'].strip().lower()
    if spin not in ('alpha', 'beta'):
      raise ValueError('Only Alpha/Beta real-orbital spin labels are supported.')
    spin = spin.capitalize()
    occupation = _float(metadata['occup'])
    if not 0 <= occupation <= 2.000001:
      raise ValueError('MO occupations must lie in [0, 2] for this analysis.')
    counters[spin] += 1
    orbitals.append(Orbital(len(orbitals)+1, counters[spin], spin,
                            _float(metadata['ene']), occupation,
                            metadata.get('sym', ''), grouped, coefficients=values))
    metadata, coefficients = {}, {}

  for line in section[1]:
    text = _clean(line)
    if not text:
      continue
    match = re.match(r'^(Sym|Ene|Spin|Occup)\s*=\s*(.*)$', text, re.IGNORECASE)
    if match:
      key, value = match[1].lower(), match[2].strip()
      if coefficients and key in ('sym', 'ene'):
        finish()
      if key in metadata:
        raise ValueError(f'Duplicate MO metadata field: {key}')
      if coefficients:
        raise ValueError('MO metadata appeared after coefficients without a new orbital header.')
      metadata[key] = value
    else:
      fields = text.split()
      if len(fields) != 2 or not fields[0].isdigit():
        raise ValueError(f'Expected AO index and one real coefficient, not: {text}')
      number = int(fields[0])
      if not 1 <= number <= basis_count or number in coefficients:
        raise ValueError(f'Invalid/duplicate AO index: {number}')
      coefficients[number] = _float(fields[1])
  finish()
  if not orbitals:
    raise ValueError('No molecular orbitals were found.')
  return orbitals, missing_count


def read_molden(filename, allow_sparse=False):
  path = resolve_molden(filename)
  raw = path.read_bytes()
  lines = raw.decode('utf-8-sig').splitlines()
  sections, dims = _sections(lines)
  atoms, unit = _atoms(sections['atoms'])
  groups, ao_group, _ = _basis(sections['gto'], atoms, dims)
  orbitals, sparse = _orbitals(sections['mo'], ao_group, len(groups), allow_sparse)
  warnings = []
  if sparse:
    warnings.append(f'{sparse} sparse MO records: omitted AO indices explicitly treated as zero.')
  if any(letter in ('d', 'f', 'g') and dims[letter] in (6, 10, 15) for _, letter in groups):
    warnings.append('Cartesian high-l SCPA shell weights can depend on orientation/basis convention.')
  warnings.append('SCPA is basis-character screening, not strict core/valence/Rydberg classification.')
  warnings.append('A small frontier-rank window cannot certify the absence of diffuse or polarization character.')
  return Model(path, hashlib.sha256(raw).hexdigest(), atoms, groups, orbitals,
                len(ao_group), dims, unit, warnings, ao_groups=ao_group)


# -----------------------------------------------------------------------------
# True spatial overlaps and population-defined Pipek-Mezey orbital diagnostics

class MetricsDependencyError(ImportError):
  """The optional integral backend is not installed, not a numerical failure."""


def _unavailable_metrics(model, reason, status='unavailable', channel_map=None):
  lookup = {mo.identifier:mo for mo in model.orbitals}
  for mo in model.orbitals:
    peer_id = f'{mo.number}{"B" if mo.spin=="Alpha" else "A"}'
    peer = lookup.get(peer_id)
    partner = None if peer is None else {
      'MO_ID':peer.identifier,'spin':peer.spin,'occupation':peer.occupation,'energy':peer.energy,
      'relative':channel_map[peer.spin]['labels'][peer.identifier] if channel_map else None,
    }
    mo.metrics = {
      'status':status,'reason':reason,'same_index_partner':partner,
      'ab_abs_overlap':None,'ab_signed_overlap':None,
      'ab_raw_overlap':None,'ab_self_ratio':None,'ab_self_ratio_signed':None,
      'spatial_self_overlap':None,'ab_peer_self_overlap':None,
      'pm2_lowdin':None,'effective_atoms_lowdin':None,'pm2_mulliken':None,
    }
  return {
    'status':status,'reason':reason,'warnings':[reason],
    'population_method':'Lowdin (symmetric orthogonalization), exponent 2',
    'matching_method':'Normalized spatial |Ca.T S Cb| for identical within-spin indices',
  }


def _source_basis_counts(model):
  from collections import Counter
  return dict(Counter(model.groups[int(group)] for group in model.ao_groups))


def _pyscf_blocks(model, loaded):
  """Normalize the documented combined-matrix or Alpha/Beta block loader output."""
  import numpy as np
  if len(loaded) != 6:
    raise ValueError('Unexpected pyscf.tools.molden.load return layout; expected six items.')
  mol, energies, coefficients, occupations, _, spins = loaded
  source = {spin:[mo for mo in model.orbitals if mo.spin==spin] for spin in ('Alpha','Beta')}
  source = {spin:rows for spin,rows in source.items() if rows}

  def spin_names(values, count):
    try:
      flat = np.asarray(values).reshape(-1).tolist()
    except (ValueError,TypeError):
      return None
    if len(flat) != count:
      return None
    names = []
    for value in flat:
      key = str(value).strip().lower()
      if key in ('alpha','0','0.0'):
        names.append('Alpha')
      elif key in ('beta','1','1.0'):
        names.append('Beta')
      else:
        return None
    return names

  if isinstance(coefficients,np.ndarray) and coefficients.ndim==2:
    matrix = np.asarray(coefficients,dtype=float)
    energy = np.asarray(energies,dtype=float).reshape(-1)
    occ = np.asarray(occupations,dtype=float).reshape(-1)
    if len(energy)!=matrix.shape[1] or len(occ)!=matrix.shape[1]:
      raise ValueError('PySCF orbital metadata length mismatch.')
    names = spin_names(spins,matrix.shape[1])
    if names is None:
      if len(source)==1:
        names = [next(iter(source))]*matrix.shape[1]
      else:
        raise ValueError('Combined PySCF MO matrix has no unambiguous spin labels.')
    blocks = {}
    for spin in source:
      mask = np.array([name==spin for name in names])
      blocks[spin] = (matrix[:,mask],energy[mask],occ[mask])
  elif isinstance(coefficients,(tuple,list,np.ndarray)) and len(coefficients) in (1,2):
    ordered = ['Alpha','Beta'] if len(coefficients)==2 else list(source)
    if set(ordered)!=set(source):
      raise ValueError('PySCF/source spin-channel layout differs.')
    if len(energies)!=len(coefficients) or len(occupations)!=len(coefficients):
      raise ValueError('Blocked PySCF coefficients require matching energy/occupation blocks.')
    blocks = {
      spin:(np.asarray(coefficients[i],dtype=float),
            np.asarray(energies[i],dtype=float).reshape(-1),
            np.asarray(occupations[i],dtype=float).reshape(-1))
      for i,spin in enumerate(ordered)
    }
  else:
    raise ValueError('Unsupported PySCF Molden MO layout; do not guess orbital ordering.')
  columns = {}
  for spin,rows in source.items():
    matrix,energy,occ = blocks[spin]
    if matrix.shape!=(model.basis_count,len(rows)) or len(energy)!=len(rows) or len(occ)!=len(rows):
      raise ValueError(f'PySCF/source MO dimensions differ in {spin}.')
    if (not np.all(np.isfinite(matrix))
        or not np.allclose(energy,[mo.energy for mo in rows],rtol=1e-9,atol=1e-8)
        or not np.allclose(occ,[mo.occupation for mo in rows],rtol=1e-9,atol=1e-8)):
      raise ValueError(f'PySCF/source orbital metadata differ in {spin}; refusing a guessed correspondence.')
    for i,mo in enumerate(rows):
      columns[mo.identifier] = matrix[:,i]
  return mol,np.column_stack([columns[mo.identifier] for mo in model.orbitals])


def _pyscf_metric_inputs(model, threads):
  """Use PySCF's Molden loader so AO reorder/Cartesian normalization are handled together."""
  try:
    import pyscf
  except ModuleNotFoundError as exc:
    if exc.name=='pyscf':
      raise MetricsDependencyError(
        'True AO-metric diagnostics need PySCF or --ao-overlap with a matching Molden-order S matrix. '
        'PySCF is not installed; no coefficient-cosine/SCPA substitutes are printed.'
      ) from exc
    raise
  import numpy as np
  from collections import Counter
  from pyscf import lib
  from pyscf.tools import molden
  previous = lib.num_threads()
  try:
    lib.num_threads(threads)
    loaded = molden.load(str(model.path),verbose=0)
    mol,coefficients = _pyscf_blocks(model,loaded)
    if mol.natm!=len(model.atoms) or mol.nao_nr()!=model.basis_count:
      raise ValueError('PySCF/source atom or AO counts differ.')
    if [int(mol.atom_charge(i)) for i in range(mol.natm)]!=[atom.z for atom in model.atoms]:
      raise ValueError('PySCF/source element order differs.')
    coords = np.asarray([atom.coordinates for atom in model.atoms],dtype=float)
    if model.coordinate_unit=='angstrom':
      coords /= lib.param.BOHR
    if not np.allclose(mol.atom_coords(),coords,rtol=0,atol=1e-7):
      raise ValueError('PySCF/source physical coordinates differ.')
    native_counts = Counter()
    for shell in range(mol.nbas):
      l = int(mol.bas_angular(shell))
      if l>=len(L_ORDER):
        raise ValueError('PySCF basis exceeds the supported angular-shell map.')
      components = (l+1)*(l+2)//2 if mol.cart else 2*l+1
      native_counts[(int(mol.bas_atom(shell))+1,L_ORDER[l])] += components*int(mol.bas_nctr(shell))
    if dict(native_counts)!=_source_basis_counts(model):
      raise ValueError('PySCF/source atom-shell dimensions differ; mixed/convention-dependent basis is unsupported.')
    overlap = np.asarray(mol.intor_symmetric('int1e_ovlp'),dtype=float)
    atom_of_ao = np.full(model.basis_count,-1,dtype=int)
    for atom,(_,_,first,last) in enumerate(mol.aoslice_by_atom()):
      atom_of_ao[int(first):int(last)] = atom
    if np.any(atom_of_ao<0):
      raise ValueError('Incomplete PySCF AO-to-atom map.')
    return overlap,coefficients,atom_of_ao,{
      'backend':'pyscf','pyscf_version':pyscf.__version__,
      'coefficient_order':'PySCF loader order, validated against source MO metadata',
      'ao_convention':'PySCF loader handles angular reorder and Cartesian coefficient normalization',
      'integrals':'int1e_ovlp; no SCF or orbital optimization performed',
      'integral_threads':threads,
    }
  finally:
    lib.num_threads(previous)


def _external_metric_inputs(model, filename):
  """Read an explicitly supplied overlap matrix in ORIGINAL Molden AO order."""
  import numpy as np
  path = Path(filename).expanduser().resolve()
  if not path.is_file():
    raise FileNotFoundError(path)
  if path.suffix.lower()=='.npy':
    overlap = np.load(path,allow_pickle=False)
  elif path.suffix.lower() in ('.txt','.dat'):
    overlap = np.loadtxt(path)
  else:
    raise ValueError('--ao-overlap supports .npy or whitespace numeric .txt/.dat matrices.')
  coefficients = np.column_stack([mo.coefficients for mo in model.orbitals])
  atom_of_ao = np.array([model.groups[int(group)][0]-1 for group in model.ao_groups])
  return overlap,coefficients,atom_of_ao,{
    'backend':'external_overlap','overlap_path':str(path),
    'overlap_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
    'coefficient_order':'Original Molden AO order',
    'ao_convention':'Caller explicitly declares the same AO order, normalization and coordinates as the source',
  }


def _evaluate_metric_inputs(model, channel_map, overlap, coefficients, atom_of_ao,
                            info, norm_tolerance=None):
  """Accept nonunit MO scales and compute explicit directed and symmetric ratios.

  norm_tolerance is accepted for call compatibility but is not applied.
  Input MO columns are not rescaled, orthogonalized, mixed or optimized.
  A valid real AO metric and positive finite self overlaps remain mandatory.
  """
  import numpy as np
  S,C = np.asarray(overlap),np.asarray(coefficients)
  if np.iscomplexobj(S) or np.iscomplexobj(C):
    raise ValueError('Only real AO overlap/coefficients are supported.')
  S,C = np.asarray(S,dtype=float),np.asarray(C,dtype=float)
  n = model.basis_count
  if S.shape!=(n,n) or C.shape!=(n,len(model.orbitals)):
    raise ValueError('Overlap/coefficient dimensions do not match the source.')
  if not np.all(np.isfinite(S)) or not np.all(np.isfinite(C)):
    raise ValueError('AO overlap and coefficients must be finite.')
  if not np.allclose(S,S.T,rtol=1e-9,atol=1e-10):
    raise ValueError('AO overlap is not symmetric.')
  S = (S+S.T)/2
  if np.any(np.diag(S)<=0):
    raise ValueError('AO overlap has nonpositive diagonal entries.')
  eigenvalues,eigenvectors = np.linalg.eigh(S)
  if eigenvalues[-1]<=0 or eigenvalues[0] < -1e-9*max(1.,float(eigenvalues[-1])):
    raise ValueError('AO overlap is not positive semidefinite.')
  clipped = int(np.count_nonzero(eigenvalues<0))
  sqrtS = (eigenvectors*np.sqrt(np.maximum(eigenvalues,0))) @ eigenvectors.T
  SC = S @ C
  norms = np.einsum('ij,ij->j',C,SC)
  if not np.all(np.isfinite(norms)) or np.any(norms<=0):
    raise ValueError('A spatial MO has zero, negative or nonfinite self overlap.')
  # No comparison with one and no unit-norm rejection are performed.
  roots = np.sqrt(norms)
  lowdin_coefficients = sqrtS @ C
  atom_of_ao = np.asarray(atom_of_ao,dtype=int)
  if atom_of_ao.shape!=(n,) or np.any(atom_of_ao<0) or np.any(atom_of_ao>=len(model.atoms)):
    raise ValueError('Invalid metric AO-to-atom assignment.')
  qL = np.zeros((len(model.atoms),len(model.orbitals)))
  qM = np.zeros_like(qL)
  for atom in range(len(model.atoms)):
    mask = atom_of_ao==atom
    qL[atom] = np.sum(lowdin_coefficients[mask]**2,axis=0)
    qM[atom] = np.sum(C[mask]*SC[mask],axis=0) / norms
  totals = qL.sum(axis=0)
  if (not np.all(np.isfinite(totals)) or np.any(totals<=0)
      or np.any(np.abs(totals/norms-1)>1e-5)):
    raise ValueError('Lowdin metric consistency failed; verify the AO overlap conditioning.')
  qL /= totals
  pmL,pmM = np.sum(qL**2,axis=0),np.sum(qM**2,axis=0)
  if not np.all(np.isfinite(qM)) or not np.all(np.isfinite(pmM)) or np.any(pmL<=0):
    raise ValueError('Nonfinite or invalid atomic populations/PM scores.')
  warnings = []
  if np.count_nonzero(qM < -1e-8):
    warnings.append('Some Mulliken orbital populations are negative; their PM inverse is not an atom-count scale.')
  column = {mo.identifier:i for i,mo in enumerate(model.orbitals)}
  for i,mo in enumerate(model.orbitals):
    mo.metrics = {
      'status':'computed','population_method':'Lowdin','pm_exponent':2,
      'pm2_lowdin':float(pmL[i]),'effective_atoms_lowdin':float(1/pmL[i]),
      'pm2_mulliken':float(pmM[i]),
      'spatial_self_overlap':float(norms[i]),
      'spatial_norm_before_normalization':float(norms[i]),
      'atomic_populations_lowdin':qL[:,i].tolist(),
      'atomic_populations_mulliken':qM[:,i].tolist(),
      'same_index_partner':None,'ab_peer_self_overlap':None,
      'ab_raw_overlap':None,'ab_self_ratio_signed':None,'ab_self_ratio':None,
      'ab_abs_overlap':None,'ab_signed_overlap':None,
    }
    peer = f'{mo.number}{"B" if mo.spin=="Alpha" else "A"}'
    if peer not in column:
      continue
    j = column[peer]
    other = model.orbitals[j]
    raw = float(np.dot(C[:,i],SC[:,j]))
    directed = raw / norms[i]
    symmetric = (raw / roots[i]) / roots[j]
    if not math.isfinite(raw) or not math.isfinite(directed) or not math.isfinite(symmetric):
      raise ValueError('Nonfinite spatial overlap or ratio.')
    if abs(symmetric)>1+1e-7:
      raise ValueError('Symmetric normalized overlap violates the AO metric bound; verify S and AO conventions.')
    mo.metrics.update({
      'same_index_partner':{
        'MO_ID':peer,'spin':other.spin,'occupation':other.occupation,'energy':other.energy,
        'relative':channel_map[other.spin]['labels'][peer],
      },
      'ab_peer_self_overlap':float(norms[j]),'ab_raw_overlap':raw,
      'ab_self_ratio_signed':float(directed),'ab_self_ratio':float(abs(directed)),
      'ab_abs_overlap':min(1.0,abs(symmetric)),
      'ab_signed_overlap':max(-1.0,min(1.0,float(symmetric))),
    })
  sample_checks = {}
  for spin,channel in channel_map.items():
    ids = [column[mo.identifier] for mo in channel['members']]
    sampled = np.unique(np.linspace(0,len(ids)-1,min(64,len(ids)),dtype=int))
    selected = [ids[k] for k in sampled]
    gram = C[:,selected].T @ SC[:,selected]
    gram = gram / roots[selected,None] / roots[None,selected]
    np.fill_diagonal(gram,0)
    maximum = float(np.max(np.abs(gram)))
    sample_checks[spin] = {'sample_size':len(selected),'maximum_off_diagonal_overlap':maximum}
    if maximum>1e-3:
      warnings.append(f'{spin}: sampled normalized same-spin overlap is {maximum:.3g}; check source orbital type/AO convention before physical interpretation.')
  return {
    **info,'status':'computed','warnings':warnings,
    'matching_method':'abs(Ca.T S Cb)/sqrt((Ca.T S Ca)*(Cb.T S Cb)); spatial part only',
    'directed_ratio_method':'abs(C_current.T S C_peer)/(C_current.T S C_current); scale-dependent, not bounded by one',
    'raw_overlap_method':'C_current.T S C_peer; signed, unnormalized spatial integral',
    'matching_pairing':'Identical original within-spin MO numbers; NOT best-match search or UCO',
    'population_method':'Lowdin symmetric AO metric; atomic fractions divided by the MO self norm',
    'pm2_lowdin_formula':'sum_A(q_A_Lowdin**2)',
    'effective_atoms_lowdin_formula':'1/PM2_Lowdin; equivalent participating atoms, not a length',
    'pm2_mulliken_formula':'sum_A(q_A_Mulliken**2)',
    'occupation_weighted':False,'localization_optimization_performed':False,
    'MO_columns_rescaled_or_rotated':False,'unit_norm_gate_enabled':False,
    'overlap_min_eigenvalue':float(eigenvalues[0]),
    'tiny_negative_eigenvalues_clipped':clipped,
    'input_self_overlap_min':float(norms.min()),'input_self_overlap_max':float(norms.max()),
    'orthogonality_sample':sample_checks,
    'interpretation_note':'Allowing column-scale differences does not certify AO ordering, normalization conventions or the physical correctness of S.',
  }


def compute_orbital_metrics(model, channel_map, args):
  """Never replace the requested true S-metric with SCPA/cosine heuristics."""
  if args.orbital_metrics=='off':
    if args.ao_overlap:
      raise ValueError('--ao-overlap cannot be combined with --orbital-metrics off.')
    info = _unavailable_metrics(model,'AO-metric diagnostics explicitly disabled.','disabled',channel_map)
  else:
    try:
      if args.ao_overlap:
        if args.orbital_metrics=='pyscf':
          raise ValueError('Choose auto with --ao-overlap, or omit the matrix for the PySCF backend.')
        data = _external_metric_inputs(model,args.ao_overlap)
      else:
        data = _pyscf_metric_inputs(model,args.metric_threads)
      overlap,coefficients,atom_of_ao,provenance = data
      info = _evaluate_metric_inputs(model,channel_map,overlap,coefficients,atom_of_ao,
                                      provenance,args.metric_norm_tol)
    except MetricsDependencyError as exc:
      if args.orbital_metrics=='pyscf':
        raise RuntimeError(str(exc)) from exc
      info = _unavailable_metrics(model,str(exc),channel_map=channel_map)
  model.metric_info = info
  return info


def _metric_record(mo):
  return mo.metrics if mo.metrics is not None else {
    'status':'unavailable','reason':'Metrics were not evaluated for this object.',
    'same_index_partner':None,'ab_abs_overlap':None,
    'ab_raw_overlap':None,'ab_self_ratio':None,'ab_self_ratio_signed':None,
    'spatial_self_overlap':None,'ab_peer_self_overlap':None,
    'pm2_lowdin':None,'effective_atoms_lowdin':None,'pm2_mulliken':None,
  }


def _metric_header(show_pm=False):
  """Distinguish the directed self-ratio from the symmetric 0..1 match."""
  return 'ABpeer   AB/self   AB|S|' + ('    PM2_L  N_eff_L' if show_pm else '')


def _metric_columns(mo, show_pm=False):
  metric = _metric_record(mo)
  peer = metric.get('same_index_partner')
  peer = peer['MO_ID'] if peer else '--'
  directed = metric.get('ab_self_ratio')
  overlap = metric.get('ab_abs_overlap')
  directed = f'{directed:.5g}' if directed is not None else 'NA'
  overlap = f'{overlap:.4f}' if overlap is not None else 'NA'
  text = f'{peer:<7} {directed:>9} {overlap:>7}'
  if show_pm:
    pm,number = metric.get('pm2_lowdin'),metric.get('effective_atoms_lowdin')
    pm = f'{pm:.4f}' if pm is not None else 'NA'
    number = f'{number:.2f}' if number is not None else 'NA'
    text += f' {pm:>8} {number:>8}'
  return text


def _metric_text(mo, show_pm=False):
  return _metric_header(show_pm)+': '+_metric_columns(mo,show_pm).strip()


# -----------------------------------------------------------------------------
# Channel ranks, selections and target-channel policies

def parse_atoms(text, atom_count):
  result = []
  for item in text.split(','):
    item = item.strip()
    match = re.fullmatch(r'(\d+)(?:\s*-\s*(\d+))?', item)
    if not match:
      raise ValueError('Use 1-based atom selections such as 1,3-8,12.')
    first, last = int(match[1]), int(match[2] or match[1])
    if first < 1 or last < first or last > atom_count:
      raise ValueError('Atom selections must be positive, ascending and in range.')
    result.extend(range(first, last+1))
  if len(result) != len(set(result)):
    raise ValueError('Duplicate atoms are not allowed in one selection.')
  return result


def window_pair(values):
  values = list(values)
  if len(values) not in (1, 2) or any(v < 0 for v in values):
    raise ValueError('--window requires one or two nonnegative integer offsets.')
  return (values[0], values[0]) if len(values) == 1 else tuple(values)


def target_channels(atom, profile):
  if atom.z in S_BLOCK:
    return 's'
  if atom.z in D_BLOCK:
    return 'd' if profile == 'dominant' else 'sd'
  if atom.z in F_REGION:
    return ''
  return 'p' if profile == 'dominant' else 'sp'


def channel_policy(model, profile, overrides):
  policy = {a.number: target_channels(a, profile) for a in model.atoms}
  changed = set()
  if overrides:
    for entry in overrides.split(';'):
      if entry.count(':') != 1:
        raise ValueError('--target-l syntax: "1:ds;2:sp" (atom selection:channels).')
      selection, letters = (part.strip() for part in entry.split(':'))
      letters = letters.lower()
      if not letters or any(c not in L_ORDER for c in letters) or len(set(letters)) != len(letters):
        raise ValueError('Target channels must be unique letters from s,p,d,f,g; ds means d+s.')
      for number in parse_atoms(selection, len(model.atoms)):
        if number in changed:
          raise ValueError('An atom has more than one --target-l override.')
        policy[number] = ''.join(c for c in L_ORDER if c in letters)
        changed.add(number)
  return policy


def channels(model, occ_tol):
  result, warnings = {}, []
  for spin in ('Alpha', 'Beta'):
    members = [mo for mo in model.orbitals if mo.spin == spin]
    if not members:
      continue
    occupied = sorted((mo for mo in members if mo.occupation > occ_tol), key=lambda mo:(mo.energy, mo.number))
    virtual = sorted((mo for mo in members if mo.occupation <= occ_tol), key=lambda mo:(mo.energy, mo.number))
    labels = {mo.identifier: 'HOMO' if i == len(occupied)-1 else f'HOMO-{len(occupied)-1-i}'
              for i, mo in enumerate(occupied)}
    labels.update({mo.identifier: 'LUMO' if i == 0 else f'LUMO+{i}' for i, mo in enumerate(virtual)})
    result[spin] = {'members': members, 'occupied': occupied, 'virtual': virtual, 'labels': labels}
    if not occupied or not virtual:
      warnings.append(f'{spin}: one frontier endpoint is unavailable in the supplied orbital set.')
    if occupied and virtual and occupied[-1].energy >= virtual[0].energy:
      warnings.append(f'{spin}: occupied/virtual energies touch or overlap; labels are rank-based, not proof of an Aufbau state.')
    if len({mo.energy for mo in members}) == 1 and len(members) > 1:
      warnings.append(f'{spin}: all energies are identical; ranks only use file order to break ties. Check noncanonical/placeholder energies.')
    max_occ = 2.0 if max(mo.occupation for mo in members) > 1+occ_tol else 1.0
    if any(occ_tol < mo.occupation < max_occ-occ_tol for mo in members):
      warnings.append(f'{spin}: partial occupations detected; occupied means Occup > {occ_tol:g}, not an assumed closed-shell classification.')
    energies = [mo.energy for mo in members]
    if energies != sorted(energies):
      warnings.append(f'{spin}: file order is not energy order. MO_ID stays in original within-spin file order; Ref uses energy ranks.')
  if len(result) == 1:
    warnings.append('Only one spin-labelled orbital set is present. No second set is fabricated; do not infer closed-shell status from that alone.')
  return result, warnings


def window_members(channel, window):
  below, above = window
  return channel['occupied'][-below-1:] + channel['virtual'][:above+1]


def _window_warnings(channel, window):
  """Report a count-window boundary that splits an energy-degenerate group."""
  result = []
  occupied, virtual = channel['occupied'], channel['virtual']
  nocc, nvir = window[0]+1, window[1]+1
  if len(occupied) > nocc and math.isclose(occupied[-nocc].energy, occupied[-nocc-1].energy,
                                          rel_tol=0, abs_tol=1e-8):
    result.append('Occupied window boundary splits a tied/near-tied energy group; widen the count window if needed.')
  if len(virtual) > nvir and math.isclose(virtual[nvir-1].energy, virtual[nvir].energy,
                                        rel_tol=0, abs_tol=1e-8):
    result.append('Virtual window boundary splits a tied/near-tied energy group; widen the count window if needed.')
  return result


def _group_sum(mo, model, atoms, policy=None):
  atom_set = set(atoms)
  return float(sum(w for (atom, letter), w in zip(model.groups, mo.weights)
                    if atom in atom_set and (policy is None or letter in policy[atom])))


def _parse_requirements(text, atoms, policy, model):
  """Parse AND constraints: atom-selection:channels:whole-MO-minimum."""
  selected = set(atoms)
  result, seen = [], set()
  for item in text.split(';'):
    parts = [part.strip() for part in item.split(':')]
    if len(parts) != 3:
      raise ValueError('--require-target syntax: "1-2:d:0.10;3:p:0.10".')
    centers = parse_atoms(parts[0],len(model.atoms))
    letters = parts[1].lower()
    if not letters or any(l not in L_ORDER for l in letters) or len(set(letters)) != len(letters):
      raise ValueError('Requirement channels must be unique letters from s,p,d,f,g.')
    letters = ''.join(l for l in L_ORDER if l in letters)
    minimum = float(parts[2])
    if not math.isfinite(minimum) or not 0 < minimum <= 1:
      raise ValueError('Each joint requirement minimum must be a fraction in (0,1].')
    if not set(centers).issubset(selected):
      raise ValueError('Joint requirement atoms must belong to the selected fragment.')
    if any(not set(letters).issubset(set(policy[n])) for n in centers):
      raise ValueError('Requirement channels must be part of --target-profile/--target-l for every listed atom.')
    key = (tuple(sorted(centers)),letters)
    if key in seen:
      raise ValueError('Duplicate joint requirement.')
    seen.add(key)
    result.append({'atoms':centers,'channels':letters,'minimum':minimum})
  return result


def _joint_requirements(model, atoms, policy, args):
  """Keep summed enrichment separate from explicitly requested co-participation."""
  if args.fragment_mode == 'sum':
    if args.require_target:
      raise ValueError('--require-target requires --fragment-mode joint.')
    return []
  if args.require_target:
    return _parse_requirements(args.require_target,atoms,policy,model)
  return [{'atoms':[n],'channels':policy[n],'minimum':args.per_atom_min} for n in atoms]


def _requirement_weight(mo, model, requirement):
  centers = set(requirement['atoms'])
  letters = set(requirement['channels'])
  return float(math.fsum(float(w) for (a,l),w in zip(model.groups,mo.weights)
                          if a in centers and l in letters))


def _target_details(mo, model, atoms, policy, requirements=()):
  """Return every requested atom/channel, including zero terms and their sums."""
  weights = {key:float(w) for key,w in zip(model.groups,mo.weights)}
  breakdown = [
    {'atom':n,'element':model.atoms[n-1].symbol,'angular_shell':letter,
      'weight':weights.get((n,letter),0.0),'basis_present':(n,letter) in weights}
    for n in atoms for letter in L_ORDER if letter in policy[n]
  ]
  target = math.fsum(row['weight'] for row in breakdown)
  fragment = _group_sum(mo,model,atoms)
  by_atom = [{'atom':n,'element':model.atoms[n-1].symbol,
              'weight':math.fsum(row['weight'] for row in breakdown if row['atom']==n)}
              for n in atoms]
  checks = []
  for requirement in requirements:
    value = _requirement_weight(mo,model,requirement)
    checks.append({**requirement,'weight':value,'passed':value+1e-12>=requirement['minimum']})
  return {
    'MO_ID':mo.identifier,'orbital_metrics':_metric_record(mo),
    'fragment_weight':fragment,'target_weight':target,
    'fragment_non_target_weight':max(0.0,fragment-target),
    'target_breakdown':breakdown,'target_by_atom':by_atom,'joint_checks':checks,
  }


def select_candidates(model, channel, atoms, policy, window, target_min, fragment_min=0,
                      *, requirements=()):
  """Apply whole-MO composition/joint gates, then the per-spin rank window."""
  composition_matches = [mo for mo in channel['members']
                          if _group_sum(mo,model,atoms,policy) >= target_min-1e-12
                          and _group_sum(mo,model,atoms) >= fragment_min-1e-12
                          and all(_requirement_weight(mo,model,r) >= r['minimum']-1e-12
                                  for r in requirements)]
  eligible = {mo.identifier for mo in window_members(channel,window)}
  selected = sorted((mo for mo in composition_matches if mo.identifier in eligible),
                    key=lambda mo:(mo.energy,mo.number))
  return composition_matches,selected


# -----------------------------------------------------------------------------
# Energy-close composition reminders (not symmetry or orbital-overlap proofs)

HARTREE_TO_EV = 27.211386245988


def _mass_overlap(first, second):
  """Similarity of normalized weights: sum(min(p,q)) = 1 - TV distance."""
  import numpy as np
  return max(0.0,min(1.0,float(np.minimum(first,second).sum())))


def _character_vectors(model, mo, type_keys):
  import numpy as np
  lookup = {key:i for i,key in enumerate(type_keys)}
  types, angular = np.zeros(len(type_keys)), np.zeros(len(L_ORDER))
  for (atom,letter),weight in zip(model.groups,mo.weights):
    types[lookup[(model.atoms[atom-1].z,letter)]] += weight
    angular[L_ORDER.index(letter)] += weight
  return types,angular


def _character_text(vector, keys):
  terms = sorted(((float(w),key) for w,key in zip(vector,keys)),key=lambda t:(-t[0],t[1]))
  shown = [term for term in terms if term[0] >= .05][:4] or terms[:1]
  return '; '.join(f'{SYMBOLS[z]}:{letter}={100*w:.1f}%' for w,(z,letter) in shown)


def near_level_report(model, channel_map, window, args):
  """Print same-spin near-level candidates in a preselected, noninteractive window."""
  record = {
    'enabled':not args.no_near,'window_offsets':list(window),
    'window_source':'--near-window; independent of subsequent fragment queries',
    'max_gap_eV':args.near_energy,'minimum_element_l_similarity':args.near_sim,
    'spin_policy':'same','occupation_policy':args.near_occupation,
    'similarity_definition':'sum(min(p_i,q_i)) for normalized whole-MO composition weights',
    'atom_pattern_threshold':args.near_sim,
    'limitations':[
      'Near-level screening itself uses SCPA character, not symmetry, spatial overlap or phase tests.',
      'Element/l similarity can hide localization on different atoms of the same element.',
    ],
    'skipped_constant_energy_channels':[],'pairs':[],'total_pairs':0,
    'display_cap_per_category':args.near_max_pairs,
  }
  if args.no_near:
    print('\n[NEAR-LEVEL REMINDERS] Disabled by --no-near.')
    return record
  print('\n[ENERGY-CLOSE, SIMILAR-CHARACTER ORBITAL REMINDERS]')
  print(f'Audit scope: unfiltered HOMO-{window[0]} through LUMO+{window[1]} per spin; '
        'independent of fragment composition thresholds and subsequent queries.')
  print('Only orbitals within the SAME spin channel are compared.')
  print(f'Require |Delta E| <= {args.near_energy:g} eV AND element+shell similarity >= {args.near_sim:.3f}.')
  print('Similarity = sum(min(p,q)); 1 means identical composition weights, NOT identical orbitals. '
        f'Occupation policy={args.near_occupation}.')
  print('O/V means occupied/virtual by the printed Occup tolerance; actual occupations are also shown.')
  print('Matches are possible near-level/degeneracy candidates only; different atom-patterns are labelled explicitly.')
  type_keys = sorted({(model.atoms[a-1].z,l) for a,l in model.groups})
  factor = HARTREE_TO_EV if args.energy_unit == 'hartree' else 1.0
  for spin,channel in channel_map.items():
    if len(channel['members']) > 1 and len({mo.energy for mo in channel['members']}) == 1:
      print(f'WARNING: {spin} has constant energies for the entire set; skipping an uninformative near-level audit.')
      record['skipped_constant_energy_channels'].append(spin)
      continue
    members = sorted(window_members(channel,window),key=lambda mo:(mo.energy,mo.number))
    vectors = {mo.identifier:_character_vectors(model,mo,type_keys) for mo in members}
    for i,first in enumerate(members):
      first_occ = first.occupation > args.occ_tol
      for second in members[i+1:]:
        difference = abs(second.energy-first.energy)
        gap = difference * factor
        if gap > args.near_energy + (1e-12 if args.near_energy > 0 else 0):
          break
        second_occ = second.occupation > args.occ_tol
        if args.near_occupation == 'same' and first_occ != second_occ:
          continue
        first_types,first_l = vectors[first.identifier]
        second_types,second_l = vectors[second.identifier]
        similarity = _mass_overlap(first_types,second_types)
        if similarity < args.near_sim-1e-12:
          continue
        atom_similarity = _mass_overlap(first.weights,second.weights)
        record['pairs'].append({
          'first':first.identifier,'second':second.identifier,
          'first_spin':spin,'second_spin':spin,
          'first_relative':channel['labels'][first.identifier],
          'second_relative':channel['labels'][second.identifier],
          'first_energy_input_unit':first.energy,'second_energy_input_unit':second.energy,
          'energy_gap_input_unit':difference,'energy_gap_eV':gap,
          'first_occupation':first.occupation,'second_occupation':second.occupation,
          'occupation_classes':('O' if first_occ else 'V')+'/'+('O' if second_occ else 'V'),
          'category':'SAME_SPIN_NEAR_LEVEL_CANDIDATE','element_l_similarity':similarity,
          'angular_l_similarity':_mass_overlap(first_l,second_l),
          'atom_l_similarity':atom_similarity,
          'atomic_pattern':'similar' if atom_similarity >= args.near_sim-1e-12 else 'different',
          'first_element_l_summary':_character_text(first_types,type_keys),
          'second_element_l_summary':_character_text(second_types,type_keys),
          'first_orbital_metrics':_metric_record(first),
          'second_orbital_metrics':_metric_record(second),
          'first_atom_l_summary':_components(first,model,args.component_min,args.top_components),
          'second_atom_l_summary':_components(second,model,args.component_min,args.top_components),
        })
  record['pairs'].sort(key=lambda pair:(pair['energy_gap_eV'],-pair['element_l_similarity'],pair['first'],pair['second']))
  record['total_pairs'] = len(record['pairs'])
  record['pairs_by_spin'] = {
    spin:sum(pair['first_spin']==spin for pair in record['pairs']) for spin in channel_map
  }
  print(f'\nSAME_SPIN_NEAR_LEVEL_CANDIDATE: {len(record["pairs"])} qualifying pairs; '
        f'displaying at most {args.near_max_pairs} across available spin channels.')
  print('Counts by spin: '+', '.join(f'{spin}={count}' for spin,count in record['pairs_by_spin'].items()))
  for pair in record['pairs'][:args.near_max_pairs]:
    print(f'  {pair["first"]} ({pair["first_relative"]}) ~ '
          f'{pair["second"]} ({pair["second_relative"]}): '
          f'DeltaE={pair["energy_gap_eV"]:.6g} eV; Occ={pair["first_occupation"]:g}/{pair["second_occupation"]:g} '
          f'[{pair["occupation_classes"]}]; '
          f'S(element,l)={pair["element_l_similarity"]:.3f}; '
          f'S(l)={pair["angular_l_similarity"]:.3f}; '
          f'S(atom,l)={pair["atom_l_similarity"]:.3f}; atom-pattern={pair["atomic_pattern"]}')
    unit = 'Eh' if args.energy_unit == 'hartree' else 'eV'
    print(f'    First E={pair["first_energy_input_unit"]:.9f} {unit}; element/l: '+pair['first_element_l_summary'])
    print('      Atom/l: '+pair['first_atom_l_summary'])
    print('      '+_metric_text(next(mo for mo in model.orbitals if mo.identifier==pair['first']),args.show_pm))
    print(f'    Second E={pair["second_energy_input_unit"]:.9f} {unit}; element/l: '+pair['second_element_l_summary'])
    print('      Atom/l: '+pair['second_atom_l_summary'])
    print('      '+_metric_text(next(mo for mo in model.orbitals if mo.identifier==pair['second']),args.show_pm))
  if not record['pairs']:
    print('  (none)')
  print('Reminder only: use orbital visualization, symmetry and/or spatial AO overlap for confirmation.')
  return record


# -----------------------------------------------------------------------------
# Reports and analysis-only arguments

def add_analysis_arguments(parser):
  """Register analysis options; return their names for strict mode separation."""
  specifications = [
    (('--frontier',), dict(type=int, default=10, help='Frontier summary offsets, default 10 on each side.')),
    (('--window',), dict(type=int, nargs='+', default=[20,20], metavar='N', help='Candidate rank window: N or occupied_N virtual_N; default 20 20.')),
    (('--fragment',), dict(help='1-based atom rows, e.g. "1,3-8,12". Suppresses automatic prompts.')),
    (('--target-profile',), dict(choices=('dominant','expanded'), default='dominant', help='dominant: s/p/d by block; expanded: s/sp/sd.')),
    (('--target-l',), dict(help='Per-atom overrides, e.g. "1:ds;2:sp". ds means d+s, not a periodic-table block.')),
    (('--target-min',), dict(type=float, default=.30, help='Minimum whole-MO fraction on selected atom+channel targets.')),
    (('--fragment-min',), dict(type=float, default=.30, help='Minimum whole-MO fraction on selected atoms, all channels.')),
    (('--show-metal',), dict(action='store_true', help='Explicitly enable the optional metal-d report; disabled by default.')),
    (('--metal-min',), dict(type=float, default=.30, help='Total d fraction cutoff for the optional --show-metal/--metals report.')),
    (('--metals',), dict(help='Explicit atom rows for a metal-d query; also enables that optional report.')),
    (('--component-min',), dict(type=float, default=.05, help='Minimum displayed atom/l term; does not alter screening.')),
    (('--top-components',), dict(type=int, default=5)),
    (('--occ-tol',), dict(type=float, default=1e-6, help='Occup > tolerance is classified as occupied for ranking.')),
    (('--energy-unit',), dict(choices=('hartree','eV'), default='hartree', help='Unit of Ene in this input; default is an explicit Hartree assumption, no guessed conversion.')),
    (('--allow-sparse',), dict(action='store_true', help='Treat omitted coefficients as zero only for verified sparse input.')),
    (('--fragment-mode',), dict(choices=('sum','joint'), default='sum', help='sum: combined target enrichment; joint: also require named groups or every selected atom.')),
    (('--per-atom-min',), dict(type=float, default=.05, help='Joint-mode per-atom target minimum when no --require-target is supplied.')),
    (('--require-target',), dict(help='Joint AND groups, e.g. "52-53:d:0.10;54:p:0.10"; whole-MO fractions.')),
    (('--near-energy',), dict(type=float, default=.05, help='Maximum absolute energy gap in eV for near-level reminders (default 0.05 eV).')),
    (('--near-sim',), dict(type=float, default=.90, help='Minimum element+angular-shell distribution similarity, in [0,1].')),
    (('--near-window',), dict(type=int, nargs='+', default=[20,20], metavar='N', help='Near-level audit offsets N or occupied_N virtual_N; default 20 20, independent of fragment --window.')),
    (('--near-spin',), dict(choices=('same',), default='same', help='Compatibility option: only same-spin comparisons are supported.')),
    (('--near-occupation',), dict(choices=('same','any'), default='any', help='Compare same occupation class only, or label all occupied/virtual combinations.')),
    (('--near-max-pairs',), dict(type=int, default=30, help='Terminal cap for all same-spin pairs; JSON retains all qualifying matches.')),
    (('--no-near',), dict(action='store_true', help='Disable near-level composition reminders.')),
    (('--show-pm',), dict(action='store_true', help='Also print PM2_L/N_eff_L columns; they are always retained in computed JSON metrics.')),
    (('--orbital-metrics',), dict(choices=('auto','pyscf','off'), default='auto', help='True spatial overlap/Lowdin-PM backend. auto uses --ao-overlap or PySCF; missing PySCF gives explicit NA, never a heuristic substitute.')),
    (('--ao-overlap',), dict(type=Path, help='Optional .npy/.txt/.dat overlap matrix in ORIGINAL Molden AO order/normalization; use with auto.')),
    (('--metric-threads',), dict(type=int, default=1, help='PySCF native integral thread limit; default 1, restored after loading.')),
    (('--metric-norm-tol',), dict(type=float, default=.02, help='Compatibility-only option; ignored. Nonunit MO self overlaps are allowed.')),
    (('--no-prompt',), dict(action='store_true', help='Print reports and exit without reading stdin.')),
    (('--json',), dict(type=Path, help='Write the complete analysis/selection record to a NEW JSON file.')),
  ]
  names = set()
  for flags, settings in specifications:
    parser.add_argument(*flags, **settings)
    names.update(flags)
  return names


def _components(mo, model, cutoff, limit, *, atoms=None, policy=None):
  allowed = set(atoms) if atoms is not None else None
  terms = []
  for (atom, letter), weight in zip(model.groups, mo.weights):
    if allowed is not None and atom not in allowed:
      continue
    if policy is not None and letter not in policy[atom]:
      continue
    terms.append((float(weight), atom, letter))
  terms.sort(key=lambda t:(-t[0], t[1], L_ORDER.index(t[2])))
  shown = [term for term in terms if term[0] >= cutoff][:limit]
  if not shown and terms:
    shown = terms[:1]
  text = '; '.join(f'{model.atoms[a-1].symbol}{a}:{l}={100*w:.1f}%' for w,a,l in shown)
  if atoms is None:
    text += ('; ' if text else '') + f'Other={100*max(0.,1-sum(t[0] for t in shown)):.1f}%'
  return text or '(no weight on these atom/channel targets)'


def _print_target_breakdown(detail, model):
  """Print the complete deciding contributions, independently of display cutoffs."""
  terms = []
  for row in detail['target_breakdown']:
    label = f'{row["element"]}{row["atom"]}:{row["angular_shell"]}'
    term = f'{label}={100*row["weight"]:.3f}%'
    if not row['basis_present']:
      term += '(no such basis shell)'
    terms.append(term)
  print('       Target breakdown (NO display cutoff): ' + '; '.join(terms))
  atom_text = '; '.join(f'{r["element"]}{r["atom"]}={100*r["weight"]:.3f}%'
                        for r in detail['target_by_atom'])
  print('       Target by atom: ' + atom_text)
  print(f'       Sum={100*detail["target_weight"]:.3f}%; '
        f'all selected atoms={100*detail["fragment_weight"]:.3f}%; '
        f'non-target on selected atoms={100*detail["fragment_non_target_weight"]:.3f}%')
  for check in detail['joint_checks']:
    label = '+'.join(f'{model.atoms[n-1].symbol}{n}[{check["channels"]}]' for n in check['atoms'])
    print(f'       Joint check: {label}={100*check["weight"]:.3f}% '
          f'>= {100*check["minimum"]:.3f}%: {"PASS" if check["passed"] else "FAIL"}')


def _table(model, channel, selected, args, *, atoms=None, policy=None, metal=False,
            requirements=()):
  unit = 'Eh' if args.energy_unit == 'hartree' else 'eV'
  extra = ' Metal-d%' if metal else ' Frag%  Target%' if atoms is not None else ''
  print(f'MO_ID  Spin   Occup    Energy({unit})       Ref         {_metric_header(args.show_pm)}{extra}  Main whole-MO atom/l terms')
  details = []
  for mo in selected:
    label = channel['labels'][mo.identifier]
    if metal:
      suffix = f' {_group_sum(mo,model,atoms,policy)*100:8.2f}'
    elif atoms is not None:
      suffix = f' {_group_sum(mo,model,atoms)*100:6.2f} {_group_sum(mo,model,atoms,policy)*100:8.2f}'
    else:
      suffix = ''
    print(f'{mo.identifier:<6} {mo.spin:<6} {mo.occupation:7.3f} {mo.energy:15.7f} {label:<11} '
          f'{_metric_columns(mo,args.show_pm)}{suffix}  {_components(mo,model,args.component_min,args.top_components)}')
    if atoms is not None:
      detail = _target_details(mo,model,atoms,policy,requirements)
      details.append(detail)
      _print_target_breakdown(detail,model)
  if not selected:
    print('(no matching orbitals; review thresholds, joint requirements, channels or rank window)')
  return details


def _selection_report(model, channel_map, atoms, policy, window, args, *, metal=False,
                      requirements=()):
  records = []
  title = 'D-BLOCK / SELECTED-CENTER d-CHARACTER CANDIDATES' if metal else 'FRAGMENT TARGET-CHANNEL CANDIDATES'
  print(f'\n[{title}]')
  print('Selected atom/channel pairs: ' + ', '.join(
    f'{model.atoms[n-1].symbol}{n}[{",".join(policy[n]) or "UNSET"}]' for n in atoms))
  threshold = args.metal_min if metal else args.target_min
  minimum = 0 if metal else args.fragment_min
  mode = 'sum' if metal else args.fragment_mode
  print(f'Selection mode: {mode}. All target percentages are fractions of the WHOLE MO.')
  print(f'Composition first: summed Target >= {threshold*100:.1f}%; '
        f'all-channel fragment >= {minimum*100:.1f}%.')
  if minimum <= threshold:
    print('The fragment-total gate adds no extra restriction here because Target is a subset of Frag.')
  if mode == 'sum':
    print('SUM means enrichment on the union; it does NOT require every selected center to participate.')
  else:
    print('JOINT adds ALL of the following whole-MO minimum-contribution conditions:')
    for requirement in requirements:
      label = '+'.join(f'{model.atoms[n-1].symbol}{n}[{requirement["channels"]}]' for n in requirement['atoms'])
      print(f'  {label} >= {100*requirement["minimum"]:.2f}%')
    print('These gates identify shared character, not proof of bonding, antibonding or a multicenter bond.')
  print(f'Then rank window: HOMO-{window[0]} through LUMO+{window[1]}, inclusive, independently per spin.')
  print('This is an orbital-count window, NOT an energy-in-eV cutoff or a certified valence classifier.')
  for spin, channel in channel_map.items():
    before, selected = select_candidates(model,channel,atoms,policy,window,threshold,minimum,
                                          requirements=requirements)
    print(f'\n{spin}: composition matches before window={len(before)}; after window={len(selected)}')
    boundary_warnings = _window_warnings(channel,window)
    for warning in boundary_warnings:
      print('WARNING:',warning)
    details = _table(model,channel,selected,args,atoms=atoms,policy=policy,
                      metal=metal,requirements=requirements)
    records.append({'spin':spin,'composition_match_count':len(before),
                    'selected_ids':[mo.identifier for mo in selected],
                    'selected_details':details,'window_warnings':boundary_warnings})
  return {'type':'metal_d' if metal else 'fragment','selection_mode':mode,'atoms':list(atoms),
          'targets':{str(n):policy[n] for n in atoms},'window_offsets':list(window),
          'target_min':threshold,'fragment_min':minimum,
          'joint_requirements':list(requirements),'channels':records}


def _header(model, channel_map, policy, args, warnings):
  print('\n=== CUB1C ORBITAL INSPECTION: NUMERICAL REPORT ONLY ===')
  print('Input:',model.path)
  print(f'Atoms={len(model.atoms)}; AO basis functions={model.basis_count}; MOs={len(model.orbitals)}')
  print('Composition method: SCPA = C_mu^2 / sum_nu(C_nu^2), aggregated by atom and angular shell.')
  print('AO-metric diagnostics:',model.metric_info['status'],model.metric_info.get('backend',''))
  print('ABpeer = opposite-spin MO with the SAME original within-spin index (not best-match/UCO).')
  print('AB/self = |<current|peer>|/<current|current>: directed, scale-dependent, and MAY EXCEED 1.')
  print('AB|S| = |<current|peer>|/sqrt(<current|current><peer|peer>): symmetric spatial match in [0,1].')
  print('Spin functions are excluded; global orbital sign is ignored in both displayed ratios.')
  print('AB|S|=1 means equal spatial shape up to sign and scale; 0 means orthogonality, not necessarily different atoms.')
  print('Nonunit MO self overlaps are allowed. Input MO columns are not rotated, rescaled or optimized.')
  print('Removing the unit-norm gate does not verify AO conventions; S must still be physically valid and consistent.')
  print('Low same-index overlap alone does not prove spin polarization; orbital ordering/rotations can lower it.')
  if args.show_pm:
    print('PM2_L = sum_A(q_A^Lowdin)^2 using self-norm-divided atomic fractions; no localization rotation is performed.')
    print('N_eff_L = 1/PM2_L: equivalent participating atoms, NOT a distance or integer atom count.')
  else:
    print('PM diagnostics are retained in JSON; add --show-pm to print their columns.')
  print('Raw signed overlap, both self overlaps, signed directed ratio and population metrics are stored in JSON.')
  print('NA means unavailable diagnostics, never zero matching or zero localization.')
  print('Basis dimensions:', ', '.join(f'{l}={model.dimensions[l]}' for l in L_ORDER))
  print(f'Energy unit: {args.energy_unit} (user/default convention for Ene; no energy conversion).')
  print('MO_ID is the original 1-based WITHIN-SPIN file order + A/B; FileMO is recorded in JSON.')
  print('Ref is an energy-rank label within the occupied or virtual set, NOT a renumbering of MO_ID.')
  print(f'Occupied rule: Occup > {args.occ_tol:g}; virtual rule: Occup <= {args.occ_tol:g}.')
  print('All percentages are WHOLE-MO SCPA fractions, never renormalized inside a fragment.')
  print('Frag% = all s/p/d/f/g weights on selected atoms; Target% = only their selected channels.')
  print('Target breakdowns include EVERY requested atom/channel, including zeros; main-term display cutoffs do not apply there.')
  print('Printed percentages are rounded; JSON stores full-precision sums and deciding contributions.')
  print('s/p/d are basis-shell types, not inferred principal quantum numbers such as 2p or 3d.')
  print(f'Display terms: >= {args.component_min*100:g}%, up to {args.top_components}; if none pass, show the largest. Other is undisplayed whole-MO weight.')
  print('Target presets: dominant=s-block:s, p-block:p, d-block:d; expanded adds p-block:s and d-block:s.')
  print('ds is a user d+s channel preset, not a separate periodic-table block. f-region targets require an override.')
  if args.show_metal or args.metals is not None:
    print('Automatic d-center ranges include group-12 elements; this is not a strict transition-metal definition.')
  print('Atom IDs are 1-based [Atoms] ROWS. Source atom numbers are retained separately in JSON.')
  labels = [f'{a.symbol}{a.number}[{policy[a.number] or "UNSET"}]' for a in model.atoms]
  for i in range(0,len(labels),8):
    print('Atoms: ' + ' '.join(labels[i:i+8]))
  for spin, channel in channel_map.items():
    occupied, virtual = channel['occupied'], channel['virtual']
    print(f'{spin} endpoints: HOMO={occupied[-1].identifier if occupied else "unavailable"}; '
          f'LUMO={virtual[0].identifier if virtual else "unavailable"}')
  for warning in warnings:
    print('WARNING:',warning)


def run_analysis(filename, args):
  """Print frontier, immediate same-spin reminders, then optional fragment reports."""
  window = window_pair(args.window)
  near_window = window_pair(args.near_window)
  if args.frontier < 0 or args.top_components < 1 or not math.isfinite(args.occ_tol) or not 0 < args.occ_tol < .1:
    raise ValueError('Require frontier >= 0, top-components >= 1 and 0 < occ-tol < 0.1.')
  for name in ('target_min','fragment_min','metal_min','component_min'):
    value = getattr(args,name)
    if not math.isfinite(value) or not 0 <= value <= 1:
      raise ValueError(f'{name.replace("_","-")} must be a fraction in [0, 1].')
  if not math.isfinite(args.per_atom_min) or not 0 < args.per_atom_min <= 1:
    raise ValueError('per-atom-min must be a fraction in (0,1].')
  if args.require_target and args.fragment_mode != 'joint':
    raise ValueError('--require-target requires --fragment-mode joint.')
  if args.near_spin != 'same':
    raise ValueError('Only same-spin near-level comparisons are supported.')
  if not math.isfinite(args.near_energy) or args.near_energy < 0:
    raise ValueError('near-energy must be a nonnegative finite value in eV.')
  if not math.isfinite(args.near_sim) or not 0 < args.near_sim <= 1 or args.near_max_pairs < 1:
    raise ValueError('Require 0 < near-sim <= 1 and near-max-pairs >= 1.')
  output = args.json.expanduser().absolute() if args.json else None
  if output is not None and (output.exists() or output.is_symlink() or not output.parent.is_dir()):
    raise FileExistsError('JSON target must be new and its parent directory must exist.')
  if args.metric_threads<1:
    raise ValueError('Require metric-threads >= 1.')
  model = read_molden(filename,args.allow_sparse)
  channel_map, warnings = channels(model,args.occ_tol)
  metric_info = compute_orbital_metrics(model,channel_map,args)
  warnings = model.warnings + warnings + metric_info['warnings']
  policy = channel_policy(model,args.target_profile,args.target_l)
  _header(model,channel_map,policy,args,warnings)
  print(f'\n[FRONTIER SUMMARY: HOMO-{args.frontier} THROUGH LUMO+{args.frontier}]')
  print('Unfiltered by fragment/target composition; endpoints are inclusive.')
  for spin, channel in channel_map.items():
    print('\n'+spin)
    for warning in _window_warnings(channel,(args.frontier,args.frontier)):
      print('WARNING:',warning)
    _table(model,channel,window_members(channel,(args.frontier,args.frontier)),args)
  # Complete the fixed near-level report before any fragment input is read.
  near = near_level_report(model,channel_map,near_window,args)
  queries = []
  metal_report_enabled = args.show_metal or args.metals is not None
  if metal_report_enabled:
    metal_atoms = parse_atoms(args.metals,len(model.atoms)) if args.metals is not None else [
      a.number for a in model.atoms if a.z in D_BLOCK
    ]
    if metal_atoms:
      metal_policy = {a.number:'d' for a in model.atoms}
      queries.append(_selection_report(model,channel_map,metal_atoms,metal_policy,window,args,metal=True))
    else:
      print('\n[METAL d SEARCH] No automatic d-block center; use --metals for explicit atom rows.')

  def fragment_report(selection, selected_window):
    atoms = parse_atoms(selection,len(model.atoms))
    unset = [n for n in atoms if not policy[n]]
    if unset:
      raise ValueError('Define --target-l for selected f-region atoms: ' + ','.join(map(str,unset)))
    requirements = _joint_requirements(model,atoms,policy,args)
    queries.append(_selection_report(model,channel_map,atoms,policy,selected_window,args,requirements=requirements))

  if args.fragment:
    fragment_report(args.fragment,window)
  elif not args.no_prompt and sys.stdin.isatty():
    while True:
      try:
        selection = input('\nFragment atom rows (e.g. 1,3-8,12; Enter/q to finish): ').strip()
        if not selection or selection.lower() in ('q','quit','exit'):
          break
        text = input(f'Candidate offsets [occupied virtual] (one N or two integers; Enter={window[0]} {window[1]}): ').strip()
        selected_window = window_pair([int(v) for v in text.replace(',',' ').split()]) if text else window
        fragment_report(selection,selected_window)
      except (EOFError,KeyboardInterrupt):
        print('\nFragment input ended.')
        break
      except ValueError as exc:
        print('Invalid selection:',exc)
  else:
    print('\nNo interactive input requested. Use --fragment to select atom rows.')
  print('\nAnalysis complete. No cube, plot, OBJ or rendering was generated.')
  print('Candidates are target-shell/enriched near-frontier orbitals, NOT certified valence/core assignments.')
  record = {
    'schema_version':4,'mode':'orbital_inspection','method':'SCPA',
    'input_molden':str(model.path),'input_sha256':model.sha256,
    'formula':'w_mu = C_mu^2 / sum_nu(C_nu^2)',
    'energy_unit':args.energy_unit,'coordinate_unit':model.coordinate_unit,
    'index_convention':'MO_ID preserves original within-spin file order; Ref uses separate occupied/virtual energy ranks',
    'occupation_tolerance':args.occ_tol,'basis_count':model.basis_count,
    'basis_dimensions':model.dimensions,'frontier_offsets':args.frontier,
    'default_candidate_window':list(window),'near_window_offsets':list(near_window),
    'metal_report_enabled':metal_report_enabled,'target_profile':args.target_profile,
    'target_policy':{str(n):p for n,p in policy.items()},
    'atoms':[vars(a) for a in model.atoms],
    'groups':[{'atom':a,'angular_shell':l} for a,l in model.groups],
    'orbitals':[{
      'MO_ID':mo.identifier,'FileMO':mo.file_number,'channel_index':mo.number,
      'spin':mo.spin,'energy':mo.energy,'occupation':mo.occupation,'symmetry':mo.symmetry,
      'relative':channel_map[mo.spin]['labels'][mo.identifier],
      'group_weights':mo.weights.tolist(),'orbital_metrics':_metric_record(mo),
    } for mo in model.orbitals],
    'screening_defaults':{'target_min':args.target_min,'fragment_min':args.fragment_min,
                          'metal_min':args.metal_min},
    'fragment_mode':args.fragment_mode,'per_atom_min':args.per_atom_min,
    'require_target_specification':args.require_target,
    'queries':queries,'near_level_reminders':near,
    'orbital_metric_method':metric_info,'warnings':warnings,
  }
  if output is not None:
    with output.open('x',encoding='utf-8') as stream:
      json.dump(record,stream,indent=2,ensure_ascii=False,allow_nan=False)
      stream.write('\n')
    print('Saved analysis JSON:',output)
  return record


def main(argv=None):
  parser = argparse.ArgumentParser(description=__doc__,allow_abbrev=False)
  parser.add_argument('molden')
  add_analysis_arguments(parser)
  args = parser.parse_args(argv)
  try:
    run_analysis(args.molden,args)
  except Exception as exc:
    parser.exit(2,f'{parser.prog}: {type(exc).__name__}: {exc}\n')


if __name__ == '__main__':
  main()
