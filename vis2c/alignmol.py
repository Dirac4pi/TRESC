#!/home/lky/miniconda3/envs/vis2c/bin/python
'''
Align an XYZ probe using manual paired atoms or an automatic MCS.
author: Dirac4pi
env: vis2c
'''

import argparse
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdFMCS


#------------------------------------------------------------------------------
def _parse_indices(text: str, label: str) -> list[int]:
  """Expand comma-separated, inclusive, 1-based atom indices and ranges.

  Parameters:
    text: Examples: "1-9,16-23,25" or "3,5,8-10".
    label: Name of the selection used in error messages.

  Returns:
    Atom numbers in the exact order supplied; ranges expand in ascending order.
  """
  indices = []
  for item in text.split(","):
    item = item.strip()
    if not item:
      raise ValueError(f"{label}: empty item in index selection.")
    if "-" in item:
      parts = item.split("-")
      if len(parts) != 2 or not all(part.strip().isdigit() for part in parts):
        raise ValueError(f"{label}: invalid range {item!r}.")
      first, last = (int(part.strip()) for part in parts)
      if first < 1 or last < first:
        raise ValueError(f"{label}: range must be positive and ascending: {item!r}.")
      indices.extend(range(first, last + 1))
    else:
      if not item.isdigit() or int(item) < 1:
        raise ValueError(f"{label}: invalid atom number {item!r}.")
      indices.append(int(item))
  if len(indices) != len(set(indices)):
    raise ValueError(f"{label}: duplicate atom numbers are not allowed.")
  return indices

#------------------------------------------------------------------------------
def _mcs_graph(mol):
  """Infer connectivity in Angstrom without assigning bond orders."""
  import qcelemental as qcel
  positions = np.asarray(mol.GetConformer().GetPositions())
  if not np.all(np.isfinite(positions)):
    raise ValueError("XYZ coordinates must be finite.")
  radii = np.array([qcel.covalentradii.get(a.GetAtomicNum(), units="angstrom")
    for a in mol.GetAtoms()], dtype=float)
  graph = Chem.RWMol(mol)
  for bond in list(graph.GetBonds()):
    graph.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
  for i in range(len(positions) - 1):
    distances = np.linalg.norm(positions[i+1:] - positions[i], axis=1)
    if np.any(distances <= 1e-12):
      raise ValueError("Coincident atoms cannot define MCS connectivity.")
    for offset in np.flatnonzero(distances <= 1.15 * (radii[i] + radii[i+1:])):
      graph.AddBond(i, i + 1 + int(offset), Chem.BondType.SINGLE)
  result = graph.GetMol()
  result.UpdatePropertyCache(strict=False)
  Chem.GetSymmSSSR(result)
  return result

#------------------------------------------------------------------------------
def _valid_fit(mol, indices):
  """Test whether selected positions support a non-collinear rigid fit."""
  points = np.asarray(mol.GetConformer().GetPositions())[list(indices)]
  if len(points) < 3:
    return False
  values = np.linalg.svd(points - points.mean(axis=0), compute_uv=False)
  return values[0] >= 1e-8 and values[1] / values[0] >= 1e-3

#------------------------------------------------------------------------------
def _mcs_indices(ref, probe, timeout=60):
  """Find an element/connectivity MCS, including explicit hydrogen atoms."""
  ref_graph, probe_graph = _mcs_graph(ref), _mcs_graph(probe)
  result = rdFMCS.FindMCS([ref_graph, probe_graph],
    maximizeBonds=False, timeout=timeout,
    atomCompare=rdFMCS.AtomCompare.CompareElements,
    bondCompare=rdFMCS.BondCompare.CompareAny,
    matchValences=False, ringMatchesRingOnly=False,
    completeRingsOnly=False, matchChiralTag=False)
  if result.canceled:
    raise ValueError("MCS search timed out; the maximum is not certified. "
      "Use manual selections or increase the search timeout.")
  if result.numAtoms < 3:
    raise ValueError("MCS has fewer than three atoms; use manual fitting atoms.")
  query = Chem.MolFromSmarts(result.smartsString)
  if query is None:
    raise ValueError("Could not parse the MCS SMARTS query.")
  # Enumerate symmetry-related embeddings, with explicit resource limits.
  limit, pair_limit = 1000, 100000
  refs = ref_graph.GetSubstructMatches(query, uniquify=False, maxMatches=limit+1)
  probes = probe_graph.GetSubstructMatches(query, uniquify=False, maxMatches=limit+1)
  truncated = len(refs) > limit or len(probes) > limit
  refs, probes = sorted(set(refs[:limit])), sorted(set(probes[:limit]))
  refs = [m for m in refs if _valid_fit(ref, m)]
  probes = [m for m in probes if _valid_fit(probe, m)]
  best, tested = None, 0
  for rm in refs:
    for pm in probes:
      if tested >= pair_limit:
        truncated = True
        break
      tested += 1
      rmsd, _ = rdMolAlign.GetAlignmentTransform(probe, ref,
        atomMap=list(zip(pm, rm)), reflect=False)
      candidate = (float(rmsd), rm, pm)
      if best is None or candidate < best:
        best = candidate
    if tested >= pair_limit:
      break
  if best is None:
    raise ValueError("No non-collinear MCS embedding supports alignment.")
  if truncated:
    print("WARNING: symmetry mappings were capped; RMSD is best among tested mappings only.")
  print(f"MCS: {result.numAtoms} atoms, {result.numBonds} bonds; mappings tested: {tested}")
  print("MCS SMARTS:", result.smartsString)
  print("MCS connectivity: d <= 1.15 * covalent-radius sum; explicit H included.")
  _, rm, pm = best
  print("REF fitting atoms:", ",".join(str(i+1) for i in rm))
  print("PROBE fitting atoms:", ",".join(str(i+1) for i in pm))
  return [i+1 for i in rm], [i+1 for i in pm]

#------------------------------------------------------------------------------
def align_xyz_fragment(ref_xyz: str, probe_xyz: str) -> np.ndarray:
  """Interactively select fragments and align probe to ref.

  Parameters:
    ref_xyz: Reference XYZ file path; its coordinates remain unchanged.
    probe_xyz: Probe XYZ file path; its coordinates remain unchanged.

  Returns:
    Probe-to-reference 4x4 rigid transform with translation in Angstrom.
    For an (N, 3) array of probe positions in Angstrom:
    new_positions = old_positions @ T[:3, :3].T + T[:3, 3]

  Notes:
    Both selections use 1-based XYZ atom numbers, e.g. "1-9,16-23,25".
    Manual entries are paired by position. Enter mcs for automatic mapping.
    Importing and calling this function prompts for input; it does not write
    a file. Running the module as a script also writes the aligned probe XYZ.
  """
  print(f"Reference XYZ (first argument): {ref_xyz}")
  print(f"Probe XYZ (second argument): {probe_xyz}")
  ref = Chem.MolFromXYZFile(str(ref_xyz))
  probe = Chem.MolFromXYZFile(str(probe_xyz))
  if ref is None or probe is None:
    raise ValueError("Failed to load an XYZ file.")

  ref_text = input(
    "Enter REF fitting atoms (1-based ranges, or mcs for automatic mapping): "
  ).strip()
  if ref_text.lower() == "mcs":
    ref_indices, probe_indices = _mcs_indices(ref, probe)
  else:
    probe_text = input(
      "Enter PROBE fitting atoms (same order as REF, or mcs): "
    ).strip()
    if probe_text.lower() == "mcs":
      print("MCS mode replaces both manual selections.")
      ref_indices, probe_indices = _mcs_indices(ref, probe)
    else:
      ref_indices = _parse_indices(ref_text, "ref_indices")
      probe_indices = _parse_indices(probe_text, "probe_indices")

  if len(ref_indices) != len(probe_indices) or len(ref_indices) < 3:
    raise ValueError("Provide at least three paired atom indices in each selection.")
  atom_map = []
  for probe_number, ref_number in zip(probe_indices, ref_indices):
    if not 1 <= probe_number <= probe.GetNumAtoms():
      raise ValueError(f"Probe atom index out of range: {probe_number}")
    if not 1 <= ref_number <= ref.GetNumAtoms():
      raise ValueError(f"Reference atom index out of range: {ref_number}")
    probe_idx = probe_number - 1
    ref_idx = ref_number - 1
    probe_element = probe.GetAtomWithIdx(probe_idx).GetSymbol()
    ref_element = ref.GetAtomWithIdx(ref_idx).GetSymbol()
    if probe_element != ref_element:
      raise ValueError(
        f"Element mismatch: probe {probe_number} ({probe_element}) "
        f"vs ref {ref_number} ({ref_element})"
      )
    atom_map.append((probe_idx, ref_idx))

  # Reject collinear or nearly collinear fitting fragments.
  for label, mol, indices in (
    ("probe", probe, [p for p, _ in atom_map]),
    ("ref", ref, [r for _, r in atom_map]),
  ):
    positions = mol.GetConformer().GetPositions()[indices]
    centered = positions - positions.mean(axis=0)
    singular_values = np.linalg.svd(centered, compute_uv=False)
    if singular_values[0] < 1e-8:
      raise ValueError(f"{label} fitting fragment has no spatial extent.")
    if singular_values[1] / singular_values[0] < 1e-3:
      raise ValueError(
        f"{label} fitting fragment is collinear or nearly collinear; "
        "add an off-axis atom."
      )

  rmsd, transform = rdMolAlign.GetAlignmentTransform(
    probe,
    ref,
    atomMap=atom_map,
    reflect=False
  )
  transform = np.asarray(transform, dtype=float)
  print(f"Fragment RMSD: {rmsd:.6f} Å")
  return transform


#==============================================================================
def main() -> None:
  """Write the aligned probe XYZ and show both geometries in one Mayavi figure.

  Returns:
    None. The reference is opaque; the aligned probe has opacity 0.5.
    Visualization is CLI-only, so imported alignment stays non-graphical.
  """
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("ref_xyz", help="Reference XYZ file path.")
  parser.add_argument("probe_xyz", help="Probe XYZ file path.")
  args = parser.parse_args()

  output = Path.cwd() / f"{Path(args.probe_xyz).stem}_aligned.xyz"
  if output.resolve() in (Path(args.ref_xyz).resolve(), Path(args.probe_xyz).resolve()):
    raise ValueError("Output filename would overwrite an input XYZ file.")

  transform = align_xyz_fragment(args.ref_xyz, args.probe_xyz)
  probe = Chem.MolFromXYZFile(args.probe_xyz)
  ref = Chem.MolFromXYZFile(args.ref_xyz)
  if probe is None or ref is None:
    raise ValueError("Failed to reload an XYZ file.")
  conformer = probe.GetConformer()
  old_positions = conformer.GetPositions()
  new_positions = old_positions @ transform[:3, :3].T + transform[:3, 3]
  for atom_index, position in enumerate(new_positions):
    conformer.SetAtomPosition(atom_index, position)
  Chem.MolToXYZFile(probe, str(output))
  print(f"Aligned probe XYZ: {output}")

  # Import graphics only when the CLI is used.
  from mayavi import mlab
  from plot3d import atoms_plot

  BOHR_TO_ANGSTROM = 0.529177210903
  ref_positions = ref.GetConformer().GetPositions() / BOHR_TO_ANGSTROM
  probe_plot_positions = new_positions / BOHR_TO_ANGSTROM
  ref_numbers = [atom.GetAtomicNum() for atom in ref.GetAtoms()]
  probe_numbers = [atom.GetAtomicNum() for atom in probe.GetAtoms()]
  figure = mlab.figure(
    figure=f"Alignment: {Path(args.ref_xyz).name} / {Path(args.probe_xyz).name}",
    size=(900, 700),
    bgcolor=(1, 1, 1),
    fgcolor=(0, 0, 0)
  )
  mlab.clf(figure=figure)
  mlab.figure(figure)
  atoms_plot(ref_numbers, ref_positions[:, 0], ref_positions[:, 1],
             ref_positions[:, 2], opacity=1.0)
  mlab.figure(figure)
  atoms_plot(probe_numbers, probe_plot_positions[:, 0], probe_plot_positions[:, 1],
             probe_plot_positions[:, 2], opacity=0.5)
  mlab.view(figure=figure)
  mlab.show()
  return None


if __name__ == "__main__":
  main()
