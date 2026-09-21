#!/usr/bin/env python3
"""
Generate the ground-state MPS of the J1-J2 Heisenberg chain using DMRG and
write it in the plain-text format read by "./run --mps <file>".

    H = J1 * sum_i S_i . S_{i+1}  +  J2 * sum_i S_i . S_{i+2}

with S_i . S_j = (X_i X_j + Y_i Y_j + Z_i Z_j) / 4, open boundary conditions.

J2/J1 = 0.5 is the Majumdar-Ghosh point, where the ground state is an exactly
known product of dimers and every correlation vanishes beyond nearest
neighbours. Away from that point there is no closed form, which is the whole
reason for using this as a benchmark. J2/J1 = 0.35 sits comfortably inside the
gapped dimerised phase (J2/J1 > 0.2411), so DMRG converges quickly.

Usage:
    python3 src/dmrg.py -n 50 --j2 0.35 -o data/j1j2_50.mps
"""

import argparse
import sys

import numpy as np

try:
    from tenpy.models.spins_nnn import SpinChainNNN2
    from tenpy.networks.mps import MPS
    from tenpy.algorithms import dmrg
except ImportError:
    sys.exit("TeNPy not found - install it with: pip install physics-tenpy")


def run_dmrg(num_sites, j1, j2, chi, sweeps, conserve):
    """Find the ground state of the J1-J2 chain, returning (energy, psi)."""

    model = SpinChainNNN2({
        "L": num_sites,
        "S": 0.5,
        "bc_MPS": "finite",
        "conserve": conserve,
        "Jx": j1, "Jy": j1, "Jz": j1,
        "Jxp": j2, "Jyp": j2, "Jzp": j2,
        "sort_charge": True,
    })

    # Start from a Neel state, which has a good overlap with the true ground
    # state in the dimerised phase and carries total Sz = 0
    product_state = ["up", "down"] * (num_sites // 2) + ["up"] * (num_sites % 2)
    psi = MPS.from_product_state(model.lat.mps_sites(), product_state, bc=model.lat.bc_MPS)

    params = {
        "mixer": True,
        "max_E_err": 1.0e-12,
        "max_S_err": 1.0e-10,
        "max_sweeps": sweeps,
        "trunc_params": {
            "chi_max": chi,
            "svd_min": 1.0e-12,
        },
        "combine": True,
    }

    eng = dmrg.TwoSiteDMRGEngine(psi, model, params)
    energy, psi = eng.run()
    psi.canonical_form()
    return energy, psi


def extract_tensors(psi):
    """
    Pull the MPS out as dense (chiL, d, chiR) arrays in exactly left-canonical
    form, which is what src/mps.h assumes.

    The 'B' (right-canonical) form is what TeNPy stores natively, so taking it
    and gauging it ourselves avoids TeNPy's 'A' form, which is built by dividing
    by singular values and loses several digits when any of them are tiny.
    """

    # TeNPy orders the local basis by increasing Sz, so its index 0 is 'down'
    # (Sz = -1/2), i.e. the qubit |1>. Flip it so that index 0 is |0> and the
    # usual Pauli matrices apply, which is the convention main.cpp works in.
    tensors = []
    for i in range(psi.L):
        b = psi.get_B(i, form="B").transpose(["vL", "p", "vR"]).to_ndarray()
        tensors.append(np.ascontiguousarray(b[:, ::-1, :]))

    # Sweep left to right with QR, which changes the gauge but not the state
    for i in range(len(tensors) - 1):
        chi_l, phys, chi_r = tensors[i].shape
        q, r = np.linalg.qr(tensors[i].reshape(chi_l * phys, chi_r))
        tensors[i] = np.ascontiguousarray(q.reshape(chi_l, phys, -1))
        tensors[i + 1] = np.ascontiguousarray(np.tensordot(r, tensors[i + 1], axes=([1], [0])))

    # The final tensor carries whatever norm is left over
    tensors[-1] = tensors[-1] / np.linalg.norm(tensors[-1])
    return tensors


def check_left_canonical(tensors, tol=1.0e-8):
    """Verify sum_s A_s^dag A_s = I, which the C++ reader relies on."""

    worst = 0.0
    for a in tensors:
        chi_r = a.shape[2]
        gram = np.einsum("lsr,lst->rt", a.conj(), a)
        worst = max(worst, np.abs(gram - np.eye(chi_r)).max())
    if worst > tol:
        print(f"  WARNING: left-canonical check failed, max deviation {worst:.3e}")
    return worst


PAULIS = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def expectation(tensors, term):
    """
    Expectation of a Pauli string from left-canonical tensors, mirroring the
    algorithm in src/mps.h so the two can be checked against each other.

    term is a list of (pauli_letter, zero_based_site).
    """

    if not term:
        return 1.0

    sites = {i: p for p, i in term}
    lo, hi = min(sites), max(sites)

    # Everything left of the support contracts to the identity, because the
    # tensors are left-canonical
    left = np.eye(tensors[lo].shape[0], dtype=complex)
    for i in range(lo, hi + 1):
        a = tensors[i]
        op = PAULIS[sites.get(i, "I")]
        # L'[r,q] = sum conj(A[l,s,r]) L[l,m] O[s,t] A[m,t,q]
        left = np.einsum("lsr,lm,st,mtq->rq", a.conj(), left, op, a, optimize=True)

    # Everything right of the support contracts to the right environment, built
    # by sweeping inwards from the trivial bond at the far right
    right = np.eye(tensors[-1].shape[2], dtype=complex)
    for i in range(len(tensors) - 1, hi, -1):
        a = tensors[i]
        right = np.einsum("lsr,rt,mst->lm", a.conj(), right, a, optimize=True)

    return float(np.real(np.sum(left * right)))


def energy_from_tensors(tensors, num_sites, j1, j2):
    """
    Rebuild the Hamiltonian expectation from the exported tensors alone, using
    the same contraction the C++ uses. This exercises X, Y and Z, so it is the
    real end-to-end check that the export is correct.
    """

    total = 0.0
    for coupling, offset in ((j1, 1), (j2, 2)):
        if coupling == 0.0:
            continue
        for i in range(num_sites - offset):
            j = i + offset
            for pauli in "XYZ":
                total += coupling * 0.25 * expectation(tensors, [(pauli, i), (pauli, j)])
    return total


def validate(tensors, psi, num_checks=8):
    """
    Cross-check the exported tensors against TeNPy's own expectation values.

    Only Sz-type terms can be asked of TeNPy directly when charge conservation
    is on, so the X and Y sectors are covered by the energy check instead.
    """

    worst = 0.0
    rng = np.random.default_rng(0)
    checks = [[("Z", i)] for i in range(min(3, psi.L))]
    for _ in range(num_checks):
        i, j = sorted(rng.choice(psi.L, size=2, replace=False))
        checks.append([("Z", int(i)), ("Z", int(j))])

    for term in checks:
        mine = expectation(tensors, term)
        # Z = -2 Sz, the sign coming from TeNPy's 'down, up' basis ordering
        # which extract_tensors flips when exporting
        tenpy_term = [("S" + p.lower(), i) for p, i in term]
        theirs = np.real(psi.expectation_value_term(tenpy_term)) * ((-2) ** len(term))
        worst = max(worst, abs(mine - theirs))

    return worst


def write_mps(filename, tensors, energy, num_sites, meta):
    """Write the plain-text MPS format consumed by src/mps.h."""

    with open(filename, "w") as f:
        f.write("# thermoNPA MPS (left-canonical)\n")
        for line in meta:
            f.write(f"# {line}\n")
        f.write("# numSites physDim energy\n")
        # .17g round-trips an IEEE double exactly, and avoids numpy's repr
        f.write(f"{num_sites} {tensors[0].shape[1]} {float(energy):.17g}\n")
        for i, a in enumerate(tensors):
            chi_l, phys, chi_r = a.shape
            f.write(f"# site {i}: chiL chiR then chiL*d*chiR (real imag) ordered (l, s, r)\n")
            f.write(f"{chi_l} {chi_r}\n")
            flat = np.ascontiguousarray(a).reshape(-1)
            parts = np.empty(2 * flat.size, dtype=float)
            parts[0::2] = flat.real
            parts[1::2] = flat.imag
            f.write(" ".join(f"{v:.17g}" for v in parts))
            f.write("\n")


def main():
    parser = argparse.ArgumentParser(description="DMRG ground state of the J1-J2 chain")
    parser.add_argument("-n", "--sites", type=int, default=50, help="number of sites")
    parser.add_argument("--j1", type=float, default=1.0, help="nearest-neighbour coupling")
    parser.add_argument("--j2", type=float, default=0.35, help="next-nearest-neighbour coupling")
    parser.add_argument("--chi", type=int, default=128, help="maximum bond dimension")
    parser.add_argument("--sweeps", type=int, default=100, help="maximum DMRG sweeps")
    parser.add_argument("--conserve", default="Sz", help="conserved quantity ('Sz' or 'None')")
    parser.add_argument("-o", "--out", default=None, help="output file")
    args = parser.parse_args()

    conserve = None if args.conserve in ("None", "none", "") else args.conserve
    out = args.out or f"data/j1j2_{args.sites}_{args.j2}.mps"

    print(f"Running DMRG: N={args.sites}, J1={args.j1}, J2={args.j2}, chi<={args.chi}")
    energy, psi = run_dmrg(args.sites, args.j1, args.j2, args.chi, args.sweeps, conserve)
    energy = float(energy)
    print(f"  Ground state energy       = {energy!r}")
    print(f"  Energy per site           = {energy / args.sites!r}")
    print(f"  Max bond dimension used   = {max(psi.chi)}")
    print(f"  Max entanglement entropy  = {max(psi.entanglement_entropy()):.6f}")

    tensors = extract_tensors(psi)
    print(f"  Left-canonical deviation  = {check_left_canonical(tensors):.3e}")
    print(f"  Validation vs TeNPy       = {validate(tensors, psi):.3e}")
    rebuilt = energy_from_tensors(tensors, args.sites, args.j1, args.j2)
    print(f"  Energy rebuilt from MPS   = {rebuilt!r}  (diff {abs(rebuilt - energy):.3e})")

    meta = [
        "H = J1 sum_i S_i.S_{i+1} + J2 sum_i S_i.S_{i+2}, open boundaries",
        f"J1 = {args.j1}, J2 = {args.j2}, N = {args.sites}, chi_max = {max(psi.chi)}",
    ]
    write_mps(out, tensors, energy, args.sites, meta)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
