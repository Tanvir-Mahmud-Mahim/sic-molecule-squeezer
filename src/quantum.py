"""Linearized quantum fluctuation analysis around a Kerr-comb mean field.

Follows the input-output treatment of Guidry et al., Nat. Photon. 16, 52 (2022)
and Optica 10, 694 (2023): linearize the intracavity Kerr Hamiltonian around
the classical comb amplitudes, solve the multimode Heisenberg-Langevin
equations in the frequency domain, and extract squeezing spectra, supermodes,
and logarithmic negativity at the collection port.

Working units: rates in units of kappa (the reference loaded linewidth).
Normalized comb amplitudes psi_mu satisfy g0|A|^2 = (kappa/2)|psi|^2, so all
nonlinear couplings appear as (1/2)*psi products in kappa units.

Mode ordering: mu = -M..M (2M+1 modes), index i = mu + M.
Doubled vector v = (b_{-M}..b_{M}, b\dagger_{-M}..b\dagger_{M}).
Each mode has two baths: intrinsic (rate kappa_i_mu) and collection
(rate kappa_c_mu); total kappa_mu = kappa_i_mu + kappa_c_mu.
"""
import numpy as np


def coupling_matrices(psi_mu_centered, M):
    """P_{qk} = sum_{m+n=q+k} psi_m psi_n ; R_{qn} = sum_j psi_j^* psi_{j+q-n}.
    psi_mu_centered: complex array indexed mu=-Mc..Mc (length 2Mc+1) of the
    classical comb; returns (2M+1)x(2M+1) matrices for fluctuation modes."""
    Mc = (len(psi_mu_centered) - 1) // 2

    def psi(m):
        return psi_mu_centered[m + Mc] if -Mc <= m <= Mc else 0.0

    n = 2 * M + 1
    P = np.zeros((n, n), complex)
    R = np.zeros((n, n), complex)
    for q in range(-M, M + 1):
        for k in range(-M, M + 1):
            s = q + k
            P[q + M, k + M] = sum(psi(m) * psi(s - m)
                                  for m in range(s - Mc, Mc + 1)
                                  if -Mc <= m <= Mc and -Mc <= s - m <= Mc)
            d = q - k
            R[q + M, k + M] = sum(np.conj(psi(j)) * psi(j + d)
                                  for j in range(-Mc, Mc + 1)
                                  if -Mc <= j + d <= Mc)
    return P, R


class CombQuantum:
    def __init__(self, psi_mu_centered, zeta_mu, kappa_i, kappa_c, M):
        """zeta_mu: 2*(delta0+Dint)/kappa for mu=-M..M (length 2M+1).
        kappa_i, kappa_c: arrays (units of kappa)."""
        self.M = M
        n = 2 * M + 1
        self.n = n
        P, R = coupling_matrices(psi_mu_centered, M)
        Delta = 0.5 * np.asarray(zeta_mu)          # detuning in kappa units
        self.kappa_i = np.asarray(kappa_i, float)
        self.kappa_c = np.asarray(kappa_c, float)
        kt = self.kappa_i + self.kappa_c
        # db_q/dt = (-kt/2 - i Delta_q) b_q + i*(1/2)[ sum_k P_qk b_k^dag
        #            + 2 sum_k R_qk b_k ] + noise    (kappa units)
        A = np.diag(-kt / 2 - 1j * Delta) + 1j * R      # R carries the *2/2
        # note: i*(1/2)*2*R = i*R
        B = 0.5j * P
        self.Mdrift = np.block([[A, B], [np.conj(B), np.conj(A)]])

    def scattering(self, omega):
        """S(omega): maps (b_in^i, b_in^c ; daggers) -> (b_out^c ; daggers).
        omega in kappa units. Returns Sc (n x 2n complex for annihilators...)
        Full doubled form: returns 2n x 4n matrix T such that
        v_out_c = T v_in, v_in = (bi, bc, bi^dag, bc^dag)."""
        n = self.n
        sqki = np.sqrt(self.kappa_i)
        sqkc = np.sqrt(self.kappa_c)
        # noise input matrix for doubled vector
        Xi = np.zeros((2 * n, 4 * n))
        Xi[:n, :n] = np.diag(sqki)
        Xi[:n, n:2 * n] = np.diag(sqkc)
        Xi[n:, 2 * n:3 * n] = np.diag(sqki)
        Xi[n:, 3 * n:] = np.diag(sqkc)
        G = np.linalg.solve(-1j * omega * np.eye(2 * n) - self.Mdrift, Xi)
        # output at collection port: b_out = sqrt(kc) b - b_in^c
        Pc = np.zeros((2 * n, 2 * n))
        Pc[:n, :n] = np.diag(sqkc)
        Pc[n:, n:] = np.diag(sqkc)
        T = Pc @ G
        T[:n, n:2 * n] -= np.eye(n)
        T[n:, 3 * n:] -= np.eye(n)
        return T

    def quad_covariance(self, omega, detect):
        """Symmetrized quadrature noise-spectral matrix of the collection-port
        output restricted to detected modes (indices into 0..n-1).
        Vacuum inputs. Returns V (2d x 2d), vacuum = identity."""
        n = self.n
        T = self.scattering(omega)
        d = len(detect)
        rows = list(detect) + [n + i for i in detect]
        Td = T[rows, :]
        # <v_in v_in^dag> for vacuum in doubled basis (b, b^dag):
        # <b b^dag>=1, others 0 => C_in = diag(I_{2n}, 0_{2n}) in ordering
        # (bi, bc, bi+, bc+): <v v^dag> = blockdiag(I_2n, 0)
        Cin = np.zeros((4 * n, 4 * n))
        Cin[:2 * n, :2 * n] = np.eye(2 * n)
        Cbb = Td @ Cin @ Td.conj().T       # <v_out v_out^dag>
        # quadratures q = L v, L = [[1,1],[-i,i]]/sqrt2 blockwise
        I = np.eye(d)
        L = np.block([[I, I], [-1j * I, 1j * I]]) / np.sqrt(2)
        Vq = L @ Cbb @ L.conj().T
        Vq = np.real(Vq + Vq.conj().T)  # 2x symmetrized: vacuum -> identity
        return Vq

    def quad_covariance_sym(self, omega, detect):
        """Two-sided (homodyne-measured) spectrum: average of +/- omega."""
        V = 0.5 * (self.quad_covariance(omega, detect)
                   + self.quad_covariance(-omega, detect))
        return V

    def squeezing_spectrum(self, omegas, detect):
        """min/max eigenvalue of V(omega) over the detected subspace."""
        smin, smax = [], []
        for om in omegas:
            w = np.linalg.eigvalsh(self.quad_covariance_sym(om, detect))
            smin.append(w[0])
            smax.append(w[-1])
        return np.array(smin), np.array(smax)

    def supermodes(self, omega, detect, k=4):
        """Eigen-decomposition of V(omega): returns eigenvalues (sorted
        ascending) and eigenvectors for the k most-squeezed quadratures."""
        V = self.quad_covariance_sym(omega, detect)
        w, U = np.linalg.eigh(V)
        return w[:k], U[:, :k]

    def steady_covariance(self):
        """Intracavity symmetrized covariance from the Lyapunov equation
        M V + V M^T + D = 0 in the quadrature basis (vacuum=identity/1?).
        Quadrature convention: V_vac = I. Returns 2n x 2n real matrix."""
        from scipy.linalg import solve_lyapunov
        n = self.n
        I = np.eye(n)
        L = np.block([[I, I], [-1j * I, 1j * I]]) / np.sqrt(2)
        Mq = np.real(L @ self.Mdrift @ np.linalg.inv(L))
        kt = self.kappa_i + self.kappa_c
        D = np.kron(np.eye(2), np.diag(kt))
        return solve_lyapunov(Mq, -D)


def log_negativity(V, i, j, n):
    """E_N of modes i,j from full 2n x 2n quadrature covariance V
    (ordering x_0..x_{n-1}, p_0..p_{n-1}; vacuum = identity)."""
    idx = [i, j, n + i, n + j]
    Vs = V[np.ix_(idx, idx)]              # (x1,x2,p1,p2)
    A = Vs[np.ix_([0, 2], [0, 2])]
    B = Vs[np.ix_([1, 3], [1, 3])]
    C = Vs[np.ix_([0, 2], [1, 3])]
    delta = np.linalg.det(A) + np.linalg.det(B) - 2 * np.linalg.det(C)
    detV = np.linalg.det(Vs)
    arg = delta ** 2 - 4 * detV
    if arg < 0:
        return 0.0
    nu_minus2 = (delta - np.sqrt(arg)) / 2
    if nu_minus2 <= 0:
        return 0.0
    nu = np.sqrt(nu_minus2)
    return max(0.0, -np.log2(nu))
