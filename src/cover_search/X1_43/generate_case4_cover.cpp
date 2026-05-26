#include <array>
#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

/**
 * Let plc1 be the list of degree 1 places of X1(43) over F2, as defined in magma code:
 *      plc1  := Places(Curve(C),  1);
 *
 * Let S = { D : D >= 0, deg D = 18, D is supported by places from plc1,
 *               D has at least 7 different places in support,
 *               plc1[1] has max multiplicity }
 *
 *
 * We will construct a set T such that:
 *  - T consists of effective divisors supported on plc1
 *  - each divisor in T has degree at most 57
 *  - T "dominates" S, i.e, for each s in S, there exists t in T such that s <= t
 *
 *
 * Here we describe the approach for constructing T.
 *
 * We write a divisor in S as D = k*P_1 + sum_{i=2}^{21} m_i*P_i, where P_i = plc1[i].
 * Since D has at least 7 support points and P_1 has maximal multiplicity, we have k <= 12.
 * (note that 12+1+1+1+1+1+1=18).
 *
 * For each k, we enumerate all "excess" divisors e = sum_{i=2}^{21} e_i*P_i, where
 * e_i = max(m_i - 2, 0). This measures by how much (k-2) * P1 + 2 * plc1sum fails to cover D.
 *
 * We then greedily join these excess divisors. Start with an uncovered excess divisor f = e. Scan
 * through remaining uncovered excess divisors e', and replace f by max(f, e') whenever
 * this maintains deg(f) <= 17-k. Here max is taken coefficientwise. Once it's no longer possible
 * to enlarge f, we add a new element in T:
 *      (k-2)*P_1 + 2*plc1sum + f
 *
 * We repeat this until every excess divisor e is dominated by some chosen f.
 *
 */

constexpr int MAX_T_DEGREE = 57;
constexpr int MIN_SUPPORT = 7;
constexpr int N_OTHER = 20;
constexpr int DEGREE_S = 18;

using Coeffs = std::array<int, N_OTHER>; // multiplicities at plc1[2], ..., plc1[21]
using CoeffsFull = std::array<int, N_OTHER+1>; // multiplicities at plc1[1], plc1[2], ..., plc1[21]

int minPossibleDegreeForExcess(int k, const Coeffs& e)
{
    int excessSupport = 0; // how many non-zero excess coeffs
    int totalExcess = 0; // sum of all excess coeffs
    for (int x : e) {
        if (x != 0) {
            excessSupport++;
            totalExcess += x;
        }
    }

    // P_1 is already in the support, so we need at least 6 more support points.
    int extraSupportNeeded = std::max(0, 6 - excessSupport);

    // Min possible degree of the divisor in S with current excess vector e
    return k + 2*excessSupport + totalExcess + extraSupportNeeded;
}

bool isExcessVectorFeasible(int k, const Coeffs &e)
{
    for (int x : e) {
        // If e_i > 0, then the actual multiplicity is 2+e_i, this can't exceed k
        if (x != 0 && 2 + x > k) {
            return false;
        }
    }

    return minPossibleDegreeForExcess(k, e) <= DEGREE_S;
}

void enumerateExcessDivisors(int k, int firstFreePos, Coeffs &current, std::vector<Coeffs> &out)
{
    if (!isExcessVectorFeasible(k, current)) {
        return;
    }

    // The current excess vector can occur, so record it.
    out.push_back(current);

    const int maxExcessAtOnePlace = k - 2;
    if (maxExcessAtOnePlace <= 0) {
        return;
    }

    // Add one new positive excess coefficient at a later position.
    for (int pos = firstFreePos; pos < N_OTHER; ++pos) {
        for (int x = 1; x <= maxExcessAtOnePlace; ++x) {
            current[pos] = x;
            enumerateExcessDivisors(k, pos + 1, current, out);
            current[pos] = 0;
        }
    }
}

// Find all possible excess divisors for a given k
std::vector<Coeffs> allExcessDivisors(int k)
{
    std::vector<Coeffs> ret;
    Coeffs current = {};
    enumerateExcessDivisors(k, 0, current, ret);

    std::cout << "Found: " << ret.size() << std::endl;

    return ret;
}

Coeffs join(const Coeffs &a, const Coeffs &b)
{
    Coeffs ret;
    for (size_t i = 0; i < ret.size(); i++) {
        ret[i] = std::max(a[i], b[i]);
    }
    return ret;
}

int sumCoeffs(const Coeffs &a)
{
    int ret = 0;
    for (int x : a) {
        ret += x;
    }
    return ret;
}

bool dominates(const Coeffs &a, const Coeffs &b)
{
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] < b[i]) {
            return false;
        }
    }

    return true;
}

std::vector<Coeffs> greedyDominatingSet(const std::vector<Coeffs>& E, int k)
{
    // We are adding divisors of form k*P_1 + 2*P_2 + ... + 2*P_21 + f to T
    const int fBudget = MAX_T_DEGREE - k - 2*N_OTHER;

    std::vector<bool> covered(E.size(), 0);
    std::vector<Coeffs> ECover;

    for (std::size_t i = 0; i < E.size(); ++i) {
        if (covered[i]) {
            continue;
        }

        Coeffs f = E[i];
        covered[i] = true;

        // Greedily enlarge f by joining it with uncovered excess divisors,
        // as long as deg(f) <= fBudget
        bool enlarged = true;
        while (enlarged && sumCoeffs(f) < fBudget) {
            enlarged = false;

            // Check if there is an element in E that we can use to enlarge f
            for (std::size_t j = 0; j < E.size(); j++) {
                if (covered[j]) {
                    continue;
                }

                Coeffs candidate = join(f, E[j]);

                if (sumCoeffs(candidate) <= fBudget) {
                    f = candidate;
                    covered[j] = true;
                    enlarged = true;
                    break;
                }
            }
        }

        // Mark all divisors covered by final f
        for (std::size_t j = 0; j < E.size(); ++j) {
            if (!covered[j] && dominates(f, E[j])) {
                covered[j] = true;
            }
        }

        ECover.push_back(f);
    }

    return ECover;
}

// Function to verify that we didn't make a mistake
bool isCoverCorrect(const std::vector<Coeffs> &E, const std::vector<Coeffs> &ECover)
{
    for (const Coeffs &e : E) {
        bool covered = false;

        for (const Coeffs &f : ECover) {
            if (dominates(f, e)) {
                covered = true;
                break;
            }
        }

        if (!covered) {
            return false;
        }
    }

    return true;
}

CoeffsFull makeTCoeffs(int k, const Coeffs &e)
{
    CoeffsFull ret;
    ret[0] = k;
    for (size_t i = 0; i < e.size(); i++) {
        ret[i + 1] = 2 + e[i];
    }
    return ret;
}

std::vector<CoeffsFull> pruneDominatedDivisors(const std::vector<CoeffsFull> &T)
{
    // Checks if A>=B
    auto dominates = [](const CoeffsFull &A, const CoeffsFull &B) {
        for (std::size_t i = 0; i < A.size(); ++i) {
            if (A[i] < B[i]) {
                return false;
            }
        }
        return true;
    };

    std::vector<CoeffsFull> pruned;

    for (std::size_t i = 0; i < T.size(); ++i) {
        bool redundant = false;

        for (std::size_t j = 0; j < T.size(); ++j) {
            if (i == j) {
                continue;
            }

            // Remove T[i] if a different divisor dominates it.
            // For equal divisors, keep only the first copy.
            if (dominates(T[j], T[i]) && (T[j] != T[i] || j < i)) {
                redundant = true;
                break;
            }
        }

        if (!redundant) {
            pruned.push_back(T[i]);
        }
    }

    return pruned;
}

void exportCoverToFile(const std::vector<CoeffsFull> &T, const std::string &filename)
{
    std::ofstream out(filename);

    if (!out) {
        std::cerr << "Could not open file " << filename << " for writing.\n";
        std::exit(1);
    }

    // One divisor per line, written as 21 multiplicities:
    // coeffs of plc1[1], plc1[2], ..., plc1[21].
    for (const CoeffsFull &t : T) {
        for (std::size_t i = 0; i < t.size(); ++i) {
            if (i != 0) {
                out << " ";
            }
            out << t[i];
        }
        out << "\n";
    }
}

int main()
{
    // final dominating set
    std::vector<CoeffsFull> T;

    for (int k = 1; k <= 12; ++k) {
        std::vector<Coeffs> E = allExcessDivisors(k);
        std::vector<Coeffs> ECover = greedyDominatingSet(E, k);
        assert(isCoverCorrect(E, ECover));

        for (const Coeffs& e : ECover) {
            T.push_back(makeTCoeffs(k, e));
        }

        std::cout << "k = " << k << ", excess divisors = " << E.size()
                  << ", cover elements = " << ECover.size() << "\n";
    }

    std::cout << "Before pruning: " << T.size() << "\n";
    T = pruneDominatedDivisors(T);
    std::cout << "After pruning: " << T.size() << "\n";

    exportCoverToFile(T, "case4_cover.txt");

    return 0;
}
