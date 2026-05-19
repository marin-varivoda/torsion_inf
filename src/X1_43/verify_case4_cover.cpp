#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


/**
 * Let plc1 be the list of degree 1 places of X1(43) over F2, as defined in magma code:
 *      plc1  := Places(Curve(C),  1);
 *
 * Let S = { D : D >= 0, deg D = 18, D is supported by places from plc1,
 *               D has at least 7 different places in support,
 *               plc1[1] has max multiplicity }
 *
 * This program will iterate over all elements in s (2317189897 of them) and check that 
 * each is covered by some element defined in `case4_cover.txt`.
 *
 */ 


constexpr int MAX_T_DEGREE = 57;
constexpr int MIN_SUPPORT = 7;
constexpr int N_PLACES = 21;
constexpr int N_OTHER = 20;
constexpr int DEGREE_S = 18;

// Safety bound for bitset implementation, the cover has ~700 elements, so 4096 is plenty
constexpr int MAX_COVER_SIZE = 4096;
constexpr int WORDS = (MAX_COVER_SIZE + 63) / 64;

using Divisor = std::array<int, N_PLACES>;
using Bits = std::array<std::uint64_t, WORDS>;

std::vector<Divisor> T;
int numWords = 0;

// atLeast[i][m] is the set of cover divisors t such that t[i] >= m.
Bits atLeast[N_PLACES][DEGREE_S + 1];

std::uint64_t totalChecked = 0;
std::uint64_t currentKChecked = 0;
constexpr std::uint64_t PROGRESS_EVERY = 10000000ULL;

auto startTime = std::chrono::steady_clock::now();

int degree(const Divisor &D)
{
    int s = 0;
    for (int x : D) {
        s += x;
    }
    return s;
}

void setBit(Bits &b, int i)
{
    b[i / 64] |= (1ULL << (i % 64));
}

bool nonempty(const Bits &b)
{
    for (int i = 0; i < numWords; ++i) {
        if (b[i] != 0) {
            return true;
        }
    }
    return false;
}

Bits intersect(const Bits &a, const Bits &b)
{
    Bits c{};
    for (int i = 0; i < numWords; ++i) {
        c[i] = a[i] & b[i];
    }
    return c;
}

void printDivisor(const Divisor &D)
{
    for (int i = 0; i < N_PLACES; ++i) {
        if (i)
            std::cerr << " ";
        std::cerr << D[i];
    }
    std::cerr << "\n";
}

std::vector<Divisor> loadCover(const std::string &filename)
{
    std::ifstream in(filename);

    if (!in) {
        std::cerr << "Could not open " << filename << "\n";
        std::exit(1);
    }

    std::vector<Divisor> cover;
    std::string line;
    int lineNumber = 0;

    while (std::getline(in, line)) {
        lineNumber++;

        if (line.empty()) {
            continue;
        }

        std::istringstream ss(line);
        Divisor D{};

        for (int i = 0; i < N_PLACES; ++i) {
            if (!(ss >> D[i])) {
                std::cerr << "Bad line " << lineNumber << ": expected 21 integers.\n";
                std::exit(1);
            }

            if (D[i] < 0) {
                std::cerr << "Bad line " << lineNumber << ": negative coefficient.\n";
                std::exit(1);
            }
        }

        std::string extra;
        if (ss >> extra) {
            std::cerr << "Bad line " << lineNumber << ": extra data after 21 integers.\n";
            std::exit(1);
        }

        if (degree(D) > MAX_T_DEGREE) {
            std::cerr << "Degree bound failed on line " << lineNumber << ".\n";
            std::cerr << "Degree = " << degree(D) << "\n";
            std::exit(1);
        }

        cover.push_back(D);
    }

    if (cover.empty()) {
        std::cerr << "Cover file is empty.\n";
        std::exit(1);
    }

    if (cover.size() > MAX_COVER_SIZE) {
        std::cerr << "Cover has " << cover.size() << " elements, but MAX_COVER_SIZE is only "
                  << MAX_COVER_SIZE << ".\n";
        std::exit(1);
    }

    return cover;
}

void buildBitsets()
{
    numWords = static_cast<int>((T.size() + 63) / 64);

    for (int j = 0; j < static_cast<int>(T.size()); ++j) {
        for (int i = 0; i < N_PLACES; ++i) {
            for (int m = 0; m <= DEGREE_S; ++m) {
                if (T[j][i] >= m) {
                    setBit(atLeast[i][m], j);
                }
            }
        }
    }
}

void printProgress(int k)
{
    auto now = std::chrono::steady_clock::now();
    double seconds = std::chrono::duration<double>(now - startTime).count();

    double rate = totalChecked / std::max(seconds, 1e-9);

    std::cout << "k = " << k << ", checked for this k = " << currentKChecked
              << ", total checked = " << totalChecked << ", rate = " << std::fixed
              << std::setprecision(0) << rate << " divisors/sec\n";
}

bool verifyRecursively(int k,
                       int pos,
                       int remainingDegree,
                       int supportAwayFromP1,
                       const Bits &candidates,
                       Divisor &current)
{
    if (remainingDegree < 0) {
        return true;
    }

    int positionsLeft = N_OTHER - pos;

    // Even putting k at every remaining position cannot reach the required degree.
    if (remainingDegree > positionsLeft * k) {
        return true;
    }

    // Even using all remaining positions, we cannot reach support size 7.
    if (supportAwayFromP1 + std::min(remainingDegree, positionsLeft) < MIN_SUPPORT - 1) {
        return true;
    }

    // Only now do we know this partial divisor can still lead to an element of S.
    if (!nonempty(candidates)) {
        std::cerr << "No cover divisor can dominate this partial divisor.\n";
        std::cerr << "Partial divisor:\n";
        printDivisor(current);
        return false;
    }

    if (pos == N_OTHER) {
        if (remainingDegree != 0) {
            return true;
        }

        if (supportAwayFromP1 < MIN_SUPPORT - 1) {
            return true;
        }

        totalChecked++;
        currentKChecked++;

        if (totalChecked % PROGRESS_EVERY == 0) {
            printProgress(k);
        }

        return true;
    }

    int coord = pos + 1;

    for (int m = 0; m <= k && m <= remainingDegree; ++m) {
        current[coord] = m;

        int newSupport = supportAwayFromP1 + (m > 0);

        if (m <= 2) {
            if (!verifyRecursively(k, pos + 1, remainingDegree - m, newSupport, candidates, current)) {
                return false;
            }
        } else {
            Bits nextCandidates = intersect(candidates, atLeast[coord][m]);

            if (!verifyRecursively(k,
                                   pos + 1,
                                   remainingDegree - m,
                                   newSupport,
                                   nextCandidates,
                                   current)) {
                return false;
            }
        }
    }

    current[coord] = 0;
    return true;
}

bool allCoverDivisorsHaveBaseTwoAwayFromP1()
{
    for (const Divisor &D : T) {
        for (int i = 1; i < N_PLACES; ++i) {
            if (D[i] < 2) {
                return false;
            }
        }
    }

    return true;
}

int main()
{
    T = loadCover("case4_cover.txt");

    std::cout << "Loaded " << T.size() << " cover divisors.\n";
    std::cout << "All have degree <= " << MAX_T_DEGREE << ".\n";

    if (!allCoverDivisorsHaveBaseTwoAwayFromP1()) {
        std::cerr << "This verifier assumes every cover divisor has coefficient >= 2 "
                  << "at P_2,...,P_21.\n";
        std::cerr << "That should be true for the cover generated by the construction.\n";
        return 1;
    }

    buildBitsets();

    for (int k = 1; k <= 12; ++k) {
        std::cout << "Starting verification for k = " << k << "...\n";

        currentKChecked = 0;

        Divisor current{};
        current[0] = k;

        Bits candidates = atLeast[0][k];

        bool ok = verifyRecursively(k, 0, DEGREE_S - k, 0, candidates, current);

        if (!ok) {
            std::cerr << "Verification failed for k = " << k << ".\n";
            return 1;
        }

        std::cout << "Finished k = " << k << ", checked " << currentKChecked << " divisors.\n";
    }

    std::cout << "SUCCESS: the cover dominates every divisor in S.\n";
    std::cout << "Total divisors checked: " << totalChecked << "\n";

    return 0;
}
