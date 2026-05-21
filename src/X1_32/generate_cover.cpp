#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

/**
 * Let S = { D : D >= 0, deg(D) = 9 }.
 *
 * We construct a covering set T such that:
 *  - every divisor in T is effective
 *  - every divisor in T has degree at most 18
 *  - T "dominates" S, i.e, for each s in S, there exists t in T such that s <= t
 *
 * We enumerate all divisors in S, randomly order them with a fixed seed. We maintain
 * a boolean array keeping track of which divisors in S are already covered. The idea
 * is to take an uncovered divisor E, we then scan through remaining uncovered divisors D.
 * If we can find D such that max(E, D) has degree <= 18, we replace E by max(E, D), mark
 * D as covered and scan again. Once it's not possible to further enlarge E, we add E to
 * T and mark every divisor D <= T as covered.
 *
 */

constexpr int TARGET_DEG = 9;
constexpr int MAX_T_DEGREE = 18;
constexpr uint64_t SEED = 1337;

// How many places of each degree X1(32) has over F3
constexpr int PLACE_COUNTS[TARGET_DEG + 1] = {
    0,
    12,  // degree 1
    4,   // degree 2
    8,   // degree 3
    9,   // degree 4
    16,  // degree 5
    140, // degree 6
    312, // degree 7
    901, // degree 8
    1976 // degree 9
};

struct Place
{
    int degree;
    int index; // valid values are [1..#places of this degree], i.e., 1-indexed
};

std::vector<Place> buildPlaces()
{
    std::vector<Place> ret;
    for (int d = 1; d <= TARGET_DEG; d++) {
        for (int idx = 1; idx <= PLACE_COUNTS[d]; idx++) {
            ret.push_back(Place{d, idx});
        }
    }
    return ret;
}

/* Global list of all possible places */
std::vector<Place> places = buildPlaces();

typedef int PlaceID; // referes to index in places
struct Divisor
{
    PlaceID places[40]; // in practice this array will have at most 18 elements, we add 40 for safety buffer
    int placesLen = 0; // how long the array actually is

    int degree = 0;

    bool operator==(const Divisor &other) const
    {
        if (this->degree != other.degree || this->placesLen != other.placesLen)
            return false;

        for (int i = 0; i < placesLen; i++) {
            if (this->places[i] != other.places[i])
                return false;
        }

        return true;
    }

    static_assert(sizeof(places)/sizeof(places[1]) >= 2*MAX_T_DEGREE);
};

// Set of all effective divisors of degree 9
std::vector<Divisor> S;

void buildSRecursive(PlaceID lastAddedPlaceID, Divisor &D)
{
    if (D.degree == TARGET_DEG) {
        S.push_back(D);
        return;
    } else if (D.degree > TARGET_DEG) {
        return;
    }

    const int origPlacesLen = D.placesLen;
    const int origDeg = D.degree;

    for (PlaceID i = lastAddedPlaceID; i < places.size(); i++) {
        int deg = places[i].degree;

        // we can stop iterating here, since places are ordered by degree
        // so for each consecutive place, we would still have D.degree + deg > TARGET_DEG
        if (D.degree + deg > TARGET_DEG)
            break;

        // otherwise try to add i to D
        D.places[D.placesLen] = i;
        D.placesLen++;
        D.degree += deg;
        buildSRecursive(i, D);

        // Discard the added place
        D.placesLen = origPlacesLen;
        D.degree = origDeg;
    }
}

// Returns the smallest C such that C >= A and C >= B (A, B, C are all effective)
Divisor join(const Divisor &A, const Divisor &B)
{
    Divisor C;

    int i = 0;
    int j = 0; // pointer to place in A and place in B

    while (i < A.placesLen || j < B.placesLen) {
        // the smaller of these will be a valid place
        const PlaceID A_ID = (i < A.placesLen) ? A.places[i] : INT_MAX;
        const PlaceID B_ID = (j < B.placesLen) ? B.places[j] : INT_MAX;

        if (A_ID < B_ID) {
            const int deg = places[A_ID].degree;
            C.places[C.placesLen] = A_ID;
            C.placesLen++;
            C.degree += deg;

            i++;
        } else if (B_ID < A_ID) {
            const int deg = places[B_ID].degree;
            C.places[C.placesLen] = B_ID;
            C.placesLen++;
            C.degree += deg;

            j++;
        } else if (A_ID == B_ID) {
            const int deg = places[A_ID].degree;
            C.places[C.placesLen] = A_ID;
            C.placesLen++;
            C.degree += deg;

            i++;
            j++;
        }
    }

    return C;
}

// Checks if A >= B
bool dominates(const Divisor &A, const Divisor &B)
{
    Divisor C = join(A, B);
    return (A == C);
}

void exportCover(const std::vector<Divisor> &T, std::filesystem::path outFile);

int main()
{
    Divisor tmp;
    buildSRecursive(0, tmp);

    std::cout << "Found total of " << S.size() << " degree 9 effective divisors." << std::endl;

    // Random shuffle S
    std::mt19937_64 rng(SEED);
    std::shuffle(S.begin(), S.end(), rng);

    // Prepare a bool vector to keep track which elements from S have been covered
    std::vector<bool> isCovered(S.size(), false);
    std::vector<Divisor> T;

    for (int i = 0; i < S.size(); i++) {
        if (isCovered[i])
            continue;

        Divisor E = S[i]; // we found an uncovered divisor
        isCovered[i] = true;

        // enlarge E as much as possible
        for (int j = i + 1; j < S.size(); j++) {
            if (isCovered[j])
                continue;

            Divisor E_candidate = join(E, S[j]);
            if (E_candidate.degree <= MAX_T_DEGREE) {
                E = E_candidate;
                isCovered[j] = true;
            }
        }

        // After this is done, add E to the final set T
        T.push_back(E);

        // Make sure to mark all the newly covered elements
        // we start iterating at i+1, since all the previous ones are already covered
        for (int k = i + 1; k < S.size(); k++) {
            isCovered[k] = isCovered[k] || dominates(E, S[k]);
        }
    }

    // Sanity checks
    for (bool b : isCovered) {
        if (!b) {
            std::cerr << "ERROR: not all divisors were covered." << std::endl;
            std::exit(1);
        }
    }
    for (const Divisor &E : T) {
        if (E.degree > MAX_T_DEGREE) {
            std::cerr << "ERROR: divisor of degree " << E.degree << " exceeds " << MAX_T_DEGREE << "." << std::endl;
            std::exit(1);
        }
    }

    std::cout << "Final T size: " << T.size() << "\n";

    exportCover(T, "deg9_divisors_cover.txt");
    std::cout << "Exported cover to deg9_divisors_cover.txt" << std::endl;

    return 0;
}

void exportCover(const std::vector<Divisor> &T, const std::filesystem::path& outFile)
{
    std::ofstream out(outFile);
    if (!out) {
        std::cerr << "ERROR: Could not open file " << outFile << " for writing.\n";
        std::exit(1);
    }

    auto placeName = [](const PlaceID id) -> std::string {
        const Place &P = places[id];
        return "P_" + std::to_string(P.degree) + "_" + std::to_string(P.index);
    };

    for (const Divisor &D : T) {
        for (int i = 0; i < D.placesLen; i++) {
            if (i != 0) {
                out << ",";
            }
            out << placeName(D.places[i]);
        }
        out << "\n";
    }
}
