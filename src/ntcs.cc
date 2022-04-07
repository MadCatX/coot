// vim: set sw=4 ts=4 sts=4 expandtab :

#include "ntcs.h"

#include "coot-utils/atom-tree.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "ligand/monomer-utils.hh"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <iterator>
#include <stdexcept>

namespace IBT {

using IndexQuad = std::array<int, 4>;
using ContactIndices = std::vector<std::vector<int>>;

class Step {
public:
    std::vector<mmdb::PAtom> atoms;
    mmdb::PResidue firstResidue;
    mmdb::PResidue secondResidue;
};

static const std::array<AtomQuad, 7> BACKBONE_ATOMQUADS{{
    {{ Atom(" C5'", Atom::FIRST),  Atom(" C4'", Atom::FIRST),  Atom(" C3'", Atom::FIRST),  Atom(" O3'", Atom::FIRST)  }}, // delta_1
    {{ Atom(" C4'", Atom::FIRST),  Atom(" C3'", Atom::FIRST),  Atom(" O3'", Atom::FIRST),  Atom(" P  ", Atom::SECOND) }}, // epsilon_1
    {{ Atom(" C3'", Atom::FIRST),  Atom(" O3'", Atom::FIRST),  Atom(" P  ", Atom::SECOND), Atom(" O5'", Atom::SECOND) }}, // zeta_1
    {{ Atom(" O3'", Atom::FIRST),  Atom(" P  ", Atom::SECOND), Atom(" O5'", Atom::SECOND), Atom(" C5'", Atom::SECOND) }}, // alpha_2
    {{ Atom(" P  ", Atom::SECOND), Atom(" O5'", Atom::SECOND), Atom(" C5'", Atom::SECOND), Atom(" C4'", Atom::SECOND) }}, // beta_2
    {{ Atom(" O5'", Atom::SECOND), Atom(" C5'", Atom::SECOND), Atom(" C4'", Atom::SECOND), Atom(" C3'", Atom::SECOND) }}, // gamma_2
    {{ Atom(" C5'", Atom::SECOND), Atom(" C4'", Atom::SECOND), Atom(" C3'", Atom::SECOND), Atom(" O3'", Atom::SECOND) }}  // delta_2
}};

static const std::array<std::array<std::string, 4>, 5> NU_TORSIONS_ATOM_NAMES{{
    {{ " C4'", " O4'", " C1'", " C2'" }}, // nu_0
    {{ " O4'", " C1'", " C2'", " C3'" }}, // nu_1
    {{ " C1'", " C2'", " C3'", " C4'" }}, // nu_2
    {{ " C2'", " C3'", " C4'", " O4'" }}, // nu_3
    {{ " C3'", " C4'", " O4'", " C1'" }}, // nu_4
}};

static
const std::array<const double NtC::*, 7> BACKBONE_TORSION_VALUES{
    &NtC::delta_1, &NtC::epsilon_1, &NtC::zeta_1,
    &NtC::alpha_2, &NtC::beta_2, &NtC::gamma_2, &NtC::delta_2
};

using TorsionAtomNames = std::array<std::string, 4>;
static const TorsionAtomNames PURINE_ATOM_NAMES      { " O4'", " C1'", " N9 ", " C4 " };
static const TorsionAtomNames PYRIMIDINE_ATOM_NAMES  { " O4'", " C1'", " N1 ", " C2 " };
static const std::array<std::string, 4> PYRAMID_ATOMS{ " P  ", " OP1", " OP2", " O5'" };

static const std::array<std::string, 9> NUCLEOTIDES{ "A", "C", "G", "U", "DA", "DC", "DG", "DT", "T" };

template <typename T>
constexpr
T _A(const T &v) {
    return T(v > 180) * (v - 360) + T(v <= 180) * v;
}

template <size_t N>
static
void BackboneTorsionsGenerator(std::vector<Torsion> &torsions, const NtC &ntc) {
    torsions[N] = Torsion(BACKBONE_ATOMQUADS[N], ntc.*BACKBONE_TORSION_VALUES[N]);
    BackboneTorsionsGenerator<N-1>(torsions, ntc);
}
template <>
void BackboneTorsionsGenerator<0>(std::vector<Torsion> &torsions, const NtC &ntc) {
    torsions[0] = Torsion(BACKBONE_ATOMQUADS[0], ntc.*BACKBONE_TORSION_VALUES[0]);
}

static
std::vector<Torsion> BackboneTorsions(const NtC &ntc) {
    std::vector<Torsion> torsions{};
    torsions.resize(7);
    BackboneTorsionsGenerator<6>(torsions, ntc);
    return torsions;
}

template <size_t N>
static
void BaseTorsionGenerator(AtomQuad &quad, const TorsionAtomNames &atoms, const Atom::WhichResidue whichResidue) {
    quad[N].name = atoms[N];
    quad[N].whichResidue = whichResidue;
    BaseTorsionGenerator<N-1>(quad, atoms, whichResidue);
}
template <>
void BaseTorsionGenerator<0>(AtomQuad &quad, const TorsionAtomNames &atoms, const Atom::WhichResidue whichResidue) {
    quad[0].name = atoms[0];
    quad[0].whichResidue = whichResidue;
}

static
AtomQuad BaseTorsion(const TorsionAtomNames &atoms, const Atom::WhichResidue whichResidue) {
    AtomQuad quad = { Atom(), Atom(), Atom(), Atom() };
    BaseTorsionGenerator<3>(quad, atoms, whichResidue);
    return quad;
}

static
std::string trim(const std::string &str) {
    auto first = str.find_first_not_of(' ');
    auto last = str.find_last_not_of(' ');;
    return str.substr(first, last + 1);
}

static
std::vector<mmdb::PAtom> filterAtoms(mmdb::PPAtom atoms, size_t nAtoms, int seqNumFirst) {
    std::vector<mmdb::PAtom> filtered{};
    for (size_t idx = 0; idx < nAtoms; idx++) {
        const auto at = atoms[idx];

        /* We don't want the upper phosphate pyramid */
        if (at->GetSeqNum() == seqNumFirst && std::find(PYRAMID_ATOMS.cbegin(), PYRAMID_ATOMS.cend(), at->GetAtomName()) != PYRAMID_ATOMS.cend())
            continue;

        filtered.push_back(at);
    }

    return filtered;
}

static
clipper::Vec3<double> getC5PrimeCoords(const mmdb::PResidue residue) {
    auto at = residue->GetAtom(" C5'");
    return { at->x, at->y, at->z };
}

static
clipper::Vec3<double> getO3PrimeCoords(const mmdb::PResidue residue) {
    auto at = residue->GetAtom(" O3'");
    return { at->x, at->y, at->z };
}

static
mmdb::PAtom getAtom(const Atom &atom, const Step &step) {
    auto tgtSeqNum = atom.whichResidue == Atom::SECOND ? step.secondResidue->GetSeqNum() : step.firstResidue->GetSeqNum();

    for (auto &at : step.atoms) {
        auto seqNum = at->GetResidue()->GetSeqNum();

        if (seqNum == tgtSeqNum && at->GetAtomName() == atom.name)
            return at;
    }

    return nullptr;
}

static
int getAtomIndex(const Atom &atom, const Step &step) {
    auto tgtSeqNum = atom.whichResidue == Atom::SECOND ? step.secondResidue->GetSeqNum() : step.firstResidue->GetSeqNum();

    for (size_t idx = 0; idx < step.atoms.size(); idx++) {
        const auto at = step.atoms[idx];
        auto seqNum = at->GetResidue()->GetSeqNum();

        if (seqNum == tgtSeqNum && at->GetAtomName() == atom.name)
            return idx;
    }

    return -1;
}

static
IndexQuad getTorsionIndices(const Torsion &torsion, const Step &step) {
    IndexQuad aIx;
    for (size_t idx = 0; idx < 4; idx++) {
        auto x = getAtomIndex(torsion.quad[idx], step);
        if (x < 0)
            throw std::runtime_error{"Atom " + torsion.quad[idx].name + " was not found"};
        aIx[idx] = x;
    }

    return aIx;
}

static
int makeSelection(mmdb::Manager *mol, const Step &step) {
    int hSel = mol->NewSelection();
    // Just select everything - we expect to get a filtered array of atoms
    for (auto &at : step.atoms)
        mol->SelectAtom(hSel, at);

    return hSel;
}

static
double torsionAngle(mmdb::PAtom a, mmdb::PAtom b, mmdb::PAtom c, mmdb::PAtom d) {
   clipper::Coord_orth p1(a->x, a->y, a->z);
   clipper::Coord_orth p2(b->x, b->y, b->z);
   clipper::Coord_orth p3(c->x, c->y, c->z);
   clipper::Coord_orth p4(d->x, d->y, d->z);

   double tors = clipper::Coord_orth::torsion(p1, p2, p3, p4);
   return clipper::Util::rad2d(tors);
}

static
double measureTorsionAngle(const AtomQuad &quad, const Step &step) {
    auto at1 = getAtom(quad[0], step);
    auto at2 = getAtom(quad[1], step);
    auto at3 = getAtom(quad[2], step);
    auto at4 = getAtom(quad[3], step);

    if (!at1 || !at2 || !at3 || !at4)
        return 0;
    return torsionAngle(at1, at2, at3, at4);
}

static
double bondLength(mmdb::PAtom a, mmdb::PAtom b) {
    auto dx = b->x - a->x;
    auto dy = b->y - a->y;
    auto dz = b->z - a->z;

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

static
double measureBondLength(const Atom &atom1, const Atom &atom2, const Step &step) {
    auto at1 = getAtom(atom1, step);
    auto at2 = getAtom(atom2, step);

    if (!at1 || !at2)
        return 0;
    return bondLength(at1, at2);
}

template <typename T>
clipper::Mat33<T> operator*(const clipper::Mat33<T> &mat, const T &scale) {
    clipper::Mat33<T> scaled = mat;
    for (size_t idx = 0; idx < 3; idx++) for (size_t jdx = 0; jdx < 3; jdx++) scaled(idx, jdx) = scaled(idx, jdx) * scale;

    return scaled;
}

template <typename T>
clipper::Mat33<T> operator*(const T &scale, const clipper::Mat33<T> &mat) {
    return mat * scale;
}

template <typename T>
static
clipper::Mat33<T> rodriguez(const clipper::Vec3<T> &axis, const T &angle) {
    static const clipper::Mat33<T> I{
        1, 0, 0,
        0, 1, 0,
        0, 0, 1
    };

    const clipper::Mat33<T> K{
               0, -axis[2],  axis[1],
         axis[2],        0, -axis[0],
        -axis[1],  axis[0],        0
    };

    return I + (std::sin(angle) * K) + (1 - std::cos(angle)) * K * K;
}

template <typename T>
static
std::ostream & operator<<(std::ostream &os, const clipper::Vec3<T> &vec) {
    os << "[ " << vec[0] << ", " << vec[1] << ", " << vec[2] << " ]";
    return os;
}

template <typename T>
static
T magnitude(const clipper::Vec3<T> &vec) {
    return std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

static
void realign(const clipper::Vec3<double> &C5Pinitial, const clipper::Vec3<double> &O3Pinitial, Step &step) {
    auto C5P = getAtom({ " C5'", Atom::FIRST }, step);
    auto O3P = getAtom({ " O3'", Atom::SECOND }, step);

    clipper::Vec3<double> C5Pfinal{ C5P->x, C5P->y, C5P->z };

    // Move to origin so we can rotate
    for (auto &at : step.atoms) {
        at->x -= C5Pfinal[0];
        at->y -= C5Pfinal[1];
        at->z -= C5Pfinal[2];
    }

    clipper::Vec3<double> O3Pfinal{ O3P->x, O3P->y, O3P->z };
    auto posVecInitial = O3Pinitial - C5Pinitial;
    auto posVecFinal = O3Pfinal; // Same as O3Pfinal because we've translated to origin

    auto magInitial = magnitude(posVecInitial);
    auto magFinal = magnitude(posVecFinal);

    auto axis = clipper::Vec3<double>::cross(posVecFinal, posVecInitial).unit();
    auto cosAngle = (posVecFinal * posVecInitial) / (magInitial * magFinal);
    auto angle = std::acos(cosAngle);

    auto RR = rodriguez(axis, angle);

    clipper::Vec3<double> rv{};
    // With his mighty matrices he shifted everything!
    for (auto &at : step.atoms) {
        rv[0] = at->x; rv[1] = at->y; rv[2] = at->z;
        // Rotate
        rv = RR * rv;
        at->x = rv[0]; at->y = rv[1]; at->z = rv[2];
        // Move back
        at->x += C5Pfinal[0]; at->y += C5Pfinal[1]; at->z += C5Pfinal[2];
        // Align C5' with the original C5' position
        at->x += (C5Pinitial[0] - C5Pfinal[0]); at->y += (C5Pinitial[1] - C5Pfinal[1]); at->z += (C5Pinitial[2] - C5Pfinal[2]);
    }
}

static
void rotateBond(coot::atom_tree_t &tree, mmdb::PPAtom atoms, const IndexQuad &aIx, double angle) {
    /* rotate_about() adds a value to the current torsion of the bond. Therefore, we need to calculate
     * the difference between the actual and target angle to get it right.
     */
    auto actualAngle = torsionAngle(atoms[aIx[0]], atoms[aIx[1]], atoms[aIx[2]], atoms[aIx[3]]);
    auto rotateBy = angle - actualAngle;
    rotateBy += (rotateBy < -180) * 360 - (rotateBy > 180) * 360;

    /* Spin me 'round 'n around */
    tree.rotate_about(aIx[1], aIx[2], rotateBy, false);
}

static
void setTorsion(const Torsion &torsion, coot::atom_tree_t &tree, Step &step) {
    IndexQuad aIx = getTorsionIndices(torsion, step);
    rotateBond(tree, step.atoms.data(), aIx, torsion.angle);
}

static
void reshapeRibose(coot::atom_tree_t &tree, Atom::WhichResidue whichResidue, const NtC &ntc, Step &step) {
    static const std::array<std::tuple<const double NtC::*, size_t>, 3> NUS_FIRST{{
        { &NtC::nu4first, 4 },
        { &NtC::nu0first, 0 },
        { &NtC::nu1first, 1 }
        //&NtC::nu2first,
        //&NtC::nu3first,
    }};

    static const std::array<std::tuple<const double NtC::*, size_t>, 3> NUS_SECOND{{
        { &NtC::nu4second, 4 },
        { &NtC::nu0second, 0 },
        { &NtC::nu1second, 1 }
        //&NtC::nu2second,
        //&NtC::nu3second,
    }};

    std::vector<AtomQuad> quads{};
    for (const auto &names : NU_TORSIONS_ATOM_NAMES)
        quads.push_back({Atom(names[0], whichResidue), Atom(names[1], whichResidue), Atom(names[2], whichResidue), Atom(names[3], whichResidue)});

    for (const auto &nu : whichResidue == Atom::SECOND ? NUS_SECOND : NUS_FIRST) {
        auto ptr = std::get<0>(nu);
        auto idx = std::get<1>(nu);
        setTorsion(Torsion(quads[idx], ntc.*ptr), tree, step);
    }
}

const std::array<std::string, 7> NTC_CLASSES{ "AA", "AB", "BA", "BB", "IC", "OP", "Z" };

const std::vector<NtC> NTCS{{
	{
		"AA00",
		_A(82.08), _A(206.27), _A(287.91), _A(293.46), _A(172.56), _A(54.93), _A(81.86),
		_A(198.67), _A(200.44),
		5.45, 4.767, _A(18.24),
		_A(1.05), _A(334.06), _A(39.44), _A(320.30), _A(24.38), _A(0.13), _A(334.75), _A(39.17), _A(319.98), _A(25.17)
	},
	{
		"AA01",
		_A(81.46), _A(197.06), _A(291.02), _A(149.0), _A(192.01), _A(182.48), _A(85.41),
		_A(204.2), _A(187.82),
		5.275, 4.716, _A(13.67),
		_A(0.88), _A(335.55), _A(37.40), _A(322.20), _A(23.29), _A(5.07), _A(333.18), _A(37.15), _A(324.81), _A(19.03)
	},
	{
		"AA02",
		_A(87.99), _A(202.37), _A(274.15), _A(293.26), _A(160.56), _A(53.88), _A(88.23),
		_A(244.84), _A(245.64),
		5.219, 4.668, _A(17.88),
		_A(342.60), _A(354.51), _A(24.25), _A(324.84), _A(33.23), _A(348.86), _A(348.57), _A(27.82), _A(324.80), _A(29.39)
	},
	{
		"AA03",
		_A(80.06), _A(223.2), _A(263.7), _A(337.16), _A(155.41), _A(27.11), _A(80.33),
		_A(191.9), _A(202.12),
		5.655, 4.946, _A(16.14),
		_A(3.79), _A(333.22), _A(38.38), _A(322.89), _A(21.00), _A(3.02), _A(334.07), _A(37.75), _A(323.04), _A(21.40)
	},
	{
		"AA04",
		_A(80.18), _A(201.93), _A(298.84), _A(259.42), _A(175.12), _A(85.38), _A(79.67),
		_A(199.19), _A(193.06),
		5.387, 4.708, _A(19.19),
		_A(1.10), _A(335.39), _A(37.34), _A(322.32), _A(23.11), _A(358.89), _A(337.02), _A(36.97), _A(321.41), _A(25.01)
	},
	{
		"AA05",
		_A(82.27), _A(213.77), _A(290.56), _A(139.74), _A(225.3), _A(182.49), _A(84.94),
		_A(208.52), _A(184.15),
		4.829, 4.327, _A(56.57),
		_A(3.67), _A(334.50), _A(36.46), _A(324.64), _A(20.03), _A(7.20), _A(332.62), _A(36.10), _A(327.05), _A(16.32)
	},
	{
		"AA06",
		_A(80.11), _A(214.2), _A(256.03), _A(135.45), _A(228.89), _A(180.07), _A(81.55),
		_A(195.99), _A(180.25),
		5.156, 4.585, _A(15.63),
		_A(1.35), _A(335.66), _A(36.78), _A(322.90), _A(22.56), _A(5.64), _A(332.85), _A(37.32), _A(324.88), _A(18.57)
	},
	{
		"AA07",
		_A(82.52), _A(243.34), _A(215.9), _A(295.72), _A(143.69), _A(52.71), _A(82.35),
		_A(203.43), _A(203.11),
		6.815, 6.805, _A(345.11),
		_A(0.92), _A(337.15), _A(34.89), _A(324.65), _A(21.74), _A(2.67), _A(335.52), _A(35.77), _A(324.81), _A(20.54)
	},
	{
		"AA08",
		_A(82.2), _A(233.46), _A(274.78), _A(305.92), _A(152.66), _A(55.37), _A(79.54),
		_A(189.38), _A(196.92),
		5.36, 4.602, _A(23.03),
		_A(5.63), _A(332.57), _A(37.69), _A(324.53), _A(18.84), _A(1.57), _A(334.63), _A(38.20), _A(321.71), _A(23.18)
	},
	{
		"AA09",
		_A(87.07), _A(231.54), _A(272.02), _A(301.66), _A(154.39), _A(52.24), _A(85.04),
		_A(216.66), _A(233.02),
		5.444, 4.705, _A(24.4),
		_A(359.68), _A(340.49), _A(30.71), _A(328.24), _A(20.32), _A(358.70), _A(340.23), _A(31.98), _A(326.31), _A(22.16)
	},
	{
		"AA10",
		_A(78.99), _A(202.77), _A(313.33), _A(210.02), _A(153.65), _A(143.25), _A(80.4),
		_A(198.59), _A(183.2),
		5.351, 4.801, _A(13.23),
		_A(1.70), _A(335.15), _A(37.35), _A(322.63), _A(22.52), _A(3.03), _A(334.85), _A(36.63), _A(324.12), _A(20.71)
	},
	{
		"AA11",
		_A(81.31), _A(260.55), _A(225.94), _A(91.23), _A(256.64), _A(192.33), _A(83.43),
		_A(189.41), _A(178.83),
		5.041, 4.605, _A(10.5),
		_A(3.60), _A(334.47), _A(36.57), _A(324.49), _A(20.13), _A(7.57), _A(332.05), _A(36.78), _A(326.62), _A(16.30)
	},
	{
		"AA12",
		_A(82.37), _A(197.73), _A(260.8), _A(298.31), _A(176.21), _A(49.65), _A(81.0),
		_A(200.62), _A(197.09),
		6.354, 6.245, _A(354.35),
		_A(2.38), _A(335.47), _A(36.15), _A(324.29), _A(21.08), _A(3.10), _A(334.02), _A(37.73), _A(323.06), _A(21.33)
	},
	{
		"AA13",
		_A(82.94), _A(190.01), _A(240.94), _A(306.57), _A(173.66), _A(50.26), _A(82.57),
		_A(199.7), _A(191.92),
		6.925, 7.277, _A(342.95),
		_A(1.49), _A(336.81), _A(34.78), _A(325.10), _A(21.14), _A(3.76), _A(334.43), _A(36.45), _A(324.74), _A(19.93)
	},
	{
		"AAS1",
		_A(78.53), _A(221.26), _A(301.25), _A(291.1), _A(171.11), _A(58.54), _A(80.63),
		_A(15.02), _A(195.93),
		5.294, 4.367, _A(44.82),
		_A(359.93), _A(336.38), _A(37.02), _A(321.96), _A(24.03), _A(3.41), _A(333.96), _A(37.55), _A(323.42), _A(20.95)
	},
	{
		"AB01",
		_A(86.28), _A(186.23), _A(281.16), _A(301.0), _A(178.53), _A(54.45), _A(141.82),
		_A(222.77), _A(255.95),
		5.282, 4.688, _A(17.72),
		_A(340.26), _A(357.78), _A(21.66), _A(326.18), _A(33.65), _A(339.20), _A(33.39), _A(327.23), _A(21.45), _A(359.51)
	},
	{
		"AB02",
		_A(93.78), _A(58.67), _A(55.69), _A(207.73), _A(188.43), _A(65.72), _A(130.51),
		_A(238.83), _A(250.42),
		4.811, 4.353, _A(23.29),
		_A(340.19), _A(3.67), _A(12.48), _A(335.81), _A(27.75), _A(330.81), _A(35.62), _A(331.58), _A(12.28), _A(10.49)
	},
	{
		"AB03",
		_A(103.02), _A(194.7), _A(254.51), _A(320.5), _A(161.6), _A(40.11), _A(136.29),
		_A(228.47), _A(255.99),
		5.131, 4.578, _A(22.53),
		_A(323.23), _A(25.00), _A(355.00), _A(343.82), _A(33.17), _A(332.51), _A(37.21), _A(327.36), _A(17.66), _A(6.02)
	},
	{
		"AB04",
		_A(86.6), _A(214.51), _A(297.16), _A(279.63), _A(208.92), _A(55.32), _A(139.05),
		_A(193.6), _A(233.41),
		6.572, 5.582, _A(7.84),
		_A(2.89), _A(335.40), _A(35.68), _A(325.06), _A(20.25), _A(337.76), _A(32.55), _A(329.84), _A(18.03), _A(2.54)
	},
	{
		"AB05",
		_A(83.04), _A(213.57), _A(287.44), _A(303.05), _A(176.29), _A(58.31), _A(145.45),
		_A(196.15), _A(238.84),
		5.73, 4.82, _A(16.68),
		_A(0.81), _A(335.71), _A(37.08), _A(322.37), _A(23.31), _A(336.93), _A(36.27), _A(325.13), _A(22.55), _A(0.16)
	},
	{
		"AB1S",
		_A(90.42), _A(213.72), _A(280.1), _A(294.93), _A(175.95), _A(55.69), _A(138.73),
		_A(238.66), _A(67.92),
		6.069, 5.433, _A(1.05),
		_A(2.84), _A(346.83), _A(27.78), _A(326.45), _A(3.75), _A(-27.80), _A(38.50), _A(317.95), _A(16.49), _A(3.15)
	},
	{
		"AB2S",
		_A(83.75), _A(220.3), _A(286.55), _A(306.81), _A(170.92), _A(52.55), _A(139.75),
		_A(202.69), _A(69.42),
		5.656, 4.782, _A(10.58),
		_A(357.84), _A(340.05), _A(33.11), _A(324.70), _A(23.69), _A(330.97), _A(38.41), _A(327.52), _A(16.69), _A(7.52)
	},
	{
		"BA01",
		_A(136.06), _A(188.69), _A(254.87), _A(299.54), _A(161.48), _A(53.28), _A(88.05),
		_A(253.93), _A(225.3),
		4.681, 4.161, _A(29.42),
		_A(336.78), _A(33.42), _A(329.35), _A(17.70), _A(3.44), _A(335.09), _A(4.01), _A(16.81), _A(327.97), _A(35.81)
	},
	{
		"BA05",
		_A(131.4), _A(184.36), _A(268.77), _A(295.91), _A(168.82), _A(52.31), _A(104.24),
		_A(250.85), _A(235.4),
		4.643, 4.106, _A(24.17),
		_A(334.32), _A(32.68), _A(332.65), _A(13.05), _A(7.85), _A(320.35), _A(28.61), _A(352.24), _A(344.87), _A(34.24)
	},
	{
		"BA08",
		_A(138.8), _A(208.0), _A(212.9), _A(300.69), _A(141.44), _A(49.0), _A(88.86),
		_A(263.31), _A(214.82),
		4.912, 4.433, _A(33.23),
		_A(326.58), _A(42.89), _A(324.39), _A(17.40), _A(9.85), _A(344.66), _A(354.18), _A(23.04), _A(327.64), _A(30.08)
	},
	{
		"BA09",
		_A(134.09), _A(199.88), _A(286.87), _A(256.41), _A(68.19), _A(171.69), _A(90.26),
		_A(264.52), _A(186.05),
		3.954, 3.807, _A(22.73),
		_A(323.59), _A(42.92), _A(327.15), _A(13.01), _A(14.43), _A(351.61), _A(349.17), _A(24.63), _A(329.70), _A(24.17)
	},
	{
		"BA10",
		_A(136.22), _A(200.19), _A(235.8), _A(95.17), _A(218.75), _A(205.06), _A(89.61),
		_A(253.85), _A(200.48),
		4.426, 3.988, _A(28.2),
		_A(333.28), _A(35.58), _A(329.32), _A(16.06), _A(6.52), _A(357.36), _A(341.69), _A(30.82), _A(326.83), _A(22.75)
	},
	{
		"BA13",
		_A(141.38), _A(219.89), _A(200.34), _A(81.88), _A(230.89), _A(196.27), _A(87.78),
		_A(264.73), _A(197.4),
		4.601, 4.154, _A(30.96),
		_A(327.56), _A(43.07), _A(323.08), _A(19.32), _A(8.02), _A(3.99), _A(335.24), _A(35.10), _A(326.17), _A(18.80)
	},
	{
		"BA16",
		_A(146.48), _A(245.83), _A(189.59), _A(61.4), _A(228.56), _A(198.59), _A(84.77),
		_A(265.92), _A(199.15),
		4.815, 4.436, _A(28.76),
		_A(328.74), _A(43.50), _A(321.50), _A(21.71), _A(5.77), _A(1.13), _A(337.28), _A(34.56), _A(325.14), _A(21.24)
	},
	{
		"BA17",
		_A(149.28), _A(253.24), _A(176.54), _A(294.9), _A(130.74), _A(44.22), _A(97.74),
		_A(271.16), _A(232.57),
		4.894, 4.548, _A(38.19),
		_A(329.50), _A(44.04), _A(319.46), _A(24.55), _A(3.45), _A(331.67), _A(13.74), _A(4.93), _A(338.48), _A(31.22)
	},
	{
		"BB00",
		_A(137.82), _A(183.11), _A(258.18), _A(303.72), _A(179.59), _A(44.23), _A(138.14),
		_A(252.55), _A(258.25),
		4.947, 4.373, _A(25.54),
		_A(337.49), _A(33.24), _A(329.00), _A(18.59), _A(2.33), _A(333.12), _A(36.96), _A(327.35), _A(17.88), _A(5.47)
	},
	{
		"BB01",
		_A(130.75), _A(180.95), _A(265.55), _A(301.4), _A(176.22), _A(48.62), _A(120.15),
		_A(247.52), _A(243.88),
		4.853, 4.303, _A(25.57),
		_A(334.91), _A(32.11), _A(332.91), _A(13.00), _A(7.48), _A(326.56), _A(33.47), _A(338.63), _A(2.29), _A(19.47)
	},
	{
		"BB02",
		_A(140.58), _A(193.95), _A(246.44), _A(30.93), _A(195.32), _A(297.1), _A(150.05),
		_A(251.88), _A(253.34),
		5.11, 4.458, _A(21.84),
		_A(338.29), _A(34.05), _A(327.09), _A(21.06), _A(0.26), _A(353.17), _A(23.26), _A(330.28), _A(26.41), _A(347.51)
	},
	{
		"BB03",
		_A(145.07), _A(175.18), _A(274.45), _A(162.78), _A(165.53), _A(174.65), _A(146.24),
		_A(241.19), _A(232.96),
		5.211, 4.512, _A(26.72),
		_A(350.40), _A(24.73), _A(330.38), _A(24.61), _A(350.40), _A(345.79), _A(29.12), _A(327.78), _A(24.65), _A(353.29)
	},
	{
		"BB04",
		_A(140.09), _A(201.29), _A(214.17), _A(314.75), _A(152.68), _A(46.11), _A(139.98),
		_A(262.54), _A(252.56),
		5.121, 4.648, _A(28.58),
		_A(332.63), _A(38.07), _A(326.17), _A(18.88), _A(5.14), _A(336.69), _A(34.91), _A(327.24), _A(19.97), _A(1.93)
	},
	{
		"BB05",
		_A(141.92), _A(219.93), _A(197.39), _A(76.4), _A(233.15), _A(213.25), _A(129.04),
		_A(266.1), _A(208.45),
		4.685, 4.189, _A(28.79),
		_A(328.96), _A(41.72), _A(324.03), _A(19.26), _A(0.24), _A(354.09), _A(22.48), _A(341.35), _A(12.02), _A(4.02)
	},
	{
		"BB07",
		_A(143.72), _A(247.34), _A(169.46), _A(296.5), _A(140.91), _A(46.12), _A(141.14),
		_A(270.68), _A(260.38),
		5.202, 4.955, _A(46.01),
		_A(329.13), _A(42.16), _A(322.92), _A(20.39), _A(6.34), _A(340.90), _A(32.30), _A(327.10), _A(22.46), _A(357.74)
	},
	{
		"BB08",
		_A(147.04), _A(248.93), _A(180.96), _A(66.3), _A(225.1), _A(208.51), _A(148.29),
		_A(270.22), _A(234.89),
		5.065, 4.691, _A(29.76),
		_A(335.23), _A(38.73), _A(322.76), _A(23.99), _A(0.29), _A(343.37), _A(31.97), _A(325.72), _A(25.44), _A(354.31)
	},
	{
		"BB10",
		_A(138.02), _A(195.81), _A(191.53), _A(21.95), _A(106.39), _A(18.94), _A(129.25),
		_A(257.24), _A(257.78),
		4.874, 4.389, _A(25.73),
		_A(329.72), _A(40.97), _A(324.24), _A(19.37), _A(6.59), _A(326.61), _A(39.97), _A(328.76), _A(12.95), _A(12.59)
	},
	{
		"BB11",
		_A(145.33), _A(199.35), _A(200.27), _A(122.97), _A(226.75), _A(187.37), _A(143.91),
		_A(256.33), _A(222.6),
		5.138, 4.648, _A(29.49),
		_A(338.40), _A(35.49), _A(324.76), _A(23.52), _A(358.61), _A(342.46), _A(31.00), _A(327.79), _A(22.58), _A(356.72)
	},
	{
		"BB12",
		_A(140.27), _A(195.62), _A(279.63), _A(256.87), _A(76.46), _A(171.35), _A(139.93),
		_A(268.93), _A(204.76),
		4.15, 3.866, _A(24.07),
		_A(332.60), _A(38.03), _A(325.91), _A(19.13), _A(4.96), _A(336.70), _A(34.52), _A(327.58), _A(19.73), _A(2.04)
	},
	{
		"BB13",
		_A(142.55), _A(187.45), _A(293.09), _A(219.28), _A(98.07), _A(161.21), _A(145.67),
		_A(253.0), _A(218.57),
		4.86, 4.359, _A(24.71),
		_A(341.60), _A(32.00), _A(327.12), _A(22.86), _A(357.03), _A(348.06), _A(27.69), _A(327.99), _A(25.73), _A(351.15)
	},
	{
		"BB14",
		_A(109.7), _A(104.08), _A(305.19), _A(219.61), _A(255.25), _A(82.61), _A(132.81),
		_A(258.66), _A(264.88),
		4.854, 4.332, _A(28.57),
		_A(-17.10), _A(28.26), _A(345.02), _A(357.41), _A(13.50), _A(-29.60), _A(38.86), _A(328.14), _A(14.97), _A(5.90)
	},
	{
		"BB15",
		_A(143.55), _A(189.14), _A(256.68), _A(344.79), _A(188.83), _A(350.19), _A(147.56),
		_A(249.9), _A(262.4),
		4.999, 4.374, _A(22.35),
		_A(342.52), _A(31.66), _A(326.91), _A(23.74), _A(356.02), _A(346.81), _A(29.03), _A(327.10), _A(25.96), _A(351.83)
	},
	{
		"BB16",
		_A(137.56), _A(220.68), _A(282.24), _A(284.37), _A(172.53), _A(47.58), _A(139.88),
		_A(204.18), _A(270.1),
		5.459, 4.897, _A(32.71),
		_A(343.98), _A(26.92), _A(332.90), _A(18.39), _A(358.35), _A(332.53), _A(38.31), _A(325.71), _A(19.34), _A(4.93)
	},
	{
		"BB17",
		_A(128.61), _A(144.73), _A(275.34), _A(230.08), _A(241.47), _A(79.29), _A(135.55),
		_A(245.49), _A(269.64),
		5.242, 4.609, _A(24.09),
		_A(334.13), _A(33.85), _A(331.45), _A(14.32), _A(7.06), _A(329.48), _A(39.86), _A(326.31), _A(17.22), _A(8.12)
	},
	{
		"BB1S",
		_A(139.29), _A(198.98), _A(281.0), _A(306.17), _A(258.05), _A(307.67), _A(151.23),
		_A(237.38), _A(66.34),
		6.588, 6.507, _A(357.41),
		_A(320.37), _A(43.67), _A(322.09), _A(20.52), _A(1.32), _A(354.78), _A(30.61), _A(323.97), _A(29.50), _A(343.77)
	},
	{
		"BB20",
		_A(142.6), _A(294.29), _A(109.51), _A(150.07), _A(198.78), _A(53.73), _A(151.67),
		_A(261.4), _A(185.24),
		5.233, 5.478, _A(276.62),
		_A(330.16), _A(41.48), _A(323.47), _A(20.45), _A(5.66), _A(352.34), _A(25.72), _A(327.18), _A(29.11), _A(346.35)
	},
	{
		"BB2S",
		_A(134.28), _A(194.19), _A(223.47), _A(46.11), _A(180.4), _A(289.97), _A(145.52),
		_A(251.23), _A(66.29),
		6.214, 6.003, _A(358.38),
		_A(326.62), _A(41.55), _A(326.19), _A(15.37), _A(11.16), _A(341.54), _A(34.00), _A(324.18), _A(26.01), _A(355.12)
	},
	{
		"BBS1",
		_A(146.33), _A(187.34), _A(274.44), _A(296.19), _A(171.66), _A(51.54), _A(134.92),
		_A(63.71), _A(259.73),
		4.549, 3.937, _A(45.02),
		_A(341.70), _A(33.72), _A(324.59), _A(25.70), _A(355.22), _A(329.28), _A(39.18), _A(327.48), _A(15.65), _A(9.28)
	},
	{
		"IC01",
		_A(82.76), _A(219.69), _A(289.59), _A(297.18), _A(222.82), _A(54.38), _A(145.18),
		_A(202.97), _A(283.4),
		7.535, 7.245, _A(10.6),
		_A(0.56), _A(336.77), _A(35.73), _A(323.61), _A(22.60), _A(335.28), _A(37.42), _A(324.81), _A(21.95), _A(1.54)
	},
	{
		"IC02",
		_A(81.89), _A(222.38), _A(279.44), _A(298.96), _A(224.78), _A(52.41), _A(143.3),
		_A(201.43), _A(243.53),
		7.847, 7.659, _A(1.88),
		_A(0.93), _A(336.62), _A(35.58), _A(323.94), _A(22.19), _A(334.99), _A(36.97), _A(325.81), _A(20.77), _A(2.47)
	},
	{
		"IC03",
		_A(81.46), _A(238.9), _A(256.55), _A(69.41), _A(178.54), _A(301.33), _A(146.02),
		_A(201.19), _A(265.54),
		7.76, 7.5, _A(359.93),
		_A(359.20), _A(337.92), _A(35.22), _A(323.32), _A(23.65), _A(339.32), _A(34.13), _A(326.40), _A(22.68), _A(358.56)
	},
	{
		"IC04",
		_A(84.76), _A(206.16), _A(286.84), _A(194.08), _A(180.92), _A(183.42), _A(149.09),
		_A(211.61), _A(244.45),
		7.199, 6.967, _A(24.68),
		_A(1.70), _A(337.83), _A(32.95), _A(327.16), _A(19.65), _A(341.24), _A(34.37), _A(324.01), _A(25.94), _A(355.37)
	},
	{
		"IC05",
		_A(140.84), _A(255.23), _A(173.84), _A(288.7), _A(178.43), _A(49.59), _A(146.97),
		_A(273.47), _A(275.37),
		7.114, 7.261, _A(45.6),
		_A(320.59), _A(42.01), _A(326.50), _A(16.50), _A(10.15), _A(349.91), _A(30.13), _A(324.67), _A(28.39), _A(350.43)
	},
	{
		"IC06",
		_A(136.51), _A(236.09), _A(280.3), _A(287.83), _A(174.16), _A(47.06), _A(141.78),
		_A(207.0), _A(269.42),
		6.232, 6.211, _A(39.52),
		_A(338.68), _A(30.45), _A(331.62), _A(16.82), _A(2.97), _A(326.34), _A(43.76), _A(322.40), _A(19.67), _A(8.60)
	},
	{
		"IC07",
		_A(84.29), _A(208.65), _A(291.48), _A(175.85), _A(127.02), _A(177.24), _A(82.99),
		_A(218.08), _A(203.57),
		7.354, 7.517, _A(333.17),
		_A(1.54), _A(335.01), _A(37.31), _A(322.49), _A(22.85), _A(0.35), _A(335.82), _A(37.24), _A(321.94), _A(23.99)
	},
	{
		"OP01",
		_A(82.63), _A(220.94), _A(122.0), _A(279.24), _A(144.5), _A(44.78), _A(81.73),
		_A(204.65), _A(193.58),
		7.39, 8.281, _A(195.33),
		_A(1.03), _A(335.79), _A(36.70), _A(322.87), _A(22.87), _A(1.56), _A(334.82), _A(37.75), _A(322.16), _A(22.96)
	},
	{
		"OP02",
		_A(83.1), _A(226.42), _A(156.11), _A(291.54), _A(159.28), _A(42.91), _A(86.33),
		_A(206.92), _A(176.35),
		8.25, 8.979, _A(259.49),
		_A(0.37), _A(337.37), _A(34.98), _A(324.20), _A(22.34), _A(8.08), _A(332.75), _A(35.01), _A(328.64), _A(14.73)
	},
	{
		"OP03",
		_A(77.75), _A(226.17), _A(300.47), _A(174.79), _A(137.98), _A(50.93), _A(83.67),
		_A(199.33), _A(194.51),
		7.458, 7.521, _A(202.12),
		_A(0.83), _A(334.88), _A(38.53), _A(320.91), _A(24.15), _A(0.97), _A(337.42), _A(34.33), _A(325.28), _A(21.34)
	},
	{
		"OP04",
		_A(80.11), _A(220.84), _A(284.73), _A(165.82), _A(171.42), _A(52.29), _A(84.85),
		_A(205.88), _A(194.39),
		7.147, 7.243, _A(205.98),
		_A(0.39), _A(335.60), _A(37.77), _A(321.48), _A(24.07), _A(3.35), _A(335.20), _A(35.62), _A(325.33), _A(19.79)
	},
	{
		"OP05",
		_A(78.42), _A(204.45), _A(50.33), _A(69.01), _A(126.12), _A(45.39), _A(84.4),
		_A(201.0), _A(186.67),
		8.922, 9.111, _A(263.8),
		_A(358.87), _A(336.58), _A(37.71), _A(320.72), _A(25.55), _A(1.66), _A(336.95), _A(34.37), _A(325.67), _A(20.67)
	},
	{
		"OP06",
		_A(81.62), _A(212.56), _A(142.53), _A(300.88), _A(172.05), _A(47.71), _A(80.7),
		_A(204.28), _A(201.06),
		8.983, 10.121, _A(243.35),
		_A(357.76), _A(339.49), _A(34.12), _A(323.64), _A(24.36), _A(1.29), _A(336.11), _A(36.14), _A(323.63), _A(22.14)
	},
	{
		"OP07",
		_A(81.99), _A(247.41), _A(193.84), _A(292.75), _A(149.02), _A(43.7), _A(80.96),
		_A(196.47), _A(187.5),
		7.506, 7.758, _A(323.87),
		_A(1.38), _A(336.34), _A(35.66), _A(324.13), _A(21.73), _A(5.20), _A(332.72), _A(37.87), _A(324.14), _A(19.36)
	},
	{
		"OP08",
		_A(79.86), _A(202.8), _A(278.28), _A(252.48), _A(83.0), _A(167.32), _A(84.08),
		_A(201.0), _A(177.04),
		5.941, 7.19, _A(317.94),
		_A(5.61), _A(331.08), _A(40.15), _A(322.14), _A(20.33), _A(9.64), _A(330.19), _A(37.69), _A(326.88), _A(14.86)
	},
	{
		"OP09",
		_A(81.92), _A(197.88), _A(268.7), _A(202.92), _A(146.65), _A(51.38), _A(147.88),
		_A(203.67), _A(249.91),
		6.722, 7.817, _A(162.57),
		_A(1.72), _A(336.35), _A(35.39), _A(324.64), _A(21.26), _A(338.04), _A(36.14), _A(324.18), _A(24.08), _A(358.42)
	},
	{
		"OP10",
		_A(146.66), _A(217.66), _A(152.28), _A(288.77), _A(174.5), _A(39.62), _A(87.05),
		_A(247.5), _A(185.99),
		7.344, 7.753, _A(36.01),
		_A(333.81), _A(39.15), _A(323.70), _A(22.15), _A(2.35), _A(10.95), _A(330.65), _A(35.70), _A(329.58), _A(12.34)
	},
	{
		"OP11",
		_A(147.03), _A(267.02), _A(300.98), _A(297.43), _A(187.0), _A(57.0), _A(83.95),
		_A(234.53), _A(197.25),
		7.801, 8.882, _A(110.16),
		_A(337.34), _A(36.66), _A(324.03), _A(23.92), _A(359.02), _A(4.95), _A(332.42), _A(38.49), _A(323.37), _A(20.00)
	},
	{
		"OP12",
		_A(141.26), _A(257.49), _A(285.43), _A(275.32), _A(183.8), _A(42.74), _A(82.74),
		_A(237.46), _A(191.26),
		7.268, 7.506, _A(72.63),
		_A(336.90), _A(34.68), _A(327.50), _A(20.15), _A(1.67), _A(3.77), _A(334.00), _A(37.09), _A(324.07), _A(20.29)
	},
	{
		"OP13",
		_A(146.31), _A(267.65), _A(248.24), _A(62.01), _A(153.45), _A(46.84), _A(85.84),
		_A(239.3), _A(186.39),
		8.914, 10.499, _A(188.5),
		_A(336.70), _A(36.82), _A(324.43), _A(23.11), _A(359.91), _A(4.58), _A(335.52), _A(33.94), _A(327.73), _A(17.51)
	},
	{
		"OP14",
		_A(147.74), _A(269.5), _A(227.34), _A(57.88), _A(196.25), _A(60.7), _A(88.21),
		_A(259.3), _A(175.99),
		8.159, 9.562, _A(229.29),
		_A(334.90), _A(38.52), _A(323.44), _A(23.07), _A(1.07), _A(6.16), _A(335.20), _A(32.87), _A(329.67), _A(15.30)
	},
	{
		"OP15",
		_A(148.61), _A(200.52), _A(151.33), _A(292.37), _A(148.55), _A(41.24), _A(85.51),
		_A(264.6), _A(188.23),
		6.611, 6.91, _A(22.31),
		_A(336.06), _A(38.86), _A(321.90), _A(25.31), _A(358.93), _A(2.75), _A(336.51), _A(34.05), _A(326.59), _A(19.36)
	},
	{
		"OP16",
		_A(147.47), _A(268.22), _A(149.69), _A(308.7), _A(177.51), _A(47.48), _A(79.98),
		_A(229.22), _A(194.72),
		7.409, 8.859, _A(100.8),
		_A(337.87), _A(36.20), _A(324.30), _A(23.89), _A(358.66), _A(2.41), _A(334.40), _A(37.82), _A(322.60), _A(22.09)
	},
	{
		"OP17",
		_A(145.2), _A(266.61), _A(293.74), _A(290.97), _A(137.7), _A(177.42), _A(83.81),
		_A(233.91), _A(194.89),
		7.875, 10.107, _A(190.64),
		_A(338.28), _A(34.78), _A(326.06), _A(22.20), _A(359.48), _A(358.25), _A(340.00), _A(32.93), _A(325.21), _A(23.14)
	},
	{
		"OP18",
		_A(149.01), _A(290.82), _A(105.84), _A(65.69), _A(198.39), _A(54.55), _A(147.22),
		_A(224.63), _A(242.31),
		6.536, 8.591, _A(230.07),
		_A(339.51), _A(35.53), _A(323.87), _A(25.27), _A(356.79), _A(337.27), _A(36.51), _A(324.35), _A(23.54), _A(359.27)
	},
	{
		"OP19",
		_A(144.89), _A(225.37), _A(64.28), _A(73.97), _A(186.21), _A(188.23), _A(125.56),
		_A(249.73), _A(254.59),
		7.348, 8.031, _A(20.38),
		_A(340.05), _A(32.49), _A(327.67), _A(21.64), _A(358.80), _A(327.43), _A(35.32), _A(335.27), _A(6.54), _A(16.19)
	},
	{
		"OP1S",
		_A(143.65), _A(206.64), _A(59.95), _A(82.45), _A(203.34), _A(191.33), _A(146.7),
		_A(242.48), _A(68.2),
		6.992, 7.55, _A(65.02),
		_A(338.91), _A(34.19), _A(326.14), _A(22.53), _A(358.96), _A(341.96), _A(34.19), _A(323.63), _A(26.87), _A(354.35)
	},
	{
		"OP20",
		_A(140.07), _A(271.03), _A(282.88), _A(297.58), _A(190.48), _A(55.65), _A(148.09),
		_A(259.78), _A(217.13),
		7.836, 8.562, _A(66.96),
		_A(334.26), _A(35.50), _A(328.48), _A(17.53), _A(4.94), _A(341.12), _A(34.46), _A(323.95), _A(26.11), _A(355.31)
	},
	{
		"OP21",
		_A(149.38), _A(242.21), _A(79.53), _A(66.73), _A(177.49), _A(62.49), _A(143.94),
		_A(228.26), _A(243.91),
		8.237, 10.261, _A(244.6),
		_A(340.38), _A(35.44), _A(323.04), _A(26.60), _A(355.46), _A(339.77), _A(33.77), _A(326.12), _A(23.05), _A(358.08)
	},
	{
		"OP22",
		_A(147.16), _A(244.71), _A(125.59), _A(285.88), _A(161.54), _A(48.75), _A(141.95),
		_A(233.13), _A(230.89),
		6.437, 6.926, _A(44.33),
		_A(339.49), _A(34.97), _A(324.44), _A(24.74), _A(357.18), _A(337.02), _A(34.61), _A(327.51), _A(20.03), _A(1.74)
	},
	{
		"OP23",
		_A(147.44), _A(260.12), _A(167.8), _A(271.73), _A(80.49), _A(175.19), _A(148.9),
		_A(220.63), _A(258.53),
		5.49, 6.978, _A(60.24),
		_A(334.98), _A(38.76), _A(323.08), _A(23.48), _A(0.78), _A(331.67), _A(43.09), _A(319.40), _A(25.28), _A(1.68)
	},
	{
		"OP24",
		_A(147.69), _A(284.48), _A(96.35), _A(81.93), _A(247.78), _A(189.67), _A(84.52),
		_A(242.06), _A(185.03),
		3.975, 5.212, _A(91.01),
		_A(336.17), _A(37.88), _A(323.35), _A(23.89), _A(359.76), _A(8.88), _A(331.23), _A(36.78), _A(327.37), _A(15.02)
	},
	{
		"OP25",
		_A(145.13), _A(260.59), _A(177.94), _A(94.65), _A(208.55), _A(63.44), _A(85.98),
		_A(283.16), _A(191.82),
		7.209, 8.071, _A(259.89),
		_A(335.61), _A(36.81), _A(325.31), _A(21.55), _A(1.62), _A(2.34), _A(336.06), _A(35.05), _A(325.35), _A(20.44)
	},
	{
		"OP26",
		_A(144.93), _A(265.48), _A(161.05), _A(210.96), _A(166.13), _A(49.76), _A(81.15),
		_A(239.99), _A(190.87),
		6.921, 6.959, _A(15.71),
		_A(339.10), _A(33.80), _A(326.87), _A(21.99), _A(359.16), _A(4.45), _A(333.22), _A(37.74), _A(323.85), _A(20.06)
	},
	{
		"OP27",
		_A(156.85), _A(260.09), _A(70.8), _A(92.54), _A(162.04), _A(170.23), _A(82.06),
		_A(211.25), _A(195.38),
		7.461, 8.093, _A(3.54),
		_A(342.86), _A(35.52), _A(320.73), _A(30.46), _A(351.43), _A(356.82), _A(341.03), _A(32.61), _A(324.63), _A(24.38)
	},
	{
		"OP28",
		_A(82.47), _A(225.5), _A(166.19), _A(291.74), _A(159.42), _A(43.17), _A(144.63),
		_A(196.24), _A(232.81),
		7.922, 8.573, _A(235.47),
		_A(0.35), _A(336.60), _A(36.34), _A(322.94), _A(23.17), _A(336.66), _A(36.09), _A(325.54), _A(21.92), _A(0.68)
	},
	{
		"OP29",
		_A(81.91), _A(242.81), _A(253.94), _A(71.62), _A(185.71), _A(58.91), _A(85.58),
		_A(206.56), _A(185.19),
		8.566, 9.564, _A(143.71),
		_A(359.58), _A(337.44), _A(35.59), _A(323.21), _A(23.53), _A(5.16), _A(334.46), _A(35.27), _A(326.81), _A(17.73)
	},
	{
		"OP30",
		_A(82.12), _A(244.32), _A(195.72), _A(75.55), _A(171.43), _A(49.37), _A(89.46),
		_A(196.97), _A(181.07),
		9.675, 10.477, _A(80.55),
		_A(0.83), _A(336.99), _A(35.26), _A(324.30), _A(21.98), _A(6.81), _A(335.45), _A(32.04), _A(330.99), _A(14.11)
	},
	{
		"OP31",
		_A(82.7), _A(214.14), _A(63.86), _A(67.06), _A(103.95), _A(185.23), _A(81.67),
		_A(210.37), _A(186.12),
		9.301, 9.618, _A(55.51),
		_A(1.06), _A(336.30), _A(36.11), _A(323.53), _A(22.34), _A(4.11), _A(333.50), _A(37.71), _A(323.68), _A(20.32)
	},
	{
		"OPS1",
		_A(146.34), _A(263.27), _A(288.78), _A(283.11), _A(184.02), _A(53.06), _A(82.36),
		_A(62.73), _A(192.18),
		7.504, 8.21, _A(93.81),
		_A(335.86), _A(37.73), _A(323.91), _A(23.16), _A(0.38), _A(4.44), _A(333.79), _A(36.85), _A(324.72), _A(19.44)
	},
	{
		"ZZ01",
		_A(81.45), _A(210.23), _A(48.84), _A(165.87), _A(149.51), _A(48.91), _A(147.29),
		_A(207.63), _A(225.61),
		6.019, 4.612, _A(305.99),
		_A(0.95), _A(336.58), _A(35.73), _A(323.84), _A(22.25), _A(337.14), _A(36.89), _A(323.92), _A(23.86), _A(359.14)
	},
	{
		"ZZ02",
		_A(144.08), _A(269.37), _A(78.27), _A(228.54), _A(174.6), _A(55.14), _A(86.37),
		_A(232.14), _A(276.86),
		6.198, 6.217, _A(60.81),
		_A(337.01), _A(35.07), _A(326.82), _A(20.89), _A(1.16), _A(4.47), _A(336.08), _A(33.08), _A(328.56), _A(17.14)
	},
	{
		"ZZ1S",
		_A(147.07), _A(262.7), _A(76.1), _A(66.26), _A(185.95), _A(178.35), _A(95.62),
		_A(206.57), _A(60.97),
		6.264, 6.389, _A(358.92),
		_A(333.66), _A(37.89), _A(325.31), _A(20.92), _A(3.43), _A(355.79), _A(347.27), _A(23.52), _A(333.38), _A(19.10)
	},
	{
		"ZZ2S",
		_A(141.24), _A(262.82), _A(71.06), _A(77.75), _A(179.9), _A(184.79), _A(146.77),
		_A(208.12), _A(76.58),
		6.238, 6.366, _A(355.87),
		_A(332.34), _A(37.79), _A(326.62), _A(18.26), _A(5.81), _A(341.71), _A(32.28), _A(326.54), _A(23.79), _A(356.56)
	},
	{
		"ZZS1",
		_A(97.14), _A(243.46), _A(292.8), _A(210.17), _A(230.97), _A(55.64), _A(144.06),
		_A(63.21), _A(205.56),
		6.782, 5.663, _A(331.88),
		_A(354.12), _A(349.38), _A(21.64), _A(334.38), _A(19.73), _A(332.17), _A(38.12), _A(326.10), _A(19.17), _A(5.35)
	},
	{
		"ZZS2",
		_A(94.87), _A(186.65), _A(64.06), _A(168.69), _A(161.84), _A(44.04), _A(142.85),
		_A(56.24), _A(212.63),
		6.737, 5.675, _A(329.45),
		_A(355.39), _A(347.50), _A(23.73), _A(333.14), _A(19.73), _A(332.62), _A(37.33), _A(326.74), _A(18.73), _A(5.47)
	}
}};

Atom::Atom() :
    name{""},
    whichResidue{FIRST}
{}

Atom::Atom(std::string name, WhichResidue whichResidue) :
    name{std::move(name)},
    whichResidue{whichResidue}
{}

const std::vector<Torsion> & NtC::backboneTorsions() const {
    return m_backboneTorsions;
}

AtomQuad NtC::chiAtomQuad(const mmdb::Residue *residue, const Atom::WhichResidue whichResidue) {
    switch (getRingType(trim(residue->name))) {
    case PURINE:
        return BaseTorsion(PURINE_ATOM_NAMES, whichResidue);
    case PYRIMIDINE:
        return BaseTorsion(PYRIMIDINE_ATOM_NAMES, whichResidue);
    }
}

Torsion NtC::chi1Torsion(const mmdb::Residue *residue) const {
    return Torsion(chiAtomQuad(residue, Atom::FIRST), chi_1);
}

Torsion NtC::chi2Torsion(const mmdb::Residue *residue) const {
    return Torsion(chiAtomQuad(residue, Atom::SECOND), chi_2);
}

Torsion::Torsion() :
    quad{{ Atom(), Atom(), Atom(), Atom() }},
    angle(0)
{}

Torsion::Torsion(AtomQuad quad, double angle) :
    quad{std::move(quad)},
    angle{angle}
{}

NtC::NtC(
    std::string name,
    double delta_1, double epsilon_1, double zeta_1,
    double alpha_2, double beta_2, double gamma_2, double delta_2,
    double chi_1, double chi_2,
    double cc, double nn, double mu,
    double nu0first, double nu1first, double nu2first, double nu3first, double nu4first,
    double nu0second, double nu1second, double nu2second, double nu3second, double nu4second
) :
    name{std::move(name)},
    delta_1{delta_1}, epsilon_1{epsilon_1}, zeta_1{zeta_1},
    alpha_2{alpha_2}, beta_2{beta_2}, gamma_2{gamma_2}, delta_2{delta_2},
    chi_1{chi_1}, chi_2{chi_2},
    cc{cc}, nn{nn}, mu{mu},
    nu0first{nu0first}, nu1first{nu1first}, nu2first{nu2first}, nu3first{nu3first}, nu4first{nu4first},
    nu0second{nu0second}, nu1second{nu1second}, nu2second{nu2second}, nu3second{nu3second}, nu4second{nu4second},
    NtCClass{nameToClass(this->name)},
    m_backboneTorsions(BackboneTorsions(*this))
{}


NtC::RingType NtC::getRingType(const std::string &compound) {
    if (compound == "A" || compound == "DA" || compound == "G" || compound == "DG")
        return PURINE;
    else if (compound == "C" || compound == "DC" || compound == "DT" || compound == "U" || compound == "T")
        return PYRIMIDINE;

    assert(false);
}

std::string NtC::nameToClass(const std::string &name) {
    auto prov = name.substr(0, 2);
    return prov == "ZZ" ? "Z" : prov;
}

void apply_NtC(mmdb::Manager *mol, coot::protein_geometry *geom, const NtC &ntc) {
    /*
     *
     * My dear friend from the future. I'm sorry for whatever happened in your life that has lead you
     * to this code. Here are a few notes to help you suffer through this.
     *
     */

    /* We expect to get the two residues that make up the NtC step and nothing else.
     * If the called gave us something else, the ball is on them.
     * Note that the step shall contain only the altpositions we are interested in
     */
    mmdb::Chain *chain = mol->GetModel(1)->GetChain(0);
    mmdb::Residue *firstResidue = chain->GetResidue(0);
    mmdb::Residue *secondResidue = chain->GetResidue(1);

    /* COMMENT THIS */
    mmdb::PPAtom allAtoms = nullptr;
    int nAllAtoms = 0;
    mol->GetAtomTable(allAtoms, nAllAtoms);
    Step step{filterAtoms(allAtoms, nAllAtoms, firstResidue->GetSeqNum()), firstResidue, secondResidue};
    int hSel = makeSelection(mol, step);

    auto C5Pinitial = getC5PrimeCoords(firstResidue);
    auto O3Pinitial = getO3PrimeCoords(secondResidue);

    mmdb::PContact contacts = nullptr;
    int nContacts = 0;
    mol->SeekContacts(
        step.atoms.data(), step.atoms.size(),
        step.atoms.data(), step.atoms.size(),
        0, 4,
        0,
        contacts, nContacts
    );

    coot::contact_info ci(step.atoms.data(), contacts, nContacts);
    auto contactIndices = ci.get_contact_indices_with_reverse_contacts();

    std::vector<std::tuple<IndexQuad, double>> backboneRotations{};
    for (const auto &tor : ntc.backboneTorsions()) {
        IndexQuad aIx = getTorsionIndices(tor, step);
        backboneRotations.emplace_back(std::move(aIx), tor.angle);
    }

    coot::atom_tree_t tree(contactIndices, 0, mol, hSel);

    /* Now rotate the backbone.
     * We cannot use set_dihedral() series of functions because they use atom names as input.
     * We can use rotate_about() which takes atom indices of the two middle atoms (the bond) to rotate.
     */
    for (const auto &rot : backboneRotations) {
        auto aIx = std::get<0>(rot);
        auto angle = std::get<1>(rot);

        rotateBond(tree, step.atoms.data(), aIx, angle);
    }

    /* We're done with the backbone.
     * Now reshape the (deoxy)ribose so that the N-C bonds point in the right direction.
     * NOTE: Trying to do this with the explicit NtC parameters - NN, CC and mu,
     * does not work because these values describe the effect, not the cause.
     * NOTE2: We can currently transition only between NtCs whose ribose
     * shape is reasonably close to the source structure. More drastic
     * transition would require us to reshape the ribose and the C1'-N bond
     * to make it match a template corresponding to the shape of the
     * ribose of the given NtC
     */
    reshapeRibose(tree, Atom::FIRST, ntc, step);
    reshapeRibose(tree, Atom::SECOND, ntc, step);

    /* Okay, so this was the backbone. Now rotate the bases. */
    auto chi1 = ntc.chi1Torsion(firstResidue);
    IndexQuad aIxChi1 = getTorsionIndices(chi1, step);
    rotateBond(tree, step.atoms.data(), aIxChi1, chi1.angle);

    auto chi2 = ntc.chi2Torsion(secondResidue);
    IndexQuad aIxChi2 = getTorsionIndices(chi2, step);
    rotateBond(tree, step.atoms.data(), aIxChi2, chi2.angle);

    std::vector<IndexQuad> aIxs;
    std::transform(backboneRotations.cbegin(), backboneRotations.cend(), std::back_inserter(aIxs), [](const std::tuple<IndexQuad, double> &v) { return std::get<0>(v); });
    aIxs.push_back(aIxChi1);
    aIxs.push_back(aIxChi2);

    auto C5Pfinal = getC5PrimeCoords(firstResidue);
    auto O3Pfinal = getO3PrimeCoords(secondResidue);

    realign(C5Pinitial, O3Pinitial, step);

    mol->DeleteSelection(hSel);
}

bool is_nucleotide(const mmdb::Residue *residue) {
    auto name = trim(residue->name);
    return std::find(NUCLEOTIDES.cbegin(), NUCLEOTIDES.cend(), name) != NUCLEOTIDES.cend();
}

double measure_NtC(mmdb::Manager *mol, NtC::Parameter param) {
    /* We expect to get the two residues that make up the NtC step and nothing else */
    mmdb::PChain chain = mol->GetModel(1)->GetChain(0);
    mmdb::PResidue firstResidue = chain->GetResidue(0);
    mmdb::PResidue secondResidue = chain->GetResidue(1);

    mmdb::PPAtom allAtoms = nullptr;
    int nAllAtoms = 0;
    mol->GetAtomTable(allAtoms, nAllAtoms);

    Step step{filterAtoms(allAtoms, nAllAtoms, firstResidue->GetSeqNum()), firstResidue, secondResidue};

    switch (param) {
    case NtC::DELTA_1:
        return measureTorsionAngle(BACKBONE_ATOMQUADS[0], step);
    case NtC::EPSILON_1:
        return measureTorsionAngle(BACKBONE_ATOMQUADS[1], step);
    case NtC::ZETA_1:
        return measureTorsionAngle(BACKBONE_ATOMQUADS[2], step);
    case NtC::ALPHA_2:
        return measureTorsionAngle(BACKBONE_ATOMQUADS[3], step);
    case NtC::BETA_2:
        return measureTorsionAngle(BACKBONE_ATOMQUADS[4], step);
    case NtC::GAMMA_2:
        return measureTorsionAngle(BACKBONE_ATOMQUADS[5], step);
    case NtC::DELTA_2:
        return measureTorsionAngle(BACKBONE_ATOMQUADS[6], step);
    case NtC::CHI_1:
    {
        auto quad = NtC::chiAtomQuad(firstResidue, Atom::FIRST);
        return measureTorsionAngle(quad, step);
    }
    case NtC::CHI_2:
    {
        auto quad = NtC::chiAtomQuad(secondResidue, Atom::SECOND);
        return measureTorsionAngle(quad, step);
    }
    case NtC::CC:
    {
        auto quad1 = NtC::chiAtomQuad(firstResidue, Atom::FIRST);
        auto quad2 = NtC::chiAtomQuad(secondResidue, Atom::SECOND);
        return measureBondLength(quad1[1], quad2[1], step);
    }
    case NtC::NN:
    {
        auto quad1 = NtC::chiAtomQuad(firstResidue, Atom::FIRST);
        auto quad2 = NtC::chiAtomQuad(secondResidue, Atom::SECOND);
        return measureBondLength(quad1[2], quad2[2], step);
    }
    case NtC::MU:
    {
        auto quad1 = NtC::chiAtomQuad(firstResidue, Atom::FIRST);
        auto quad2 = NtC::chiAtomQuad(secondResidue, Atom::SECOND);
        AtomQuad muQuad{quad1[2], quad1[1], quad2[1], quad2[2]};
        return measureTorsionAngle(muQuad, step);
    }
    }

    return 0;
}

} // namespace IBT
