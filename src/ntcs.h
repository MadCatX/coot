// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef NTCS_H_
#define NTCS_H_

#include "geometry/protein-geometry.hh"

#include <mmdb2/mmdb_manager.h>

#include <array>
#include <string>
#include <vector>

namespace IBT {

class Atom {
public:
    enum WhichResidue {
        FIRST,
        SECOND
    };

    Atom();
    Atom(std::string name, WhichResidue whichResidue);

    std::string name;
    WhichResidue whichResidue;
};

using AtomQuad = std::array<Atom, 4>;

class Torsion {
public:
    Torsion();
    Torsion(AtomQuad quad, double angle);

    AtomQuad quad;
    double angle;
};

class NtC {
public:
    enum Parameter {
        DELTA_1,
        EPSILON_1,
        ZETA_1,
        ALPHA_2,
        BETA_2,
        GAMMA_2,
        DELTA_2,
        CHI_1,
        CHI_2,
        CC,
        NN,
        MU
    };

    NtC() = delete;
    NtC(
        std::string name,
        double delta_1, double epsilon_1, double zeta_1,
        double alpha_2, double beta_2, double gamma_2, double delta_2,
        double chi_1, double chi_2,
        double cc, double nn, double mu,
        double nu0first, double nu1first, double nu2first, double nu3first, double nu4first,
        double nu0second, double nu1second, double nu2second, double nu3second, double nu4second
    );

    const std::string name;

    const double delta_1;
    const double epsilon_1;
    const double zeta_1;
    const double alpha_2;
    const double beta_2;
    const double gamma_2;
    const double delta_2;

    const double chi_1;
    const double chi_2;

    const double cc;
    const double nn;
    const double mu;

    const double nu0first;
    const double nu1first;
    const double nu2first;
    const double nu3first;
    const double nu4first;
    const double nu0second;
    const double nu1second;
    const double nu2second;
    const double nu3second;
    const double nu4second;

    const std::string NtCClass;

    const std::vector<Torsion> & backboneTorsions() const;

    Torsion chi1Torsion(const mmdb::Residue *residue) const;
    Torsion chi2Torsion(const mmdb::Residue *residue) const;

    static
    AtomQuad chiAtomQuad(const mmdb::Residue *residue, Atom::WhichResidue whichResidue);

    static
    std::string nameToClass(const std::string &name);

private:
    enum RingType {
        PURINE,
        PYRIMIDINE,
    };

    const std::vector<Torsion> m_backboneTorsions;

    static
    RingType getRingType(const std::string &compound);
};

extern const std::array<std::string, 7> NTC_CLASSES;
extern const std::vector<NtC> NTCS;

void apply_NtC(mmdb::Manager *mol, coot::protein_geometry *geom, const NtC &ntc);
double measure_NtC(mmdb::Manager *mol, NtC::Parameter param);
bool is_nucleotide(const mmdb::Residue *residue);

} // namespace IBT

#endif // NTCS_H_
