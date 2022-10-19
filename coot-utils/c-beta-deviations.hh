
#include <map>
#include <clipper/core/coords.h>

// We need to put this here before anything from mmdb2 gets included
#include "compat/coot-inhibit-win32-api.h"
#include <mmdb2/mmdb_manager.h>
#include "mini-mol/atom-quads.hh"

namespace coot {

   class c_beta_deviation_t {

      // there can be several c_beta deviations per residue (alt confs)

   public:
      c_beta_deviation_t(mmdb::Atom *at_in, const clipper::Coord_orth &pos_ideal_in,  double dist_in) :
	 at(at_in), pos_ideal(pos_ideal_in), dist(dist_in) { }
      c_beta_deviation_t() : at(nullptr), pos_ideal(clipper::Coord_orth(-1,-1,-1)), dist(0.0) {} // needed for map
      mmdb::Atom *at;
      clipper::Coord_orth pos_ideal;
      double dist;
   };


   std::map<mmdb::Residue *, std::map<std::string, c_beta_deviation_t> >
   get_c_beta_deviations(mmdb::Manager *mol);
   // maybe other, more fine-grained arguments needed?

   std::map<std::string, c_beta_deviation_t> get_c_beta_deviations(mmdb::Residue *residue_p);

   clipper::Coord_orth make_CB_ideal_pos(const atom_quad &q, const std::string &res_name);
   
}
