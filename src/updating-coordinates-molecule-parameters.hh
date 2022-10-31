
#ifndef UPDATING_COORDINATES_MOLECULE_PARAMETERS_T
#define UPDATING_COORDINATES_MOLECULE_PARAMETERS_T

#include <cstddef>
#include<string>

// This class is for reading the output of Refmac
//
class updating_coordinates_molecule_parameters_t {
public:
   int imol;
   std::string pdb_file_name;
   uint_least64_t ctime;

   updating_coordinates_molecule_parameters_t() {
      imol = -1;
      ctime = 0;
   }
   explicit updating_coordinates_molecule_parameters_t(const std::string &file_name) : pdb_file_name(file_name) {
      ctime = 0;
      imol = -1; // is this used?
   }
   void update_time(uint_least64_t time) {
      ctime = time;
   }
};

// This class is for updating the difference map when the model has changed.
//
class updating_model_molecule_parameters_t {
public:
   int imol_coords;
   int imol_map_with_data_attached;
   int imol_2fofc_map;  // sigmaA weighted, that is
   int imol_fofc_map; // ditto
   updating_model_molecule_parameters_t() {
      imol_coords = -1;
      imol_map_with_data_attached = -1;
      imol_2fofc_map = -1;
      imol_fofc_map = -1;
   }
   updating_model_molecule_parameters_t(int imol_coords_in, int imol_d, int imol_map_in) :
      imol_coords(imol_coords_in), imol_map_with_data_attached(imol_d), imol_fofc_map(imol_map_in) { imol_2fofc_map = -1;
   }
   updating_model_molecule_parameters_t(int imol_coords_in, int imol_data, int imol_map_2fofc_in, int imol_map_fofc_in) :
      imol_coords(imol_coords_in), imol_map_with_data_attached(imol_data), imol_2fofc_map(imol_map_2fofc_in), imol_fofc_map(imol_map_fofc_in) {}
};

#endif
