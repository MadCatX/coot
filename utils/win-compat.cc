#include "win-compat.hh"

std::string
coot::uri_to_file_name(const std::string &uri) {

   std::string r = uri;

   // Why this? https://en.wikipedia.org/wiki/File_URI_scheme
   
#ifdef WINDOWS_MINGW
	 r = uri.substr(8);
#else
	 r = uri.substr(7);
#endif

   return r;

}
