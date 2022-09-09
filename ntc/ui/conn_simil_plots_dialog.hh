// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_UI_CONN_SIMIL_DIALOG_HH
#define _NTC_UI_CONN_SIMIL_DIALOG_HH

#include "../types.hh"

#include <functional>
#include <string>
#include <vector>

struct NtCConnSimilPlotsDialog;

using OnSimilaritySelected = std::function<void (NtCSimilarity)>;
using OnWidgetResized = std::function<void (int, int)>;

void ntc_csp_dialog_destroy(NtCConnSimilPlotsDialog *dlg);
bool ntc_csp_dialog_is_valid(NtCConnSimilPlotsDialog *dlg);
NtCConnSimilPlotsDialog * ntc_csp_dialog_make(OnSimilaritySelected onSimilaritySelected = nullptr, OnWidgetResized onWidgetResized = nullptr);
void ntc_csp_dialog_show(NtCConnSimilPlotsDialog *dlg, int width = 0, int height = 0);
void ntc_csp_dialog_update_connectivities(NtCConnSimilPlotsDialog *dlg, const std::vector<LLKA_Connectivity> &connectivities);
void ntc_csp_dialog_update_similarities(NtCConnSimilPlotsDialog *dlg, const std::vector<NtCSimilarity> &similarities);

#endif // _NTC_UI_CONN_SIMIL_DIALOG_HH
