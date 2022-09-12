// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_UI_NTC_DIALOG_HH
#define _NTC_UI_NTC_DIALOG_HH

#include "../types.hh"

#include <LLKA/llka_connectivity_similarity.h>
#include <LLKA/llka_ntc.h>

#include <functional>
#include <string>
#include <vector>

struct LLKA_ClassifiedStep;
struct LLKA_StepMetrics;

struct NtCDialog;
using OnDisplayedNtCChanged = std::function<void(NtCDialog *, LLKA_NtC)>;
using OnNtCDialogAccepted = std::function<void(NtCDialog *, LLKA_NtC)>;
using OnNtCDialogRejected = std::function<void(NtCDialog *, LLKA_NtC)>;

struct NtCDialogOptions {
    NtCDialogOptions();

    OnDisplayedNtCChanged onDisplayedNtCChanged;
    OnNtCDialogAccepted onAccepted;
    OnNtCDialogRejected onRejected;

    int connSimilDlgWidth;
    int connSimilDlgHeight;
};

void ntc_dialog_destroy(NtCDialog *dlg);
void ntc_dialog_display_classification(NtCDialog *dlg, const LLKA_ClassifiedStep &classified);
void ntc_dialog_display_differences(NtCDialog *dlg, const LLKA_StepMetrics &differences);
void ntc_dialog_display_rmsd(NtCDialog *dlg, double rmsd);
NtCDialogOptions ntc_dialog_get_options(NtCDialog *dlg);
bool ntc_dialog_is_valid(NtCDialog *dlg);
NtCDialog * ntc_dialog_make(const NtCDialogOptions &options);
void ntc_dialog_show(NtCDialog *dlg);
void ntc_dialog_update_similarities(NtCDialog *dlg, std::vector<NtCSimilarity> similarities);

#endif //_NTC_UI_NTC_DIALOG_HH
