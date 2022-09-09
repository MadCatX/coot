// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_UI_NOTIFICATIONS_HH
#define _NTC_UI_NOTIFICATIONS_HH

#include <string>

void ntc_notify_error(const std::string &text);
void ntc_notify_info(const std::string &text);
void ntc_notify_warning(const std::string &text);

#endif // _NTC_UI_NOTIFICATIONS_HH
