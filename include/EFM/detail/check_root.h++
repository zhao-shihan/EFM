#pragma once

#if not __has_include("TFile.h")
#error "EFM requires ROOT::RIO to function normally"
#endif

#if not __has_include("TNtuple.h")
#error "EFM requires ROOT::Tree to function normally"
#endif
