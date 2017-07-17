#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class std::string+;
#pragma link C++ class std::set<std::string>+;
#pragma link C++ class std::vector<std::string>+;
#pragma link C++ class std::vector<double>+;

#pragma link C++ namespace sp+;
#pragma link C++ class sp::LoadSP+;
#pragma link C++ class sp::Parameter+;
#pragma link C++ class sp::Response+;
#pragma link C++ class std::vector<sp::Parameter>+;
#pragma link C++ class std::vector<sp::Response>+;
#pragma link C++ class sp::SPIO+;
#pragma link C++ class sp::SPManager+;
#pragma link C++ class sp::UnfoldAlgoBase+;
#pragma link C++ class sp::UnfoldAlgoDAgnostini+;
#

#endif
