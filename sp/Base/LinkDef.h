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
#pragma link C++ namespace sp::msg+;
#pragma link C++ enum sp::msg::Level_t+;
#pragma link C++ class sp::logger+;
#pragma link C++ class sp::sp_base+;
#pragma link C++ class sp::sperr+;

#endif
