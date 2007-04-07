#include <iostream>

int ARIADNE_TEST_FAILURES=0;

// This needs to be a function since we do not want to evaluate the result twice,
// and can't store it in a variable since we don't know it's type.
template<class R, class ER> 
bool 
_ariadne_test(std::ostream& os, const R& r, const ER& er) {
  os << r; return (r==er);
}

#define ARIADNE_ASSERT(command) \
{ \
  std::string command_string=#command; \
  std::cout << command_string << ": " << std::flush; \
  bool result = (command); \
  if(result) { \
    std::cout << "true\n" << std::endl; \
  } else { \
    ++ARIADNE_TEST_FAILURES; \
    std::cout << "false\n" << std::endl; \
    std::cerr << __FILE__ << ":" << __LINE__ << ": " << __PRETTY_FUNCTION__ << ": Assertion`" << command_string << "' failed." << std::endl; \
  } \
}\


#define ARIADNE_TEST(command) \
{ \
  std::string command_string=#command; \
  std::cout << command_string << ": " << std::flush; \
  std::cout << command << "\n" << std::endl;          \
}\

#define ARIADNE_TEST_ASSERT(command,expected) \
{ \
  std::cout << #command << ": " << std::flush; \
  bool ok = ariadne_test(std::cout,command,expected); \
  if(ok) { \
    std::cout << "\n" << std::endl; \
  } else { \
    ++ARIADNE_TEST_FAILURES; \
    std::cout << "  (expected: " << #expected << ")\n" << std::endl;        \
    std::cerr << __FILE__ << ":" << __LINE__ << ": " << __PRETTY_FUNCTION__ << ": Assertion`" << #command << " == " << #expected << "' failed." << std::endl; \
  } \
}\

