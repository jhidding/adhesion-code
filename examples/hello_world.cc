// ~\~ language=C++ filename=examples/hello_world.cc
// ~\~ begin <<adhesion_example.md|examples/hello_world.cc>>[init]
#include <cstdlib>
#include <iostream>

// ~\~ begin <<adhesion_example.md|example-main-function>>[init]
int main(int argc, char **argv) {
  // ~\~ begin <<adhesion_example.md|hello-world>>[init]
  std::cout << "Hello, World!" << std::endl;
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|hello-world>>[1]
  return EXIT_SUCCESS;
  // ~\~ end
}
// ~\~ end
// ~\~ end
