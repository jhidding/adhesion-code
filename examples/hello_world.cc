// ~\~ language=C++ filename=examples/hello_world.cc
// ~\~ begin <<adhesion_example.md|examples/hello_world.cc>>[0]
#include <cstdlib>
#include <iostream>

// ~\~ begin <<adhesion_example.md|example-main-function>>[0]
int main(int argc, char **argv) {
  // ~\~ begin <<adhesion_example.md|hello-world>>[0]
  std::cout << "Hello, World!" << std::endl;
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|hello-world>>[1]
  return EXIT_SUCCESS;
  // ~\~ end
}
// ~\~ end
// ~\~ end
